#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Element/Timoshenko2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Slender-beam limit: as L/h -> infinity the Timoshenko shear contribution
// becomes negligible and Timoshenko2D must agree with EB2D.
//
// Setup: a very slender rectangular cantilever, L/h = 1000, well-refined mesh
// (40 elements). For this geometry the shear-deflection fraction is
//   3 EI / (kappa GA L^2)  ~  0.65/kappa * (h/L)^2  ~  10^{-7}
// so Timoshenko and EB should agree to relative ~1e-6 (slack to allow for the
// finite shear correction itself).
int main()
{
    using namespace beamlib;

    const double E   = 200e9;
    const double nu  = 0.3;
    const double G   = E / (2.0 * (1.0 + nu));
    const double rho = 7800.0;

    const double Lbeam   = 1.0;
    const double h       = Lbeam / 1000.0;          // very slender
    const double b       = 0.01;
    const double A       = b * h;
    const double Iz      = b * h * h * h / 12.0;
    const double kappa_z = 5.0 / 6.0;

    const int    nElems = 40;
    const double dx     = Lbeam / nElems;
    const double Fz     = -1.0;                     // small load (very flexible beam)

    auto buildModel = [&](auto& model) {
        model.props.E       = E;
        model.props.G       = G;
        model.props.rho     = rho;
        model.props.A       = A;
        model.props.Iz      = Iz;
        model.props.kappa_z = kappa_z;
        model.nodes.resize(nElems + 1);
        for (int i = 0; i <= nElems; ++i) {
            model.nodes[i].x0 = Vec3(i * dx, 0.0, 0.0);
        }
        model.nodes[0].fixAll();
        model.nodes[nElems].load[1] = Fz;
        for (int e = 0; e < nElems; ++e) {
            ElementConn c;
            c.nodeA = e;
            c.nodeB = e + 1;
            model.elements.push_back(c);
        }
        model.buildDofMap();
    };

    BeamModel<EulerBernoulli2D> mEB;
    buildModel(mEB);
    NRResult rEB =
        NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(mEB, mEB.getExternalForceVector());
    if (!rEB.converged) { std::printf("FAIL: EB NR failed\n"); return 1; }
    const double uz_EB   = mEB.nodes[nElems].dof[1];
    const double th_EB   = mEB.nodes[nElems].dof[2];

    BeamModel<Timoshenko2D> mTimo;
    buildModel(mTimo);
    NRResult rTimo =
        NewtonRaphsonSolver<Timoshenko2D>::solveOneStep(mTimo, mTimo.getExternalForceVector());
    if (!rTimo.converged) { std::printf("FAIL: Timo NR failed\n"); return 1; }
    const double uz_Timo = mTimo.nodes[nElems].dof[1];
    const double th_Timo = mTimo.nodes[nElems].dof[2];

    const double tol     = 1e-6;
    const double err_uz  = std::fabs(uz_Timo - uz_EB) / std::fabs(uz_EB);
    const double err_th  = std::fabs(th_Timo - th_EB) / std::fabs(th_EB);

    if (err_uz > tol) {
        std::printf("FAIL: slender u_z: Timo=%g EB=%g rel err=%g (tol %g)\n",
                    uz_Timo, uz_EB, err_uz, tol);
        return 1;
    }
    if (err_th > tol) {
        std::printf("FAIL: slender theta_y: Timo=%g EB=%g rel err=%g (tol %g)\n",
                    th_Timo, th_EB, err_th, tol);
        return 1;
    }

    // Sanity: shear contribution magnitude should be ~ 1e-7 of total deflection
    // for this geometry, much smaller than the test tolerance 1e-6.
    const double EI    = E * Iz;
    const double GAk   = kappa_z * G * A;
    const double frac  = (3.0 * EI / (GAk * Lbeam * Lbeam));
    std::printf("PASS test_timo2d_slender_limit: L/h=%.0f, u_z_EB=%g, u_z_Timo=%g, "
                "rel diff=%.2e (shear frac %.2e, tol %g)\n",
                Lbeam / h, uz_EB, uz_Timo, err_uz, frac, tol);
    return 0;
}
