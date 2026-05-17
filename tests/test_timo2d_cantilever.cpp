#include <BeamLib/Element/Timoshenko2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Timoshenko 2D cantilever under tip transverse load.
//
// Analytical (single closed-form, valid for any L/h):
//   u_z(L)    = F_z * L^3 / (3 E I)  +  F_z * L / (kappa * G * A)
//   theta_y(L)= F_z * L^2 / (2 E I)
//
// Geometry / material are chosen so that the shear term is NOT negligible
// (~ 3% of the bending term) so a relative tolerance of 1e-9 actually
// exercises both contributions. This is the main analytical correctness test.
int main()
{
    using namespace beamlib;

    const double E       = 200e9;
    const double G       = 80e9;
    const double rho     = 7800.0;
    const double A       = 0.01;
    const double Iz      = 8.333e-6;
    const double kappa_z = 5.0 / 6.0;

    const double Lbeam = 1.0;
    const int    nElems = 10;
    const double dx     = Lbeam / nElems;
    const double Fz     = -1000.0;

    BeamModel<Timoshenko2D> model;
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
    VecX Fext = model.getExternalForceVector();
    NRConfig cfg;
    cfg.tol     = 1e-7;
    cfg.maxIter = 5;
    NRResult res =
        NewtonRaphsonSolver<Timoshenko2D>::solveOneStep(model, Fext, cfg);

    if (!res.converged) {
        std::printf("FAIL: NR did not converge (iters=%d res=%g)\n",
                    res.iterations, res.finalResidual);
        return 1;
    }
    if (res.iterations != 1) {
        std::printf("FAIL: expected 1 NR correction step for linear problem, got %d\n",
                    res.iterations);
        return 1;
    }

    const double EI  = E * Iz;
    const double GAk = kappa_z * G * A;

    const double uz_bend   = Fz * Lbeam * Lbeam * Lbeam / (3.0 * EI);
    const double uz_shear  = Fz * Lbeam / GAk;
    const double uz_expect = uz_bend + uz_shear;
    const double ty_expect = Fz * Lbeam * Lbeam / (2.0 * EI);

    const double uz = model.nodes[nElems].dof[1];
    const double ty = model.nodes[nElems].dof[2];
    const double ux = model.nodes[nElems].dof[0];

    const double tol     = 1e-9;
    const double err_uz  = std::fabs(uz - uz_expect) / std::fabs(uz_expect);
    const double err_ty  = std::fabs(ty - ty_expect) / std::fabs(ty_expect);

    if (err_uz > tol) {
        std::printf("FAIL: tip u_z=%g, expected %g (bend=%g + shear=%g), rel %g\n",
                    uz, uz_expect, uz_bend, uz_shear, err_uz);
        return 1;
    }
    if (err_ty > tol) {
        std::printf("FAIL: tip theta_y=%g, expected %g, rel %g\n",
                    ty, ty_expect, err_ty);
        return 1;
    }
    if (uz >= 0.0 || ty >= 0.0) {
        std::printf("FAIL: expected negative u_z and theta_y for downward F_z; got u_z=%g theta_y=%g\n",
                    uz, ty);
        return 1;
    }
    if (std::fabs(ux) > 1e-14) {
        std::printf("FAIL: tip u_x should be ~0 (no axial load), got %g\n", ux);
        return 1;
    }

    // Sanity-print the shear contribution magnitude to make sure the test is
    // exercising the Timoshenko-specific term, not coincidentally passing on
    // EB physics alone.
    const double shear_frac = std::fabs(uz_shear) / std::fabs(uz_expect);
    std::printf("PASS test_timo2d_cantilever: u_z=%g (bend=%g, shear=%g, %.2f%% of total), "
                "theta_y=%g, iters=%d\n",
                uz, uz_bend, uz_shear, 100.0 * shear_frac, ty, res.iterations);
    return 0;
}
