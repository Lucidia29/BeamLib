#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Element/Timoshenko2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Shear-locking check: very slender beam (L/h = 1000) on a COARSE mesh
// (only 4 elements). A locking element would over-stiffen and produce a
// tip deflection orders of magnitude smaller than the analytical value.
//
// The Phi-corrected closed-form Timoshenko element must match the EB
// analytical solution to relative ~1e-6 even with 4 elements -- the same
// accuracy as the well-refined slender-limit test but with 10x coarser mesh.
// This is the canonical justification that the chosen formulation is
// locking-free; PROJECT_SPEC sec 13.4 specifies this test family.
int main()
{
    using namespace beamlib;

    const double E   = 200e9;
    const double nu  = 0.3;
    const double G   = E / (2.0 * (1.0 + nu));
    const double rho = 7800.0;

    const double Lbeam   = 1.0;
    const double h       = Lbeam / 1000.0;       // extreme slenderness
    const double b       = 0.01;
    const double A       = b * h;
    const double Iz      = b * h * h * h / 12.0;
    const double kappa_z = 5.0 / 6.0;

    const int    nElems = 4;                     // intentionally coarse
    const double dx     = Lbeam / nElems;
    const double Fz     = -1.0;

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
    NRResult res =
        NewtonRaphsonSolver<Timoshenko2D>::solveOneStep(model, model.getExternalForceVector());
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL: NR failed (iters=%d res=%g)\n",
                    res.iterations, res.finalResidual);
        return 1;
    }

    // EB cantilever analytical solution (with the small Timoshenko shear
    // correction for completeness, since L/h is finite).
    const double EI   = E * Iz;
    const double GAk  = kappa_z * G * A;
    const double uz_an = Fz * Lbeam * Lbeam * Lbeam / (3.0 * EI)
                       + Fz * Lbeam / GAk;
    const double th_an = Fz * Lbeam * Lbeam / (2.0 * EI);
    const double uz    = model.nodes[nElems].dof[1];
    const double th    = model.nodes[nElems].dof[2];

    // Sanity: a locked element would give a tip deflection orders of
    // magnitude smaller than analytical. We REQUIRE the result to match
    // analytical to relative 1e-6 (same accuracy as the refined 40-element
    // slender-limit test). If the formulation locked, this test would fail
    // by many orders of magnitude.
    const double tol     = 1e-6;
    const double err_uz  = std::fabs(uz - uz_an) / std::fabs(uz_an);
    const double err_th  = std::fabs(th - th_an) / std::fabs(th_an);

    if (err_uz > tol) {
        std::printf("FAIL [LOCKING SUSPECTED]: u_z=%g, analytical %g, rel err %g (tol %g)\n",
                    uz, uz_an, err_uz, tol);
        std::printf("   Ratio FEM/analytical = %.6f (1.0 = no locking; <<1 = locking)\n",
                    uz / uz_an);
        return 1;
    }
    if (err_th > tol) {
        std::printf("FAIL [LOCKING SUSPECTED]: theta_y=%g, analytical %g, rel err %g\n",
                    th, th_an, err_th);
        return 1;
    }

    std::printf("PASS test_timo2d_shear_locking: L/h=%.0f, %d elements, "
                "u_z=%g vs analytical %g (rel %.2e). Locking-free verified.\n",
                Lbeam / h, nElems, uz, uz_an, err_uz);
    return 0;
}
