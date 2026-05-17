#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Element/Timoshenko2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Deep-beam check: Timoshenko deflection must exceed EB deflection in a
// physically correct, quantitatively predictable way.
//
// For a tip-loaded fixed-free cantilever:
//   u_z_Timo = F_z L^3 / (3 EI) + F_z L / (kappa GA)        ("bend + shear")
//   u_z_EB   = F_z L^3 / (3 EI)                              ("bend only")
//   ratio R  = u_z_Timo / u_z_EB = 1 + 3 EI / (kappa GA L^2)
//
// We choose a deep rectangular cross-section (h = L/4, so L/h = 4):
//   A    = b * h
//   Iz   = b * h^3 / 12
//   kappa = 5/6 (standard rectangular)
//   G    = E / (2(1+nu)), with nu = 0.3 -> G = E / 2.6
// giving 3 EI/(kappa GA L^2) = 3 * E * (b h^3/12) / (kappa * E/2.6 * b h * L^2)
//                             = 3 * 2.6 / (12 * kappa) * (h/L)^2
//                             = 0.65 / kappa * (h/L)^2
// With kappa = 5/6 and h/L = 1/4: ratio - 1 = 0.65 * (6/5) * (1/16)
//                                            = 0.0488 (~ 5% shear contribution)
// Large enough to be unambiguous, small enough that the EB result is still a
// meaningful comparison point.
int main()
{
    using namespace beamlib;

    const double E   = 200e9;
    const double nu  = 0.3;
    const double G   = E / (2.0 * (1.0 + nu));    // = E / 2.6
    const double rho = 7800.0;

    const double Lbeam   = 1.0;
    const double h       = Lbeam / 4.0;            // deep beam: L/h = 4
    const double b       = 0.05;
    const double A       = b * h;
    const double Iz      = b * h * h * h / 12.0;
    const double kappa_z = 5.0 / 6.0;

    const int    nElems = 10;
    const double dx     = Lbeam / nElems;
    const double Fz     = -1000.0;

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

    // EB solve
    BeamModel<EulerBernoulli2D> mEB;
    buildModel(mEB);
    NRResult rEB =
        NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(mEB, mEB.getExternalForceVector());
    if (!rEB.converged) { std::printf("FAIL: EB NR failed\n"); return 1; }
    const double uz_EB = mEB.nodes[nElems].dof[1];

    // Timoshenko solve
    BeamModel<Timoshenko2D> mTimo;
    buildModel(mTimo);
    NRResult rTimo =
        NewtonRaphsonSolver<Timoshenko2D>::solveOneStep(mTimo, mTimo.getExternalForceVector());
    if (!rTimo.converged) { std::printf("FAIL: Timo NR failed\n"); return 1; }
    const double uz_Timo = mTimo.nodes[nElems].dof[1];

    // Sanity: Timoshenko must give a strictly larger deflection magnitude
    // (more negative for Fz < 0).
    if (std::fabs(uz_Timo) <= std::fabs(uz_EB)) {
        std::printf("FAIL: deep-beam Timoshenko deflection magnitude (%g) "
                    "must exceed EB magnitude (%g)\n", uz_Timo, uz_EB);
        return 1;
    }

    // Quantitative: ratio = 1 + 3 EI / (kappa GA L^2)
    const double EI  = E * Iz;
    const double GAk = kappa_z * G * A;
    const double ratio_an = 1.0 + 3.0 * EI / (GAk * Lbeam * Lbeam);
    const double ratio    = uz_Timo / uz_EB;
    const double err      = std::fabs(ratio - ratio_an) / ratio_an;

    if (err > 1e-9) {
        std::printf("FAIL: deflection ratio %g vs analytical %g (rel err %g)\n",
                    ratio, ratio_an, err);
        return 1;
    }

    std::printf("PASS test_timo2d_deep_beam: L/h=%.1f, u_z_EB=%g, u_z_Timo=%g, "
                "ratio=%.6f (analytical %.6f, %.2f%% shear contribution)\n",
                Lbeam / h, uz_EB, uz_Timo, ratio, ratio_an,
                100.0 * (ratio - 1.0) / ratio);
    return 0;
}
