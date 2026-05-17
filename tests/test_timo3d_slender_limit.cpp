#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Element/Timoshenko3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Timoshenko 3D -> EB3D consistency in the slender limit, exercised in BOTH
// bending planes simultaneously. A slender beam (L/h = 1000) loaded with
// independent F_y and F_z must produce Timoshenko-3D tip displacements that
// agree with the direct EB3D solve to relative ~1e-6 in u_y and u_z, and to
// the same accuracy in theta_y and theta_z.
int main()
{
    using namespace beamlib;

    const double E       = 200e9;
    const double nu      = 0.3;
    const double G       = E / (2.0 * (1.0 + nu));
    const double rho     = 7800.0;

    const double Lbeam   = 1.0;
    // Square cross-section so BOTH bending planes are equally slender
    // (L/h = L/b = 1000). With an asymmetric section the two L-over-section-
    // dimension ratios differ, and one bending plane would not actually be in
    // the slender regime, masking the Timoshenko->EB convergence claim.
    const double h       = Lbeam / 1000.0;
    const double b       = h;
    const double A       = b * h;
    const double Iy      = b * h * h * h / 12.0;
    const double Iz      = h * b * b * b / 12.0;
    const double Ix      = Iy + Iz;
    const double kappa_y = 5.0 / 6.0;
    const double kappa_z = 5.0 / 6.0;

    const int    nElems = 40;
    const double dx     = Lbeam / nElems;
    const double Fy     = -1e-6;       // small load: tiny but slender section
    const double Fz     = -1e-6;

    auto buildEB = [&]() {
        BeamModel<EulerBernoulli3D> m;
        m.props.E = E; m.props.G = G; m.props.rho = rho;
        m.props.A = A; m.props.Iy = Iy; m.props.Iz = Iz; m.props.Ix = Ix;
        m.nodes.resize(nElems + 1);
        for (int i = 0; i <= nElems; ++i) m.nodes[i].x0 = Vec3(i * dx, 0, 0);
        m.nodes[0].fixAll();
        m.nodes[nElems].load[1] = Fy;
        m.nodes[nElems].load[2] = Fz;
        for (int e = 0; e < nElems; ++e) {
            ElementConn c; c.nodeA = e; c.nodeB = e + 1;
            m.elements.push_back(c);
        }
        m.buildDofMap();
        return m;
    };
    auto buildTimo = [&]() {
        BeamModel<Timoshenko3D> m;
        m.props.E = E; m.props.G = G; m.props.rho = rho;
        m.props.A = A; m.props.Iy = Iy; m.props.Iz = Iz; m.props.Ix = Ix;
        m.props.kappa_y = kappa_y; m.props.kappa_z = kappa_z;
        m.nodes.resize(nElems + 1);
        for (int i = 0; i <= nElems; ++i) m.nodes[i].x0 = Vec3(i * dx, 0, 0);
        m.nodes[0].fixAll();
        m.nodes[nElems].load[1] = Fy;
        m.nodes[nElems].load[2] = Fz;
        for (int e = 0; e < nElems; ++e) {
            ElementConn c; c.nodeA = e; c.nodeB = e + 1;
            m.elements.push_back(c);
        }
        m.buildDofMap();
        return m;
    };

    auto mEB = buildEB();
    NRResult rEB =
        NewtonRaphsonSolver<EulerBernoulli3D>::solveOneStep(mEB, mEB.getExternalForceVector());
    if (!rEB.converged) { std::printf("FAIL EB3D NR\n"); return 1; }

    auto mT = buildTimo();
    NRResult rT =
        NewtonRaphsonSolver<Timoshenko3D>::solveOneStep(mT, mT.getExternalForceVector());
    if (!rT.converged) { std::printf("FAIL Timo3D NR\n"); return 1; }

    const auto& tipEB = mEB.nodes[nElems].dof;
    const auto& tipT  = mT .nodes[nElems].dof;

    auto check_rel = [](double a, double b, double tol, const char* tag) {
        const double denom = std::fabs(b) > 1e-30 ? std::fabs(b) : 1.0;
        const double err   = std::fabs(a - b) / denom;
        if (err > tol) {
            std::printf("FAIL [%s]: Timo %g vs EB %g (rel %g, tol %g)\n",
                        tag, a, b, err, tol);
            return 1;
        }
        return 0;
    };

    const double tol = 1e-6;
    if (check_rel(tipT[1], tipEB[1], tol, "u_y slender")) return 1;
    if (check_rel(tipT[2], tipEB[2], tol, "u_z slender")) return 1;
    if (check_rel(tipT[4], tipEB[4], tol, "theta_y slender")) return 1;
    if (check_rel(tipT[5], tipEB[5], tol, "theta_z slender")) return 1;

    std::printf("PASS test_timo3d_slender_limit: L/h=%.0f, "
                "u_y_T/EB=%.6e/%.6e (rel %.2e), u_z_T/EB=%.6e/%.6e (rel %.2e)\n",
                Lbeam / h,
                tipT[1], tipEB[1], std::fabs(tipT[1] - tipEB[1]) / std::fabs(tipEB[1]),
                tipT[2], tipEB[2], std::fabs(tipT[2] - tipEB[2]) / std::fabs(tipEB[2]));
    return 0;
}
