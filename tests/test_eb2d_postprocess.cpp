#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

int main() {
    using namespace beamlib;

    const int    nElems = 10;
    const double Lbeam  = 1.0;
    const double dx     = Lbeam / nElems;

    BeamModel<EulerBernoulli2D> model;
    model.props.E  = 200e9;
    model.props.A  = 0.01;
    model.props.Iz = 8.333e-6;

    model.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        model.nodes[i].x0 = Vec3(i * dx, 0.0, 0.0);
    }
    model.nodes[0].fixAll();

    const double Fz = -1000.0;
    model.nodes[nElems].load[1] = Fz;

    for (int e = 0; e < nElems; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        model.elements.push_back(c);
    }

    model.buildDofMap();
    VecX Fext = model.getExternalForceVector();
    NRResult res = NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(model, Fext);
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL: NR not converged in 1 step (iters=%d, residual=%g)\n",
                    res.iterations, res.finalResidual);
        return 1;
    }

    // Reactions at the fixed root (node 0). External F_z applied at tip,
    // so equilibrium gives R_z(0) = +1000 and M_y(0) = +1000.
    auto reactions = model.computeReactions();
    if (reactions.size() != static_cast<size_t>(nElems + 1)) {
        std::printf("FAIL: expected %d reaction entries, got %zu\n",
                    nElems + 1, reactions.size());
        return 1;
    }

    const double R_x_root = reactions[0].R[0];
    const double R_z_root = reactions[0].R[1];
    const double M_y_root = reactions[0].R[2];

    const double tol_rel = 1e-9;
    const double tol_abs = 1e-9;

    if (std::fabs(R_x_root) > tol_abs) {
        std::printf("FAIL: R_x(0) should be ~0, got %g\n", R_x_root);
        return 1;
    }
    if (std::fabs(R_z_root - 1000.0) / 1000.0 > tol_rel) {
        std::printf("FAIL: R_z(0) expected 1000, got %g\n", R_z_root);
        return 1;
    }
    if (std::fabs(M_y_root - 1000.0) / 1000.0 > tol_rel) {
        std::printf("FAIL: M_y(0) expected 1000, got %g\n", M_y_root);
        return 1;
    }

    // Free DOFs should report zero reaction.
    for (int i = 1; i <= nElems; ++i) {
        for (int d = 0; d < 3; ++d) {
            if (std::fabs(reactions[i].R[d]) > tol_abs) {
                std::printf("FAIL: reaction at free DOF (node %d, dof %d) = %g\n",
                            i, d, reactions[i].R[d]);
                return 1;
            }
        }
    }

    // Element internal forces. Cantilever with constant +1000 shear; bending
    // moment M(x) = 1000 * (L - x) -> M_y_A of root element = 1000,
    // M_y_B of tip element = 0.
    auto ifs = model.computeInternalForces();
    if (ifs.size() != static_cast<size_t>(nElems)) {
        std::printf("FAIL: expected %d element internal-force entries, got %zu\n",
                    nElems, ifs.size());
        return 1;
    }

    // Check every element: V_z = +1000 at both ends, M_y linear in (L - x).
    for (int e = 0; e < nElems; ++e) {
        const double xA = e * dx;
        const double xB = (e + 1) * dx;
        const double M_A_expected = 1000.0 * (Lbeam - xA);
        const double M_B_expected = 1000.0 * (Lbeam - xB);

        if (std::fabs(ifs[e].V_z_A - 1000.0) / 1000.0 > tol_rel) {
            std::printf("FAIL: V_z_A at element %d expected 1000, got %g\n",
                        e, ifs[e].V_z_A);
            return 1;
        }
        if (std::fabs(ifs[e].V_z_B - 1000.0) / 1000.0 > tol_rel) {
            std::printf("FAIL: V_z_B at element %d expected 1000, got %g\n",
                        e, ifs[e].V_z_B);
            return 1;
        }
        const double M_A_denom = std::fabs(M_A_expected) > 1.0 ? std::fabs(M_A_expected) : 1.0;
        const double M_B_denom = std::fabs(M_B_expected) > 1.0 ? std::fabs(M_B_expected) : 1.0;
        if (std::fabs(ifs[e].M_y_A - M_A_expected) / M_A_denom > tol_rel) {
            std::printf("FAIL: M_y_A at element %d expected %g, got %g\n",
                        e, M_A_expected, ifs[e].M_y_A);
            return 1;
        }
        if (std::fabs(ifs[e].M_y_B - M_B_expected) / M_B_denom > tol_rel) {
            std::printf("FAIL: M_y_B at element %d expected %g, got %g\n",
                        e, M_B_expected, ifs[e].M_y_B);
            return 1;
        }
        if (std::fabs(ifs[e].N_A) > tol_abs || std::fabs(ifs[e].N_B) > tol_abs) {
            std::printf("FAIL: axial force should be ~0 in element %d, got N_A=%g N_B=%g\n",
                        e, ifs[e].N_A, ifs[e].N_B);
            return 1;
        }
    }

    // Specific endpoint sanity (also asserted above implicitly):
    if (std::fabs(ifs[0].V_z_A - 1000.0) > tol_abs ||
        std::fabs(ifs[0].M_y_A - 1000.0) > tol_abs) {
        std::printf("FAIL: root element A end V_z/M_y\n");
        return 1;
    }
    if (std::fabs(ifs[nElems - 1].V_z_B - 1000.0) > tol_abs ||
        std::fabs(ifs[nElems - 1].M_y_B) > tol_abs) {
        std::printf("FAIL: tip element B end V_z/M_y\n");
        return 1;
    }

    std::printf("PASS test_eb2d_postprocess: R_z(0)=%g M_y(0)=%g, "
                "root V_z=%g M_y=%g, tip V_z=%g M_y=%g\n",
                R_z_root, M_y_root,
                ifs[0].V_z_A, ifs[0].M_y_A,
                ifs[nElems - 1].V_z_B, ifs[nElems - 1].M_y_B);
    return 0;
}
