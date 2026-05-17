#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

namespace {

int check_near_zero(double v, double tol, const char* tag)
{
    if (std::fabs(v) > tol) {
        std::printf("FAIL [%s]: expected ~0 got %g (tol %g)\n", tag, v, tol);
        return 1;
    }
    return 0;
}

} // namespace

// Pure torsion check for EB3D.
//
//   - cantilever along global x, all 6 DOFs at node 0 fixed
//   - tip torque T_x = +500 at node N
//   - analytical: theta_x(L) = T_x * L / (G * Ix), all other DOFs ~0
//   - reaction at root: -T_x along theta_x (residual minus zero ext load)
int main()
{
    using namespace beamlib;

    const double E   = 200e9;
    const double G   = 80e9;
    const double A   = 0.01;
    const double Iy  = 8.333e-6;
    const double Iz  = 8.333e-6;
    const double Ix  = 16.667e-6;

    const double Lbeam = 1.0;
    const int    nElems = 10;
    const double dx     = Lbeam / nElems;
    const double Tx     = 500.0;

    BeamModel<EulerBernoulli3D> model;
    model.props.E  = E;
    model.props.G  = G;
    model.props.A  = A;
    model.props.Iy = Iy;
    model.props.Iz = Iz;
    model.props.Ix = Ix;

    model.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        model.nodes[i].x0 = Vec3(i * dx, 0.0, 0.0);
    }
    model.nodes[0].fixAll();
    model.nodes[nElems].load[3] = Tx;  // index 3 = theta_x torque

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
        NewtonRaphsonSolver<EulerBernoulli3D>::solveOneStep(model, Fext, cfg);

    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL torsion NR: converged=%d iters=%d res=%g\n",
                    int(res.converged), res.iterations, res.finalResidual);
        return 1;
    }

    const double tx_tip_expected = Tx * Lbeam / (G * Ix);
    const auto& tip = model.nodes[nElems].dof;
    const double ux = tip[0], uy = tip[1], uz = tip[2];
    const double tx = tip[3], ty = tip[4], tz = tip[5];

    const double err_tx = std::fabs(tx - tx_tip_expected)
                          / std::fabs(tx_tip_expected);
    if (err_tx > 1e-9) {
        std::printf("FAIL: tip theta_x = %g, expected %g, rel %g\n",
                    tx, tx_tip_expected, err_tx);
        return 1;
    }

    if (check_near_zero(ux, 1e-14, "tip u_x"))     return 1;
    if (check_near_zero(uy, 1e-14, "tip u_y"))     return 1;
    if (check_near_zero(uz, 1e-14, "tip u_z"))     return 1;
    if (check_near_zero(ty, 1e-14, "tip theta_y")) return 1;
    if (check_near_zero(tz, 1e-14, "tip theta_z")) return 1;

    // Linearly varying theta_x along the span: theta_x(x_k) = T_x * x_k / (G Ix)
    for (int k = 0; k <= nElems; ++k) {
        const double xk = k * dx;
        const double tx_expected = Tx * xk / (G * Ix);
        const double tx_got = model.nodes[k].dof[3];
        if (k == 0) {
            if (std::fabs(tx_got) > 1e-14) {
                std::printf("FAIL: theta_x at root must be 0 (fixed), got %g\n",
                            tx_got);
                return 1;
            }
        } else {
            const double err =
                std::fabs(tx_got - tx_expected) / std::fabs(tx_expected);
            if (err > 1e-9) {
                std::printf("FAIL: theta_x at node %d = %g, expected %g (rel %g)\n",
                            k, tx_got, tx_expected, err);
                return 1;
            }
        }
    }

    // Reaction recovery: at the root, the only nonzero reaction component is
    // theta_x with value -T_x (FE residual at constrained DOF = R_int - F_ext,
    // and no external load is applied there).
    auto reactions = model.computeReactions();
    const auto& Rroot = reactions[0].R;
    if (check_near_zero(Rroot[0], 1e-9, "reaction R_x")) return 1;
    if (check_near_zero(Rroot[1], 1e-9, "reaction R_y")) return 1;
    if (check_near_zero(Rroot[2], 1e-9, "reaction R_z")) return 1;
    if (check_near_zero(Rroot[4], 1e-9, "reaction M_y")) return 1;
    if (check_near_zero(Rroot[5], 1e-9, "reaction M_z")) return 1;
    const double Rtx = Rroot[3];
    if (std::fabs(Rtx - (-Tx)) / std::fabs(Tx) > 1e-9) {
        std::printf("FAIL: reaction theta_x at root = %g, expected %g\n",
                    Rtx, -Tx);
        return 1;
    }

    std::printf("PASS test_eb3d_torsion: theta_x(L)=%g, root reaction=%g, iters=%d\n",
                tx, Rtx, res.iterations);
    return 0;
}
