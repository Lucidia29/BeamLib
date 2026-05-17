#include <BeamLib/Element/Timoshenko3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Pure torsion check for Timoshenko 3D. St. Venant torsion is unchanged from
// EB3D (no shear strain in the torsion block) so the result must match the
// EB3D torsion behavior exactly.
//
//   tip torque T_x = +500 -> theta_x(L) = T_x L / (G I_x)
//   all other DOFs ~ 0
//   reaction at root theta_x = -T_x
int main()
{
    using namespace beamlib;

    const double E   = 200e9;
    const double G   = 80e9;
    const double rho = 7800.0;
    const double A   = 0.01;
    const double Iy  = 8.0e-6;
    const double Iz  = 1.6e-5;
    const double Ix  = 2.4e-5;

    const double Lbeam = 1.0;
    const int    nElems = 10;
    const double dx     = Lbeam / nElems;
    const double Tx     = 500.0;

    BeamModel<Timoshenko3D> model;
    model.props.E       = E;
    model.props.G       = G;
    model.props.rho     = rho;
    model.props.A       = A;
    model.props.Iy      = Iy;
    model.props.Iz      = Iz;
    model.props.Ix      = Ix;
    model.props.kappa_y = 5.0 / 6.0;
    model.props.kappa_z = 5.0 / 6.0;

    model.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        model.nodes[i].x0 = Vec3(i * dx, 0.0, 0.0);
    }
    model.nodes[0].fixAll();
    model.nodes[nElems].load[3] = Tx;
    for (int e = 0; e < nElems; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        model.elements.push_back(c);
    }

    model.buildDofMap();
    NRConfig cfg; cfg.tol = 1e-7; cfg.maxIter = 5;
    NRResult res =
        NewtonRaphsonSolver<Timoshenko3D>::solveOneStep(model, model.getExternalForceVector(), cfg);
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL torsion NR: converged=%d iters=%d\n",
                    int(res.converged), res.iterations);
        return 1;
    }

    const double tx_expect = Tx * Lbeam / (G * Ix);
    const auto& tip = model.nodes[nElems].dof;
    if (std::fabs(tip[3] - tx_expect) / std::fabs(tx_expect) > 1e-9) {
        std::printf("FAIL tip theta_x: got %g expected %g\n", tip[3], tx_expect);
        return 1;
    }
    // Negligible bending and axial response
    auto near0 = [](double v, double tol, const char* tag) {
        if (std::fabs(v) > tol) {
            std::printf("FAIL [%s]: expected ~0 got %g (tol %g)\n", tag, v, tol);
            return 1;
        }
        return 0;
    };
    if (near0(tip[0], 1e-14, "tip u_x")) return 1;
    if (near0(tip[1], 1e-14, "tip u_y")) return 1;
    if (near0(tip[2], 1e-14, "tip u_z")) return 1;
    if (near0(tip[4], 1e-14, "tip theta_y")) return 1;
    if (near0(tip[5], 1e-14, "tip theta_z")) return 1;

    // Linear theta_x distribution
    for (int k = 1; k <= nElems; ++k) {
        const double xk = k * dx;
        const double tx_an = Tx * xk / (G * Ix);
        const double err = std::fabs(model.nodes[k].dof[3] - tx_an) / std::fabs(tx_an);
        if (err > 1e-9) {
            std::printf("FAIL theta_x at node %d: rel err %g\n", k, err);
            return 1;
        }
    }

    // Root reaction along theta_x = -T_x
    auto rxs = model.computeReactions();
    const double Rtx = rxs[0].R[3];
    if (std::fabs(Rtx - (-Tx)) / Tx > 1e-9) {
        std::printf("FAIL root reaction: %g vs expected %g\n", Rtx, -Tx);
        return 1;
    }

    std::printf("PASS test_timo3d_torsion: theta_x(L)=%g, root reaction=%g, iters=%d\n",
                tip[3], Rtx, res.iterations);
    return 0;
}
