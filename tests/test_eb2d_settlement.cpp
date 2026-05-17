#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Propped cantilever (fix-pinned) with prescribed vertical settlement at the
// roller end, no distributed/concentrated load. This exercises the elimination
// + RHS correction path for nonzero prescribed displacement.
//
// BCs:
//   node 0     : u_x = 0, u_z = 0, theta_y = 0   (fixed)
//   node nElems: u_x = 0, u_z = delta (prescribed), theta_y free
//
// Analytical EB solution (EI w'''' = 0 with w(0)=w'(0)=0, w(L)=delta, M(L)=0):
//   w(x)  = (3 delta / (2 L^2)) x^2 - (delta / (2 L^3)) x^3
//   M(x)  = EI w''(x) = (3 EI delta / L^2) (1 - x/L)
//   V(x)  = dM/dx     = -3 EI delta / L^3   (constant)
//
// Reactions (in the BeamLib structural sign convention, theta_y = dw/dx):
//   R_A_z =  V_internal(0+)   = -3 EI delta / L^3
//   M_A_y = -M_internal(0+)   = -3 EI delta / L^2
//   R_B_z = -V_internal(L-)   = +3 EI delta / L^3
// (For delta > 0 the right support is lifted, the beam wants to rise everywhere,
//  so the left support must pull DOWN at A and the moment reaction must oppose
//  the positive-slope tendency at A. The "R = -M_internal" / "R = V_internal"
//  pattern is the same one used by test_eb2d_postprocess for the cantilever:
//  applying it to the cantilever there gives R_z = +1000 and M_y = +1000 from
//  V_internal = +1000 and M_internal = -1000.)
int main() {
    using namespace beamlib;

    const int    nElems = 10;
    const double Lbeam  = 2.0;
    const double dx     = Lbeam / nElems;
    const double delta  = 1.0e-3;   // 1 mm settlement at the right support

    BeamModel<EulerBernoulli2D> model;
    model.props.E  = 200e9;
    model.props.A  = 0.01;
    model.props.Iz = 8.333e-6;

    model.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        model.nodes[i].x0 = Vec3(i * dx, 0.0, 0.0);
    }
    model.nodes[0].fixAll();

    // Roller at right end: vertical position prescribed = delta, axial fixed
    // to remove the rigid axial mode, rotation free.
    model.nodes[nElems].fixed[0]      = true;   // u_x = 0
    model.nodes[nElems].fixed[1]      = true;   // u_z = delta (prescribed below)
    model.nodes[nElems].prescribed[1] = delta;

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

    const double EI = model.props.E * model.props.Iz;
    const double L  = Lbeam;
    const double L2 = L * L;
    const double L3 = L * L * L;

    auto w_analytic = [&](double x) {
        return (3.0 * delta / (2.0 * L2)) * x * x
             - (delta / (2.0 * L3)) * x * x * x;
    };

    const double tol_rel = 1e-9;
    const double tol_abs = 1e-12;

    for (int i = 0; i <= nElems; ++i) {
        const double x = i * dx;
        const double w_an = w_analytic(x);
        const double w    = model.nodes[i].totalDof(1);
        const double denom = std::fabs(w_an) > 1e-12 ? std::fabs(w_an) : 1.0;
        const double err   = std::fabs(w - w_an) / denom;
        if (std::fabs(w_an) < 1e-12) {
            if (std::fabs(w) > tol_abs) {
                std::printf("FAIL: node %d expected w=0, got %g\n", i, w);
                return 1;
            }
        } else if (err > tol_rel) {
            std::printf("FAIL: w(x=%g) expected %g, got %g (rel err %g)\n",
                        x, w_an, w, err);
            return 1;
        }
    }

    // Reactions.
    auto reactions = model.computeReactions();
    const double R_A_z = reactions[0].R[1];
    const double M_A_y = reactions[0].R[2];
    const double R_B_z = reactions[nElems].R[1];

    const double R_A_z_expected = -3.0 * EI * delta / L3;
    const double M_A_y_expected = -3.0 * EI * delta / L2;
    const double R_B_z_expected = +3.0 * EI * delta / L3;

    auto checkRel = [&](double got, double expected, const char* tag) {
        const double denom = std::fabs(expected);
        const double err   = std::fabs(got - expected) / denom;
        if (err > tol_rel) {
            std::printf("FAIL [%s]: expected %g, got %g (rel err %g)\n",
                        tag, expected, got, err);
            return 1;
        }
        return 0;
    };

    if (checkRel(R_A_z, R_A_z_expected, "R_A_z")) return 1;
    if (checkRel(M_A_y, M_A_y_expected, "M_A_y")) return 1;
    if (checkRel(R_B_z, R_B_z_expected, "R_B_z")) return 1;

    // Equilibrium sanity: sum of vertical reactions = 0 (no external load).
    if (std::fabs(R_A_z + R_B_z) / std::fabs(R_A_z) > 1e-9) {
        std::printf("FAIL: vertical reaction imbalance %g vs %g\n", R_A_z, R_B_z);
        return 1;
    }

    std::printf("PASS test_eb2d_settlement: R_A_z=%g M_A_y=%g R_B_z=%g, "
                "midspan w=%g (analytic %g)\n",
                R_A_z, M_A_y, R_B_z,
                model.nodes[nElems / 2].totalDof(1),
                w_analytic((nElems / 2) * dx));
    return 0;
}
