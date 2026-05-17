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

    const double My = 500.0;
    model.nodes[nElems].load[2] = My;

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
    NRResult res = NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(model, Fext, cfg);

    if (!res.converged) {
        std::printf("FAIL: NR did not converge\n");
        return 1;
    }
    if (res.iterations != 1) {
        std::printf("FAIL: linear problem should converge in 1 correction step, got %d\n",
                    res.iterations);
        return 1;
    }

    const double EI            = model.props.E * model.props.Iz;
    const double uz_tip_exp    = My * Lbeam * Lbeam / (2.0 * EI);
    const double thy_tip_exp   = My * Lbeam / EI;
    const double uz_tip        = model.nodes[nElems].dof[1];
    const double thy_tip       = model.nodes[nElems].dof[2];

    const double tol = 1e-9;
    const double err_uz_tip  = std::fabs(uz_tip  - uz_tip_exp)  / std::fabs(uz_tip_exp);
    const double err_thy_tip = std::fabs(thy_tip - thy_tip_exp) / std::fabs(thy_tip_exp);
    if (err_uz_tip > tol) {
        std::printf("FAIL u_z tip: %g vs %g (rel %g)\n", uz_tip, uz_tip_exp, err_uz_tip);
        return 1;
    }
    if (err_thy_tip > tol) {
        std::printf("FAIL theta_y tip: %g vs %g (rel %g)\n", thy_tip, thy_tip_exp, err_thy_tip);
        return 1;
    }
    if (uz_tip <= 0 || thy_tip <= 0) {
        std::printf("FAIL signs: positive M_y should give positive u_z, theta_y. got u_z=%g theta_y=%g\n",
                    uz_tip, thy_tip);
        return 1;
    }

    // Verify quadratic deflection and linear rotation along the span
    for (int i = 1; i <= nElems; ++i) {
        const double xi             = i * dx;
        const double uz_expected_i  = My * xi * xi / (2.0 * EI);
        const double thy_expected_i = My * xi / EI;
        const double uz_i           = model.nodes[i].dof[1];
        const double thy_i          = model.nodes[i].dof[2];
        const double err_uz   = std::fabs(uz_i  - uz_expected_i)  / std::fabs(uz_tip_exp);
        const double err_thy  = std::fabs(thy_i - thy_expected_i) / std::fabs(thy_tip_exp);
        if (err_uz > tol || err_thy > tol) {
            std::printf("FAIL at node %d: u_z=%g (exp %g, rel %g), theta_y=%g (exp %g, rel %g)\n",
                        i, uz_i, uz_expected_i, err_uz,
                        thy_i, thy_expected_i, err_thy);
            return 1;
        }
    }

    // Axial DOFs should remain zero
    for (int i = 0; i <= nElems; ++i) {
        const double ux_i = model.nodes[i].dof[0];
        if (std::fabs(ux_i) > 1e-14) {
            std::printf("FAIL axial leak at node %d: u_x=%g\n", i, ux_i);
            return 1;
        }
    }

    std::printf("PASS test_eb2d_moment: u_z tip=%g, theta_y tip=%g, iters=%d\n",
                uz_tip, thy_tip, res.iterations);
    return 0;
}
