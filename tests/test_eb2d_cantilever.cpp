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

    NRConfig cfg;
    cfg.tol     = 1e-7;
    cfg.maxIter = 5;
    NRResult res = NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(model, Fext, cfg);

    if (!res.converged) {
        std::printf("FAIL: NR did not converge. iters=%d residual=%g\n",
                    res.iterations, res.finalResidual);
        return 1;
    }
    if (res.iterations != 1) {
        std::printf("FAIL: linear problem should converge in exactly 1 correction step, got %d\n",
                    res.iterations);
        return 1;
    }

    const double EIz          = model.props.E * model.props.Iz;
    const double uz_expected  = Fz * Lbeam * Lbeam * Lbeam / (3.0 * EIz);
    const double thy_expected = Fz * Lbeam * Lbeam / (2.0 * EIz);

    const double uz  = model.nodes[nElems].dof[1];
    const double thy = model.nodes[nElems].dof[2];

    const double tol    = 1e-9;
    const double err_uz  = std::fabs(uz  - uz_expected)  / std::fabs(uz_expected);
    const double err_thy = std::fabs(thy - thy_expected) / std::fabs(thy_expected);

    if (err_uz > tol) {
        std::printf("FAIL: tip u_z relative error %g (got %g, expected %g)\n",
                    err_uz, uz, uz_expected);
        return 1;
    }
    if (err_thy > tol) {
        std::printf("FAIL: tip theta_y relative error %g (got %g, expected %g)\n",
                    err_thy, thy, thy_expected);
        return 1;
    }
    if (uz >= 0.0 || thy >= 0.0) {
        std::printf("FAIL: expected both u_z and theta_y negative, got u_z=%g theta_y=%g\n",
                    uz, thy);
        return 1;
    }

    // Tip u_x should be zero (no axial load)
    const double ux_tip = model.nodes[nElems].dof[0];
    if (std::fabs(ux_tip) > 1e-14) {
        std::printf("FAIL: tip u_x should be ~0, got %g\n", ux_tip);
        return 1;
    }

    std::printf("PASS test_eb2d_cantilever: u_z=%g (rel err %g), theta_y=%g (rel err %g), iters=%d\n",
                uz, err_uz, thy, err_thy, res.iterations);
    return 0;
}
