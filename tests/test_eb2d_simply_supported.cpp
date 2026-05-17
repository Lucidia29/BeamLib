#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Simply supported beam under uniform load q.
//
// Boundary conditions:
//   node 0     : u_x = 0, u_z = 0,           theta_y free
//   node nElems: u_x free, u_z = 0,          theta_y free
// (a single u_x constraint at the left support kills the axial rigid-body mode
//  without overconstraining; the right support is a roller.)
//
// Equivalent consistent nodal loads for uniform q on a single Hermite EB
// element of length L_e are
//   F_z at each end = q L_e / 2
//   M_y at each end = +- q L_e^2 / 12   (positive at left, negative at right)
// Assembled over the chain, interior nodes get full q L_e (the two halves sum)
// and the interior moment contributions cancel between adjacent elements,
// leaving only F_z = q L_e / 2 at the two supports and a single moment of
// +q L_e^2 / 12 / -q L_e^2 / 12 at the very ends (which are absorbed by the
// rotational reactions of a fix-end; here the supports have theta_y free so
// the end moments DO NOT cancel and must be applied as external moments).
//
// Manual assembly here matches: F_z = q L_e at every interior node,
// F_z = q L_e / 2 at each end node, plus M_y = +q L_e^2 / 12 at node 0 and
// M_y = -q L_e^2 / 12 at the last node.
int main() {
    using namespace beamlib;

    const int    nElems = 20;
    const double Lbeam  = 2.0;
    const double Le     = Lbeam / nElems;
    // Downward uniform distributed load: p(x) = -q with q > 0. The
    // equivalent nodal loads below carry the explicit minus sign so the
    // physical interpretation is the textbook "gravity load" case.
    const double q      = 1000.0;

    BeamModel<EulerBernoulli2D> model;
    model.props.E  = 200e9;
    model.props.A  = 0.01;
    model.props.Iz = 8.333e-6;

    model.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        model.nodes[i].x0 = Vec3(i * Le, 0.0, 0.0);
    }
    model.nodes[0].fixed[0]      = true;
    model.nodes[0].fixed[1]      = true;
    model.nodes[nElems].fixed[1] = true;

    for (int e = 0; e < nElems; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        model.elements.push_back(c);
    }

    // Manually assemble consistent equivalent nodal loads for p = -q
    // (downward). Per-element contributions: F_z = -q L_e/2 at each end,
    // M_y = -q L_e^2/12 at A, +q L_e^2/12 at B.
    const double half = -0.5 * q * Le;
    const double mom  = -q * Le * Le / 12.0;
    for (int e = 0; e < nElems; ++e) {
        model.nodes[e].load[1]     += half;
        model.nodes[e].load[2]     += mom;
        model.nodes[e + 1].load[1] += half;
        model.nodes[e + 1].load[2] -= mom;
    }

    model.buildDofMap();
    VecX Fext = model.getExternalForceVector();
    NRResult res = NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(model, Fext);
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL: NR not converged in 1 step (iters=%d, residual=%g)\n",
                    res.iterations, res.finalResidual);
        return 1;
    }

    // Midspan deflection. With nElems = 20, midspan is node nElems/2 at x = L/2.
    const double EI       = model.props.E * model.props.Iz;
    const double w_max_an = -5.0 * q * std::pow(Lbeam, 4) / (384.0 * EI);
    const double w_mid    = model.nodes[nElems / 2].dof[1];

    const double tol_rel = 1e-8;
    if (std::fabs(w_mid - w_max_an) / std::fabs(w_max_an) > tol_rel) {
        std::printf("FAIL: midspan deflection expected %g, got %g (rel err %g)\n",
                    w_max_an, w_mid,
                    std::fabs(w_mid - w_max_an) / std::fabs(w_max_an));
        return 1;
    }

    // Reactions: each support carries q L / 2 in the upward (+z) direction
    // to balance the downward distributed load.
    auto reactions = model.computeReactions();
    const double R_A_z = reactions[0].R[1];
    const double R_B_z = reactions[nElems].R[1];
    const double R_expected = 0.5 * q * Lbeam;   // +z, upward

    if (std::fabs(R_A_z - R_expected) / R_expected > 1e-9) {
        std::printf("FAIL: R_A_z expected %g, got %g\n", R_expected, R_A_z);
        return 1;
    }
    if (std::fabs(R_B_z - R_expected) / R_expected > 1e-9) {
        std::printf("FAIL: R_B_z expected %g, got %g\n", R_expected, R_B_z);
        return 1;
    }

    // Axial reaction at A should be zero (no axial load applied).
    if (std::fabs(reactions[0].R[0]) > 1e-8) {
        std::printf("FAIL: R_A_x should be ~0, got %g\n", reactions[0].R[0]);
        return 1;
    }

    std::printf("PASS test_eb2d_simply_supported: w_mid=%g (rel err %g), "
                "R_A=%g R_B=%g\n",
                w_mid,
                std::fabs(w_mid - w_max_an) / std::fabs(w_max_an),
                R_A_z, R_B_z);
    return 0;
}
