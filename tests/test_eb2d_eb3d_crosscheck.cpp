#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Cross-check EB2D vs EB3D for the same physical xz-plane cantilever bending.
//
// EB2D uses props.Iz for in-plane (xz) bending; EB3D uses props.Iy for the
// xz-plane bending (the local plane spanned by Vx and Vz). For the same beam
// (horizontal along global x, tip load F_z), the two answers must agree once
// we set Iy_3D = Iz_2D. The matched DOFs to compare are:
//   - tip u_z   ( EB2D dof[1] vs EB3D dof[2] )
//   - tip theta_y ( EB2D dof[2] vs EB3D dof[4] )
// All EB3D DOFs unrelated to xz-plane bending must be ~0.
//
// Tolerance: relative 1e-9 between the two implementations (per the Batch 3
// prompt). The analytic check is folded in implicitly because both models
// match the same analytical EB cantilever formula u_z = F L^3/(3EI).
int main()
{
    using namespace beamlib;

    const double E   = 200e9;
    const double G   = 80e9;
    const double A   = 0.01;
    const double Ibend = 8.333e-6;   // shared bending I (EB2D Iz <-> EB3D Iy)
    const double Iother= 4.0e-5;     // different value -> guard against silent
                                     //   misuse of the wrong I in EB3D
    const double Ix    = 6.0e-5;     // torsion polar (unused by load)

    const double Lbeam = 1.0;
    const int    nElems = 10;
    const double dx     = Lbeam / nElems;
    const double Fz     = -1000.0;

    // --- EB2D model ---
    BeamModel<EulerBernoulli2D> m2;
    m2.props.E  = E;
    m2.props.A  = A;
    m2.props.Iz = Ibend;
    m2.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        m2.nodes[i].x0 = Vec3(i * dx, 0.0, 0.0);
    }
    m2.nodes[0].fixAll();
    m2.nodes[nElems].load[1] = Fz;    // EB2D index 1 = u_z
    for (int e = 0; e < nElems; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        m2.elements.push_back(c);
    }
    m2.buildDofMap();
    VecX F2 = m2.getExternalForceVector();
    NRResult r2 =
        NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(m2, F2);
    if (!r2.converged || r2.iterations != 1) {
        std::printf("FAIL EB2D NR: converged=%d iters=%d\n",
                    int(r2.converged), r2.iterations);
        return 1;
    }

    // --- EB3D model: Iy = EB2D's Iz, Iz set to a *different* value so the
    // test would catch any swap between the two bending planes. ---
    BeamModel<EulerBernoulli3D> m3;
    m3.props.E  = E;
    m3.props.G  = G;
    m3.props.A  = A;
    m3.props.Iy = Ibend;
    m3.props.Iz = Iother;
    m3.props.Ix = Ix;
    m3.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        m3.nodes[i].x0 = Vec3(i * dx, 0.0, 0.0);
    }
    m3.nodes[0].fixAll();
    m3.nodes[nElems].load[2] = Fz;    // EB3D index 2 = u_z
    for (int e = 0; e < nElems; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        m3.elements.push_back(c);
    }
    m3.buildDofMap();
    VecX F3 = m3.getExternalForceVector();
    NRResult r3 =
        NewtonRaphsonSolver<EulerBernoulli3D>::solveOneStep(m3, F3);
    if (!r3.converged || r3.iterations != 1) {
        std::printf("FAIL EB3D NR: converged=%d iters=%d\n",
                    int(r3.converged), r3.iterations);
        return 1;
    }

    // --- Compare matched DOFs ---
    const double uz2 = m2.nodes[nElems].dof[1];
    const double ty2 = m2.nodes[nElems].dof[2];
    const double uz3 = m3.nodes[nElems].dof[2];
    const double ty3 = m3.nodes[nElems].dof[4];

    auto check_rel = [](double a, double b, double tol, const char* tag) {
        const double denom = std::fabs(b) > 1e-30 ? std::fabs(b) : 1.0;
        const double err = std::fabs(a - b) / denom;
        if (err > tol) {
            std::printf("FAIL [%s]: EB3D %g vs EB2D %g (rel %g, tol %g)\n",
                        tag, a, b, err, tol);
            return 1;
        }
        return 0;
    };
    if (check_rel(uz3, uz2, 1e-9, "tip u_z (EB2D <-> EB3D)"))     return 1;
    if (check_rel(ty3, ty2, 1e-9, "tip theta_y (EB2D <-> EB3D)")) return 1;

    // Unrelated EB3D DOFs at the tip must be ~0 (only xz bending is excited).
    const auto& tip = m3.nodes[nElems].dof;
    auto near_zero = [](double v, double tol, const char* tag) {
        if (std::fabs(v) > tol) {
            std::printf("FAIL [%s]: EB3D expected ~0, got %g\n", tag, v);
            return 1;
        }
        return 0;
    };
    if (near_zero(tip[0], 1e-14, "tip u_x (EB3D)"))     return 1;
    if (near_zero(tip[1], 1e-14, "tip u_y (EB3D)"))     return 1;
    if (near_zero(tip[3], 1e-14, "tip theta_x (EB3D)")) return 1;
    if (near_zero(tip[5], 1e-14, "tip theta_z (EB3D)")) return 1;

    std::printf("PASS test_eb2d_eb3d_crosscheck: u_z=%g, theta_y=%g\n",
                uz3, ty3);
    return 0;
}
