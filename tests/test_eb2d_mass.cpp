#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <cmath>
#include <cstdio>

namespace {

using namespace beamlib;

int check_symmetry(const MatX& M, const char* tag) {
    const int n = static_cast<int>(M.rows());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double d = std::fabs(M(i, j) - M(j, i));
            const double s = 1.0 + std::fabs(M(i, j));
            if (d > 1e-12 * s) {
                std::printf("FAIL [%s] symmetry M(%d,%d)=%g vs M(%d,%d)=%g\n",
                            tag, i, j, M(i, j), j, i, M(j, i));
                return 1;
            }
        }
    }
    return 0;
}

int check_relative(double got, double expected, const char* tag) {
    const double denom = std::fabs(expected) > 1e-30 ? std::fabs(expected) : 1.0;
    const double err = std::fabs(got - expected) / denom;
    if (err > 1e-12) {
        std::printf("FAIL [%s]: got %g expected %g (rel %g)\n",
                    tag, got, expected, err);
        return 1;
    }
    return 0;
}

// Single horizontal element, all DOFs free. The assembled global mass equals
// the element mass directly (no transformation, no DOF reduction).
int test_single_horizontal() {
    BeamModel<EulerBernoulli2D> model;
    model.props.E   = 200e9;
    model.props.A   = 0.01;
    model.props.Iz  = 8.333e-6;
    model.props.rho = 7800.0;

    const double Le = 1.0;
    model.nodes.resize(2);
    model.nodes[0].x0 = Vec3(0, 0, 0);
    model.nodes[1].x0 = Vec3(Le, 0, 0);

    ElementConn c;
    c.nodeA = 0;
    c.nodeB = 1;
    model.elements.push_back(c);

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);

    if (M.rows() != 6 || M.cols() != 6) {
        std::printf("FAIL single-horizontal size: %lld x %lld\n",
                    static_cast<long long>(M.rows()), static_cast<long long>(M.cols()));
        return 1;
    }

    MatX Md = MatX(M);
    if (check_symmetry(Md, "single-horizontal")) return 1;
    if (Md.norm() < 1e-30) {
        std::printf("FAIL single-horizontal: mass matrix is zero\n");
        return 1;
    }

    const double rhoAL = model.props.rho * model.props.A * Le;
    // Axial diagonal: 2 * rho*A*L / 6
    if (check_relative(Md(0, 0), 2.0 * rhoAL / 6.0, "M(0,0) axial diag"))   return 1;
    if (check_relative(Md(3, 3), 2.0 * rhoAL / 6.0, "M(3,3) axial diag"))   return 1;
    if (check_relative(Md(0, 3), 1.0 * rhoAL / 6.0, "M(0,3) axial off"))    return 1;
    // Bending diagonal: 156 * rho*A*L / 420
    if (check_relative(Md(1, 1), 156.0 * rhoAL / 420.0, "M(1,1) bending"))   return 1;
    if (check_relative(Md(4, 4), 156.0 * rhoAL / 420.0, "M(4,4) bending"))   return 1;
    // Coupled bending-rotation: 22*L * rho*A*L / 420 at (1,2) and (4,5) sign
    if (check_relative(Md(1, 2),  22.0 * Le * rhoAL / 420.0, "M(1,2)"))      return 1;
    if (check_relative(Md(4, 5), -22.0 * Le * rhoAL / 420.0, "M(4,5)"))      return 1;
    // Pure rotation diagonal: 4*L^2 * rho*A*L / 420
    if (check_relative(Md(2, 2), 4.0 * Le * Le * rhoAL / 420.0, "M(2,2)"))   return 1;

    // Total mass via rigid axial translation u = [1,0,0, 1,0,0]:
    //   u^T M u = sum of axial block entries = rho*A*L
    VecX u = VecX::Zero(6);
    u[0] = 1.0;
    u[3] = 1.0;
    const double rigid_mass = u.dot(Md * u);
    if (check_relative(rigid_mass, rhoAL, "rigid axial translation"))         return 1;

    return 0;
}

// Same element rotated 90 degrees (vertical, +z). The local axial mass must
// appear at the global u_z DOF; local bending mass at the global u_x DOF.
int test_single_vertical() {
    BeamModel<EulerBernoulli2D> model;
    model.props.E   = 200e9;
    model.props.A   = 0.01;
    model.props.Iz  = 8.333e-6;
    model.props.rho = 7800.0;

    const double Le = 1.0;
    model.nodes.resize(2);
    model.nodes[0].x0 = Vec3(0, 0, 0);
    model.nodes[1].x0 = Vec3(0, 0, Le);   // vertical, +z direction

    ElementConn c;
    c.nodeA = 0;
    c.nodeB = 1;
    model.elements.push_back(c);

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);

    MatX Md = MatX(M);
    if (Md.rows() != 6) {
        std::printf("FAIL vertical size: %lld\n", static_cast<long long>(Md.rows()));
        return 1;
    }
    if (check_symmetry(Md, "vertical")) return 1;

    const double rhoAL = model.props.rho * model.props.A * Le;
    // Vertical: T_n = [[0,-1,0],[1,0,0],[0,0,1]]. After T^T M_l T:
    //   M_g(0,0) = local bending diagonal 156*rhoAL/420 (global u_x at node A)
    //   M_g(1,1) = local axial diagonal  2*rhoAL/6     (global u_z at node A)
    //   M_g(2,2) = unchanged 4 L^2 rhoAL/420            (theta_y is invariant)
    if (check_relative(Md(0, 0), 156.0 * rhoAL / 420.0, "vertical M_g(0,0)"))   return 1;
    if (check_relative(Md(1, 1),   2.0 * rhoAL / 6.0,   "vertical M_g(1,1)"))   return 1;
    if (check_relative(Md(2, 2), 4.0 * Le * Le * rhoAL / 420.0,
                       "vertical M_g(2,2)"))                                    return 1;
    // u_x and u_z should not couple within a node block for axial+bending split
    if (std::fabs(Md(0, 1)) > 1e-12 * (1.0 + Md.norm())) {
        std::printf("FAIL vertical M_g(0,1) should be ~0, got %g\n", Md(0, 1));
        return 1;
    }
    return 0;
}

// Two horizontal elements, node 0 fully fixed: shared interior node 1 gets a
// contribution from both elements; size after DOF reduction is 6 (= 2 free
// nodes * 3 DOFs).
int test_two_elements_with_fixed_root() {
    BeamModel<EulerBernoulli2D> model;
    model.props.E   = 200e9;
    model.props.A   = 0.01;
    model.props.Iz  = 8.333e-6;
    model.props.rho = 7800.0;

    const double Le = 1.0;
    model.nodes.resize(3);
    for (int i = 0; i < 3; ++i) model.nodes[i].x0 = Vec3(i * Le, 0, 0);
    model.nodes[0].fixAll();

    for (int e = 0; e < 2; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        model.elements.push_back(c);
    }

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);

    if (M.rows() != 6) {
        std::printf("FAIL two-element size: %lld\n",
                    static_cast<long long>(M.rows()));
        return 1;
    }
    MatX Md = MatX(M);
    if (check_symmetry(Md, "two-element")) return 1;

    const double rhoAL = model.props.rho * model.props.A * Le;
    // Free DOF layout: node1 = idx 0..2, node2 = idx 3..5
    // Node 1 u_x diagonal: B-end of element 0 (M_l(3,3) = 2 rhoAL/6) +
    //                      A-end of element 1 (M_l(0,0) = 2 rhoAL/6)
    //                    = 4 rhoAL/6 = 2 rhoAL/3
    if (check_relative(Md(0, 0), 4.0 * rhoAL / 6.0, "node1 axial diag"))         return 1;
    // Node 2 u_x diagonal: only B-end of element 1 contributes
    if (check_relative(Md(3, 3), 2.0 * rhoAL / 6.0, "node2 axial diag"))         return 1;
    // Node 1 bending diagonal: B-end (156/420) of elem 0 + A-end (156/420) of elem 1
    if (check_relative(Md(1, 1), 2.0 * 156.0 * rhoAL / 420.0,
                       "node1 bending diag"))                                    return 1;
    // Node 2 bending diagonal: only one contribution
    if (check_relative(Md(4, 4), 156.0 * rhoAL / 420.0,
                       "node2 bending diag"))                                    return 1;
    return 0;
}

} // namespace

int main() {
    if (test_single_horizontal())          return 1;
    if (test_single_vertical())            return 1;
    if (test_two_elements_with_fixed_root()) return 1;
    std::printf("PASS test_eb2d_mass (all sub-tests)\n");
    return 0;
}
