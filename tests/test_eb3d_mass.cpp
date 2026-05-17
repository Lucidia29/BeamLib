#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Math/BeamMath3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <cmath>
#include <cstdio>

namespace {

using namespace beamlib;

int check_symmetry(const MatX& M, const char* tag)
{
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

int check_rel(double got, double expected, double tol, const char* tag)
{
    const double denom = std::fabs(expected) > 1e-30 ? std::fabs(expected) : 1.0;
    const double err = std::fabs(got - expected) / denom;
    if (err > tol) {
        std::printf("FAIL [%s]: got %g expected %g (rel %g, tol %g)\n",
                    tag, got, expected, err, tol);
        return 1;
    }
    return 0;
}

struct Sec {
    double E   = 200e9;
    double G   = 80e9;
    double A   = 0.01;
    double Iy  = 8.333e-6;
    double Iz  = 1.667e-5;     // intentionally Iy != Iz
    double Ix  = 2.5e-5;
    double rho = 7800.0;
};

void setProps(BeamModel<EulerBernoulli3D>& m, const Sec& s)
{
    m.props.E   = s.E;
    m.props.G   = s.G;
    m.props.A   = s.A;
    m.props.Iy  = s.Iy;
    m.props.Iz  = s.Iz;
    m.props.Ix  = s.Ix;
    m.props.rho = s.rho;
}

// Single horizontal element along +x with refVector=(0,1,0) -> lambda = I, so
// the assembled global mass matrix equals the local element mass directly.
// Verify all 4 blocks of the EB3D local mass:
//   axial (0,6), torsion (3,9), xy-bending (1,5,7,11), xz-bending (2,4,8,10).
int test_single_horizontal()
{
    Sec s;
    BeamModel<EulerBernoulli3D> model;
    setProps(model, s);

    const double Le = 1.0;
    model.nodes.resize(2);
    model.nodes[0].x0 = Vec3(0, 0, 0);
    model.nodes[1].x0 = Vec3(Le, 0, 0);

    ElementConn c;
    c.nodeA = 0;
    c.nodeB = 1;
    c.refVector = Vec3(0, 1, 0);
    model.elements.push_back(c);

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);

    if (M.rows() != 12 || M.cols() != 12) {
        std::printf("FAIL horizontal size: %lld x %lld (expected 12 x 12)\n",
                    (long long)M.rows(), (long long)M.cols());
        return 1;
    }
    MatX Md = MatX(M);
    if (check_symmetry(Md, "single-horizontal")) return 1;
    if (Md.norm() < 1e-30) {
        std::printf("FAIL horizontal: mass matrix is zero\n");
        return 1;
    }

    const double rhoAL  = s.rho * s.A  * Le;
    const double rhoIxL = s.rho * s.Ix * Le;
    const double L2     = Le * Le;
    const double bend   = rhoAL / 420.0;

    // --- Axial block at (0, 6) ---
    if (check_rel(Md(0, 0), 2.0 * rhoAL / 6.0, 1e-12, "axial M(0,0)")) return 1;
    if (check_rel(Md(6, 6), 2.0 * rhoAL / 6.0, 1e-12, "axial M(6,6)")) return 1;
    if (check_rel(Md(0, 6), 1.0 * rhoAL / 6.0, 1e-12, "axial M(0,6)")) return 1;

    // --- Torsional rotary inertia block at (3, 9) ---
    if (check_rel(Md(3, 3), 2.0 * rhoIxL / 6.0, 1e-12, "torsion M(3,3)")) return 1;
    if (check_rel(Md(9, 9), 2.0 * rhoIxL / 6.0, 1e-12, "torsion M(9,9)")) return 1;
    if (check_rel(Md(3, 9), 1.0 * rhoIxL / 6.0, 1e-12, "torsion M(3,9)")) return 1;

    // --- xy bending block at indices (1, 5, 7, 11) ---
    if (check_rel(Md(1, 1),  156.0 * bend,      1e-12, "xy M(1,1)"))   return 1;
    if (check_rel(Md(7, 7),  156.0 * bend,      1e-12, "xy M(7,7)"))   return 1;
    if (check_rel(Md(1, 7),   54.0 * bend,      1e-12, "xy M(1,7)"))   return 1;
    if (check_rel(Md(1, 5),   22.0 * Le * bend, 1e-12, "xy M(1,5)"))   return 1;
    if (check_rel(Md(1, 11), -13.0 * Le * bend, 1e-12, "xy M(1,11)"))  return 1;
    if (check_rel(Md(5, 7),   13.0 * Le * bend, 1e-12, "xy M(5,7)"))   return 1;
    if (check_rel(Md(7, 11), -22.0 * Le * bend, 1e-12, "xy M(7,11)"))  return 1;
    if (check_rel(Md(5, 5),    4.0 * L2 * bend, 1e-12, "xy M(5,5)"))   return 1;
    if (check_rel(Md(5, 11),  -3.0 * L2 * bend, 1e-12, "xy M(5,11)"))  return 1;
    if (check_rel(Md(11, 11),  4.0 * L2 * bend, 1e-12, "xy M(11,11)")) return 1;

    // --- xz bending block at indices (2, 4, 8, 10) ---
    if (check_rel(Md(2, 2),  156.0 * bend,      1e-12, "xz M(2,2)"))   return 1;
    if (check_rel(Md(8, 8),  156.0 * bend,      1e-12, "xz M(8,8)"))   return 1;
    if (check_rel(Md(2, 8),   54.0 * bend,      1e-12, "xz M(2,8)"))   return 1;
    if (check_rel(Md(2, 4),   22.0 * Le * bend, 1e-12, "xz M(2,4)"))   return 1;
    if (check_rel(Md(2, 10), -13.0 * Le * bend, 1e-12, "xz M(2,10)"))  return 1;
    if (check_rel(Md(4, 8),   13.0 * Le * bend, 1e-12, "xz M(4,8)"))   return 1;
    if (check_rel(Md(8, 10), -22.0 * Le * bend, 1e-12, "xz M(8,10)"))  return 1;
    if (check_rel(Md(4, 4),    4.0 * L2 * bend, 1e-12, "xz M(4,4)"))   return 1;
    if (check_rel(Md(4, 10),  -3.0 * L2 * bend, 1e-12, "xz M(4,10)"))  return 1;
    if (check_rel(Md(10, 10),  4.0 * L2 * bend, 1e-12, "xz M(10,10)")) return 1;

    // --- Rigid-translation total mass: u_x at both nodes = 1 -> rho*A*L ---
    VecX u = VecX::Zero(12);
    u[0] = 1.0;
    u[6] = 1.0;
    const double rigid_axial = u.dot(Md * u);
    if (check_rel(rigid_axial, rhoAL, 1e-12, "rigid axial translation"))
        return 1;
    // Rigid u_y both = 1
    u.setZero();
    u[1] = 1.0;
    u[7] = 1.0;
    const double rigid_uy = u.dot(Md * u);
    if (check_rel(rigid_uy, rhoAL, 1e-12, "rigid u_y translation"))
        return 1;
    // Rigid u_z both = 1
    u.setZero();
    u[2] = 1.0;
    u[8] = 1.0;
    const double rigid_uz = u.dot(Md * u);
    if (check_rel(rigid_uz, rhoAL, 1e-12, "rigid u_z translation"))
        return 1;
    // Rigid theta_x both = 1 -> rho*Ix*L (polar inertia)
    u.setZero();
    u[3] = 1.0;
    u[9] = 1.0;
    const double rigid_thx = u.dot(Md * u);
    if (check_rel(rigid_thx, rhoIxL, 1e-12, "rigid theta_x translation"))
        return 1;

    return 0;
}

// Element rotated 90 degrees so the local x is the global +y direction.
// With refVector = (0, 0, 1), the local frame is:
//   Vx = (0, 1, 0), Vy = (0, 0, 1), Vz = (1, 0, 0)
//   lambda = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
// Translation block uses lambda. So T^T M_l T maps the local axial mass to
// global u_y (since local x = global y), local xy-bending mass to global u_z
// (since local y = global z), local xz-bending mass to global u_x (since
// local z = global x). We check these three diagonals at node A directly.
int test_rotated_along_y()
{
    Sec s;
    BeamModel<EulerBernoulli3D> model;
    setProps(model, s);

    const double Le = 1.0;
    model.nodes.resize(2);
    model.nodes[0].x0 = Vec3(0, 0, 0);
    model.nodes[1].x0 = Vec3(0, Le, 0);

    ElementConn c;
    c.nodeA     = 0;
    c.nodeB     = 1;
    c.refVector = Vec3(0, 0, 1);
    model.elements.push_back(c);

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);
    MatX Md = MatX(M);
    if (check_symmetry(Md, "rotated-along-y")) return 1;

    const double rhoAL = s.rho * s.A  * Le;
    const double bend  = rhoAL / 420.0;

    // Node A is at indices 0..5 of the assembled (12-DOF) matrix.
    //
    // Local axial mass M_l(0,0) = 2 rhoAL/6 -> after T^T M T mapping,
    // appears at global u_y diagonal Md(1,1).
    if (check_rel(Md(1, 1), 2.0 * rhoAL / 6.0, 1e-12,
                  "rotated: axial mass on global u_y")) return 1;
    // Local xy-bending u_y diagonal M_l(1,1) = 156 bend (translational,
    // along local y = global z) -> global u_z diagonal Md(2,2).
    if (check_rel(Md(2, 2), 156.0 * bend, 1e-12,
                  "rotated: xy bending mass on global u_z")) return 1;
    // Local xz-bending u_z diagonal M_l(2,2) = 156 bend (along local z =
    // global x) -> global u_x diagonal Md(0,0).
    if (check_rel(Md(0, 0), 156.0 * bend, 1e-12,
                  "rotated: xz bending mass on global u_x")) return 1;

    return 0;
}

// Two horizontal elements with node 0 fully fixed.
// Free-DOF reduction: 12 free DOFs (2 free nodes x 6 each). Shared interior
// node accumulates contributions from both elements.
int test_two_elements_fixed_root()
{
    Sec s;
    BeamModel<EulerBernoulli3D> model;
    setProps(model, s);

    const double Le = 1.0;
    model.nodes.resize(3);
    for (int i = 0; i < 3; ++i) model.nodes[i].x0 = Vec3(i * Le, 0, 0);
    model.nodes[0].fixAll();

    for (int e = 0; e < 2; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        c.refVector = Vec3(0, 1, 0);
        model.elements.push_back(c);
    }

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);

    if (M.rows() != 12) {
        std::printf("FAIL two-element size: %lld (expected 12)\n",
                    (long long)M.rows());
        return 1;
    }
    MatX Md = MatX(M);
    if (check_symmetry(Md, "two-element")) return 1;

    const double rhoAL  = s.rho * s.A  * Le;
    const double rhoIxL = s.rho * s.Ix * Le;
    const double bend   = rhoAL / 420.0;

    // Free-DOF layout (Node A is fixed):
    //   node 1 (interior): free DOFs 0..5  = [u_x, u_y, u_z, th_x, th_y, th_z]
    //   node 2 (tip)     : free DOFs 6..11 = same per-node order
    //
    // Node 1 u_x diagonal: B-end of elem 0 (2 rhoAL/6) + A-end of elem 1 (2 rhoAL/6).
    if (check_rel(Md(0, 0), 4.0 * rhoAL / 6.0, 1e-12,
                  "node1 axial diag")) return 1;
    // Node 2 u_x diagonal: only B-end of elem 1.
    if (check_rel(Md(6, 6), 2.0 * rhoAL / 6.0, 1e-12,
                  "node2 axial diag")) return 1;
    // Node 1 theta_x diagonal: B-end of elem 0 + A-end of elem 1 (torsional).
    if (check_rel(Md(3, 3), 4.0 * rhoIxL / 6.0, 1e-12,
                  "node1 torsion diag")) return 1;
    // Node 2 theta_x diagonal: only one element contributes.
    if (check_rel(Md(9, 9), 2.0 * rhoIxL / 6.0, 1e-12,
                  "node2 torsion diag")) return 1;
    // Node 1 u_y diagonal: B-end of elem 0 + A-end of elem 1 (xy bending).
    if (check_rel(Md(1, 1), 2.0 * 156.0 * bend, 1e-12,
                  "node1 u_y bending diag")) return 1;
    if (check_rel(Md(7, 7), 156.0 * bend, 1e-12,
                  "node2 u_y bending diag")) return 1;
    // Node 1 u_z diagonal: B-end of elem 0 + A-end of elem 1 (xz bending).
    if (check_rel(Md(2, 2), 2.0 * 156.0 * bend, 1e-12,
                  "node1 u_z bending diag")) return 1;
    if (check_rel(Md(8, 8), 156.0 * bend, 1e-12,
                  "node2 u_z bending diag")) return 1;

    return 0;
}

// Generic rotated element: verify T^T M_l T transformation is correctly
// applied for the rotation block under the S * lambda * S adapter. We solve
// using the same assemble path and check that the assembled matrix
// reproduces M_g = T^T M_l T with the T from BeamMath3D.
int test_generic_rotated_consistency()
{
    Sec s;
    BeamModel<EulerBernoulli3D> model;
    setProps(model, s);

    const double L = 1.5;
    // Generic non-axis-aligned beam: from arbitrary xA along an arbitrary
    // direction, with a non-default refVector. The geometry is chosen so
    // that the local frame's lambda does not commute with S = diag(1,-1,1),
    // exercising the structural-convention adapter in the rotation blocks.
    const Vec3 xA(0.4, -0.3, 0.6);
    const Vec3 dir = Vec3(2.0, 1.0, -1.0).normalized();
    const Vec3 xB  = xA + L * dir;
    const Vec3 refVec(0.1, 1.0, 0.3);   // generic, not aligned with axes

    model.nodes.resize(2);
    model.nodes[0].x0 = xA;
    model.nodes[1].x0 = xB;

    ElementConn c;
    c.nodeA     = 0;
    c.nodeB     = 1;
    c.refVector = refVec;
    model.elements.push_back(c);

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);

    if (M.rows() != 12) return 1;
    MatX Md = MatX(M);
    if (check_symmetry(Md, "generic-rotated")) return 1;

    // Re-derive expected M_g = T^T M_l T independently, using the same
    // BeamMath3D utilities, and compare entry-by-entry.
    MatMN<12, 12> T = buildTransformation3D(xA, xB, refVec);
    ElementMassResult<6> mr =
        EulerBernoulli3D::computeMass(Vec3::Zero(), Vec3(L, 0, 0), model.props);
    MatMN<12, 12> Mg_expected = T.transpose() * mr.me * T;

    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            const double got = Md(i, j);
            const double exp = Mg_expected(i, j);
            const double tol = 1e-10 * (1.0 + std::fabs(exp));
            if (std::fabs(got - exp) > tol) {
                std::printf("FAIL [generic-rotated] M(%d,%d)=%g vs expected %g\n",
                            i, j, got, exp);
                return 1;
            }
        }
    }
    return 0;
}

} // namespace

int main()
{
    if (test_single_horizontal())          return 1;
    if (test_rotated_along_y())            return 1;
    if (test_two_elements_fixed_root())    return 1;
    if (test_generic_rotated_consistency()) return 1;
    std::printf("PASS test_eb3d_mass (all sub-tests)\n");
    return 0;
}
