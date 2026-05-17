#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Math/Rotation2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

namespace {

using namespace beamlib;

struct CantileverResult {
    double u_xB;
    double u_zB;
    double th_yB;
};

CantileverResult solveSingleElement(const Vec3& xA, const Vec3& xB,
                                    const VecN<6>& tipLoadGlobal) {
    BeamModel<EulerBernoulli2D> model;
    model.props.E  = 200e9;
    model.props.A  = 0.01;
    model.props.Iz = 8.333e-6;

    model.nodes.resize(2);
    model.nodes[0].x0 = xA;
    model.nodes[1].x0 = xB;
    model.nodes[0].fixAll();
    for (int d = 0; d < 3; ++d) {
        model.nodes[1].load[d] = tipLoadGlobal[3 + d];
    }
    ElementConn c;
    c.nodeA = 0;
    c.nodeB = 1;
    model.elements.push_back(c);

    model.buildDofMap();
    VecX Fext = model.getExternalForceVector();
    NRResult res = NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(model, Fext);
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL: NR not converged in single-element solve\n");
        std::exit(1);
    }
    return CantileverResult{
        model.nodes[1].dof[0],
        model.nodes[1].dof[1],
        model.nodes[1].dof[2],
    };
}

int subtest_local_transverse_force() {
    const double L = 1.0;
    const double F_local_z = -1000.0;  // local transverse force, downward in local frame

    // Reference: horizontal beam, force directly in global (since T = I).
    VecN<6> loadHoriz = VecN<6>::Zero();
    loadHoriz[3 + 1] = F_local_z;
    CantileverResult ref =
        solveSingleElement(Vec3(0, 0, 0), Vec3(L, 0, 0), loadHoriz);

    // 45-degree beam, same parameters. Transform local force into global, solve,
    // transform global tip displacement back to local, compare with ref.
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    const Vec3 xA(0, 0, 0);
    const Vec3 xB(L * inv_sqrt2, 0, L * inv_sqrt2);

    // Global -> local for displacement: dispLocal = T * dispGlobal.
    // For forces: contravariant, so forceGlobal = T^T * forceLocal.
    MatMN<6, 6> T = Rotation2D::compute(xA, xB);
    VecN<6> loadLocal = VecN<6>::Zero();
    loadLocal[3 + 1] = F_local_z;
    VecN<6> loadGlobal = T.transpose() * loadLocal;

    CantileverResult rot = solveSingleElement(xA, xB, loadGlobal);

    // Transform global tip displacement back to local frame using the same T.
    VecN<6> dispGlobal = VecN<6>::Zero();
    dispGlobal[3 + 0] = rot.u_xB;
    dispGlobal[3 + 1] = rot.u_zB;
    dispGlobal[3 + 2] = rot.th_yB;
    VecN<6> dispLocal = T * dispGlobal;

    const double tol = 1e-9;
    auto cmp = [&](double a, double b, const char* tag) {
        const double denom = std::fabs(b) > 1e-30 ? std::fabs(b) : 1.0;
        const double err = std::fabs(a - b) / denom;
        if (err > tol) {
            std::printf("FAIL [%s]: local %g vs ref %g (rel %g)\n",
                        tag, a, b, err);
            return 1;
        }
        return 0;
    };

    if (cmp(dispLocal[3 + 0], ref.u_xB,  "u_x_local"))  return 1;
    if (cmp(dispLocal[3 + 1], ref.u_zB,  "u_z_local"))  return 1;
    if (cmp(dispLocal[3 + 2], ref.th_yB, "theta_y"))    return 1;
    return 0;
}

int subtest_global_vertical_force() {
    // Same 45-degree beam, but apply pure global F_z at tip. With the local
    // axis tilted at 45 degrees, this global force has BOTH a local axial and
    // a local transverse component, so both global u_x and global u_z must be
    // non-zero. Cross-check via the local-frame projection.
    const double L = 1.0;
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    const Vec3 xA(0, 0, 0);
    const Vec3 xB(L * inv_sqrt2, 0, L * inv_sqrt2);

    const double Fz_global = -1000.0;
    VecN<6> loadGlobal = VecN<6>::Zero();
    loadGlobal[3 + 1] = Fz_global;
    CantileverResult res = solveSingleElement(xA, xB, loadGlobal);

    // Both translational global DOFs must couple in.
    if (std::fabs(res.u_xB) < 1e-12) {
        std::printf("FAIL: global u_x at tip should be non-zero for "
                    "45-deg beam under global F_z, got %g\n", res.u_xB);
        return 1;
    }
    if (std::fabs(res.u_zB) < 1e-12) {
        std::printf("FAIL: global u_z at tip should be non-zero, got %g\n",
                    res.u_zB);
        return 1;
    }

    // Decompose the global force into local components using Rotation2D:
    // F_local = T * F_global (force vector transforms the same way as a
    // co-vector here because T is orthogonal).
    MatMN<6, 6> T = Rotation2D::compute(xA, xB);
    VecN<6> loadLocal = T * loadGlobal;
    const double F_local_axial = loadLocal[3 + 0];
    const double F_local_trans = loadLocal[3 + 1];

    // For F_z_global = -1000 with c = s = 1/sqrt(2), Rotation2D yields:
    //   F_local_axial = c * 0 - s * (-1000) = +1000 / sqrt(2)
    //   F_local_trans = s * 0 + c * (-1000) = -1000 / sqrt(2)
    const double expected_axial = -Fz_global * inv_sqrt2;
    const double expected_trans =  Fz_global * inv_sqrt2;
    if (std::fabs(F_local_axial - expected_axial) > 1e-9 ||
        std::fabs(F_local_trans - expected_trans) > 1e-9) {
        std::printf("FAIL: local force decomposition mismatch "
                    "(axial=%g vs %g, trans=%g vs %g)\n",
                    F_local_axial, expected_axial,
                    F_local_trans, expected_trans);
        return 1;
    }

    // Solve a horizontal beam with this same local force and compare local
    // tip displacement; this confirms the axial/transverse coupling matches
    // an independent direct local solve.
    VecN<6> loadHoriz = VecN<6>::Zero();
    loadHoriz[3 + 0] = F_local_axial;
    loadHoriz[3 + 1] = F_local_trans;
    CantileverResult ref =
        solveSingleElement(Vec3(0, 0, 0), Vec3(L, 0, 0), loadHoriz);

    VecN<6> dispGlobal = VecN<6>::Zero();
    dispGlobal[3 + 0] = res.u_xB;
    dispGlobal[3 + 1] = res.u_zB;
    dispGlobal[3 + 2] = res.th_yB;
    VecN<6> dispLocal = T * dispGlobal;

    const double tol = 1e-9;
    auto cmp = [&](double a, double b, const char* tag) {
        const double denom = std::fabs(b) > 1e-30 ? std::fabs(b) : 1.0;
        const double err = std::fabs(a - b) / denom;
        if (err > tol) {
            std::printf("FAIL [%s]: rotated %g vs reference %g (rel %g)\n",
                        tag, a, b, err);
            return 1;
        }
        return 0;
    };
    if (cmp(dispLocal[3 + 0], ref.u_xB,  "coupled u_x_local")) return 1;
    if (cmp(dispLocal[3 + 1], ref.u_zB,  "coupled u_z_local")) return 1;
    if (cmp(dispLocal[3 + 2], ref.th_yB, "coupled theta_y"))   return 1;
    return 0;
}

} // namespace

int main() {
    if (subtest_local_transverse_force()) return 1;
    if (subtest_global_vertical_force())  return 1;
    std::printf("PASS test_eb2d_rotated_beam (both subtests)\n");
    return 0;
}
