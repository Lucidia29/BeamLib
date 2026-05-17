#include <BeamLib/Element/Timoshenko3D.h>
#include <BeamLib/Math/BeamMath3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// 3D Timoshenko transformation tests. Same coverage as the EB3D transform
// test in spirit: confirm that the local response of a rotated element is
// recovered by transforming global tip displacements back to local using
// the EB3D S*lambda*S adapter; confirm that refVector reorients which
// bending plane engages (and thus which I and kappa apply); and confirm
// that rigid-body motion produces near-zero residual.
//
// Because Timoshenko3D has directional shear stiffness, an axis-swap bug
// is easy to miss in pure-bending tests; the refVector-reorient sub-test
// here explicitly distinguishes (I_y, kappa_z) from (I_z, kappa_y).

namespace {

using namespace beamlib;

struct Sec {
    double E       = 200e9;
    double G       = 80e9;
    double rho     = 7800.0;
    double A       = 0.01;
    double Iy      = 8.0e-6;
    double Iz      = 2.4e-5;       // 3x Iy so the two planes are distinguishable
    double Ix      = 3.2e-5;
    double kappa_y = 5.0 / 6.0;
    double kappa_z = 0.85;          // != kappa_y
};

VecN<6> solveSingleElement(const Vec3& xA, const Vec3& xB,
                           const Vec3& refVector,
                           const VecN<6>& tipLoad,
                           const Sec& sec)
{
    BeamModel<Timoshenko3D> model;
    model.props.E       = sec.E;
    model.props.G       = sec.G;
    model.props.rho     = sec.rho;
    model.props.A       = sec.A;
    model.props.Iy      = sec.Iy;
    model.props.Iz      = sec.Iz;
    model.props.Ix      = sec.Ix;
    model.props.kappa_y = sec.kappa_y;
    model.props.kappa_z = sec.kappa_z;

    model.nodes.resize(2);
    model.nodes[0].x0 = xA;
    model.nodes[1].x0 = xB;
    model.nodes[0].fixAll();
    for (int d = 0; d < 6; ++d) model.nodes[1].load[d] = tipLoad[d];

    ElementConn c; c.nodeA = 0; c.nodeB = 1; c.refVector = refVector;
    model.elements.push_back(c);

    model.buildDofMap();
    NRResult res =
        NewtonRaphsonSolver<Timoshenko3D>::solveOneStep(model, model.getExternalForceVector());
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL single-element solve\n");
        std::exit(1);
    }
    return model.nodes[1].dof;
}

int check_rel(double a, double b, double tol, const char* tag)
{
    const double denom = std::fabs(b) > 1e-30 ? std::fabs(b) : 1.0;
    const double err = std::fabs(a - b) / denom;
    if (err > tol) {
        std::printf("FAIL [%s]: %g vs %g (rel %g, tol %g)\n", tag, a, b, err, tol);
        return 1;
    }
    return 0;
}

// Subtest A: rotated 45-deg-in-xz beam, refVector = (0, 1, 0). Apply a local
// +z transverse force, transform to global, solve, then transform tip back
// to local. Both translations and rotations (with the S*lambda*S adapter)
// must reproduce the horizontal-reference result.
int subtest_local_force_rotation()
{
    Sec sec;
    const double L = 1.0;
    const double F = -1000.0;

    // Reference: horizontal beam, lambda = I, so local = global.
    VecN<6> loadHoriz = VecN<6>::Zero();
    loadHoriz[2] = F;
    VecN<6> refTip = solveSingleElement(
        Vec3(0,0,0), Vec3(L,0,0), Vec3(0,1,0), loadHoriz, sec);

    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    const Vec3 xA(0, 0, 0);
    const Vec3 xB(L * inv_sqrt2, 0.0, L * inv_sqrt2);
    const Vec3 ref(0, 1, 0);

    MatMN<12, 12> T = buildTransformation3D(xA, xB, ref);
    Mat3 lambda = T.block<3, 3>(0, 0);

    VecN<6> loadRot = VecN<6>::Zero();
    Vec3 fLocal(0, 0, F);
    Vec3 fGlobal = lambda.transpose() * fLocal;
    loadRot[0] = fGlobal[0]; loadRot[1] = fGlobal[1]; loadRot[2] = fGlobal[2];
    VecN<6> rotTip = solveSingleElement(xA, xB, ref, loadRot, sec);

    // Use the same S*lambda*S adapter for rotations (same as test_eb3d_transform).
    Mat3 lambda_rot = lambda;
    lambda_rot.row(1) *= -1.0;
    lambda_rot.col(1) *= -1.0;

    Vec3 uG(rotTip[0], rotTip[1], rotTip[2]);
    Vec3 tG(rotTip[3], rotTip[4], rotTip[5]);
    Vec3 uL = lambda * uG;
    Vec3 tL = lambda_rot * tG;

    if (check_rel(uL[0], refTip[0], 1e-9, "local u_x")) return 1;
    if (check_rel(uL[1], refTip[1], 1e-9, "local u_y")) return 1;
    if (check_rel(uL[2], refTip[2], 1e-9, "local u_z")) return 1;
    if (check_rel(tL[0], refTip[3], 1e-9, "local theta_x")) return 1;
    if (check_rel(tL[1], refTip[4], 1e-9, "local theta_y")) return 1;
    if (check_rel(tL[2], refTip[5], 1e-9, "local theta_z")) return 1;

    std::printf("PASS subtest_local_force_rotation\n");
    return 0;
}

// Subtest B: same beam axis, two different refVectors. With refVector =
// (0,1,0), the local xy plane is the global xy plane and a global F_y
// engages local xy bending (uses E I_z, kappa_y G A). With refVector =
// (0,0,1), local xy = global xz and a global F_y engages local xz bending
// (uses E I_y, kappa_z G A). Asymmetric Iy != Iz and kappa_y != kappa_z
// make these distinguishable to high precision.
int subtest_refvector_reorient_plane()
{
    Sec sec;
    const double L = 1.0;
    const double Fy = -1000.0;
    VecN<6> load = VecN<6>::Zero();
    load[1] = Fy;

    VecN<6> tipA = solveSingleElement(
        Vec3(0,0,0), Vec3(L,0,0), Vec3(0,1,0), load, sec);
    const double uy_A_expect =
        Fy * L * L * L / (3.0 * sec.E * sec.Iz)
        + Fy * L / (sec.kappa_y * sec.G * sec.A);
    if (check_rel(tipA[1], uy_A_expect, 1e-9,
                  "u_y refVector=(0,1,0) (xy bending: EIz, kappa_y)")) return 1;

    VecN<6> tipB = solveSingleElement(
        Vec3(0,0,0), Vec3(L,0,0), Vec3(0,0,1), load, sec);
    // With refVector = (0,0,1) and beam along +x:
    //   lambda = [[1,0,0],[0,0,1],[0,-1,0]]
    //   global F_y -> local force = lambda * (0,F_y,0) = (0, 0, -F_y)
    //   so local F_z = -F_y, engaging local xz bending (EIy, kappa_z)
    //   local u_z_loc = (-Fy) * L^3 / (3 EIy) + (-Fy) * L / (kappa_z GA)
    //   global u_y    = (lambda^T * (0, 0, u_z_loc))[1]
    //                 = (lambda row 2)[1] * u_z_loc = -u_z_loc
    //   so global u_y = -u_z_loc = +Fy * (L^3/(3 EIy) + L/(kappa_z GA))
    // (Same algebraic form as the EB3D refVector test, with the extra
    //  shear-flexibility term.)
    const double uy_expected_B =
        Fy * (L * L * L / (3.0 * sec.E * sec.Iy)
              + L / (sec.kappa_z * sec.G * sec.A));
    if (check_rel(tipB[1], uy_expected_B, 1e-9,
                  "u_y refVector=(0,0,1) (xz bending: EIy, kappa_z)")) return 1;

    // The two answers must differ since (Iz, kappa_y) != (Iy, kappa_z).
    if (std::fabs(tipA[1] - tipB[1]) < 1e-9) {
        std::printf("FAIL: refVector must reorient bending plane; got same u_y "
                    "(%g, %g)\n", tipA[1], tipB[1]);
        return 1;
    }

    std::printf("PASS subtest_refvector_reorient_plane: "
                "u_y(refY)=%g (Iz/kappa_y plane), u_y(refZ)=%g (Iy/kappa_z plane)\n",
                tipA[1], tipB[1]);
    return 0;
}

// Subtest C: rigid-body patch test. Replicates the EB3D rigid-body test for
// Timoshenko 3D, ensuring the inherited S*lambda*S transformation still
// produces zero internal residual under a true rigid-body motion when the
// beam orientation is one where lambda does NOT commute with S.
int subtest_rigid_body()
{
    Sec sec;
    const double L = 1.7;
    const double alpha = 0.7;   // xy-plane rotation, lambda does not commute with S
    const Vec3 xA(0.4, -0.3, 0.6);
    const Vec3 xB = xA + L * Vec3(std::cos(alpha), std::sin(alpha), 0.0);
    const Vec3 refVector(0, 0, 1);

    BeamModel<Timoshenko3D> model;
    model.props.E = sec.E; model.props.G = sec.G; model.props.rho = sec.rho;
    model.props.A = sec.A; model.props.Iy = sec.Iy; model.props.Iz = sec.Iz; model.props.Ix = sec.Ix;
    model.props.kappa_y = sec.kappa_y; model.props.kappa_z = sec.kappa_z;
    model.nodes.resize(2);
    model.nodes[0].x0 = xA;
    model.nodes[1].x0 = xB;
    ElementConn c; c.nodeA = 0; c.nodeB = 1; c.refVector = refVector;
    model.elements.push_back(c);

    const Vec3 t(0.002, -0.001, 0.003);
    const Vec3 Omega(0.0015, -0.0008, 0.0012);
    Vec3 thetaStruct(Omega.x(), -Omega.y(), Omega.z());

    auto setNode = [&](int i, const Vec3& X) {
        const Vec3 u = t + Omega.cross(X);
        model.nodes[i].dof[0] = u.x();
        model.nodes[i].dof[1] = u.y();
        model.nodes[i].dof[2] = u.z();
        model.nodes[i].dof[3] = thetaStruct.x();
        model.nodes[i].dof[4] = thetaStruct.y();
        model.nodes[i].dof[5] = thetaStruct.z();
    };
    setNode(0, xA);
    setNode(1, xB);

    model.buildDofMap();
    VecX R; SpMat K;
    model.assemble(R, K);

    const double Rnorm = R.norm();
    if (Rnorm > 1e-9) {
        std::printf("FAIL rigid-body: ||R||=%g exceeds tol 1e-9\n", Rnorm);
        return 1;
    }
    std::printf("PASS subtest_rigid_body: ||R||=%g\n", Rnorm);
    return 0;
}

} // namespace

int main()
{
    if (subtest_local_force_rotation())      return 1;
    if (subtest_refvector_reorient_plane())  return 1;
    if (subtest_rigid_body())                return 1;
    std::printf("PASS test_timo3d_transform (all sub-tests)\n");
    return 0;
}
