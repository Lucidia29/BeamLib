#pragma once
#include "../Core/Types.h"
#include <cmath>

namespace beamlib {

// Local right-handed orthonormal frame for a 3D beam element.
//   Vx — unit vector along the deformed-free beam axis, xB - xA normalized.
//   Vy — unit vector in the plane normal to Vx, derived from the user-supplied
//        refVector via Gram-Schmidt projection. Fallback rule for the
//        degenerate refVector || Vx case is documented in buildLocalFrame3D().
//   Vz — Vx x Vy (right-handed by construction).
//
// `lambda` is the 3x3 matrix with rows (Vx, Vy, Vz). It maps global vectors to
// local components: v_local = lambda * v_global. (Equivalently, the columns of
// lambda are the global basis vectors written in local coordinates.) This
// matches the existing 2D Rotation2D convention and the convention required by
// BeamModel::assemble / assembleMass (dispLocal = T * dispGlobal,
// K_global = T^T K_local T).
struct LocalFrame3D {
    Vec3 Vx;
    Vec3 Vy;
    Vec3 Vz;
    Mat3 lambda;
};

// Construct LocalFrame3D from end positions and a per-element refVector.
//
// Algorithm:
//   1. Vx = (xB - xA) / |xB - xA|. Requires |xB - xA| > 0; degenerate (zero-
//      length) elements are not supported and would be a modeling error.
//   2. Vy_raw = refVector - (refVector . Vx) Vx  (Gram-Schmidt projection of
//      refVector onto the plane normal to Vx).
//   3. If |Vy_raw| < tol (refVector ~|| Vx), substitute (0, 0, 1) and retry.
//   4. If still |Vy_raw| < tol (beam was along z and refVector was along z),
//      substitute (0, 1, 0). One of (0,0,1) or (0,1,0) is guaranteed not to be
//      parallel to Vx, so the fallback chain always terminates.
//   5. Vy = Vy_raw / |Vy_raw|, Vz = Vx x Vy.
//
// Note: the fallback ladder is (0,0,1) then (0,1,0), in that order, regardless
// of the original refVector. This makes the fallback frame deterministic and
// independent of which refVector the user supplied for a beam that ended up
// degenerate. The PROJECT_SPEC and EB3D theory doc both specify this order.
inline LocalFrame3D buildLocalFrame3D(const Vec3& xA, const Vec3& xB,
                                      const Vec3& refVector)
{
    const Vec3 dx = xB - xA;
    const double L = dx.norm();
    const Vec3 Vx = dx / L;

    constexpr double parallelTol = 1e-12;

    auto project = [&](const Vec3& ref) {
        return ref - ref.dot(Vx) * Vx;
    };

    Vec3 Vy_raw = project(refVector);
    if (Vy_raw.norm() < parallelTol) {
        Vy_raw = project(Vec3(0.0, 0.0, 1.0));
        if (Vy_raw.norm() < parallelTol) {
            Vy_raw = project(Vec3(0.0, 1.0, 0.0));
        }
    }
    const Vec3 Vy = Vy_raw.normalized();
    const Vec3 Vz = Vx.cross(Vy);

    LocalFrame3D f;
    f.Vx = Vx;
    f.Vy = Vy;
    f.Vz = Vz;
    f.lambda.row(0) = Vx.transpose();
    f.lambda.row(1) = Vy.transpose();
    f.lambda.row(2) = Vz.transpose();
    return f;
}

// 12x12 element transformation matrix for a 3D beam with DOF order
// [u_x, u_y, u_z, theta_x, theta_y, theta_z] at each node.
//
// The translation 3x3 blocks at (0,0) and (6,6) use the plain orthogonal
// lambda from buildLocalFrame3D: translations transform as a standard
// 3-vector (v_local = lambda * v_global).
//
// The rotation 3x3 blocks at (3,3) and (9,9) use S * lambda * S, where
// S = diag(1, -1, 1). This is the "structural-convention" adapter:
//
//   BeamLib stores rotation DOFs in the convention
//       theta_x = +Omega_x   (RH about local +x)
//       theta_y = -Omega_y   (structural; equals du_z/dx, opposite of RH
//                             about +y -- inherited from EB2D Chapter 1)
//       theta_z = +Omega_z   (matches RH about +z because theta_z = du_y/dx)
//   for the RH-rule physical rotation pseudovector Omega. Equivalently,
//   theta = S * Omega in any orthonormal frame.
//
//   For a physical rotation Omega_local = lambda * Omega_global. Therefore
//       theta_local = S * Omega_local
//                   = S * lambda * Omega_global
//                   = S * lambda * S * (S * Omega_global)
//                   = (S * lambda * S) * theta_global.
//
//   So the correct global->local transformation for BeamLib's stored
//   rotation DOFs is S * lambda * S, not lambda. The plain lambda block is
//   correct only when lambda commutes with S (e.g., beams in the xz plane);
//   for beams whose lambda mixes the y-axis with another axis the difference
//   is the commutator [lambda, S] != 0, and a rigid-body patch test fails
//   (see tests/test_eb3d_rigid_body.cpp).
//
// This adapter is specific to BeamLib's EB / Timoshenko mixed sign
// convention. A future element that uses pure RH-rule pseudovector rotations
// (e.g., GE3D with explicit Rodrigues parameterization) should build its
// own transformation rather than reuse this function.
inline MatMN<12, 12> buildTransformation3D(const Vec3& xA, const Vec3& xB,
                                           const Vec3& refVector)
{
    const LocalFrame3D f = buildLocalFrame3D(xA, xB, refVector);

    // Rotation block: lambda_rot(i,j) = s(i) * lambda(i,j) * s(j) with
    // s = (1, -1, 1). Equivalent to S * lambda * S.
    Mat3 lambda_rot = f.lambda;
    lambda_rot.row(1) *= -1.0;
    lambda_rot.col(1) *= -1.0;

    MatMN<12, 12> T = MatMN<12, 12>::Zero();
    T.block<3, 3>(0, 0) = f.lambda;     // node A translations
    T.block<3, 3>(3, 3) = lambda_rot;   // node A rotations
    T.block<3, 3>(6, 6) = f.lambda;     // node B translations
    T.block<3, 3>(9, 9) = lambda_rot;   // node B rotations
    return T;
}

} // namespace beamlib
