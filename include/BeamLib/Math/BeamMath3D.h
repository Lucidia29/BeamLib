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
// [u_x, u_y, u_z, theta_x, theta_y, theta_z] at node A then at node B.
// T = block-diag(lambda, lambda, lambda, lambda), so each 3-vector DOF group
// (uA, thetaA, uB, thetaB) transforms with the same 3x3 lambda. Both
// translations and small rotations are treated as 3-component vectors here;
// this is exact for the translations and consistent for linear EB/Timoshenko
// rotation DOFs (the conventions for theta_y/theta_z in the local frame match
// EB2D's structural convention, see docs/theory/02_euler_bernoulli_3d.tex).
inline MatMN<12, 12> buildTransformation3D(const Vec3& xA, const Vec3& xB,
                                           const Vec3& refVector)
{
    const LocalFrame3D f = buildLocalFrame3D(xA, xB, refVector);
    MatMN<12, 12> T = MatMN<12, 12>::Zero();
    for (int blk = 0; blk < 4; ++blk) {
        T.block<3, 3>(blk * 3, blk * 3) = f.lambda;
    }
    return T;
}

} // namespace beamlib
