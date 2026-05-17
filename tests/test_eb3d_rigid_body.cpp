#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Math/BeamMath3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <cmath>
#include <cstdio>

// 3D rigid-body patch test for EB3D.
//
// Premise:
//   For any arbitrarily oriented element, a small rigid-body translation
//   t + Omega x X applied to both nodes produces zero internal strain and
//   therefore zero element residual. This must hold *independently of the
//   local frame chosen by refVector* and *independently of which axis the
//   beam is along*.
//
// Subtlety this test catches (Codex Batch 3 review, issue 1):
//   BeamLib stores rotation DOFs in the "structural" convention:
//       theta_x = +Omega_x    (RH about local +x; matches RH-rule)
//       theta_y = -Omega_y    (structural; opposite of RH about +y because
//                              theta_y is defined as du_z/dx)
//       theta_z = +Omega_z    (matches RH about +z because theta_z = du_y/dx)
//   So in either the global or local frame, theta = S * Omega where
//   S = diag(1, -1, 1).
//
//   For a physical rotation pseudovector, Omega_local = lambda * Omega_global.
//   Therefore the correct local-from-global transformation for the rotation
//   DOFs stored in this mixed convention is
//       theta_local = S * lambda * S * theta_global,
//   *not* theta_local = lambda * theta_global. The plain-lambda transform is
//   correct only when lambda happens to commute with S (e.g., beams in the
//   xz plane). For beams whose local frame is rotated out of the xz plane,
//   the plain-lambda transform produces a non-physical local rotation field
//   and a non-zero residual under a true rigid-body motion.
//
// This test uses a beam in the xy plane plus a generic refVector, ensuring
// lambda does NOT commute with S. It then imposes a generic rigid-body
// (t, Omega) and asserts that the internal residual norm is at machine zero.
int main()
{
    using namespace beamlib;

    const double E   = 200e9;
    const double G   = 80e9;
    const double A   = 0.01;
    const double Iy  = 8.333e-6;
    const double Iz  = 1.6e-5;     // intentionally Iy != Iz
    const double Ix  = 2.0e-5;

    const double L = 1.7;          // arbitrary, non-unit, non-special
    // Beam in the xy plane at an angle alpha (rotation about z-axis).
    const double alpha = 0.7;      // ~40 degrees
    const Vec3 xA(0.4, -0.3, 0.6); // arbitrary origin to also test
                                   //   translation-invariance of the frame
    const Vec3 xB = xA + L * Vec3(std::cos(alpha), std::sin(alpha), 0.0);
    const Vec3 refVector(0.0, 0.0, 1.0);   // perpendicular to beam axis

    // Sanity: confirm this geometry produces a lambda that does NOT commute
    // with S. If lambda commuted with S, the test would not catch the bug.
    LocalFrame3D frame = buildLocalFrame3D(xA, xB, refVector);
    {
        Mat3 S = Mat3::Zero();
        S(0, 0) =  1.0;
        S(1, 1) = -1.0;
        S(2, 2) =  1.0;
        Mat3 commutator = frame.lambda * S - S * frame.lambda;
        if (commutator.norm() < 1e-12) {
            std::printf("FAIL setup: lambda commutes with S=diag(1,-1,1); "
                        "this geometry will not exercise the bug. "
                        "Adjust beam orientation / refVector.\n");
            return 1;
        }
    }

    // Build a 2-node model with NO constraints. We will impose nodal DOFs
    // directly (matching a rigid-body field) and check that the assembled
    // residual is essentially zero.
    BeamModel<EulerBernoulli3D> model;
    model.props.E  = E;
    model.props.G  = G;
    model.props.A  = A;
    model.props.Iy = Iy;
    model.props.Iz = Iz;
    model.props.Ix = Ix;

    model.nodes.resize(2);
    model.nodes[0].x0 = xA;
    model.nodes[1].x0 = xB;
    // No DOFs fixed -> all 12 DOFs are "free" in the DofMap. The residual
    // will be a 12-vector; we check it is ~ 0 element-wise.

    ElementConn conn;
    conn.nodeA     = 0;
    conn.nodeB     = 1;
    conn.refVector = refVector;
    model.elements.push_back(conn);

    // Rigid body: small translation t + small rotation Omega (RH pseudovector,
    // GLOBAL frame). Magnitudes chosen small enough that the linear theory
    // model is in its regime of validity, but big enough that any per-DOF
    // residual much above machine eps signals a real consistency failure.
    const Vec3 t(0.002, -0.001, 0.003);
    const Vec3 Omega(0.0015, -0.0008, 0.0012);

    // BeamLib's structural rotation representation: theta = S * Omega.
    Vec3 thetaStruct(Omega.x(), -Omega.y(), Omega.z());

    auto setNode = [&](int i, const Vec3& X) {
        // u = t + Omega x X  (true rigid-body translation field, RH cross)
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

    // Build the DofMap (all DOFs free) and assemble the internal residual.
    model.buildDofMap();
    VecX R;
    SpMat K;
    model.assemble(R, K);

    // Acceptance: the rigid-body field is reproduced exactly by linear shape
    // functions (axial + torsion) and exactly by the cubic Hermite functions
    // (bending), so the integration / assembly arithmetic is exact in real
    // arithmetic. Floating-point round-off at this scale should be well
    // below ~1e-9 in the residual norm.
    const double tol = 1e-9;
    const double Rnorm = R.norm();
    if (Rnorm > tol) {
        std::printf("FAIL: rigid-body residual norm %g exceeds tol %g\n",
                    Rnorm, tol);
        for (int i = 0; i < R.size(); ++i) {
            std::printf("   R[%d] = %.6g\n", i, R[i]);
        }
        return 1;
    }

    std::printf("PASS test_eb3d_rigid_body: ||R|| = %g (tol %g)\n", Rnorm, tol);
    return 0;
}
