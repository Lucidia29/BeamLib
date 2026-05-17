#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Math/BeamMath3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

namespace {

using namespace beamlib;

struct SectionInputs {
    double E   = 200e9;
    double G   = 80e9;
    double A   = 0.01;
    double Iy  = 2.0e-5;     // intentionally Iy != Iz so the bending plane
    double Iz  = 5.0e-6;     // choice is verifiable
    double Ix  = 4.0e-5;
};

// Single-element cantilever along an arbitrary direction with custom
// refVector. Returns the tip displacement DOF vector after one NR step.
VecN<6> solveSingleElement(const Vec3& xA, const Vec3& xB,
                           const Vec3& refVector,
                           const VecN<6>& tipLoad,
                           const SectionInputs& sec)
{
    BeamModel<EulerBernoulli3D> model;
    model.props.E  = sec.E;
    model.props.G  = sec.G;
    model.props.A  = sec.A;
    model.props.Iy = sec.Iy;
    model.props.Iz = sec.Iz;
    model.props.Ix = sec.Ix;

    model.nodes.resize(2);
    model.nodes[0].x0 = xA;
    model.nodes[1].x0 = xB;
    model.nodes[0].fixAll();
    for (int d = 0; d < 6; ++d) model.nodes[1].load[d] = tipLoad[d];

    ElementConn c;
    c.nodeA     = 0;
    c.nodeB     = 1;
    c.refVector = refVector;
    model.elements.push_back(c);

    model.buildDofMap();
    VecX Fext = model.getExternalForceVector();
    NRResult res =
        NewtonRaphsonSolver<EulerBernoulli3D>::solveOneStep(model, Fext);
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL: NR did not converge in a single-element solve\n");
        std::exit(1);
    }
    return model.nodes[1].dof;
}

int check_rel(double a, double b, double tol, const char* tag)
{
    const double denom = std::fabs(b) > 1e-30 ? std::fabs(b) : 1.0;
    const double err = std::fabs(a - b) / denom;
    if (err > tol) {
        std::printf("FAIL [%s]: %g vs %g (rel %g, tol %g)\n",
                    tag, a, b, err, tol);
        return 1;
    }
    return 0;
}

// Subtest A: rotated beam with a local transverse force.
//   - reference: horizontal beam along +x, refVector=(0,1,0), apply local +z
//     tip force F_local_z = -1000 (in BeamLib's identity local frame, this is
//     just a global F_z load).
//   - rotated:  same length, but the beam axis is set so that local x makes
//     a 45-degree angle with global x in the xz plane (refVector=(0,1,0),
//     so the local y axis stays aligned with global y). Apply the SAME local
//     force, transformed into global coordinates.
//   - compare:  transform the rotated-beam tip displacement back to local and
//     verify it equals the reference tip displacement.
int subtest_local_force_rotation()
{
    SectionInputs sec;
    const double L = 1.0;
    const double F = -1000.0;          // local +z load

    // Reference: horizontal beam (lambda = I -> local = global).
    VecN<6> loadHoriz = VecN<6>::Zero();
    loadHoriz[2] = F;                  // global = local u_z load
    VecN<6> refTip = solveSingleElement(
        Vec3(0, 0, 0), Vec3(L, 0, 0), Vec3(0, 1, 0), loadHoriz, sec);

    // Rotated 45-deg in xz plane. Local x = (cos45, 0, sin45), refVector kept
    // at (0,1,0) so local y stays at (0,1,0); local z is then x cross y.
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    const Vec3 xA(0, 0, 0);
    const Vec3 xB(L * inv_sqrt2, 0.0, L * inv_sqrt2);
    const Vec3 ref(0, 1, 0);

    MatMN<12, 12> T = buildTransformation3D(xA, xB, ref);
    Mat3 lambda = T.block<3, 3>(0, 0);

    // Force vectors: with T orthogonal (lambda is a rotation), force as a
    // column vector transforms the same as a displacement: f_global = T^T f_local.
    // But for the second-node block at (3-element offset of 6 in the node-B
    // translation slot 6..8), it suffices to write the same lambda^T mapping:
    VecN<6> loadRot = VecN<6>::Zero();
    Vec3 fLocal(0, 0, F);              // local +z load at tip translations
    Vec3 fGlobal = lambda.transpose() * fLocal;
    loadRot[0] = fGlobal[0];
    loadRot[1] = fGlobal[1];
    loadRot[2] = fGlobal[2];

    VecN<6> rotTip = solveSingleElement(xA, xB, ref, loadRot, sec);

    // Convert rotated tip global displacement back to local.
    Vec3 uG(rotTip[0], rotTip[1], rotTip[2]);
    Vec3 tG(rotTip[3], rotTip[4], rotTip[5]);
    Vec3 uL = lambda * uG;
    Vec3 tL = lambda * tG;

    // Compare with reference (which IS in local frame because lambda = I).
    if (check_rel(uL[0], refTip[0], 1e-9, "local u_x")) return 1;
    if (check_rel(uL[1], refTip[1], 1e-9, "local u_y")) return 1;
    if (check_rel(uL[2], refTip[2], 1e-9, "local u_z")) return 1;
    if (check_rel(tL[0], refTip[3], 1e-9, "local theta_x")) return 1;
    if (check_rel(tL[1], refTip[4], 1e-9, "local theta_y")) return 1;
    if (check_rel(tL[2], refTip[5], 1e-9, "local theta_z")) return 1;

    std::printf("PASS subtest_local_force_rotation\n");
    return 0;
}

// Subtest B: changing refVector reorients the bending plane.
//
// Beam axis along global x, length L. Apply global F_y at the tip.
//
//   refVector = (0,1,0): local y = global y, local z = global z. Global F_y
//   appears as local F_y (load in local xy plane), so it engages local
//   xy-plane bending which uses EIz. Expected tip global u_y = F_y L^3/(3 E Iz).
//
//   refVector = (0,0,1): local y = global z, local z = -global y. Global F_y
//   appears in the local frame as a -F_y load in local z (so xz-plane bending
//   engages, with EIy). Expected tip global u_y = F_y L^3/(3 E Iy).
//
// With Iy != Iz the two answers must differ, and each must agree with the
// closed-form formula corresponding to its bending plane.
int subtest_refvector_changes_bending_plane()
{
    SectionInputs sec;
    const double L  = 1.0;
    const double Fy = -1000.0;

    VecN<6> load = VecN<6>::Zero();
    load[1] = Fy;

    // refVector = (0,1,0): xy plane bending, EIz.
    VecN<6> tipA = solveSingleElement(
        Vec3(0, 0, 0), Vec3(L, 0, 0), Vec3(0, 1, 0), load, sec);
    const double uy_expected_A = Fy * L * L * L / (3.0 * sec.E * sec.Iz);
    if (check_rel(tipA[1], uy_expected_A, 1e-9,
                  "u_y with refVector=(0,1,0) (xy bending, EIz)")) return 1;
    // Other translational DOFs should be ~0.
    if (std::fabs(tipA[0]) > 1e-14 || std::fabs(tipA[2]) > 1e-14) {
        std::printf("FAIL: refVector=(0,1,0) gave nonzero u_x or u_z (%g, %g)\n",
                    tipA[0], tipA[2]);
        return 1;
    }

    // refVector = (0,0,1): xz plane bending, EIy.
    VecN<6> tipB = solveSingleElement(
        Vec3(0, 0, 0), Vec3(L, 0, 0), Vec3(0, 0, 1), load, sec);
    const double uy_expected_B = Fy * L * L * L / (3.0 * sec.E * sec.Iy);
    if (check_rel(tipB[1], uy_expected_B, 1e-9,
                  "u_y with refVector=(0,0,1) (xz bending, EIy)")) return 1;
    if (std::fabs(tipB[0]) > 1e-14 || std::fabs(tipB[2]) > 1e-14) {
        std::printf("FAIL: refVector=(0,0,1) gave nonzero u_x or u_z (%g, %g)\n",
                    tipB[0], tipB[2]);
        return 1;
    }

    // Sanity: with Iy != Iz the two deflections must differ; this is the
    // whole point of the test, so guard against a silent fallback that
    // ignores refVector.
    if (std::fabs(tipA[1] - tipB[1]) < 1e-9) {
        std::printf("FAIL: tip u_y must depend on refVector when Iy != Iz "
                    "(got %g vs %g)\n", tipA[1], tipB[1]);
        return 1;
    }

    std::printf("PASS subtest_refvector_changes_bending_plane: "
                "u_y(refY)=%g, u_y(refZ)=%g\n", tipA[1], tipB[1]);
    return 0;
}

// Subtest C: fallback when refVector is parallel to the beam axis.
//
// Two cases exercise the fallback chain:
//   C1: beam along (0, 1, 0) with refVector = (0, 1, 0) -> first fallback to
//       (0, 0, 1) kicks in.
//   C2: beam along (0, 0, 1) with refVector = (0, 0, 1) -> first fallback
//       (0, 0, 1) is ALSO parallel, second fallback to (0, 1, 0) kicks in.
//
// For each, verify that the constructed lambda is orthogonal with positive
// determinant, and that a transverse global force gives the analytical
// EB deflection along the corresponding bending plane.
int subtest_fallback()
{
    SectionInputs sec;
    const double L = 1.0;

    auto check_orthogonal = [&](const Mat3& lambda, const char* tag) {
        Mat3 I3 = lambda * lambda.transpose();
        Mat3 dev = I3 - Mat3::Identity();
        if (dev.norm() > 1e-12) {
            std::printf("FAIL [%s]: lambda not orthogonal (||L L^T - I||=%g)\n",
                        tag, dev.norm());
            return 1;
        }
        if (lambda.determinant() < 0.0) {
            std::printf("FAIL [%s]: det(lambda) negative (%g)\n",
                        tag, lambda.determinant());
            return 1;
        }
        return 0;
    };

    // --- C1: beam along +y, refVector = (0,1,0) ---
    {
        const Vec3 xA(0, 0, 0);
        const Vec3 xB(0, L, 0);
        const Vec3 ref(0, 1, 0);
        LocalFrame3D f = buildLocalFrame3D(xA, xB, ref);
        if (check_orthogonal(f.lambda, "C1 fallback frame")) return 1;
        // The first fallback (0,0,1) is perpendicular to +y, so local y must
        // be (0,0,1) and local z must be Vx x Vy = (+y) x (+z) = (+x).
        Vec3 expVy(0, 0, 1);
        Vec3 expVz(1, 0, 0);
        if ((f.Vy - expVy).norm() > 1e-12 || (f.Vz - expVz).norm() > 1e-12) {
            std::printf("FAIL C1: fallback frame mismatch Vy=(%g,%g,%g) "
                        "Vz=(%g,%g,%g)\n",
                        f.Vy[0], f.Vy[1], f.Vy[2],
                        f.Vz[0], f.Vz[1], f.Vz[2]);
            return 1;
        }

        // Apply global F_z at tip. Local force = lambda * (0,0,F_z) =
        // (lambda col 2)*F_z = (Vx_z, Vy_z, Vz_z)*F_z = (0, 1, 0)*F_z.
        // So local F_y = F_z, engaging local xy-plane bending with EIz.
        // Global u_z = (lambda^T (0, u_y_loc, 0))[2] = (row 1 of lambda)[2]
        //            = Vy[2] = 1 -> global u_z = u_y_loc.
        const double Fz = -1000.0;
        VecN<6> load = VecN<6>::Zero();
        load[2] = Fz;
        VecN<6> tip = solveSingleElement(xA, xB, ref, load, sec);
        const double uz_expected = Fz * L * L * L / (3.0 * sec.E * sec.Iz);
        if (check_rel(tip[2], uz_expected, 1e-9,
                      "C1 global u_z under global F_z")) return 1;
        // u_x should be ~0 (no axial), u_y ~0 (beam along y, axial direction)
        // is the axial DOF and unloaded.
        if (std::fabs(tip[0]) > 1e-14 || std::fabs(tip[1]) > 1e-14) {
            std::printf("FAIL C1: orthogonal DOFs nonzero u_x=%g u_y=%g\n",
                        tip[0], tip[1]);
            return 1;
        }
    }

    // --- C2: beam along +z, refVector = (0,0,1) ---
    {
        const Vec3 xA(0, 0, 0);
        const Vec3 xB(0, 0, L);
        const Vec3 ref(0, 0, 1);
        LocalFrame3D f = buildLocalFrame3D(xA, xB, ref);
        if (check_orthogonal(f.lambda, "C2 second-fallback frame")) return 1;
        // First fallback (0,0,1) is parallel to +z, so the second fallback
        // (0,1,0) is used: Vy = (0,1,0), Vz = (+z) x (+y) = (-x,0,0).
        Vec3 expVy(0, 1, 0);
        Vec3 expVz(-1, 0, 0);
        if ((f.Vy - expVy).norm() > 1e-12 || (f.Vz - expVz).norm() > 1e-12) {
            std::printf("FAIL C2: second-fallback frame mismatch Vy=(%g,%g,%g) "
                        "Vz=(%g,%g,%g)\n",
                        f.Vy[0], f.Vy[1], f.Vy[2],
                        f.Vz[0], f.Vz[1], f.Vz[2]);
            return 1;
        }

        // Apply global F_x at tip (transverse to beam axis +z). Local force =
        // lambda * (F_x, 0, 0) = (Vx[0], Vy[0], Vz[0]) * F_x = (0, 0, -F_x)
        // -> local F_z = -F_x, engaging local xz-plane bending with EIy.
        // Global u_x = (lambda^T (0,0,u_z_loc))[0] = Vz[0] * u_z_loc = -u_z_loc
        //           = -(-F_x L^3/(3 E Iy)) = F_x L^3/(3 E Iy).
        const double Fx = -1000.0;
        VecN<6> load = VecN<6>::Zero();
        load[0] = Fx;
        VecN<6> tip = solveSingleElement(xA, xB, ref, load, sec);
        const double ux_expected = Fx * L * L * L / (3.0 * sec.E * sec.Iy);
        if (check_rel(tip[0], ux_expected, 1e-9,
                      "C2 global u_x under global F_x")) return 1;
        if (std::fabs(tip[1]) > 1e-14 || std::fabs(tip[2]) > 1e-14) {
            std::printf("FAIL C2: orthogonal DOFs nonzero u_y=%g u_z=%g\n",
                        tip[1], tip[2]);
            return 1;
        }
    }

    std::printf("PASS subtest_fallback (C1 and C2)\n");
    return 0;
}

} // namespace

int main()
{
    if (subtest_local_force_rotation())            return 1;
    if (subtest_refvector_changes_bending_plane()) return 1;
    if (subtest_fallback())                        return 1;
    std::printf("PASS test_eb3d_transform (all subtests)\n");
    return 0;
}
