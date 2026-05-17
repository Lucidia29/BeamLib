#pragma once
#include "ElementBase.h"
#include "../Core/SectionProperties.h"
#include "../Math/BeamMath3D.h"
#include "../PostProcess/InternalForces.h"

namespace beamlib {

// 3D Euler-Bernoulli beam element, 6 DOFs per node:
//   [u_x, u_y, u_z, theta_x, theta_y, theta_z]
// Element DOF layout (12 entries): [node A DOFs, node B DOFs].
//
// Local frame convention (built by BeamMath3D::buildLocalFrame3D):
//   x_loc — along beam axis (xB - xA normalized)
//   y_loc — Gram-Schmidt projection of refVector onto plane normal to x_loc
//   z_loc — x_loc x y_loc (right-handed)
//
// Local bending rotation conventions (structural, matching EB2D):
//   theta_y_loc = d u_z_loc / d x_loc   (xz-plane bending, stiffness E*Iy)
//   theta_z_loc = d u_y_loc / d x_loc   (xy-plane bending, stiffness E*Iz)
// Both planes use the standard 4x4 Hermite block. The xy plane is also right-
// hand-rule consistent under this definition; the xz plane uses the structural
// convention deliberately (see docs/theory/02_euler_bernoulli_3d.tex).
//
// Torsion theta_x uses a 2-node linear interpolation, stiffness G*Ix/L.
//
// DOF index map (local element indices):
//   axial   (u_x_A, u_x_B)                             -> (0, 6)
//   torsion (theta_x_A, theta_x_B)                     -> (3, 9)
//   xy plane (u_y_A, theta_z_A, u_y_B, theta_z_B)      -> (1, 5, 7, 11)  uses Iz
//   xz plane (u_z_A, theta_y_A, u_z_B, theta_y_B)      -> (2, 4, 8, 10)  uses Iy
struct EulerBernoulli3D {
    static constexpr int nDofsPerNode = 6;
    static constexpr bool hasTransformation = true;
    static constexpr bool isLinear = true;

    static ElementResult<6> computeElement(
        const Vec3& xA, const Vec3& xB,
        const VecN<12>& dispVec,
        const SectionProperties& props)
    {
        const double L  = (xB - xA).norm();
        const double L2 = L * L;
        const double L3 = L2 * L;

        const double EA  = props.E * props.A;
        const double GIx = props.G * props.Ix;
        const double EIy = props.E * props.Iy;
        const double EIz = props.E * props.Iz;

        ElementResult<6> result;
        auto& K = result.ke;

        // --- Axial: indices (0, 6) ---
        K(0, 0) = K(6, 6) =  EA / L;
        K(0, 6) = K(6, 0) = -EA / L;

        // --- Torsion: indices (3, 9) ---
        K(3, 3) = K(9, 9) =  GIx / L;
        K(3, 9) = K(9, 3) = -GIx / L;

        // --- xy-plane bending: indices (1, 5, 7, 11), uses EIz ---
        // Hermite block in DOF order (u_yA, theta_zA, u_yB, theta_zB):
        //   [ 12     6L    -12     6L ]
        //   [ 6L    4L^2   -6L    2L^2]  * EIz/L^3
        //   [-12    -6L     12    -6L ]
        //   [ 6L    2L^2   -6L    4L^2]
        {
            const double c_uu_d  =  12.0 * EIz / L3;       // diagonal u-u
            const double c_uu_o  = -12.0 * EIz / L3;       // u_A-u_B
            const double c_ut_p  =   6.0 * EIz / L2;       // positive u-theta couplings
            const double c_ut_n  =  -6.0 * EIz / L2;       // negative u-theta couplings
            const double c_tt_d  =   4.0 * EIz / L;        // diagonal theta-theta
            const double c_tt_o  =   2.0 * EIz / L;        // theta_A-theta_B

            K(1, 1)  = K(7, 7)   =  c_uu_d;
            K(1, 7)  = K(7, 1)   =  c_uu_o;

            K(1, 5)  = K(5, 1)   =  c_ut_p;   // u_yA  - theta_zA
            K(1, 11) = K(11, 1)  =  c_ut_p;   // u_yA  - theta_zB
            K(5, 7)  = K(7, 5)   =  c_ut_n;   // theta_zA - u_yB
            K(7, 11) = K(11, 7)  =  c_ut_n;   // u_yB  - theta_zB

            K(5, 5)  = K(11, 11) =  c_tt_d;
            K(5, 11) = K(11, 5)  =  c_tt_o;
        }

        // --- xz-plane bending: indices (2, 4, 8, 10), uses EIy ---
        // Hermite block in DOF order (u_zA, theta_yA, u_zB, theta_yB),
        // identical to EB2D bending block (with EIz replaced by EIy):
        //   [ 12     6L    -12     6L ]
        //   [ 6L    4L^2   -6L    2L^2]  * EIy/L^3
        //   [-12    -6L     12    -6L ]
        //   [ 6L    2L^2   -6L    4L^2]
        {
            const double c_uu_d  =  12.0 * EIy / L3;
            const double c_uu_o  = -12.0 * EIy / L3;
            const double c_ut_p  =   6.0 * EIy / L2;
            const double c_ut_n  =  -6.0 * EIy / L2;
            const double c_tt_d  =   4.0 * EIy / L;
            const double c_tt_o  =   2.0 * EIy / L;

            K(2, 2)  = K(8, 8)   =  c_uu_d;
            K(2, 8)  = K(8, 2)   =  c_uu_o;

            K(2, 4)  = K(4, 2)   =  c_ut_p;   // u_zA  - theta_yA
            K(2, 10) = K(10, 2)  =  c_ut_p;   // u_zA  - theta_yB
            K(4, 8)  = K(8, 4)   =  c_ut_n;   // theta_yA - u_zB
            K(8, 10) = K(10, 8)  =  c_ut_n;   // u_zB  - theta_yB

            K(4, 4)  = K(10, 10) =  c_tt_d;
            K(4, 10) = K(10, 4)  =  c_tt_o;
        }

        result.re = K * dispVec;
        return result;
    }

    static MatMN<12, 12> computeTransformation(
        const Vec3& xA, const Vec3& xB,
        const Vec3& refVector)
    {
        return buildTransformation3D(xA, xB, refVector);
    }

    // Consistent mass matrix for EB3D.
    //
    // Components (all in the local element frame; see theory doc):
    //   - Axial translational mass (u_x):     m = rho*A*L/6 * [[2,1],[1,2]]
    //   - Torsional rotary inertia (theta_x): m = rho*Ix*L/6 * [[2,1],[1,2]]
    //   - xy bending plane (u_y, theta_z):    standard 4x4 Hermite mass block
    //                                          with rho*A*L/420 scaling
    //   - xz bending plane (u_z, theta_y):    same 4x4 Hermite block
    //
    // The torsional rotary inertia rho*Ix is *required* for EB3D: without it,
    // theta_x is mass-less and the torsional mode is missing from any modal /
    // dynamic analysis. PROJECT_SPEC section 3.4 documents this. Bending rotary
    // inertia (extra terms in the theta_y / theta_z diagonals) is *excluded*
    // for EB; the rotation DOFs still have nonzero mass entries because the
    // Hermite shape functions used for translational kinetic energy depend on
    // them, but no separate rho*Iy / rho*Iz term is added. (Timoshenko adds
    // these in Batch 4/5.)
    static ElementMassResult<6> computeMass(
        const Vec3& xA, const Vec3& xB,
        const SectionProperties& props)
    {
        const double L     = (xB - xA).norm();
        const double L2    = L * L;
        const double rhoAL = props.rho * props.A  * L;
        const double rhoIxL= props.rho * props.Ix * L;
        const double axial = rhoAL  / 6.0;
        const double tors  = rhoIxL / 6.0;
        const double bend  = rhoAL  / 420.0;

        ElementMassResult<6> result;
        auto& M = result.me;

        // --- Axial translational mass: indices (0, 6) ---
        M(0, 0) = M(6, 6) = 2.0 * axial;
        M(0, 6) = M(6, 0) = 1.0 * axial;

        // --- Torsional rotary inertia: indices (3, 9) ---
        M(3, 3) = M(9, 9) = 2.0 * tors;
        M(3, 9) = M(9, 3) = 1.0 * tors;

        // --- xy bending plane: indices (1, 5, 7, 11) ---
        // Standard Hermite mass block in order (u_yA, theta_zA, u_yB, theta_zB):
        //   [156   22L    54   -13L  ]
        //   [22L   4L^2   13L  -3L^2 ] * rho*A*L/420
        //   [54    13L    156  -22L  ]
        //   [-13L  -3L^2  -22L  4L^2 ]
        M(1, 1)  = M(7, 7)   = 156.0 * bend;
        M(1, 7)  = M(7, 1)   =  54.0 * bend;
        M(1, 5)  = M(5, 1)   =  22.0 * L  * bend;
        M(1, 11) = M(11, 1)  = -13.0 * L  * bend;
        M(5, 7)  = M(7, 5)   =  13.0 * L  * bend;
        M(7, 11) = M(11, 7)  = -22.0 * L  * bend;
        M(5, 5)  = M(11, 11) =   4.0 * L2 * bend;
        M(5, 11) = M(11, 5)  =  -3.0 * L2 * bend;

        // --- xz bending plane: indices (2, 4, 8, 10) ---
        // Same Hermite mass block in order (u_zA, theta_yA, u_zB, theta_yB).
        M(2, 2)  = M(8, 8)   = 156.0 * bend;
        M(2, 8)  = M(8, 2)   =  54.0 * bend;
        M(2, 4)  = M(4, 2)   =  22.0 * L  * bend;
        M(2, 10) = M(10, 2)  = -13.0 * L  * bend;
        M(4, 8)  = M(8, 4)   =  13.0 * L  * bend;
        M(8, 10) = M(10, 8)  = -22.0 * L  * bend;
        M(4, 4)  = M(10, 10) =   4.0 * L2 * bend;
        M(4, 10) = M(10, 4)  =  -3.0 * L2 * bend;

        return result;
    }

    // Internal force postprocessing for EB3D is DEFERRED in Batch 3.
    //
    // The Phase 1 ElementInternalForces struct currently exposes only axial N,
    // xz-plane shear V_z, and xz-plane moment M_y -- a 2D-shaped output. A
    // future batch should extend it (or add an EB3D-specific variant) to also
    // carry the xy-plane shear V_y, the torsional moment T, and the xy-plane
    // moment M_z. To keep the element-trait contract satisfied so that any
    // BeamModel<EulerBernoulli3D> instantiation compiles cleanly, this stub
    // populates the available fields from the local element nodal forces
    // (axial + xz plane only) and leaves V_y, M_z, T out of scope. Tests in
    // this batch do not exercise this path; reaction recovery, which is the
    // only consumer of internal forces beyond direct test access, lives at the
    // BeamModel level and uses the full element residual, not these fields.
    //
    // The xz-plane (V_z, M_y) sign mapping is identical to EB2D's because the
    // xz-plane bending block matches EB2D's bending block.
    static ElementInternalForces computeInternalForces(
        const Vec3& xA, const Vec3& xB,
        const VecN<12>& dispVec,
        const SectionProperties& props)
    {
        ElementResult<6> r = computeElement(xA, xB, dispVec, props);
        ElementInternalForces f;
        f.N_A   = -r.re[0];
        f.V_z_A = +r.re[2];
        f.M_y_A = +r.re[4];
        f.N_B   = +r.re[6];
        f.V_z_B = -r.re[8];
        f.M_y_B = -r.re[10];
        return f;
    }
};

} // namespace beamlib
