#pragma once
#include "ElementBase.h"
#include "../Core/SectionProperties.h"
#include "../Math/BeamMath3D.h"
#include "../PostProcess/InternalForces.h"

namespace beamlib {

// 3D Timoshenko beam element, 6 DOFs per node:
//   [u_x, u_y, u_z, theta_x, theta_y, theta_z]
// Element DOF layout (12 entries): [node A DOFs, node B DOFs].
//
// Composition: axial + St. Venant torsion + two uncoupled Timoshenko
// bending/shear planes, each carrying its own Phi correction.
//
//   xy plane  (DOFs 1, 5, 7, 11):
//       deflection u_y, rotation theta_z
//       bending stiffness     EIz   (I about local +z)
//       shear   stiffness     kappa_y * G * A
//       Phi_y = 12 EIz / (kappa_y * G * A * L^2)
//
//   xz plane  (DOFs 2, 4, 8, 10):
//       deflection u_z, rotation theta_y
//       bending stiffness     EIy   (I about local +y)
//       shear   stiffness     kappa_z * G * A
//       Phi_z = 12 EIy / (kappa_z * G * A * L^2)
//
// kappa_y and kappa_z stay independent; the element uses them separately
// for the two planes. The blockwise construction is valid because in linear
// theory the two bending/shear planes are uncoupled in the local frame
// (the axial, torsion, and two bending modes all carry independent strain
// measures and constitutive matrices); the cross-section transformation
// from local to global then mixes them through lambda.
//
// EB3D limit: both Phi -> 0 reduces the 12x12 ke to the EB3D ke entry-by-
// entry. See test_timo3d_phi_limit.
//
// Sign / convention notes:
//   - theta_y and theta_z follow the BeamLib structural convention used in
//     EB3D (Chapter 2 of the theory manual). In the EB limit
//       theta_y = du_z/dx,  theta_z = du_y/dx.
//   - In Timoshenko, theta_y and theta_z are independent rotational DOFs;
//     the shear strains are
//       gamma_xy = du_y/dx - theta_z   (xy plane)
//       gamma_xz = du_z/dx - theta_y   (xz plane)
//     so that the EB limit recovers the structural slopes.
//   - The rotation transformation (translations use lambda, rotations use
//     S*lambda*S, S = diag(1,-1,1)) is inherited from BeamMath3D and is the
//     same one EB3D uses. This is required for correct multi-element 3D
//     frame behavior; see EB3D theory document section
//     "Element-level transformation".
struct Timoshenko3D {
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

        const double EA   = props.E * props.A;
        const double GIx  = props.G * props.Ix;
        const double EIy  = props.E * props.Iy;
        const double EIz  = props.E * props.Iz;
        const double GAky = props.kappa_y * props.G * props.A;
        const double GAkz = props.kappa_z * props.G * props.A;

        // Phi_y belongs to the xy bending plane (u_y / theta_z, stiffness EIz,
        // shear stiffness kappa_y * G * A). Phi_z belongs to the xz plane
        // (u_z / theta_y, stiffness EIy, shear stiffness kappa_z * G * A).
        const double Phi_y = 12.0 * EIz / (GAky * L2);
        const double Phi_z = 12.0 * EIy / (GAkz * L2);
        const double f_y   = 1.0 + Phi_y;
        const double f_z   = 1.0 + Phi_z;

        ElementResult<6> result;
        auto& K = result.ke;

        // --- Axial: indices (0, 6) ---
        K(0, 0) = K(6, 6) =  EA / L;
        K(0, 6) = K(6, 0) = -EA / L;

        // --- Torsion: indices (3, 9) ---
        K(3, 3) = K(9, 9) =  GIx / L;
        K(3, 9) = K(9, 3) = -GIx / L;

        // --- xy-plane bending/shear: indices (1, 5, 7, 11), EIz, Phi_y ---
        // Phi_y-corrected closed-form Timoshenko block in DOF order
        // (u_yA, theta_zA, u_yB, theta_zB):
        //   K_b = EIz / (L^3 (1+Phi_y)) *
        //     [ 12       6L           -12       6L          ]
        //     [ 6L     (4+Phi_y)L^2   -6L     (2-Phi_y)L^2  ]
        //     [-12     -6L             12     -6L           ]
        //     [ 6L     (2-Phi_y)L^2   -6L     (4+Phi_y)L^2  ]
        {
            const double s = EIz / (L3 * f_y);
            const double k_uu_d =  12.0 * s;
            const double k_uu_o = -12.0 * s;
            const double k_ut_p =   6.0 * L * s;
            const double k_ut_n =  -6.0 * L * s;
            const double k_tt_d = (4.0 + Phi_y) * L2 * s;
            const double k_tt_o = (2.0 - Phi_y) * L2 * s;

            K(1, 1)  = K(7, 7)   =  k_uu_d;
            K(1, 7)  = K(7, 1)   =  k_uu_o;

            K(1, 5)  = K(5, 1)   =  k_ut_p;    // u_yA  - theta_zA
            K(1, 11) = K(11, 1)  =  k_ut_p;    // u_yA  - theta_zB
            K(5, 7)  = K(7, 5)   =  k_ut_n;    // theta_zA - u_yB
            K(7, 11) = K(11, 7)  =  k_ut_n;    // u_yB  - theta_zB

            K(5, 5)  = K(11, 11) =  k_tt_d;
            K(5, 11) = K(11, 5)  =  k_tt_o;
        }

        // --- xz-plane bending/shear: indices (2, 4, 8, 10), EIy, Phi_z ---
        // Same block form, with EIy and Phi_z, on the (u_z, theta_y) DOFs.
        {
            const double s = EIy / (L3 * f_z);
            const double k_uu_d =  12.0 * s;
            const double k_uu_o = -12.0 * s;
            const double k_ut_p =   6.0 * L * s;
            const double k_ut_n =  -6.0 * L * s;
            const double k_tt_d = (4.0 + Phi_z) * L2 * s;
            const double k_tt_o = (2.0 - Phi_z) * L2 * s;

            K(2, 2)  = K(8, 8)   =  k_uu_d;
            K(2, 8)  = K(8, 2)   =  k_uu_o;

            K(2, 4)  = K(4, 2)   =  k_ut_p;    // u_zA  - theta_yA
            K(2, 10) = K(10, 2)  =  k_ut_p;    // u_zA  - theta_yB
            K(4, 8)  = K(8, 4)   =  k_ut_n;    // theta_yA - u_zB
            K(8, 10) = K(10, 8)  =  k_ut_n;    // u_zB  - theta_yB

            K(4, 4)  = K(10, 10) =  k_tt_d;
            K(4, 10) = K(10, 4)  =  k_tt_o;
        }

        result.re = K * dispVec;
        return result;
    }

    static MatMN<12, 12> computeTransformation(
        const Vec3& xA, const Vec3& xB,
        const Vec3& refVector)
    {
        // Inherit the EB3D / Batch 3 convention directly: translations use
        // lambda, rotations use S*lambda*S (S = diag(1, -1, 1)). The shear
        // strain definitions above (gamma = u' - theta with the BeamLib
        // structural-convention theta_y / theta_z) are compatible with EB3D's
        // convention, so no extra adapter is needed beyond what is already in
        // BeamMath3D::buildTransformation3D.
        return buildTransformation3D(xA, xB, refVector);
    }

    // Consistent mass matrix for Timoshenko 3D.
    //
    // Components (all in the local element frame):
    //   - Axial translational mass on (0, 6): rho*A*L/6 * [[2,1],[1,2]]
    //   - Torsional rotary inertia on (3, 9): rho*Ix*L/6 * [[2,1],[1,2]]
    //   - xy-plane translational Hermite mass on (1, 5, 7, 11): standard
    //     rho*A*L/420 block (identical form to EB3D)
    //   - xz-plane translational Hermite mass on (2, 4, 8, 10): same form
    //   - Bending rotary inertia on theta_y (added to (4, 10)):
    //       rho*Iy*L/6 * [[2,1],[1,2]]
    //   - Bending rotary inertia on theta_z (added to (5, 11)):
    //       rho*Iz*L/6 * [[2,1],[1,2]]
    //
    // The Timoshenko 3D mass differs from EB3D in two places:
    //   - rho*Iy*L on theta_y (cross-section rotation about local +y)
    //   - rho*Iz*L on theta_z (cross-section rotation about local +z)
    // Both are required for any Timoshenko dynamic / modal analysis
    // (PROJECT_SPEC section 3.4). The torsional rho*Ix term is inherited
    // unchanged from EB3D, where it is also required.
    //
    // What is *not* included: Phi-dependent corrections to the translational
    // Hermite mass (Friedman-Kosmatka closed-form mass). For slender to
    // moderately thick beams the dominant correction over EB3D is the
    // bending rotary inertia above; the Phi-dependent translational terms
    // are a small refinement at very high L/h regimes that are not in
    // Phase 1's verification scope.
    static ElementMassResult<6> computeMass(
        const Vec3& xA, const Vec3& xB,
        const SectionProperties& props)
    {
        const double L      = (xB - xA).norm();
        const double L2     = L * L;
        const double rhoAL  = props.rho * props.A  * L;
        const double rhoIxL = props.rho * props.Ix * L;
        const double rhoIyL = props.rho * props.Iy * L;
        const double rhoIzL = props.rho * props.Iz * L;
        const double axial  = rhoAL  / 6.0;
        const double tors   = rhoIxL / 6.0;
        const double bend   = rhoAL  / 420.0;
        const double rot_y  = rhoIyL / 6.0;
        const double rot_z  = rhoIzL / 6.0;

        ElementMassResult<6> result;
        auto& M = result.me;

        // --- Axial translational mass: indices (0, 6) ---
        M(0, 0) = M(6, 6) = 2.0 * axial;
        M(0, 6) = M(6, 0) = 1.0 * axial;

        // --- Torsional rotary inertia: indices (3, 9) ---
        M(3, 3) = M(9, 9) = 2.0 * tors;
        M(3, 9) = M(9, 3) = 1.0 * tors;

        // --- xy bending Hermite translational mass: indices (1, 5, 7, 11) ---
        M(1, 1)  = M(7, 7)   = 156.0 * bend;
        M(1, 7)  = M(7, 1)   =  54.0 * bend;
        M(1, 5)  = M(5, 1)   =  22.0 * L  * bend;
        M(1, 11) = M(11, 1)  = -13.0 * L  * bend;
        M(5, 7)  = M(7, 5)   =  13.0 * L  * bend;
        M(7, 11) = M(11, 7)  = -22.0 * L  * bend;
        M(5, 5)  = M(11, 11) =   4.0 * L2 * bend;
        M(5, 11) = M(11, 5)  =  -3.0 * L2 * bend;

        // --- xz bending Hermite translational mass: indices (2, 4, 8, 10) ---
        M(2, 2)  = M(8, 8)   = 156.0 * bend;
        M(2, 8)  = M(8, 2)   =  54.0 * bend;
        M(2, 4)  = M(4, 2)   =  22.0 * L  * bend;
        M(2, 10) = M(10, 2)  = -13.0 * L  * bend;
        M(4, 8)  = M(8, 4)   =  13.0 * L  * bend;
        M(8, 10) = M(10, 8)  = -22.0 * L  * bend;
        M(4, 4)  = M(10, 10) =   4.0 * L2 * bend;
        M(4, 10) = M(10, 4)  =  -3.0 * L2 * bend;

        // --- Bending rotary inertia (Timoshenko-specific) ---
        // theta_y rotary inertia (rho*Iy) on indices (4, 10):
        M(4, 4)   += 2.0 * rot_y;
        M(10, 10) += 2.0 * rot_y;
        M(4, 10)  += 1.0 * rot_y;
        M(10, 4)  += 1.0 * rot_y;
        // theta_z rotary inertia (rho*Iz) on indices (5, 11):
        M(5, 5)   += 2.0 * rot_z;
        M(11, 11) += 2.0 * rot_z;
        M(5, 11)  += 1.0 * rot_z;
        M(11, 5)  += 1.0 * rot_z;

        return result;
    }

    // Full 3D internal forces at the two element ends, cross-section
    // ("beam diagram") convention.
    //
    // dispVec must be the LOCAL nodal displacement vector; the assembly layer
    // handles the global->local transformation when hasTransformation == true.
    //
    //   end A (section just right of node A):
    //     N_A   = -f_e[0]            V_y_A = +f_e[1]
    //     V_z_A = +f_e[2]            T_x_A = +f_e[3]
    //     M_y_A = +f_e[4]            M_z_A = +f_e[5]
    //   end B (section just left of node B):
    //     N_B   = +f_e[6]            V_y_B = -f_e[7]
    //     V_z_B = -f_e[8]            T_x_B = -f_e[9]
    //     M_y_B = -f_e[10]           M_z_B = -f_e[11]
    //
    // The pattern extends the EB2D rule unchanged: at A every component is
    // taken as the "left material applies on right material" force/moment,
    // with the axial sign flipped because "tension positive" is opposite to
    // "left-on-right axial" at A. At B the section's left material is the
    // element interior, so every non-axial component flips while N_B keeps
    // its sign.
    static ElementInternalForces computeInternalForces(
        const Vec3& xA, const Vec3& xB,
        const VecN<12>& dispVec,
        const SectionProperties& props)
    {
        ElementResult<6> r = computeElement(xA, xB, dispVec, props);
        ElementInternalForces f;
        // End A
        f.N_A   = -r.re[0];
        f.V_y_A = +r.re[1];
        f.V_z_A = +r.re[2];
        f.T_x_A = +r.re[3];
        f.M_y_A = +r.re[4];
        f.M_z_A = +r.re[5];
        // End B
        f.N_B   = +r.re[6];
        f.V_y_B = -r.re[7];
        f.V_z_B = -r.re[8];
        f.T_x_B = -r.re[9];
        f.M_y_B = -r.re[10];
        f.M_z_B = -r.re[11];
        return f;
    }
};

} // namespace beamlib
