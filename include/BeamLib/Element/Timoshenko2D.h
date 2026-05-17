#pragma once
#include "ElementBase.h"
#include "../Core/SectionProperties.h"
#include "../Math/Rotation2D.h"
#include "../PostProcess/InternalForces.h"

namespace beamlib {

// 2D Timoshenko beam element for the BeamLib xz-plane convention.
//   DOFs per node: [u_x, u_z, theta_y]
//   Element DOF order: [u_xA, u_zA, theta_yA, u_xB, u_zB, theta_yB]
//
// Theory: shear-deformable Euler--Bernoulli generalization.
//   - axial strain:   eps   = du_x / dx
//   - bending strain: kappa = d(theta_y) / dx
//   - shear strain:   gamma = du_z / dx - theta_y
// EB enforces gamma = 0 (theta_y == du_z/dx). Timoshenko allows nonzero gamma
// and resists it with constitutive stiffness kappa*G*A. This makes Timoshenko
// correctly capture shear flexibility of deep beams, while still reducing to
// EB in the slender limit where the shear contribution vanishes.
//
// Stiffness: 2-node locking-free closed-form "Phi-corrected" element. Defining
//   Phi = 12 E I / (kappa G A L^2)
// the bending/shear block on indices (1, 2, 4, 5) is, in DOF order
// (u_zA, theta_yA, u_zB, theta_yB),
//
//   K_b = EI / (L^3 (1 + Phi)) *
//     [ 12       6L         -12       6L        ]
//     [ 6L     (4+Phi)L^2   -6L     (2-Phi)L^2  ]
//     [-12     -6L           12     -6L         ]
//     [ 6L     (2-Phi)L^2   -6L     (4+Phi)L^2  ]
//
// Properties:
//   - Phi -> 0 (slender beam): K_b reduces exactly to the standard EB Hermite
//     bending block (4 EI/L on the (theta, theta) diagonal, 2 EI/L on the
//     off-diagonal). Verified by test_timo2d_phi_limit.
//   - Phi large (deep beam): the (1 + Phi) denominator softens the bending
//     stiffness, capturing shear flexibility. Verified by test_timo2d_deep_beam
//     against the analytical Timoshenko tip deflection.
//   - The closed form is free of shear locking even on coarse meshes with
//     very large L/h, because the (1 + Phi) denominator scales correctly:
//     in the slender limit (Phi very small) the shear flexibility term in
//     the analytical deflection FL/(kappa G A) becomes negligible compared
//     to the bending term FL^3/(3EI), and the element matches EB. Verified
//     by test_timo2d_shear_locking.
//
// Sign convention: BeamLib structural convention theta_y = du_z/dx (same as
// EB2D, Chapter 1). The Hermite-style block has the same signs as EB2D
// because in the Phi -> 0 limit they must coincide. theta_y is the
// independent section rotation; the shear strain is then
//   gamma = du_z/dx - theta_y,
// equal to zero in the EB limit.
//
// Mass: consistent translational Hermite mass (same form as EB2D) plus an
// explicit rotary inertia contribution on the theta_y DOFs using linear shape
// functions. See computeMass() for the breakdown. The rotary inertia term is
// rho*I*L/6 * [[2,1],[1,2]] on the (theta_yA, theta_yB) DOFs and is the only
// difference from the EB2D mass.
struct Timoshenko2D {
    static constexpr int nDofsPerNode = 3;
    static constexpr bool hasTransformation = true;
    static constexpr bool isLinear = true;

    static ElementResult<3> computeElement(
        const Vec3& xA, const Vec3& xB,
        const VecN<6>& dispVec,
        const SectionProperties& props)
    {
        const double L  = (xB - xA).norm();
        const double L2 = L * L;
        const double L3 = L2 * L;

        const double EA  = props.E * props.A;
        const double EI  = props.E * props.Iz;
        const double GAk = props.kappa() * props.G * props.A;   // kappa_z * G * A
        const double Phi = 12.0 * EI / (GAk * L2);
        const double f   = 1.0 + Phi;                            // (1 + Phi) factor

        ElementResult<3> result;
        auto& K = result.ke;

        // --- Axial: indices (0, 3) ---
        K(0, 0) = K(3, 3) =  EA / L;
        K(0, 3) = K(3, 0) = -EA / L;

        // --- Phi-corrected bending/shear block on indices (1, 2, 4, 5) ---
        // DOF subset is (u_zA, theta_yA, u_zB, theta_yB).
        const double s   = EI / (L3 * f);
        const double k_uu_d = 12.0 * s;                          // u-u diagonal
        const double k_uu_o = -12.0 * s;                         // u_A - u_B
        const double k_ut_p =  6.0 * L * s;                      // positive u-theta coupling
        const double k_ut_n = -6.0 * L * s;                      // negative u-theta coupling
        const double k_tt_d = (4.0 + Phi) * L2 * s;              // theta-theta diagonal
        const double k_tt_o = (2.0 - Phi) * L2 * s;              // theta_A - theta_B

        K(1, 1) = K(4, 4) =  k_uu_d;
        K(1, 4) = K(4, 1) =  k_uu_o;

        K(1, 2) = K(2, 1) =  k_ut_p;          // u_zA - theta_yA
        K(1, 5) = K(5, 1) =  k_ut_p;          // u_zA - theta_yB
        K(2, 4) = K(4, 2) =  k_ut_n;          // theta_yA - u_zB
        K(4, 5) = K(5, 4) =  k_ut_n;          // u_zB - theta_yB

        K(2, 2) = K(5, 5) =  k_tt_d;
        K(2, 5) = K(5, 2) =  k_tt_o;

        result.re = K * dispVec;
        return result;
    }

    static MatMN<6, 6> computeTransformation(
        const Vec3& xA, const Vec3& xB,
        const Vec3& /*refVector*/)
    {
        return Rotation2D::compute(xA, xB);
    }

    // Consistent mass for Timoshenko 2D = consistent translational Hermite
    // mass + linear-rotary inertia.
    //
    // Translational part (ρA, on u_z deflection field) uses the standard
    // 4x4 Hermite mass block embedded on indices (1, 2, 4, 5), exactly the
    // same form as EB2D::computeMass. This captures kinetic energy
    // (1/2) ρA (du_z/dt)^2 with the cubic interpolation that matches the
    // Φ→0 EB limit.
    //
    // Rotary inertia part (ρI, on θ_y section-rotation field) uses linear
    // shape functions for θ_y. This adds:
    //   ρI*L/6 * [[2, 1], [1, 2]] on indices (2, 5).
    // This is the only mass-matrix difference between Timoshenko and EB2D.
    // It is required for any meaningful Timoshenko dynamic / modal analysis
    // (Batches 8/9): without it the rotational kinetic energy is missing
    // and high-frequency modes shift. PROJECT_SPEC §3.4 lists rotary inertia
    // as a required Timoshenko mass term.
    //
    // Axial part is identical to EB2D: ρA*L/6 * [[2,1],[1,2]] on (0, 3).
    //
    // The translational Hermite block does not depend on Φ. The "fully
    // consistent" Friedman--Kosmatka Φ-dependent mass is more elaborate and
    // is intentionally deferred to a later batch if higher-frequency mode
    // accuracy is required for thick beams. For slender to moderate-thick
    // beams in Phase 1, the present formulation is the standard textbook
    // teaching choice.
    static ElementMassResult<3> computeMass(
        const Vec3& xA, const Vec3& xB,
        const SectionProperties& props)
    {
        const double L     = (xB - xA).norm();
        const double L2    = L * L;
        const double rhoAL = props.rho * props.A  * L;
        const double rhoIL = props.rho * props.Iz * L;
        const double axial = rhoAL / 6.0;
        const double bend  = rhoAL / 420.0;
        const double rot   = rhoIL / 6.0;

        ElementMassResult<3> result;
        auto& M = result.me;

        // --- Axial translational mass on (0, 3) ---
        M(0, 0) = M(3, 3) = 2.0 * axial;
        M(0, 3) = M(3, 0) = 1.0 * axial;

        // --- Translational Hermite mass on (1, 2, 4, 5), same as EB2D ---
        M(1, 1) = M(4, 4) = 156.0 * bend;
        M(1, 4) = M(4, 1) =  54.0 * bend;
        M(1, 2) = M(2, 1) =  22.0 * L * bend;
        M(1, 5) = M(5, 1) = -13.0 * L * bend;
        M(2, 4) = M(4, 2) =  13.0 * L * bend;
        M(4, 5) = M(5, 4) = -22.0 * L * bend;
        M(2, 2) = M(5, 5) =   4.0 * L2 * bend;
        M(2, 5) = M(5, 2) =  -3.0 * L2 * bend;

        // --- Rotary inertia addition on (2, 5) for theta_y ---
        // Linear shape functions for theta_y give the standard 2x2 mass
        // pattern rho*I*L/6 * [[2,1],[1,2]], which we ADD on top of the
        // Hermite translational entries above.
        M(2, 2) += 2.0 * rot;
        M(5, 5) += 2.0 * rot;
        M(2, 5) += 1.0 * rot;
        M(5, 2) += 1.0 * rot;

        return result;
    }

    // Internal forces in the cross-section "beam diagram" convention.
    // Sign mapping mirrors EB2D::computeInternalForces exactly because both
    // elements use the same DOF order and the same structural sign convention
    // for theta_y. dispVec must be the LOCAL nodal displacement vector; the
    // assembly layer handles the global->local transformation when
    // hasTransformation == true.
    //
    //   N_A   = -f_e[0],   V_z_A = +f_e[1],   M_y_A = +f_e[2]
    //   N_B   = +f_e[3],   V_z_B = -f_e[4],   M_y_B = -f_e[5]
    //
    // See 01_euler_bernoulli_2d.tex §"Internal Forces" for the derivation of
    // the sign rule; the Timoshenko element shares the rule because the
    // residual semantics are identical.
    static ElementInternalForces computeInternalForces(
        const Vec3& xA, const Vec3& xB,
        const VecN<6>& dispVec,
        const SectionProperties& props)
    {
        ElementResult<3> r = computeElement(xA, xB, dispVec, props);
        ElementInternalForces f;
        f.N_A   = -r.re[0];
        f.V_z_A = +r.re[1];
        f.M_y_A = +r.re[2];
        f.N_B   = +r.re[3];
        f.V_z_B = -r.re[4];
        f.M_y_B = -r.re[5];
        return f;
    }
};

} // namespace beamlib
