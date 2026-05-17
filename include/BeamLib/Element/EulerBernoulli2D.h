#pragma once
#include "ElementBase.h"
#include "../Core/SectionProperties.h"
#include "../Math/Rotation2D.h"

namespace beamlib {

struct EulerBernoulli2D {
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
        const double EIz = props.E * props.Iz;

        ElementResult<3> result;
        auto& K = result.ke;

        // Axial DOFs at element indices (0, 3)
        K(0, 0) = K(3, 3) =  EA / L;
        K(0, 3) = K(3, 0) = -EA / L;

        // Bending DOFs at element indices (1, 2, 4, 5)
        K(1, 1) = K(4, 4) =  12.0 * EIz / L3;
        K(1, 4) = K(4, 1) = -12.0 * EIz / L3;
        K(1, 2) = K(2, 1) =   6.0 * EIz / L2;
        K(1, 5) = K(5, 1) =   6.0 * EIz / L2;
        K(2, 4) = K(4, 2) =  -6.0 * EIz / L2;
        K(4, 5) = K(5, 4) =  -6.0 * EIz / L2;
        K(2, 2) = K(5, 5) =   4.0 * EIz / L;
        K(2, 5) = K(5, 2) =   2.0 * EIz / L;

        result.re = K * dispVec;
        return result;
    }

    static MatMN<6, 6> computeTransformation(
        const Vec3& xA, const Vec3& xB,
        const Vec3& /*refVector*/)
    {
        return Rotation2D::compute(xA, xB);
    }

    static ElementMassResult<3> computeMass(
        const Vec3& xA, const Vec3& xB,
        const SectionProperties& props)
    {
        const double L     = (xB - xA).norm();
        const double L2    = L * L;
        const double rhoAL = props.rho * props.A * L;
        const double axial = rhoAL / 6.0;
        const double bend  = rhoAL / 420.0;

        ElementMassResult<3> result;
        auto& M = result.me;

        M(0, 0) = M(3, 3) =  2.0 * axial;
        M(0, 3) = M(3, 0) =  1.0 * axial;

        M(1, 1) = M(4, 4) = 156.0 * bend;
        M(1, 4) = M(4, 1) =  54.0 * bend;
        M(1, 2) = M(2, 1) =  22.0 * L * bend;
        M(1, 5) = M(5, 1) = -13.0 * L * bend;
        M(2, 4) = M(4, 2) =  13.0 * L * bend;
        M(4, 5) = M(5, 4) = -22.0 * L * bend;
        M(2, 2) = M(5, 5) =   4.0 * L2 * bend;
        M(2, 5) = M(5, 2) =  -3.0 * L2 * bend;

        return result;
    }
};

} // namespace beamlib
