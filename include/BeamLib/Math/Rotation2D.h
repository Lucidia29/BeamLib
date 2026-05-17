#pragma once
#include "../Core/Types.h"

namespace beamlib {

struct Rotation2D {
    static MatMN<6, 6> compute(const Vec3& xA, const Vec3& xB) {
        const Vec3 dx = xB - xA;
        const double L = dx.norm();
        const double c = dx[0] / L;
        const double s = dx[2] / L;

        MatMN<6, 6> T = MatMN<6, 6>::Zero();
        for (int blk = 0; blk < 2; ++blk) {
            const int o = blk * 3;
            T(o + 0, o + 0) =  c;
            T(o + 0, o + 1) = -s;
            T(o + 1, o + 0) =  s;
            T(o + 1, o + 1) =  c;
            T(o + 2, o + 2) =  1.0;
        }
        return T;
    }
};

} // namespace beamlib
