#pragma once
#include "Types.h"

namespace beamlib {

struct SectionProperties {
    double E = 0.0;
    double G = 0.0;
    double rho = 0.0;
    double A = 0.0;
    double Iy = 0.0;
    double Iz = 0.0;
    double Ix = 0.0;
    double kappa_y = 5.0 / 6.0;
    double kappa_z = 5.0 / 6.0;

    Mat3 C_axialShear() const {
        Mat3 C = Mat3::Zero();
        C(0, 0) = E * A;
        C(1, 1) = kappa_y * G * A;
        C(2, 2) = kappa_z * G * A;
        return C;
    }

    Mat3 D_bendTorsion() const {
        Mat3 D = Mat3::Zero();
        D(0, 0) = G * Ix;
        D(1, 1) = E * Iy;
        D(2, 2) = E * Iz;
        return D;
    }

    double kappa() const { return kappa_z; }
};

} // namespace beamlib
