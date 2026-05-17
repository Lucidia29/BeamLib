#pragma once
#include "../Core/Types.h"
#include <vector>

namespace beamlib {

// 2D EB beam internal forces at each end of an element, cross-section
// ("beam diagram") convention. With this convention:
//   N > 0   axial tension (member stretched along its axis).
//   V_z > 0 the left portion of the cross section exerts +z force on the right
//           portion (so for a constant-shear loading V_z has the same sign at
//           both ends of every element along the beam).
//   M_y > 0 left portion exerts a +y moment on the right portion at the section
//           (work-conjugate to the BeamLib structural rotation theta_y = dw/dx).
//
// Conversion from element nodal forces f_e = K_e * d_local (the FE residual at
// the two nodes) to cross-section values is documented in
// docs/theory/01_euler_bernoulli_2d.tex (section "Internal Forces").
struct ElementInternalForces {
    int elementId = -1;
    double N_A = 0.0;
    double V_z_A = 0.0;
    double M_y_A = 0.0;
    double N_B = 0.0;
    double V_z_B = 0.0;
    double M_y_B = 0.0;
};

} // namespace beamlib
