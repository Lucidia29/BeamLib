#pragma once
#include "../Core/Types.h"
#include <vector>

namespace beamlib {

// Element-end internal forces in the cross-section ("beam diagram") convention.
//
// Common to 2D and 3D elements. 2D elements (EB2D, Timoshenko 2D) populate
// only the {N, V_z, M_y} subset; the 3D-only members V_y, T_x, M_z stay at
// their default zero value. 3D elements (Timoshenko 3D, plus EB3D once its
// post-processing path is filled in) populate all six.
//
// Sign convention (consistent across 2D and 3D):
//   N > 0   axial tension (member stretched along its local x axis).
//   V_y > 0 the left portion of the cross section exerts +y force on the
//           right portion (constant-shear loading therefore gives the same
//           V_y sign at both ends of every element along the beam).
//   V_z > 0 same idea, force component along local +z.
//   T_x > 0 left portion exerts a +x moment on the right portion (St. Venant
//           torsion, right-hand rule about local +x).
//   M_y > 0 left portion exerts a +y moment on the right portion (work-
//           conjugate to the BeamLib structural rotation theta_y = du_z/dx).
//   M_z > 0 same idea for the moment about local +z (work-conjugate to
//           theta_z = du_y/dx).
//
// Conversion from local element nodal forces f_e = K_l * d_local to cross-
// section values follows the same A/B endpoint rule for all elements:
//   end A (section just right of A): "left = node A boundary" so f_e
//                                     entries map directly except axial
//                                     which flips because tension is the
//                                     opposite of "left-on-right axial".
//   end B (section just left of B):  "left = element interior" so all
//                                     entries flip, with the same axial
//                                     exception leaving N_B without a flip.
// Detailed derivation: docs/theory/01_euler_bernoulli_2d.tex section
// "Internal Forces"; the 3D extension simply pairs the new V_y, T_x, M_z
// fields with the corresponding local DOF indices.
struct ElementInternalForces {
    int elementId = -1;

    // Forces / moments at end A (cross-section just right of node A)
    double N_A   = 0.0;
    double V_y_A = 0.0;
    double V_z_A = 0.0;
    double T_x_A = 0.0;
    double M_y_A = 0.0;
    double M_z_A = 0.0;

    // Forces / moments at end B (cross-section just left of node B)
    double N_B   = 0.0;
    double V_y_B = 0.0;
    double V_z_B = 0.0;
    double T_x_B = 0.0;
    double M_y_B = 0.0;
    double M_z_B = 0.0;
};

} // namespace beamlib
