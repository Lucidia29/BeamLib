#pragma once
#include "../Core/Types.h"
#include <vector>

namespace beamlib {

// Reaction force at a node, expressed per DOF in the global frame.
// Entries corresponding to free DOFs are zero. Sized by the model's
// nDofsPerNode (3 for 2D, 6 for 3D).
struct ReactionForce {
    int nodeId = -1;
    VecX R;
};

// AnalysisResult itself is intentionally NOT defined in Batch 2B.
// PROJECT_SPEC section 5.2 keeps it as a forward contract for the Phase 1
// API boundary (Batch 10). PostProcess helpers in this batch return the
// individual containers (ReactionForce, ElementInternalForces) directly so
// they can be composed into AnalysisResult later without churn.

} // namespace beamlib
