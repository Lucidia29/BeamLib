#pragma once
#include "Types.h"

namespace beamlib {

template <int NDofsPerNode>
struct Node {
    static constexpr int nDofs = NDofsPerNode;
    Vec3 x0 = Vec3::Zero();
    VecN<NDofsPerNode> dof = VecN<NDofsPerNode>::Zero();
    std::array<bool, NDofsPerNode> fixed = {};
    VecN<NDofsPerNode> load = VecN<NDofsPerNode>::Zero();
    // Prescribed displacement for fixed DOFs. Used when fixed[d] == true.
    // Zero gives homogeneous Dirichlet; nonzero triggers the elimination + RHS
    // correction path in BeamModel::assemble (Batch 2B).
    VecN<NDofsPerNode> prescribed = VecN<NDofsPerNode>::Zero();

    Node() = default;
    explicit Node(const Vec3& x0_) : x0(x0_) {}
    void fixAll() { fixed.fill(true); }

    double totalDof(int d) const {
        return fixed[d] ? prescribed[d] : dof[d];
    }
};

using Node2D = Node<3>;
using Node3D = Node<6>;

} // namespace beamlib
