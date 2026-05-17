#pragma once
#include "../Core/Node.h"
#include <vector>

namespace beamlib {

template <int NDofsPerNode>
class DofMap {
public:
    void build(const std::vector<Node<NDofsPerNode>>& nodes) {
        const int N = static_cast<int>(nodes.size());
        idx_.assign(static_cast<size_t>(N) * NDofsPerNode, -1);
        int next = 0;
        for (int i = 0; i < N; ++i) {
            for (int d = 0; d < NDofsPerNode; ++d) {
                if (!nodes[i].fixed[d]) {
                    idx_[static_cast<size_t>(i) * NDofsPerNode + d] = next++;
                }
            }
        }
        numFree_ = next;
    }

    int numFree() const { return numFree_; }

    int freeDofIndex(int nodeIdx, int localDof) const {
        return idx_[static_cast<size_t>(nodeIdx) * NDofsPerNode + localDof];
    }

private:
    std::vector<int> idx_;
    int numFree_ = 0;
};

} // namespace beamlib
