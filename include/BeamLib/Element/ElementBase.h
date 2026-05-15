#pragma once
#include "../Core/Types.h"

namespace beamlib {

template <int NDofsPerNode>
struct ElementResult {
    static constexpr int elemDofs = 2 * NDofsPerNode;
    VecN<elemDofs> re = VecN<elemDofs>::Zero();
    MatMN<elemDofs, elemDofs> ke = MatMN<elemDofs, elemDofs>::Zero();
};

template <int NDofsPerNode>
struct ElementMassResult {
    static constexpr int elemDofs = 2 * NDofsPerNode;
    MatMN<elemDofs, elemDofs> me = MatMN<elemDofs, elemDofs>::Zero();
};

struct ElementConn {
    int nodeA = 0;
    int nodeB = 0;
    Vec3 refVector = Vec3(0, 1, 0);
};

} // namespace beamlib
