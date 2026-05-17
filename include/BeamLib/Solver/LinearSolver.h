#pragma once
#include "../Core/Types.h"
#include <Eigen/SparseLU>

namespace beamlib {

struct SolveResult {
    VecX x;
    bool success = false;
};

struct LinearSolver {
    static SolveResult solve(const SpMat& K, const VecX& rhs) {
        SolveResult r;
        Eigen::SparseLU<SpMat> lu;
        lu.compute(K);
        if (lu.info() != Eigen::Success) {
            r.success = false;
            return r;
        }
        r.x = lu.solve(rhs);
        if (lu.info() != Eigen::Success) {
            r.success = false;
            return r;
        }
        r.success = true;
        return r;
    }
};

} // namespace beamlib
