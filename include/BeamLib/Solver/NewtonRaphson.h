#pragma once
#include "../Model/BeamModel.h"
#include "LinearSolver.h"

namespace beamlib {

struct NRConfig {
    double tol     = 1e-7;
    int    maxIter = 20;
};

struct NRResult {
    bool   converged     = false;
    int    iterations    = 0;
    double finalResidual = 0.0;
};

template <typename ElemType>
struct NewtonRaphsonSolver {
    static NRResult solveOneStep(
        BeamModel<ElemType>& model,
        const VecX&          Fext,
        const NRConfig&      config = {})
    {
        NRResult result;
        VecX  R_int;
        SpMat K;

        for (int it = 0; it < config.maxIter; ++it) {
            model.assemble(R_int, K);
            VecX R = R_int - Fext;
            const double normR = R.norm();
            result.finalResidual = normR;
            if (normR < config.tol) {
                result.converged  = true;
                result.iterations = it;
                return result;
            }
            SolveResult ls = LinearSolver::solve(K, -R);
            if (!ls.success) {
                result.converged  = false;
                result.iterations = it;
                return result;
            }
            model.scatterFreeDofs(ls.x);
            result.iterations = it + 1;
        }

        model.assemble(R_int, K);
        VecX R = R_int - Fext;
        result.finalResidual = R.norm();
        result.converged     = (result.finalResidual < config.tol);
        return result;
    }
};

} // namespace beamlib
