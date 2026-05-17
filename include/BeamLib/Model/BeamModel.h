#pragma once
#include "../Core/Node.h"
#include "../Core/SectionProperties.h"
#include "../Element/ElementBase.h"
#include "DofMap.h"
#include <array>
#include <vector>

namespace beamlib {

template <typename ElemType>
class BeamModel {
public:
    static constexpr int nDPN     = ElemType::nDofsPerNode;
    static constexpr int elemDofs = 2 * nDPN;

    std::vector<Node<nDPN>> nodes;
    std::vector<ElementConn> elements;
    SectionProperties props;

    void buildDofMap() { dofMap_.build(nodes); }

    int numFreeDofs() const { return dofMap_.numFree(); }

    int freeDofIndex(int nodeIdx, int localDof) const {
        return dofMap_.freeDofIndex(nodeIdx, localDof);
    }

    void assemble(VecX& R, SpMat& K) const {
        const int n = dofMap_.numFree();
        R = VecX::Zero(n);

        std::vector<Triplet> triplets;
        triplets.reserve(elements.size() * static_cast<size_t>(elemDofs) * elemDofs);

        for (const auto& elem : elements) {
            const auto& nA = nodes[elem.nodeA];
            const auto& nB = nodes[elem.nodeB];

            VecN<elemDofs> dispGlobal;
            dispGlobal.template segment<nDPN>(0)    = nA.dof;
            dispGlobal.template segment<nDPN>(nDPN) = nB.dof;

            VecN<elemDofs>              reElem;
            MatMN<elemDofs, elemDofs>   keElem;

            if constexpr (ElemType::hasTransformation) {
                MatMN<elemDofs, elemDofs> T =
                    ElemType::computeTransformation(nA.x0, nB.x0, elem.refVector);
                VecN<elemDofs> dispLocal = T * dispGlobal;
                const double L = (nB.x0 - nA.x0).norm();
                const Vec3 xA_loc = Vec3::Zero();
                const Vec3 xB_loc(L, 0.0, 0.0);
                ElementResult<nDPN> r =
                    ElemType::computeElement(xA_loc, xB_loc, dispLocal, props);
                reElem = T.transpose() * r.re;
                keElem = T.transpose() * r.ke * T;
            } else {
                ElementResult<nDPN> r =
                    ElemType::computeElement(nA.x0, nB.x0, dispGlobal, props);
                reElem = r.re;
                keElem = r.ke;
            }

            std::array<int, elemDofs> loc;
            for (int d = 0; d < nDPN; ++d) {
                loc[d]        = dofMap_.freeDofIndex(elem.nodeA, d);
                loc[d + nDPN] = dofMap_.freeDofIndex(elem.nodeB, d);
            }

            for (int i = 0; i < elemDofs; ++i) {
                const int li = loc[i];
                if (li < 0) continue;
                R[li] += reElem[i];
                for (int j = 0; j < elemDofs; ++j) {
                    const int lj = loc[j];
                    if (lj < 0) continue;
                    triplets.emplace_back(li, lj, keElem(i, j));
                }
            }
        }

        K.resize(n, n);
        K.setFromTriplets(triplets.begin(), triplets.end());
    }

    VecX getExternalForceVector() const {
        const int n = dofMap_.numFree();
        VecX F = VecX::Zero(n);
        for (size_t i = 0; i < nodes.size(); ++i) {
            for (int d = 0; d < nDPN; ++d) {
                int idx = dofMap_.freeDofIndex(static_cast<int>(i), d);
                if (idx >= 0) F[idx] = nodes[i].load[d];
            }
        }
        return F;
    }

    void scatterFreeDofs(const VecX& delta) {
        for (size_t i = 0; i < nodes.size(); ++i) {
            for (int d = 0; d < nDPN; ++d) {
                int idx = dofMap_.freeDofIndex(static_cast<int>(i), d);
                if (idx >= 0) nodes[i].dof[d] += delta[idx];
            }
        }
    }

    VecX gatherFreeDofs() const {
        const int n = dofMap_.numFree();
        VecX u = VecX::Zero(n);
        for (size_t i = 0; i < nodes.size(); ++i) {
            for (int d = 0; d < nDPN; ++d) {
                int idx = dofMap_.freeDofIndex(static_cast<int>(i), d);
                if (idx >= 0) u[idx] = nodes[i].dof[d];
            }
        }
        return u;
    }

private:
    DofMap<nDPN> dofMap_;
};

} // namespace beamlib
