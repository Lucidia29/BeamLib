#pragma once
#include "../Core/Node.h"
#include "../Core/SectionProperties.h"
#include "../Element/ElementBase.h"
#include "../PostProcess/AnalysisResult.h"
#include "../PostProcess/InternalForces.h"
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

    // Compute one element's global-frame residual r_g and tangent k_g.
    // dispGlobal is collected so that fixed DOFs carry their prescribed value
    // (zero by default) and free DOFs carry the current dof iterate. This
    // means r_g already contains the K_fr * u_r RHS-correction contribution
    // for nonzero prescribed displacement (Batch 2B requirement).
    void computeElementContribution_(
        const ElementConn& elem,
        VecN<elemDofs>& reElem,
        MatMN<elemDofs, elemDofs>& keElem) const
    {
        const auto& nA = nodes[elem.nodeA];
        const auto& nB = nodes[elem.nodeB];

        VecN<elemDofs> dispGlobal;
        for (int d = 0; d < nDPN; ++d) {
            dispGlobal[d]        = nA.totalDof(d);
            dispGlobal[d + nDPN] = nB.totalDof(d);
        }

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
    }

    // Free-DOF assembly: R is the internal residual restricted to free DOFs
    // (size = numFreeDofs), K is the corresponding tangent. For nonzero
    // prescribed displacement at fixed DOFs, R[i] for a free DOF i already
    // contains the K_fr * u_r contribution (because dispGlobal in the element
    // loop carries the prescribed values).
    void assemble(VecX& R, SpMat& K) const {
        const int n = dofMap_.numFree();
        R = VecX::Zero(n);

        std::vector<Triplet> triplets;
        triplets.reserve(elements.size() * static_cast<size_t>(elemDofs) * elemDofs);

        for (const auto& elem : elements) {
            VecN<elemDofs>            reElem;
            MatMN<elemDofs, elemDofs> keElem;
            computeElementContribution_(elem, reElem, keElem);

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

    // Full-system internal residual, indexed by global DOF (numNodes * nDPN).
    // Used for reaction recovery: at fixed DOF i, reaction = R_full[i] - F_ext[i].
    VecX assembleFullResidual() const {
        const int N = static_cast<int>(nodes.size());
        VecX Rfull = VecX::Zero(static_cast<Eigen::Index>(N) * nDPN);

        for (const auto& elem : elements) {
            VecN<elemDofs>            reElem;
            MatMN<elemDofs, elemDofs> keElem;
            computeElementContribution_(elem, reElem, keElem);

            const int baseA = elem.nodeA * nDPN;
            const int baseB = elem.nodeB * nDPN;
            for (int d = 0; d < nDPN; ++d) {
                Rfull[baseA + d] += reElem[d];
                Rfull[baseB + d] += reElem[d + nDPN];
            }
        }
        return Rfull;
    }

    // Consistent mass matrix on free DOFs. Mirrors assemble(), but uses
    // ElemType::computeMass and applies the same T^T M_l T transformation
    // for elements with hasTransformation == true.
    void assembleMass(SpMat& M) const {
        const int n = dofMap_.numFree();
        std::vector<Triplet> triplets;
        triplets.reserve(elements.size() * static_cast<size_t>(elemDofs) * elemDofs);

        for (const auto& elem : elements) {
            const auto& nA = nodes[elem.nodeA];
            const auto& nB = nodes[elem.nodeB];

            MatMN<elemDofs, elemDofs> meElem;
            if constexpr (ElemType::hasTransformation) {
                MatMN<elemDofs, elemDofs> T =
                    ElemType::computeTransformation(nA.x0, nB.x0, elem.refVector);
                const double L = (nB.x0 - nA.x0).norm();
                const Vec3 xA_loc = Vec3::Zero();
                const Vec3 xB_loc(L, 0.0, 0.0);
                ElementMassResult<nDPN> mr =
                    ElemType::computeMass(xA_loc, xB_loc, props);
                meElem = T.transpose() * mr.me * T;
            } else {
                ElementMassResult<nDPN> mr =
                    ElemType::computeMass(nA.x0, nB.x0, props);
                meElem = mr.me;
            }

            std::array<int, elemDofs> loc;
            for (int d = 0; d < nDPN; ++d) {
                loc[d]        = dofMap_.freeDofIndex(elem.nodeA, d);
                loc[d + nDPN] = dofMap_.freeDofIndex(elem.nodeB, d);
            }

            for (int i = 0; i < elemDofs; ++i) {
                const int li = loc[i];
                if (li < 0) continue;
                for (int j = 0; j < elemDofs; ++j) {
                    const int lj = loc[j];
                    if (lj < 0) continue;
                    triplets.emplace_back(li, lj, meElem(i, j));
                }
            }
        }

        M.resize(n, n);
        M.setFromTriplets(triplets.begin(), triplets.end());
    }

    // Per-element internal forces (raw two-end values from the element static
    // method, in cross-section convention). Element-local dispVec is built
    // using the same prescribed/free split as the residual assembly, so this
    // is consistent with the solved state.
    std::vector<ElementInternalForces> computeInternalForces() const {
        std::vector<ElementInternalForces> result;
        result.reserve(elements.size());

        for (size_t e = 0; e < elements.size(); ++e) {
            const auto& elem = elements[e];
            const auto& nA = nodes[elem.nodeA];
            const auto& nB = nodes[elem.nodeB];

            VecN<elemDofs> dispGlobal;
            for (int d = 0; d < nDPN; ++d) {
                dispGlobal[d]        = nA.totalDof(d);
                dispGlobal[d + nDPN] = nB.totalDof(d);
            }

            ElementInternalForces ifs;
            if constexpr (ElemType::hasTransformation) {
                MatMN<elemDofs, elemDofs> T =
                    ElemType::computeTransformation(nA.x0, nB.x0, elem.refVector);
                VecN<elemDofs> dispLocal = T * dispGlobal;
                const double L = (nB.x0 - nA.x0).norm();
                const Vec3 xA_loc = Vec3::Zero();
                const Vec3 xB_loc(L, 0.0, 0.0);
                ifs = ElemType::computeInternalForces(xA_loc, xB_loc, dispLocal, props);
            } else {
                ifs = ElemType::computeInternalForces(nA.x0, nB.x0, dispGlobal, props);
            }
            ifs.elementId = static_cast<int>(e);
            result.push_back(ifs);
        }
        return result;
    }

    // Reaction at every node, per global DOF, in the global frame.
    // For each fixed DOF: R[d] = R_int_full[d] - F_ext[d]. F_ext at a fixed
    // DOF is taken from node.load[d] (loads applied at constrained DOFs are
    // valid input and must be subtracted from the residual). Free-DOF
    // entries are left at zero.
    std::vector<ReactionForce> computeReactions() const {
        VecX Rfull = assembleFullResidual();
        std::vector<ReactionForce> reactions;
        reactions.reserve(nodes.size());

        for (size_t i = 0; i < nodes.size(); ++i) {
            const auto& n = nodes[i];
            ReactionForce rf;
            rf.nodeId = static_cast<int>(i);
            rf.R = VecX::Zero(nDPN);
            const int base = static_cast<int>(i) * nDPN;
            for (int d = 0; d < nDPN; ++d) {
                if (n.fixed[d]) {
                    rf.R[d] = Rfull[base + d] - n.load[d];
                }
            }
            reactions.push_back(rf);
        }
        return reactions;
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
