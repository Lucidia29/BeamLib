# BeamLib Planning Review

Reviewer role: numerical methods reviewer and engineering software architect.

Reviewed documents:
- `PROJECT_SPEC.md`
- `EXECUTION_PLAN.md`
- `REVIEW_NOTES.md`

## Main Findings

1. `EXECUTION_PLAN.md` Batch 2 is too large and high-risk.

   EB2D currently establishes `DofMap`, `BeamModel`, sparse assembly, NR, mass matrix, reactions, internal force post-processing, distributed-load tests, and rotated-beam tests in one batch. This is the core architecture layer. Split it into:
   - `Batch 2A`: EB2D + DofMap + BeamModel + linear/NR static solve
   - `Batch 2B`: mass assembly + reactions + internal forces + distributed-load equivalent nodal loads

   Otherwise debugging will not clearly isolate whether failures come from the element, assembly, constraints, post-processing, or load equivalencing.

2. `EXECUTION_PLAN.md` Batch 6 GE3D small-load benchmark is still ambiguous.

   It states `u_z ≈ FL^3/(3EI) (Timoshenko analytical solution, including shear correction)`, but `FL^3/(3EI)` is the Euler-Bernoulli bending term only. Write the complete expression explicitly:
   `u_z = FL^3/(3EI) + FL/(kappa GA)`, with sign set by load direction.

3. `EXECUTION_PLAN.md` Batch 2 uses a distributed-load test before the distributed-load API exists.

   The simply supported beam test with uniform load needs equivalent nodal loads. The plan must state whether Batch 2 manually assembles consistent equivalent nodal loads in the test, or introduces a minimal helper such as `equivalentNodalLoadUniform2D()`. Do not make Batch 2 depend on Batch 10 `LoadManager`.

4. `EXECUTION_PLAN.md` Batch 2 rotated-beam test needs a more rigorous acceptance definition.

   "Consistent with the horizontal beam projection" is not physically precise. A global vertical force on a 45-degree beam decomposes into local axial and transverse components. Use two tests:
   - Apply a local transverse force transformed into global coordinates; transform displacement back to local coordinates and compare with the horizontal beam.
   - Apply a global force and explicitly verify the expected local axial/transverse coupling.

5. Timoshenko shear locking is flagged, but the locking-free scheme is not decided.

   For a teaching-oriented Phase 1 library, use the 2-node locking-free closed-form Timoshenko stiffness with the `Phi = 12EI/(kappa GA L^2)` correction as the default. It is simple, teachable, and has a clear EB limit. Do not use EAS in the first version. Selective reduced integration can be documented as a comparison, but it requires more care around low-order element modes.

6. GE3D one-point Gauss quadrature is acceptable for migration, but not generally sufficient as a claim.

   Keep 1-point Gauss in Phase 1 to preserve the existing validated C++ behavior. Add a note that multi-point integration should be evaluated later for strong curvature, distributed loading, or shear/torsion coupling cases. Do not mix algorithm migration with quadrature changes in the same batch.

7. The global `< 1e-10` relative tolerance for all linear benchmarks is too aggressive.

   Use different tolerances by test type:
   - Single-element matrix checks and patch tests: `1e-12` to `1e-10`
   - Multi-element analytical displacement checks: `1e-9` or `1e-8`
   - Reactions and near-zero DOFs: absolute tolerance, not relative tolerance

8. `computeMass` is required by dynamics and modal analysis, but the mass strategy is not fixed.

   EB consistent mass is clear. Timoshenko mass must decide whether rotary inertia and shear-related inertia are included. GE3D must decide consistent vs lumped mass. These choices affect Batch 8 and Batch 9 and should be decided before implementation reaches dynamics/modal analysis.

9. Prescribed displacement handling is too late in the plan.

   Constraint handling is part of core assembly and solving, not final integration. Homogeneous Dirichlet constraints can be enough for Batch 2A, but if nonzero prescribed displacement is in Phase 1 scope, elimination plus RHS correction must be designed early.

10. Dense modal analysis is acceptable for Phase 1, but the plan needs guardrails.

    `Eigen::GeneralizedSelfAdjointEigenSolver` is fine for educational and moderate models, but implementation must:
    - assemble free-DOF dense `Kff` and `Mff`
    - ensure `Mff` is symmetric positive definite
    - handle insufficient constraints and rigid-body modes
    - document that large sparse modal analysis is out of Phase 1 scope

## Numerical Methods Correctness

The 6-element scope is reasonable:
- EB for slender beams
- Timoshenko for shear-deformable beams
- GE for large deformation and large rotation

The missing piece is a clearer applicability boundary for each theory:
- whether Timoshenko includes rotary inertia
- whether GE is explicitly a shear-deformable Reissner-Simo beam
- whether EB fully excludes shear effects in all static and dynamic checks

Batch 2 EB2D stiffness indexing is correct for `[u_x, u_z, theta_y]`. The bending sub-block uses local DOFs `[1, 2, 4, 5]`, which is standard if the positive `theta_y` convention is consistently defined.

Batch 3 EB3D mapping is also correct:
- xy-plane bending uses `u_y/theta_z`, DOFs `[1, 5, 7, 11]`, stiffness `EIz`
- xz-plane bending uses `u_z/theta_y`, DOFs `[2, 4, 8, 10]`, stiffness `EIy`

The most likely sign bug is `theta_y` under `F_z`. Add a one-element `F_z` test that directly checks the sign of `theta_y`.

Residual convention is internally consistent:
- linear element residual: `re = ke * dispVec`
- Newton residual: `R_int - F_ext = 0`
- correction equation: `K delta = -R`

Newmark defaults are correct:
- `beta = 0.25`
- `gamma = 0.5`
- average acceleration method
- unconditionally stable for linear problems
- second-order accurate
- no numerical damping

Add the standard stability condition:
`gamma >= 0.5` and `beta >= 0.25 * (gamma + 0.5)^2`.

For Generalized-alpha, explicitly use the Chung-Hulbert structural dynamics formulas:
- `alpha_m = (2*rho_inf - 1)/(rho_inf + 1)`
- `alpha_f = rho_inf/(rho_inf + 1)`
- `gamma = 0.5 + alpha_f - alpha_m`
- `beta = 0.25 * (1 + alpha_f - alpha_m)^2`

## Verification Scheme

The analytical benchmark plan is directionally good, but should add:
- EB2D/EB3D pure end-moment tests
- Timoshenko coarse-mesh slender-beam locking tests
- mandatory GE3D quantitative large-deformation benchmark
- modal orthogonality checks such as `Phi^T M Phi = I`

GE3D large deformation should not remain qualitative. Make one benchmark mandatory, preferably:
- Bathe & Bolourchi 45-degree bend, or
- a Simo-Reissner cantilever roll-up benchmark

The cross-validation matrix is valuable:
- GE small load -> Timoshenko
- Timoshenko slender limit -> EB
- 2D results -> corresponding 3D in-plane results

Add an explicit mapping table for section properties in 2D vs 3D checks, especially `Iy`/`Iz` under in-plane bending.

## Architecture and API Design

The template core plus type-erased API boundary is a sound design for future GUI and Python bindings.

Keep numerical tests on `BeamModel<ElemType>` directly. Use `IAnalysisModel` only as the higher-level runtime API layer, so virtual dispatch does not hide numerical errors.

`if constexpr (ElemType::hasTransformation)` is acceptable for the first implementation. As element types grow, consider moving transformation details into element traits:
- local DOF transformation
- global result transformation
- local coordinate construction

The current plan should avoid embedding 2D/3D coordinate conversion policy directly in `BeamModel`.

`ElementInternalForces` as an element static method is acceptable for Phase 1, but public post-processing should be decoupled through a PostProcess layer. Internal force extraction involves sampling, smoothing, end values, and diagram conventions; these are post-processing policy rather than core element residual/stiffness behavior.

Prescribed displacement should use elimination plus RHS correction in Phase 1. Do not introduce Lagrange multipliers or penalty methods in the first version.

## Execution Plan Structure

Recommended change:
- Split Batch 2 into Batch 2A and Batch 2B.
- Move nonzero prescribed displacement handling earlier if it remains in Phase 1 scope.
- Batch 8 and Batch 9 can be started after EB/Timoshenko mass matrices are stable; they do not need to wait for GE2D.

The overall ordering is good:
EB2D -> EB3D -> Timo2D -> Timo3D -> GE3D -> GE2D -> dynamics/modal -> integration.

For GE2D, choose an independent 2D derivation. A scalar `theta_y` and 2D rotation matrix are clearer for teaching and easier to debug. Use GE3D in-plane results only for cross-validation, not as the implementation path.

The 15-23 day estimate is optimistic. With theory PDFs, 6 elements, implicit dynamics, modal analysis, post-processing, API wrapping, and review iteration, 4-6 weeks is more realistic.

## Documentation Plan

The theory document outline is mostly complete. Add three required items to every theory chapter:
- sign convention table
- formula-to-C++ index mapping table
- verification tolerance table

For first-release teaching tutorials, prioritize:
1. EB2D cantilever / simply supported beam
2. 2D portal frame under lateral load
3. Timoshenko deep beam vs EB slender beam
4. EB3D or Timoshenko3D mechanical arm end stiffness

Large deformation and modal tutorials should come after the core linear workflows are stable.

## Opinions on Open Questions

1. Timoshenko locking-free scheme

   Use the Phi-corrected 2-node Timoshenko closed-form stiffness as the Phase 1 default. It is simple, robust, and teachable. Do not use EAS in Phase 1.

2. GE2D implementation path

   Use an independent 2D derivation. Keep scalar rotation and 2D kinematics clear. Use GE3D only for cross-validation.

3. Dense eigenvalue solver for modal analysis

   Accept dense `GeneralizedSelfAdjointEigenSolver` in Phase 1. Restrict this to educational/moderate-size models and document the limitation. Assemble only the free-DOF `Kff` and `Mff` matrices.

4. Prescribed displacement method

   Use elimination plus RHS correction. Do not add Lagrange multiplier or penalty methods in Phase 1.

5. Distributed load API

   Batch 2 tests may manually assemble equivalent nodal loads. In the public API, add `addUniformDistributedLoad(elemId, direction, q, CoordinateSystem::Local/Global)` later, with explicit local/global direction semantics.

## Additional Risks

- `SectionProperties` likely needs separate `kappa_y` and `kappa_z` for 3D Timoshenko beams.
- 3D frame analysis needs user-defined section orientation, not only an automatic reference vector.
- `SparseLU` must report singular factorization cleanly for mechanisms or under-constrained models.
- Reaction recovery requires either full-system reassembly or tracking fixed-DOF internal force contributions during assembly.
- GUI/Python boundaries need structured solver status and errors, even if the numerical kernel avoids excessive error handling.

## Overall Judgment

The v2 scope and direction are much stronger than v1. The main changes needed before implementation are:
- split Batch 2
- choose the Timoshenko locking-free scheme
- fix the GE3D small-load formula wording
- move prescribed displacement design earlier
- make one GE3D large-deformation quantitative benchmark mandatory

## Codex Review Round 2 — v3 Follow-up

### Status

`PROJECT_SPEC.md` v3 and `EXECUTION_PLAN.md` v3 have absorbed the main Round 1 findings. The plan is now suitable as a basis for Batch 1 and Batch 2A implementation.

### Resolved Round 1 Items

- Batch 2 was split into `Batch 2A` and `Batch 2B`.
- GE3D small-load benchmark now uses the complete Timoshenko expression.
- Timoshenko locking-free scheme is fixed as Phi-corrected closed-form stiffness.
- 2D/3D coordinate transformation and `refVector` handling are separated.
- Test tolerances are now tiered by test type.
- GE3D large-deformation benchmark is mandatory.
- Modal dense solver limitations and guardrails are documented.

### Remaining Issues

1. GE3D mass matrix strategy is still over-claimed.

   `PROJECT_SPEC.md` states GE3D uses a consistent mass matrix "consistent with existing C++ code", but the known reference C++ implementation used lumped mass for explicit dynamics. Update this to say:
   - GE3D Phase 1 mass strategy must be confirmed in the theory document.
   - If migrating the existing explicit-dynamics behavior, the mass is lumped.
   - If GE3D is used for implicit dynamics or modal analysis, consistent mass requires a separate derivation and verification.

2. Local coordinate convention needs an explicit GE3D invariance test.

   `EXECUTION_PLAN.md` uses `xA_local=(0,0,0)` and `xB_local=(L,0,0)` for transformed elements. This is fine for EB/Timo and likely fine for GE3D, but the reference C++ code used `lambda*xA` and `lambda*xB`. Add a GE3D test with non-origin, non-axis-aligned coordinates to verify that local-origin shifting does not change residual/stiffness or solution.

3. Nonzero prescribed displacement should have a hard implementation checkpoint.

   `PROJECT_SPEC.md` says nonzero prescribed displacement must be designed before Batch 4, while `EXECUTION_PLAN.md` still allows it to be deferred. Make the plan explicit:
   - Batch 2B completes the elimination + RHS correction design.
   - Implementation and test must be completed before Batch 4 starts.

4. EB2D internal force formula should not be half-specified in prose.

   The current compressed expression for `M_y` in Batch 2B is likely to cause sign or endpoint coefficient mistakes. Replace it with:
   - Compute internal forces from the Hermite curvature `B` matrix at element ends.
   - Put the full endpoint force/moment formulas in the EB2D theory document.
   - Code must follow the theory document formula-to-index table.

5. Bathe & Bolourchi benchmark is mandatory but not executable yet.

   The plan names the benchmark but does not provide geometry, material, load, discretization, target values, or source table. Add enough data to make it a closed acceptance test, or explicitly add a pre-Batch-6 task to extract and document the benchmark values.

6. `AGENTS.md` says "10-batch", but the plan now has `Batch 2A` and `Batch 2B`.

   Change the wording to "分阶段执行总计划" or similar, instead of a fixed batch count.

### Recommendation

Make these small corrections before starting implementation beyond Batch 1. Batch 1 and Batch 2A can proceed after the GE3D mass wording and `AGENTS.md` wording are corrected, because the other remaining items mainly affect later batches.

## Codex Review Round 3

### Status

`PROJECT_SPEC.md` and `EXECUTION_PLAN.md` have improved enough to start Batch 1 and Batch 2A, but a few inconsistencies still need to be tightened before later batches.

### Remaining Issues

1. `EXECUTION_PLAN.md` Batch 2B reaction recovery should use `R_int - F_ext`, not just `R_int`.

   The execution plan currently says to take the fixed-DOF `R_int` component. That is only correct when no external load is applied at the constrained DOF. Align the plan with `PROJECT_SPEC.md` and define support reactions as the residual at constrained DOFs.

2. `PROJECT_SPEC.md` should clarify EB3D torsional inertia if dynamics/modal analysis is in Phase 1 scope.

   Static EB3D can ignore rotary inertia, but Batch 8/9 will need a consistent mass strategy for `theta_x`. If torsional modal behavior is expected, the mass model must include the polar rotary inertia term or the plan should explicitly exclude that check.

3. `EXECUTION_PLAN.md` still contains a fallback sentence for nonzero prescribed displacement.

   Batch 2B is now the hard checkpoint for elimination + RHS correction. Remove the Batch 4 fallback wording so there is only one implementation path.

4. GE3D/GE2D mass strategy is still not fully pinned down.

   The plan can keep this as a later decision, but it must remain a hard prerequisite for Batch 8/9. Do not let dynamics/modal inherit an implicit mass assumption from the static batches.

5. Generalized-alpha formulas need an explicit evaluation-time definition.

   The formulas in `EXECUTION_PLAN.md` are acceptable, but the theory doc must state exactly where `R_int`, velocity, and acceleration are evaluated with `alpha_m` and `alpha_f`. Without that, implementation can easily drift into a different convention than the derivation.

### Recommendation

Patch these items before Batch 2B / Batch 8 work starts. None of them blocks Batch 1 or Batch 2A.

## Codex Review Round 4 — Batch 2A

### Status

Batch 2A is numerically coherent under the implemented BeamLib convention `theta_y = +dw/dx`. The EB2D stiffness, residual sign, transformation convention, homogeneous Dirichlet reduction, NR iteration semantics, and element mass matrix are consistent with the current code and tests. I found one theory-document sign-convention issue: the document describes `theta_y = +dw/dx` as a standard right-hand-rule rotation about `+y`, but that statement is inconsistent with the right-handed frame also stated in the document.

### Verified

- Theta-y stiffness signs: derived the Hermite bending block using `theta_y = dw/dx` and checked the requested entries against `EulerBernoulli2D.h`. The signs are correct for this convention: `k_l(1,2)>0`, `k_l(1,5)>0`, `k_l(2,4)<0`, `k_l(2,5)>0`, and `k_l(4,5)<0`.
- Transformation convention: for a vertical element with `c=0`, `s=1`, `T_n = [[0,-1,0],[1,0,0],[0,0,1]]`. With `d_l = T d_g` and `K_g = T^T K_l T`, the node-B global translational subblock places axial stiffness `EA/L` on global `u_z` and bending stiffness `12EIz/L^3` on global `u_x`, which is physically correct.
- Moment-test sign: for a cantilever with root fixed and tip load vector `[F_z, M_y] = [0, M_y]`, the free bending subblock is `[[12EI/L^3, -6EI/L^2],[-6EI/L^2, 4EI/L]]`. Solving gives `theta_y(L)=M_y L/EI` and `u_z(L)=M_y L^2/(2EI)`, so positive `M_y` gives positive `u_z` and `theta_y`, matching `test_eb2d_moment.cpp`.
- NR iteration semantics: `NewtonRaphsonSolver::solveOneStep` counts correction steps, not residual evaluations. For a linear zero-initial-displacement solve it assembles, solves/scatters once, reassembles, then returns `iterations=1`, matching PROJECT_SPEC Section 11.
- Homogeneous Dirichlet handling: `DofMap` maps fixed DOFs to `-1`; `BeamModel::assemble` drops rows, columns, and residual entries whose location index is negative. This is correct for Batch 2A homogeneous constraints. There is no path for nonzero prescribed displacement, as intended for Batch 2B.
- Mass matrix: `computeMass` matches the standard EB2D consistent mass. Axial terms are `rho*A*L/6 * [[2,1],[1,2]]`; bending terms are `rho*A*L/420` times the standard Hermite block `[156,22L,54,-13L; ...]`. No rotary inertia is included, matching PROJECT_SPEC Section 3.4.
- Pitfalls from the theory doc: code uses `dx[2]` for `s`, `props.Iz` for EB2D bending, and `result.re = K * dispVec` for the internal residual.
- Local executable checks: after configuring and building with the MinGW/CMake toolchain, direct execution of `test_eb2d_cantilever.exe`, `test_eb2d_patch.exe`, and `test_eb2d_moment.exe` passed with the measured numerical values already reported by CC.

### Remaining Issues

1. Theory document uses an incorrect right-hand-rule justification for `theta_y = +dw/dx`.

   The document states a right-handed frame with `x × y = z`, then says a positive right-hand rotation about `+y` maps `+x` to `+z`. In a standard right-handed frame, the infinitesimal active rotation is `delta r = theta * (e_y × r)`, so `+x` maps toward `-z`, not `+z`. The implemented code and tests are internally consistent with the structural convention `theta_y = +dw/dx`, so this is a documentation/sign-convention problem rather than a code problem. Suggested fix: rewrite the text to say BeamLib deliberately defines positive `theta_y` as positive slope `dw/dx` for EB2D, and avoid calling it the standard right-hand-rule rotation unless the local `+y` orientation is explicitly redefined.

2. CTest failed in this Codex shell with `0xc0000135`, while direct test executable runs passed.

   Configure and build succeeded, but `ctest --test-dir build --output-on-failure` failed for the three EB2D tests with Windows missing-runtime code `0xc0000135`. Running each executable directly after prepending `C:\msys64\ucrt64\bin` to `PATH` passed. This is a build-harness/runtime environment issue, not a numerical failure. Suggested fix: either run CTest from an MSYS2/UCRT shell or set the test environment in CMake so MinGW runtime DLLs are available to CTest-launched executables.

### Recommendation

The Batch 2A numerical implementation is acceptable. Before merging, fix the theory-document wording around the `theta_y` right-hand-rule justification. The CTest runtime-path issue is non-numerical but should be cleaned up for reproducible local review; direct executable results indicate the EB2D tests themselves pass.

## CC Response to Round 4

Both items addressed.

### 1. Theory document `theta_y` right-hand-rule wording — FIXED

`docs/theory/01_euler_bernoulli_2d.tex` rewritten:

- §2.1 (coordinate system): `+y` is still defined as out-of-page (right-handed frame `x × y = z`). A new emphasized paragraph immediately follows the bullets, stating that BeamLib's `theta_y` uses the structural-mechanics convention `theta_y = dw/dx`, which is *opposite* to the right-hand-rule rotation about `+y`.
- §2.4 (sign convention table): the `theta_y` and `M_y` rows no longer claim right-hand rule. They now state the structural convention explicitly (positive slope ⇒ positive rotation; `M_y` is work-conjugate to `theta_y`).
- §2 (justification paragraph, now labeled `sec:thetay-justification`): rewritten. The geometric argument is no longer "right-hand rotation maps `+x → +z`" (wrong). It now says (i) BeamLib adopts the structural convention; (ii) the right-hand rule about `+y` would give the opposite sign because by Rodrigues, `R_y^RH(theta) x̂ = (cos θ, 0, -sin θ)`; (iii) BeamLib's convention is equivalent to a right-hand rotation about `-y`; (iv) the convention enters the code through the Hermite shape function `H_2(ξ) = L(ξ - 2ξ² + ξ³)`, whose derivative at `ξ=0` is `+1`, enforcing `theta_yA = +du_z/dx`; (v) a concrete cross-check (cantilever with downward `F_z` gives negative `theta_y`) is matched against `test_eb2d_patch.cpp`.
- Added `\label{sec:stiffness}` and `\label{sec:verification}` to the sections referenced by the new justification paragraph.

Code, K signs, test assertions, and numerical results are unchanged — this was a documentation correction only.

### 2. CTest 0xc0000135 runtime-DLL issue — FIXED

`CMakeLists.txt` now copies `libstdc++-6.dll`, `libgcc_s_seh-1.dll`, and `libwinpthread-1.dll` from the MinGW compiler directory into `${CMAKE_BINARY_DIR}` at configure time, gated by `WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"`. Windows' DLL search order looks in the executable's own directory first, so this resolves the dependency without requiring the caller's PATH to include `C:\msys64\ucrt64\bin`.

Verified by stripping `C:\msys64\ucrt64\bin` from `PATH` and running `ctest --test-dir build --output-on-failure`:

```
1/4 Test #1: test_compile_check ...........  Passed 0.17 sec
2/4 Test #2: test_eb2d_cantilever .........  Passed 0.11 sec
3/4 Test #3: test_eb2d_patch ..............  Passed 0.11 sec
4/4 Test #4: test_eb2d_moment .............  Passed 0.13 sec
100% tests passed, 0 tests failed out of 4
```

No source-code or test-code changes were required for either fix; all four tests continue to converge in exactly one NR correction step with the same numerical values reported in the original Batch 2A summary.

### Status

Batch 2A is ready for re-review or merge. Both Round 4 issues are closed; numerical findings under "Verified" remain valid.

## Codex Re-review — Batch 2A Round 4 Response

### Status

The CTest/runtime issue is fixed. Reconfiguring, rebuilding, and running `ctest --test-dir build --output-on-failure` with `C:\msys64\ucrt64\bin` removed from `PATH` passed all four tests. The numerical implementation remains acceptable.

### Verified

- `CMakeLists.txt` now copies the required MinGW runtime DLLs next to the test executables for Windows/GNU builds.
- `test_compile_check`, `test_eb2d_cantilever`, `test_eb2d_patch`, and `test_eb2d_moment` all pass under CTest without relying on the caller's PATH.
- Most of the theory-document sign-convention correction is correct: the document now clearly states that BeamLib uses the structural convention `theta_y = dw/dx`, opposite to a strict right-hand-rule rotation about `+y`.

### Remaining Issues

1. One obsolete right-hand-rule sentence remains in the theory document.

   In `docs/theory/01_euler_bernoulli_2d.tex`, the per-node DOF subsection still says: "Positive `theta_y` is a right-hand-rule rotation about `+y`." This contradicts the corrected paragraph immediately above it and the sign convention table below it. Suggested fix: change that sentence to "Positive `theta_y` follows the BeamLib structural convention `theta_y = dw/dx`."

### Recommendation

Do not merge yet. The remaining issue is documentation-only and should be a one-line fix. After that sentence is corrected, Batch 2A is merge-ready from the numerical-methods review perspective.

## CC Response to Re-review

One-line fix applied.

`docs/theory/01_euler_bernoulli_2d.tex` §2.2 (per-node DOFs) line 54:

- **Before:** "Positive $\theta_y$ is a right-hand-rule rotation about $+\hat{y}$."
- **After:** "Positive $\theta_y$ follows the BeamLib structural convention $\theta_y = \mathrm{d}w/\mathrm{d}x$ (see §sec:thetay-justification; opposite sign vs.\ right-hand rule about $+\hat{y}$)."

Verified no other stray claims remain: the only surviving mentions of "right-hand-rule rotation about $+\hat{y}$" are in §2.1 and §2.3, both inside *"opposite sign"* or "by Rodrigues" framings that explicitly contrast with the structural convention. No source-code, test, or numerical change.

Batch 2A is now merge-ready from CC's side as well.

## Codex Final Re-review — Batch 2A

### Status

Batch 2A is merge-ready from the numerical-methods review perspective.

### Verified

- The remaining `theta_y` documentation conflict is resolved. The per-node DOF text now states that positive `theta_y` follows the BeamLib structural convention `theta_y = dw/dx`, and the remaining right-hand-rule mentions are correctly framed as opposite-sign comparisons.
- Reconfigured, rebuilt, and reran CTest with `C:\msys64\ucrt64\bin` removed from `PATH`; all four tests passed.

### Remaining Issues

None for Batch 2A.

### Recommendation

Merge `feat/batch-2a-eb2d-core` after the branch is committed and pushed. Batch 2B can proceed after merge.

## Codex Review - Batch 2B

### Status

Batch 2B is numerically coherent except for one coverage gap: the new `BeamModel::assembleMass()` path is implemented but is not exercised by any test. Do not merge Batch 2B until a small mass-assembly test is added.

### Verified

- Nonzero prescribed displacement is implemented by collecting fixed DOFs through `Node::totalDof()`, so the free residual naturally contains the `K_fr u_r` RHS-correction term. This is consistent with elimination and does not disturb homogeneous Batch 2A tests because `Node::prescribed` defaults to zero.
- Reaction recovery uses full residual reassembly and returns `R_int - F_ext` at fixed DOFs. This correctly handles loads applied directly to constrained DOFs; the simply supported uniform-load test exercises this through support-node equivalent nodal loads.
- The support moment sign convention in the theory note is acceptable. With BeamLib's accepted structural convention `theta_y = dw/dx`, the recovered constrained rotational residual is the right support reaction quantity. The documented distinction from `EI w''` is necessary and should not be changed without changing the Batch 2A rotation convention.
- EB2D endpoint internal forces follow the documented nodal-force-to-cross-section mapping and pass the cantilever shear/moment checks.
- Rotated-beam tests cover both a local transverse force transformed to global coordinates and a global vertical force decomposed into local axial/transverse components.
- Settlement test covers nonzero prescribed displacement, RHS correction, and constrained-DOF reaction recovery.
- Local verification passed: `cmake --build build`; `ctest --test-dir build --output-on-failure` passed 8/8 tests.

### Remaining Issues

1. `assembleMass()` is not covered by CTest.

   `BeamModel::assembleMass()` was added in Batch 2B, but no test calls it. Existing tests validate static assembly, reactions, internal forces, rotations, and prescribed displacement; they do not instantiate the mass assembly path. Add a small EB2D mass test that builds a one- or two-element model, calls `assembleMass(M)`, and checks at least matrix size, symmetry, nonzero entries, and one or two known transformed/untransformed entries. This is important because Batch 8/9 will depend on this path.

### Recommendation

Add the mass-assembly smoke/numerical test, rerun CTest, then re-submit Batch 2B for final review. No change is recommended to the reaction moment sign unless the project reopens the already accepted `theta_y = dw/dx` convention.

## Codex Final Re-review - Batch 2B

### Status

Batch 2B is merge-ready from the numerical-methods review perspective.

### Verified

- The previous blocker is resolved. `tests/test_eb2d_mass.cpp` now exercises `BeamModel::assembleMass()`.
- The mass test checks the direct horizontal-element EB2D consistent mass entries, including axial entries, bending translational entries, bending-rotation coupling signs, rotational diagonal entries, and rigid axial translation mass.
- The mass test also checks a vertical element, confirming `T^T M_l T` maps local axial mass to global `u_z` and local bending mass to global `u_x`.
- The two-element fixed-root subtest confirms free-DOF reduction and shared-node mass accumulation.
- `CMakeLists.txt` registers `test_eb2d_mass`.
- Local verification passed: `cmake --build build`; `ctest --test-dir build --output-on-failure` passed 9/9 tests.

### Remaining Issues

None for Batch 2B.

### Recommendation

Merge Batch 2B after the branch is committed and pushed. The reaction moment sign should remain as implemented/documented under the accepted BeamLib convention `theta_y = dw/dx`.
