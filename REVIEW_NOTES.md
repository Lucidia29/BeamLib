# BeamLib Review Notes

These notes capture the review changes that must be reflected in `CODEX_TASKS.md` before implementation continues.

## Required Fixes

1. GE3D cantilever acceptance must not use a pure Euler-Bernoulli small-deflection baseline.
   - The current GE3D formulation includes shear strain, so the benchmark in Batch 2 should be aligned with the existing golden reference (`sanity_cantilever_ref.csv`) or explicitly stated as a Timoshenko-type small-deflection baseline.
   - Do not keep the current EB-only values as the verification target.

2. Linear EB element residual sign must match the solver convention.
   - `BeamSolver` uses `R_int - F_ext`.
   - Therefore the linear EB elements in Batch 4 and Batch 5 should use `re = K * dispVec`, not `re = -K * dispVec`.

3. 2D and 3D transformation handling must be separated clearly.
   - The current `BeamModel` transformation sketch is safe for 3D GE3D, but not a general 2D solution.
   - Batch 4 needs a dedicated 2D x-z transformation path, and a non-axial / rotated beam test.
   - Keep 3D transformation in `GeomExact3D::computeTransformation`.

4. Batch 3 large-deformation verification should be framed as a reference-based or stress-style test.
   - The current `F_z = 5e5` / `~0.5L` setup is too assumption-heavy for first-pass acceptance.
   - Use the MATLAB / reference C++ benchmark, or downgrade this case to a stress test with qualitative checks.

5. The explicit integrator should be scoped to GE3D first.
   - The lumped-mass rule in Batch 3 is 6DOF GE3D-specific.
   - Keep the first version as `ExplicitCentralDifference<GeomExact3D>` and extend later if EB elements need mass support.

6. Batch 1 should not rely on an empty static library source list unless the CMake behavior is explicitly known.
   - Prefer an `INTERFACE` target for the first batch, or a clearly documented placeholder `.cpp`.

7. Batch 2 should include element-level tests, not only an end-to-end cantilever.
   - Add zero-displacement, pure axial, and pure torsion checks to isolate element, assembly, and solver issues.

8. The meaning of “1 NR iteration” must be explicit.
   - For linear problems, define this as one correction step after the initial residual evaluation.
   - Keep the wording consistent between Batch 4 and Batch 5.

9. `LoadManager` needs a complete prescribed-displacement contract.
   - Define whether `addPrescribedDisplacement` also writes the displacement value into the node state.
   - Define whether it must be called before `buildDofMap()` or whether rebuilding is required afterward.

10. Git workflow should include a baseline commit.
    - Create the baseline commit first.
    - Then work on `feat/batch-*` branches.
    - Keep review notes and plan updates in git so Claude Code can read them directly.

## Recommendation

Update `CODEX_TASKS.md` to incorporate these fixes before starting the next Codex batch.
