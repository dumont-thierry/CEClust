# CEClust 0.3.0.9004

- Adaptive diagrams now default to connected-component traversal after freezing
  the pilot. Candidates are ordered by distance to the area-weighted barycentre
  in (log C, lambda), restricted by the measurement policy, and stop at the
  first point passing both diagnostic thresholds. Witness values may describe
  a region; direct measurements, raw diagnostics and independent winners remain
  identifiable. Use adaptive_control=list(diagnostic_strategy="cells") for
  exhaustive cell-by-cell coverage.
- Every measured point in the numeric grid engines checks the reference rule
  before filtering. When the best original execution disagrees with the grid
  representative below the algorithmic threshold and does not improve its
  score, fix the representative and use exactly one projected-label CEM descent
  per bootstrap. Algorithmic stability is then one by construction; RSI is
  the mean agreement of the frozen-reference bootstrap projections. Ordinary
  exhaustive diagnostics and the two trigger conditions are retained.
- Both optimizers and penalties, fixed grids, adaptive grids and EM rescoring
  use this rule. reference_bootstrap=FALSE restores independent-only diagnostics;
  adaptive_control$reference_bootstrap can set it for adaptive diagrams.
- Added per-reference bootstrap checkpoints, provenance in selection/plots,
  component geometry, ordering, trigger, single-descent and restart tests.

# CEClust 0.3.0.9003

- Replaced adaptive CEM orchestration with a separate original-data pilot,
  followed by selective independent original/bootstrap campaigns. Versions
  0.3.0.9001/9002 fitted those campaigns during pilot construction.
- Restored parent starts for pilot children, local propagation after campaign
  improvements, border/coverage rechecks and explicit pending-frontier status.
  Independent campaign outputs and their diagnostics remain frozen.
- Adaptive CEM now defaults to the saturation band (0.10); explicit policies,
  fixed-grid defaults and EM-bank defaults retain their documented meaning.
- Added per-pilot, per-campaign-execution and harmonization checkpoints, a
  phase/seed execution register, and compact selective measurement records.
  Old adaptive CEM caches are not treated as the corrected protocol.
- Fixed floating-point rejection of valid mesh ceilings and allowed interior
  EM refinements when other proposed children would exceed Cmax. Ceiling and
  budget restrictions remain explicit search limits.
- Added behavior tests for execution staging, whole-campaign gating, independent
  seeds/projection, parent starts, frozen campaigns, quadrant exemptions and
  numerical resume. Updated reference help and migration guidance.

# CEClust 0.3.0.9002

- Consolidated package published on `main`; historical `archive/*` branches
  remain available. The README explains installation from GitHub and the
  consolidation inventory maps the successive additions to their checks.
- The initial CEM harmonization visits the whole mesh; subsequent passes begin
  at new cells and propagate through affected neighbours until stable.
- Refreshed installation, compatibility and partial-checkpoint documentation.

# CEClust 0.3.0.9001

- Integrated adaptive windows: four directions, row-wise right stop after two
  unsaturated extensions with stable order, and two-sided boundary refinement.
- Both CEM and EM candidate banks, with or without BIC. Adaptive CEM retains
  neighbour harmonization and a separate agreement control. EM needs none.
- Exhaustive diagnostics, optional targeted measurement and controlled trivial
  quadrant exemptions; unmeasured values remain missing and are not retained.
- Explicit EM order lists, per-start EM and per-repetition fitting checkpoints.
- Adaptive geometry preserved in summaries, plots, extraction, compact objects
  and the explorer. Rescoring EM banks rebuilds the adaptive diagram.
- Limits and stopping evidence are retained in the result. Updated help and
  executable README cover the complete workflow.

# CEClust 0.3.0.9000

* Rewrite the complete user documentation: 19 dedicated function pages,
  shared optimizer/diagnostic/result references, executable guides, migration
  notes, and an installed quick-start script. Clarify actual defaults,
  measurement policies, complete-bank checkpoints and compatibility behavior.

* Integrate determinant-only Gaussian profiling and guarded reassignment from
  local core 0.2.0.9004, retaining categorical/mixed-data CEM.
* Separate `optimization = "cem" | "em_bank"` from `penalty = "none" | "bic"`.
  EM candidates are ranked after hard-parameter refitting under the target C.
* Compile constrained soft EM through the package Rcpp toolchain, without a
  precompiled external DLL or private Rcpp calls.
* Integrate numeric penalized CEM and deletion proposals from v075, and independent
  EM banks from v078. Grid, bootstrap, plots, presets, checkpoints and result
  metadata honor explicit optimizer choices.
* Use exhaustive original-original and original-bootstrap stability by default.
  Legacy best-reference and bootstrap-bootstrap definitions remain explicitly named.
  A single original repetition no longer reports perfect algorithmic stability.
* Add configurable RSI measurement policies and `CECrescoreGrid()` for paired
  with/without-BIC evaluation on the same banks. Skipped diagnostics stay missing.
* Preserve the historical fixed-grid interface when no new optimizer options
  are supplied. Adaptive experimental mesh runners remain separate; this
  release does not claim to consolidate their control flow.

# CEClust 0.2.0

- Promoted the current development code to a release-ready package version with
  refreshed metadata, NEWS, and package-level documentation.
- Consolidated the public workflow around fixed-lambda fitting, linked-lambda
  diagnostics, best-partition extraction, trajectory repair, and visualization.
- Renamed the main lambda-grid workflow to `CECfitLambdaGrid()`, while keeping
  `CECdiagnose_lambda_grid_linked()` as a backward-compatible alias.
- Added a minimal high-level grid API: `CECfitBoundGrid()`, `CECfitPreset()`,
  `CECselectStableLambdas()`, `CECsummariseGrid()`, `CECplotGrid()`,
  `CECplotPartition()`, `CECplotPath()`, and `CECexplore()`.
- Reworked the user-facing documentation for exported functions, including more
  explicit return values, workflow-oriented examples, and clearer descriptions
  of runtime and stability diagnostics.
- Expanded the README into a practical tutorial and aligned the vignette with
  the recommended analysis workflow.
- Cleaned the source package layout by excluding build artefacts and adding
  lightweight regression tests for the main public entry points.

