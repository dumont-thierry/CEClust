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
