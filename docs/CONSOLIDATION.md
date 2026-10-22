# Consolidated development version 0.3.0.9002

The `main` branch is the current package. The four `archive/*` branches retain
earlier implementations; they are not dependencies of this version. Numerical
engines, adaptive geometry and compiled EM are included in the package source.

This inventory was checked against the successive research implementations and
the development discussion on 22 October 2026. The maintained user references
are the function help pages and the three installed vignettes.

| Requirement | Current implementation | Verification |
|:--|:--|:--|
| Gaussian density ceiling | A determinant constraint, without an imposed eigenvalue floor; explicit singular-support profiles | `test-determinant-only`, `test-singular-affine` |
| Hard CEM and its monotone objective | Reassignment and component-deletion proposals accepted only when the refitted objective improves | `test-optimization-methods`, `test-reassignment-floor` |
| Choice of optimizer independent of penalty | `optimization="cem"` or `"em_bank"`, crossed with `penalty="none"` or `"bic"` | All four combinations tested in classification and grid entry points |
| Constrained soft EM | EM enforces the density ceiling `Cmax`; valid likelihood traces are monotone; stop reasons distinguish convergence and budgets or numerical limits | Native EM checks and `test-optimization-methods` |
| EM followed by MAP and hard refitting | Empty MAP groups are removed; hard masses and Gaussian profiles are refitted under each target C before scoring | Direct likelihood and hard-score checks |
| Common bank for the entire grid | One bank per independent execution; all orders and starts contribute candidates to every cell | Fixed and adaptive bank tests |
| Large and spaced order searches | `orders`, or every integer through `Rmax`; effective orders and failed starts are retained | `test-adaptive-grid` |
| Reuse with another penalty or grid | `CECrescoreGrid` reuses the fitted banks, reranks each execution and recomputes diagnostics | Paired fixed/adaptive tests; archived banknote regression |
| Independent algorithmic and bootstrap outputs | Each execution selects on its own training sample; bootstrap outputs use frozen learned parameters on the original rows | Winner-by-winner and diagnostic checks |
| Exhaustive stability | L(L-1)/2 original-original pairs and L B original-bootstrap pairs; standalone `CECstability` also accepts another method's labels | Independent agreement checks and `test-optimization-methods` |
| Targeted diagnostic measurement | `all`, `unsaturated`, `nontrivial_unsaturated`, `saturation_band`; skipped diagnostics remain missing and cannot pass the filter | Policy and missing-value tests |
| Three separate saturation controls | Final `sat_threshold`, frontier `adaptive_control$frontier_saturation`, and measurement `rsi_saturation_max` | Documented in `CECdiagnostics` and `CECadaptiveDiagram` |
| Initial window and four extensions | `grid_strategy="adaptive"`; supplied initial centres or a window scaled from the data | Exact geometry import and regression against the research runner |
| Row-wise right stopping | Each row is trivial/saturated, or remains unsaturated at unchanged REO over two extensions; changed comparison intervals reset its history | Dedicated mixed-row test |
| Refinement on both sides of a frontier | Partition or saturation changes, including different partitions with equal REO; local patience/displacement controls | Geometry, nonoverlap and same-order frontier tests |
| CEM neighbour harmonization | First pass on the whole mesh, then new cells and their neighbours with propagation until stable; separate `A_harm` control; independent diagnostic outputs remain frozen | Monotonicity, local propagation, final full-mesh audit and independent-winner tests |
| No extra EM-bank harmonization | All candidates are already available throughout the diagram | Bank rescoring and candidate-minimum checks |
| Optional trivial quadrants | Trivial anchors, stability checks and a persistent shared audit budget across overlapping quadrants; exempt cells remain unmeasured and excluded | Quadrant coverage and audit-budget tests |
| Parallel execution and interrupted fitting | Independent repetitions; per-repetition and per-EM-start checkpoints, with seed and input validation | Serial/parallel and partial-resume regressions |
| Correct irregular-cell displays | Actual adaptive rectangles, extraction inside a cell and ordered lambda paths | Plot/extraction tests and visual inspection |
| Compact saved results | Banks, independent winners, labels, data and adaptive geometry remain available for inspection and rescoring | Compact-object and data-fingerprint tests |
| Pedagogical documentation | Executed README, 26 help topics, three executable vignettes and an installed quick-start script | Fresh source build, help examples, vignette checks and README assertions |

## Successive designs

The early per-cell EM/MAP experiment was superseded by the shared-bank design.
The current CEM route still searches hard allocations directly. Earlier
economical and best-reference stability experiments were superseded by the
exhaustive definition for the consolidated engines. Historical diagnostic
conventions remain explicitly named on the compatibility route; their values
must not be relabelled as exhaustive RSI.

Fixed grids remain the default for compatibility. Request the adaptive strategy
explicitly. Numeric Gaussian implementations are provided for the new engines;
the historical unpenalized CEM route retains categorical and mixed inputs.

Both optimizers are approximate. EM increases a soft likelihood to generate
candidates and subsequently ranks hard partitions; it does not optimize the
hard CEC criterion during its EM iterations. The bank minimum is a minimum over
the available candidates. The diamond audit demonstrated why that distinction
matters; no new hybrid optimizer is claimed here.

The package filter uses the documented strict stability thresholds and an
inclusive saturation threshold. The inclusive .85/.90 filter used in the final
all-diamond report was a selection setting of that study; it did not change
the package defaults or the numerical fitting protocol.

## Validation before publication

The source, Windows archive, numerical tests, historical regressions and
documentation are checked together. `R CMD check --no-manual` includes help
examples and the three vignettes. The README is rendered by executing its
examples and verifying the stated partitions, diagnostics and BIC comparison.
See [CONTRIBUTING.md](../CONTRIBUTING.md) for the maintained build procedure.
