# Consolidated development version 0.3.0.9007

This local source is the repaired development version; publication to `main`
is a separate step. The four `archive/*` branches retain
earlier implementations; they are not dependencies of this version. Numerical
engines, adaptive geometry and compiled EM are included in the package source.

## Independent stability gate (23 October 2026)

The conditional reference bootstrap requires all three conditions: raw independent
Stab_algo and raw independent RSI strictly above their respective thresholds;
agreement of the best original execution
with the grid representative strictly below that threshold; and no improvement
of the representative's score. Equality with the stability threshold does not
pass. If either independent diagnostic fails, preserve the raw diagnostics and
continue the component traversal. No reference bootstrap is launched for that
point. The check uses raw stability even when an older object contains a value
of one from a conditional calculation. It applies before filtering at every
measured point of the numeric CEM/EM, fixed/adaptive and rescoring routes.

## Compiled numeric CEM (23 October 2026)

`src/cem_fused.cpp` now shares Gaussian moments across the hard score, covariance
fitting and raw density calculation. Reassignment and frozen deletion use the
same accepted moves and stopping budgets as the retained `q_descent_R` path.
The default `optimizer_control$cem_engine="compiled"` can be switched to
`"reference"`. Ambiguous ranks and compiled errors fall back per start, with
their reasons recorded. The R statistics path is preallocated; exhaustive
agreements use exact memoization and direct contingency counting.

In version 0.3.0.9005, matched 20-worker campaigns with the same 40 independent searches and all 590
agreements, measured gains of x14.6 and x11.2 including saving. The complete
fresh pilot reproduces all 4796 cell partitions and 6373 searches in 94.69 s
(maximum score difference 4.44e-16). The barycentric diagnostic phase then uses
14 points in 11.10 s, including output preparation. The 105.79 s equivalent
total sums separately timed fresh pilot and diagnostic phases; replay reading
is excluded. Early stopping under the new reference rule is a separate gain
from the matched numerical-kernel improvement. Validation is archived in
`Validation/speed_20261023`; no scientific budget was reduced.
These full-procedure timings used the earlier two-condition reference rule.
The additional gates in 0.3.0.9006/9007 can require more candidate points per component;
it does not change the compiled fitting kernel.

## Component traversal and reference bootstrap (23 October 2026)

The adaptive default is now a frozen-pilot component traversal. Exact partitions
and edge connectivity define components independently of saturation. Test points
are ordered by distance to the full area-weighted barycentre in raw (log C,
lambda), restricted by the measurement policy, and stop at first diagnostic
success. Regional values identify their measured witness; direct measurements
and raw independent winners remain separate. The previous measurement-improvement
loop below is retained as `diagnostic_strategy="cells"`, not the default.

Before filtering every measured point in numeric fixed/adaptive grids, both
optimizers and penalties check the reference rule. When both raw independent
Stab_algo and RSI strictly exceed their thresholds, a sufficiently different
best original execution that does not improve the representative's score
triggers one projected-reference CEM descent per bootstrap. SA=1 is fixed by
construction; RSI is agreement with that fixed reference. Raw independent
diagnostics, comparisons and outputs remain archived. Repeated checks of the
same reference reuse completed conditional descents, including without a disk
checkpoint. `CECstability` alone remains the raw label-only arithmetic helper.

Verification: 1773 test expectations, eight optimizer/penalty/strategy parallel
and partial-resume comparisons, archived 4796-cell geometry regression, and
the real five-group example with 20 single-descents (RSI .9817). See
`test-component-diagnostics` and `test-reference-bootstrap`, plus the local
`Validation/components_20261023` delivery checks. Fixed-grid defaults remain
fixed; adaptive component traversal is selected only with the adaptive engine.

This inventory was reviewed against the research implementations and the
archived development discussion on 22 October 2026. The earlier inventory
missed an orchestration defect: versions 0.3.0.9001/9002 performed repeated
CEM fitting during adaptive pilot construction. Version 0.3.0.9003 replaces
that path with the staged driver described below. Score agreement alone was
insufficient to verify this requirement. The maintained user references
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
| Deferred CEM campaigns | One multistart original-data pilot per new cell; no L+B campaign before pilot calibration/refinement ends | Fitting-call phase/count tests in `test-staged-cem` |
| Parent starts in the CEM pilot | A new child receives the parent partition; independent campaigns always use fresh randomizations | Warm-start isolation and archived-parent tests |
| Measurement improvements and coverage | Inject better campaign partitions, harmonize locally, recheck borders and eligibility; keep completed campaigns frozen | Controlled improvement/coverage test, no duplicate campaign executions |
| Post-measurement refinement | Reactivated frontiers or exhausted budgets prevent completed adaptation and have explicit final status | Pending-frontier test |
| Targeted diagnostic measurement | `all`, `unsaturated`, `nontrivial_unsaturated`, `saturation_band`; adaptive CEM gates whole campaigns (default saturation band .10), while EM/fixed routes gate comparisons; skipped diagnostics remain missing and excluded | Policy and missing-value tests |
| Three separate saturation controls | Final `sat_threshold`, frontier `adaptive_control$frontier_saturation`, and measurement `rsi_saturation_max` | Documented in `CECdiagnostics` and `CECadaptiveDiagram` |
| Initial window and four extensions | `grid_strategy="adaptive"`; supplied initial centres or a window scaled from the data | Exact geometry import and regression against the research runner |
| Row-wise right stopping | Each row is trivial/saturated, or remains unsaturated at unchanged REO over two extensions; changed comparison intervals reset its history | Dedicated mixed-row test |
| Refinement on both sides of a frontier | Partition or saturation changes, including different partitions with equal REO; local patience/displacement controls | Geometry, nonoverlap and same-order frontier tests |
| CEM neighbour harmonization | First pass on the whole mesh, then new cells and their neighbours with propagation until stable; separate `A_harm` control; independent diagnostic outputs remain frozen | Monotonicity, local propagation, final full-mesh audit and independent-winner tests |
| No extra EM-bank harmonization | All candidates are already available throughout the diagram | Bank rescoring and candidate-minimum checks |
| Optional trivial quadrants | Trivial anchors, stability checks and a persistent shared audit budget across overlapping quadrants; exempt cells remain unmeasured and excluded | Quadrant coverage and audit-budget tests |
| Parallel execution and interrupted fitting | Independent repetitions; per-repetition/per-EM-start caches, and staged-CEM per-pilot, per-campaign-execution and harmonization caches, with seed and input validation | Serial/parallel and partial-resume regressions |
| Correct irregular-cell displays | Actual adaptive rectangles, extraction inside a cell and ordered lambda paths | Plot/extraction tests and visual inspection |
| Compact saved results | Banks, independent winners, labels, data and adaptive geometry remain available for inspection and rescoring | Compact-object and data-fingerprint tests |
| Pedagogical documentation | Executed README, 26 help topics, three executable vignettes and an installed quick-start script | Fresh source build, help examples, vignette checks and README assertions |

## Staged adaptive CEM

`adaptive_cem.R` owns the CEM pilot, deferred campaigns and final coverage
loop. `adaptive_grid.R` dispatches to it and retains the EM-bank mesh driver.
Shared geometry is in `adaptive_geometry.R`; propagation is in
`adaptive_harmonization.R`. The fixed-grid and historical compatibility
routes remain separate. The maintained numerical kernels are unchanged by
this orchestration correction.

Adaptive CEM now defaults to `saturation_band` (0.10). Explicit policies remain
available. EM-bank and fixed-grid defaults stay `all`. The three saturation
controls, optimization choices and penalty choices remain independent.

Old adaptive CEM caches belong to a different protocol. They are preserved
but are not loaded as staged campaigns. Returned execution records identify
pilot versus measurement, seed, parent-start use and checkpoint reuse.

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
