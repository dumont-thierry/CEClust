# Maintaining CEClust

This development source combines the adaptive-grid and candidate-band engines
behind `CECdiagram`. Historical CEM and EM-bank entry points remain available.
Keep numerical changes and documentation changes
separate in review. The public interface is listed in `NAMESPACE`.

## Documentation sources

- `man/*.Rd` is the maintained reference. Each exported function has its own
  page; the conceptual pages explain shared scoring, diagnostics and objects.
- `README.md` is maintained directly for the common interface.
- `vignettes/unified-methods.Rmd` describes both methods and their defaults.
- `vignettes/intro.Rmd` is the executed introductory analysis.
- `vignettes/minimal-api.Rmd` covers controls and reproducibility.
- `vignettes/migration.Rmd` explains historical behavior.
- `inst/examples/quick_start.R` is an installed, runnable script.

The package contains historical roxygen comments, but documentation and
`NAMESPACE` are maintained explicitly. Do not run roxygen over the old source
as a routine documentation build: it can export research helpers and overwrite
the reviewed reference pages. Public source comments point to their maintained
help topics.

## Build and check

From a working R development environment with a C++ toolchain, RcppArmadillo,
testthat, knitr, rmarkdown and Pandoc available:

```sh
R CMD build CEClust
R CMD check --no-manual CEClust_0.4.0.9000.tar.gz
```

Do not disable vignette building for a distribution. The introductory and
advanced guides execute small examples; parallel and interactive examples are
explicitly marked as such. Help examples are run by R CMD check. Inspect the
rendered guides as well as the check log. Build a Windows binary with
`R CMD INSTALL --build` on the matching R platform and test that archive in a
separate library before distributing it.

When changing an argument, update its help, shared control reference,
README/vignettes if affected, returned metadata and tests. Keep the optimizer,
penalty and diagnostic measurement policy independent. Never replace the
independent winners by the pooled best partition when reporting stability.

The adaptive engine is packaged and must not source research files at runtime.
Its imported geometric rules are tested against the archived patches. Test all
four optimizer/penalty combinations, independent winners, irregular-cell plots,
checkpoint reuse and budget-stop reporting when changing it. Fixed-grid and
legacy compatibility behaviour must remain available.

## Engine boundaries and protocol tests

`diagram.R` owns the common interface, its budget profiles and the automatic
initial grid. `CECinitialGrid` uses the full-sample Gaussian ML peak as Cmin
and centres `2^(1:10)*Cmin`. The common adaptive interface disables left
extension. Direct historical interfaces keep their defaults. `bands_*.R`
contains the candidate-band engine and its separate optional cleanup. Do not
silently transfer one engine's diagnostic rules to the other. `tests/unified.R`
covers dispatch, initialization, checkpoint use and parallel execution.

`minimal_api.R` selects the public route. `adaptive_grid.R` validates the
adaptive request and dispatches CEM to `adaptive_cem.R`. The remaining adaptive
bank driver fits EM banks once and only evaluates them as the mesh changes.
`adaptive_geometry.R` owns extension, frontier and screening rules;
`adaptive_harmonization.R` owns CEM neighbour propagation. Numerical objectives
and descents are in `solver_numerics.R` and `optimization.R`. `solver_grid.R`
retains the fixed-grid engine and common output construction. The historical
linked-path engine remains for compatibility; it is not used by these drivers.

For staged CEM, enforce these invariants in `test-staged-cem.R`. The improvement
and renewed-coverage loop below applies to `diagnostic_strategy="cells"`;
the default component strategy additionally follows the rules in the next section:

- one complete original-data pilot execution per new cell, before any L+B work;
- parent starts only in the pilot, never in independent campaigns;
- whole campaigns gated by the current measurement policy on active leaves;
- frozen original/bootstrap outputs after improvement and harmonization;
- renewed border/coverage checks, with pending refinement explicitly reported;
- per-execution checkpoints and identical randomizations across worker counts;
- missing direct diagnostics and winners on unmeasured cells; regional values
  require an identified passing witness, never automatic keep.

A test of score equality or a `diagnostics=FALSE` flag alone cannot establish
these invariants. Preserve the distinction between an execution's internal
starts, independent executions, and the final representative. Test the number,
phase and input data of actual fitting calls, as well as their outputs.


## Regional and reference-bootstrap branches

`component_diagnostics.R` owns geometry, barycentric ordering and first-success
stopping. By default, raw independent diagnostics must pass and the best of the
L original-data outputs must agree with the frozen reference by at least the
algorithmic threshold. Otherwise continue the traversal without replacing the
raw diagnostic values. Test both disagreement directions and score equality.
The gate is optional and applies to selection and plotting, including failed
direct points of exhausted components. `reference_bootstrap.R` owns this shared
gate and the optional conditional alternative with exactly one CEM descent per
projected reference bootstrap. The conditional check runs before filtering
at every measured point only when explicitly enabled. These are explicit
changes to the diagnostics; preserve raw independent outputs and comparisons.
The first condition requires both raw Stab_algo and raw RSI to strictly exceed
their respective configured thresholds. Forward rsi_threshold from every caller;
never test a previously corrected value or use a fixed default instead of the
user's cutoff. Missing values or equality with either cutoff do not pass.
The component strategy freezes the pilot. It must not inject campaign winners
or substitute another score condition for the configured reference agreement.
Saturation still gates candidate
measurement and masks final display. Pointwise and regional provenance must
survive selection, extraction, compaction and plotting.

Test equal-order distinct partitions, disconnected regions, unequal cell
areas, raw axis distances, saturation gating, early stop, all three trigger
conditions (including boundaries), penalized scores, one descent/no random
starts, frozen projection, and serial/parallel/partial resume equivalence.
Keep the cell strategy, gate-disabled independent option and conditional option
as regression baselines. The shared agreement cache must preserve every pair
array, metric, ordering, random draw and threshold decision. Cache keys must
identify actual partitions across cells, not cell-local winner indices.


## Numeric CEM kernel

`src/cem_fused.cpp` shares moments, covariance fitting and raw density evaluation.
`q_descent` dispatches to it; `q_descent_R` is the retained reference. Both are
cloned into the penalty-specific environment. Keep score acceptance, proposal
flooring and singular-support conventions separate. A failed or numerically
ambiguous compiled fit must fall back without discarding its start or consuming
randomness. Preserve engine/fallback records and timings.

`test-compiled-cem` compares profiles, raw densities, starts, accepted moves,
labels and exhaustive pair vectors. Performance validation must keep initial
orders, internal starts, L/B, seeds and workers unchanged. Report whole-operation
times alongside isolated search timings; archive the tested source and package.
