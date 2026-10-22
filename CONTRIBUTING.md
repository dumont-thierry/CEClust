# Maintaining CEClust

This development source combines the corrected Gaussian core with explicit
CEM and EM-bank entry points. Keep numerical changes and documentation changes
separate in review. The public interface is listed in `NAMESPACE`.

## Documentation sources

- `man/*.Rd` is the maintained reference. Each exported function has its own
  page; the conceptual pages explain shared scoring, diagnostics and objects.
- `README.Rmd` generates `README.md` with `rmarkdown::render()`.
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
R CMD check --no-manual CEClust_0.3.0.9003.tar.gz
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

`minimal_api.R` selects the public route. `adaptive_grid.R` validates the
adaptive request and dispatches CEM to `adaptive_cem.R`. The remaining adaptive
bank driver fits EM banks once and only evaluates them as the mesh changes.
`adaptive_geometry.R` owns extension, frontier and screening rules;
`adaptive_harmonization.R` owns CEM neighbour propagation. Numerical objectives
and descents are in `solver_numerics.R` and `optimization.R`. `solver_grid.R`
retains the fixed-grid engine and common output construction. The historical
linked-path engine remains for compatibility; it is not used by these drivers.

For staged CEM, enforce these invariants in `test-staged-cem.R`:

- one complete original-data pilot execution per new cell, before any L+B work;
- parent starts only in the pilot, never in independent campaigns;
- whole campaigns gated by the current measurement policy on active leaves;
- frozen original/bootstrap outputs after improvement and harmonization;
- renewed border/coverage checks, with pending refinement explicitly reported;
- per-execution checkpoints and identical randomizations across worker counts;
- missing diagnostics and winners on unmeasured cells, never automatic keep.

A test of score equality or a `diagnostics=FALSE` flag alone cannot establish
these invariants. Preserve the distinction between an execution's internal
starts, independent executions, and the final representative. Test the number,
phase and input data of actual fitting calls, as well as their outputs.
