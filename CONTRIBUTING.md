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
R CMD check --no-manual CEClust_0.3.0.9002.tar.gz
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
