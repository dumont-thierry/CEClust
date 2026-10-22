% Maintained by hand; see CONTRIBUTING.md.
\name{CECplotGrid}
\alias{CECplotGrid}
\title{Plot the grid and its diagnostic filtering}
\description{
Display the occupied group count across C and lambda, with diagnostic rejection overlays and optional boundaries between partition regimes.
}
\usage{
CECplotGrid(x, ...)
}
\arguments{
\item{x}{A grid returned by \code{CECfitBoundGrid()} or \code{CECfitPreset()}. For a consolidated one-C diagnostic, use its \code{grid} component.}
\item{...}{Plot controls described below; no fitting is performed.}
}
\value{
The plotted summary data, returned invisibly. The stored grid is unchanged.
}
\details{
The axes use the supplied coordinates on their numerical scale. A larger
group count does not by itself imply stability. Inspect the legend and summary
for rejection reasons, especially when diagnostics were deliberately not measured.
Two cells of the same REO can describe different partitions. The function uses
base R graphics; wrap it in \code{pdf()} or \code{png()} to save a figure.
}
\section{Plot controls}{
\describe{
\item{\code{stab_algo_threshold, rsi_threshold, sat_threshold}}{Defaults from grid metadata. They change the displayed filter, not stored results.}
\item{\code{main}}{Plot title.}
\item{\code{selected_C, selected_lambda}}{Optional highlighted cell.}
\item{\code{reo_col}}{Vector of colors for occupied group counts; use contrasting colors.}
\item{\code{show_legend}}{Display the legend; default TRUE.}
\item{\code{draw_zones}}{Draw partition-regime boundaries; default TRUE.}
\item{\code{sensibility}}{Agreement threshold for display boundaries, default 0.95. It is not the RSI threshold and does not change admissibility.}
\item{\code{bad_fill, bad_hatch_col}}{Colors for rejected areas and hatching.}
\item{\code{grid_summary}}{Optional already prepared summary for an advanced custom view; otherwise thresholds are applied internally.}
}
}
\seealso{\link{CECdiagnostics}, \link{CECplotPartition}, \link{CECplotPath}, \link{CECexplore}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
CECplotGrid(g, main = "CEC without BIC", reo_col = c("#0072B2", "#D55E00", "#009E73"))
}
