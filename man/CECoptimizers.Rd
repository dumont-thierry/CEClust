% Maintained by hand; see CONTRIBUTING.md.
\name{CECoptimizers}
\alias{CECoptimizers}
\title{Optimization, hard scoring and the BIC order penalty}
\description{
Common definitions and all optimizer_control fields for the numeric Gaussian engines.
}
\section{The criterion}{
For a hard partition with occupied groups j=1,...,r, mass
\eqn{p_j=n_j/n} and reference density fitted under C, the empirical score is
\deqn{H_\lambda=\lambda\left(-\sum_j p_j\log p_j\right)
       -\frac1n\sum_i\log g_{\widehat\theta_{z_i,C}}(X_i).}
With \code{penalty="bic"}, the reported objective adds \eqn{\beta(r-1)},
where \eqn{\beta=(d+1)(d+2)\log(n)/(4n)}. This equals the full Gaussian
mixture dimension penalty up to a constant independent of r. The penalty is
not multiplied by lambda. This is a penalty on the hard composite criterion,
not ordinary mixture-likelihood BIC. \code{Hphi=objective/lambda} preserves
the historical score scale; compare scores only for the same data, C and lambda.
}
\section{CEM and EM banks}{
CEM proposes hard reassignment steps and accepts improvements of the
refitted hard score. Numeric penalized CEM also considers component-deletion
proposals. A numerical floor used to propose moves is not added to the score.
The multi-start search can increase its initialization order; it remains local.

An EM bank fits soft Gaussian mixtures at orders 1,...,Rmax, or the supplied \code{orders}, under Cmax, with
Ninit starts per order. Initial allocations are softened to positive weights.
MAP produces hard partitions; empty MAP groups are removed. A hard trivial
candidate is always retained. Equivalent candidates are deduplicated, but
start outcomes and failures remain in the bank. At every target C, hard masses
and Gaussian profiles are estimated anew before scoring each candidate.
The same bank can be used for every lambda and either penalty. CEM decreases
the hard score; EM increases the soft likelihood used to generate candidates.
These are different objectives, and neither certifies a global hard optimum.
}
\section{Control fields}{
\describe{
\item{\code{orders=NULL}}{Optional distinct integer soft orders, between 1 and n, for EM. By default all integers 1,...,Rmax are tried. When supplied, only these orders are fitted; an explicit Rmax must cover them. A trivial hard baseline is retained even when order 1 is omitted. Adaptive CEM uses Rmax and subdivision instead.}
\item{\code{Rmax}}{For EM, largest soft mixture order, capped at n; MAP can leave fewer groups. For numeric CEM, initial search order, not a final cap. Default r0, or min(20,n) in direct classification when r0 is NULL.}
\item{\code{Ninit}}{EM starts per order; numeric CEM starts per search batch. Default Nshots in direct classification, Nshots_fresh on a grid.}
\item{\code{Cmax}}{EM density ceiling, default target C or maximum C_grid on a fixed grid, and 64 times the largest initial C on an adaptive grid. Must be at least every target C. Numeric CEM uses the target C directly.}
\item{\code{max_iter}}{EM iteration limit, or CEM reassignment limit per descent phase; defaults to Nloop.}
\item{\code{tolerance}}{Nonnegative numerical stopping/improvement tolerance; default 1e-8.}
\item{\code{soften}}{EM initialization mixing fraction strictly between 0 and 1; default 0.02. Not a hard-score parameter.}
\item{\code{seed}}{Optional nonnegative integer seed. In the full-grid interface it overrides algo_seed; in the one-C interface it overrides seed.}
\item{\code{rinit_max}}{Numeric CEM initialization-order ceiling, default 320, raised to at least Rmax. It is not an EM order control.}
\item{\code{min_child}}{Minimum size used in CEM subdivision starts, default 3. This is an initialization constraint, not a minimum size imposed on every final cluster.}
\item{\code{confirmations}}{Numeric CEM search confirmations, default 2.}
\item{\code{max_batches}}{Maximum CEM search batches, default 8.}
\item{\code{saturation_stop}}{Numeric CEM search stopping saturation, default 0.05. This differs from both the RSI measurement policy and the diagram's final sat_threshold.}
}
Unknown fields are rejected. Method-specific controls have no effect on the
other optimizer. Inspect saved controls and convergence records when comparing budgets.
}
\section{Input and bounds}{
New engines accept finite numeric data and full-covariance Gaussian
references. A Gaussian peak bound C constrains the covariance determinant,
not each eigenvalue separately. Saturation marks activation of this bound;
singular groups use the constrained limiting profile. It is not a minimum
cluster-size rule. Numeric data are not standardized automatically.
Scale or transform them once before fitting if that is the intended analysis,
and retain that transformation for reproducibility. C depends on those units.
Categorical and mixed-data unpenalized CEM remain on the historical route.
}
\seealso{\link{CECclassif}, \link{CECfitBoundGrid}, \link{CECdiagnostics}, \link{CECmigration}.}
