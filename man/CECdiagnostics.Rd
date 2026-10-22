% Maintained by hand; see CONTRIBUTING.md.
\name{CECdiagnostics}
\alias{CECdiagnostics}
\title{Interpret exhaustive stability, RSI and saturation}
\description{
Definitions used by the consolidated grid workflow, and the distinction between measuring a diagnostic and retaining a cell.
}
\section{Exhaustive comparisons}{
Let P_l be the output of original execution l and Q_b the projected
output of bootstrap execution b, both labeling the same n original observations.
For an agreement a, the averages are
\deqn{\mathrm{Stab}_{algo}=\frac{2}{L(L-1)}\sum_{l<l'} a(P_l,P_{l'}),
\qquad \mathrm{Stab}_{boot}=\frac{1}{LB}\sum_{l,b}a(P_l,Q_b).}
RSI is the second divided by the first, when the first is positive.
Repetitions retain their own selected winners; the pooled best fit is not
inserted into these averages. Bootstrap indices and optimizer randomizations
are independent of the original runs. A one-group output in every repetition
has agreement one, but a single original repetition cannot establish this.

RSI is not a confidence probability or a proof that a partition is scientifically
meaningful. It compares reproducibility under resampling with agreement under
optimization randomization. It can exceed one and is not clipped. High matching
agreement can be dominated by a large group; inspect sizes and alternative
partitions, especially when small outlier groups are present.
}
\section{Saturation}{
Saturation is the fraction of observations assigned to groups whose
Gaussian density constraint is active, allowing a numerical tolerance. It is
stored as a fraction between zero and one. A singular group is counted as
saturated in the numeric profile. Saturation alone is not evidence that the
data lack structure; it reports dependence on the concentration ceiling.
}
\section{Measurement policy}{
For fixed grids and adaptive EM banks, policies are applied to the final
representative of each cell. Adaptive CEM first applies the policy to its
completed pilot, then rechecks coverage after measurement improvements and
harmonization. Existing campaigns stay frozen if the representative changes.
\describe{
\item{\code{all}}{Measure every cell, including trivial and saturated cells.}
\item{\code{unsaturated}}{Measure only cells with saturation zero; trivial cells are included.}
\item{\code{nontrivial_unsaturated}}{Measure only cells with more than one occupied group and saturation zero.}
\item{\code{saturation_band}}{Measure cells whose saturation is at most rsi_saturation_max, including trivial cells.}
}
The policy controls both stability averages as well as their ratio; it does
not change the score. In adaptive CEM it gates the entire independent fitting
campaign, after the pilot; elsewhere it gates comparisons of already fitted
outputs. Adaptive CEM defaults to \code{saturation_band}; fixed grids and
adaptive EM default to \code{all}. Explicit policies override these defaults.
Unmeasured cells contain NA and are not retained. The measured flag requires finite algorithmic and bootstrap
averages. For EM grids, a broader policy can be evaluated by rescoring the banks.
}
\section{Filtering and interpretation}{
Consolidated defaults retain cells with Stab_algo > 0.85, RSI > 0.95
and saturation <= 0. These are user-set requirements, not calibrated confidence
levels. The RSI measurement-band limit is separate from the final saturation
threshold, and both differ from the CEM search stopping control.

The output is a multiresolution collection. Deduplicate retained cells by
partition_id within the same grid, then inspect each partition and the region
where it occurs. Partitions need not be nested, and equal group counts need
not describe the same organization. No theorem ensures that filtering will
leave only a few partitions or select a scientifically correct one.
}
\seealso{\link{CECstability}, \link{CECselectStableLambdas}, \link{CECresults}.}
