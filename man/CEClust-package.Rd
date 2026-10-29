% Maintained by hand; see CONTRIBUTING.md.
\name{CEClust-package}
\alias{CEClust-package}
\alias{CEClust}
\title{CEClust: composite entropy clustering and multiresolution diagnostics}
\description{
Fit hard clusterings, explore density-bound and entropy-weight grids, and assess algorithmic and resampling reproducibility. Numeric Gaussian analyses support CEM or EM candidate banks, with or without an order penalty.
}
\section{Start here}{
Read \code{vignette("intro",package="CEClust")} for a complete runnable
analysis. The advanced guide is \code{vignette("minimal-api",package="CEClust")};
compatibility is covered by \code{vignette("migration",package="CEClust")}.
\enumerate{
\item Fit one partition with \link{CECclassif}, or a fixed or adaptive grid with \link{CECfitBoundGrid}; see \link{CECadaptiveDiagram}.
\item Inspect \link{CECsummariseGrid} and \link{CECplotGrid}.
\item Filter a C slice with \link{CECselectStableLambdas}; inspect retained partitions with \link{CECextractPartition} and \link{CECplotPartition}.
\item Compare penalties on the same EM banks with \link{CECrescoreGrid}.
}
The diagram is a multiresolution representation, not a mandatory selector of
one partition. Explicit optimizer choices are recommended for new analyses.
}
\section{Concepts and reference}{
\link{CECoptimizers} defines the criterion, penalty, bounds and budgets.
\link{CECdiagnostics} defines exhaustive stability, RSI and saturation.
\link{CECresults} describes saved objects and independent winners.
\link{CECmigration} distinguishes current and historical workflows.
Theoretical reference families can be broader than Gaussian densities; the
new optimization implementations documented here are numeric Gaussian engines.
}
\section{Version and source}{
This help describes local development version 0.3.0.9009. The published main
branch currently contains 0.3.0.9008. The archive branches preserve historical implementations.
Check \code{packageVersion("CEClust")} after installation.
Source: \url{https://github.com/dumont-thierry/CEClust}.
Maintainer: Thierry Dumont, \email{tdumont@parisnanterre.fr}.
}
\seealso{\link{CECclassif}, \link{CECfitBoundGrid}, \link{CECfitPreset}, \link{CECstability}.}
