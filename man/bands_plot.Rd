\name{bands_plot}
\alias{bands_plot}
\title{Plot the Complete or Filtered CEC Diagram}
\usage{bands_plot(fit, filtered = TRUE, main = NULL)}
\arguments{
\item{fit}{A \code{cec_bands} result.}
\item{filtered}{If true, excluded cells are grey. If false, display all final
reference winners.}
\item{main}{Optional plot title.}
}
\description{
Draw the final diagram in \eqn{(\log C,\lambda)} coordinates. Colours indicate
the number of groups and crosses identify cells where diagnostics were measured.
}
\details{
The horizontal dashed line is \eqn{\lambda=1}. Different partitions with the
same group count have the same colour; use the candidate keys in
\code{fit$tiles} and \code{\link{bands_partitions}} to distinguish them.
In component mode, accepted regions may contain cells without a measured
diagnostic. Crosses mark measurement locations, not necessarily barycentres or
the displayed region's geometrical centre.
}
\value{The input object, returned invisibly.}
\seealso{\code{\link{bands_partitions}}, \code{\link{bands_control}}}
