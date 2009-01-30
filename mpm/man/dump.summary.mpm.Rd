\name{dump.summary.mpm}
\alias{dump.summary.mpm}
\title{Dump the summary output for an mpm object to a text file}
\description{
  Output the mpm summary to a tab-demimited file for processing by other programs (Excel, Spotfire...)
  If the filename is empty, return the data instead of writing to file (useful for web services).
}
\usage{
  \method{dump}{summary.mpm}(x, filename = "")
}
\arguments{
  \item{x}{object of class \code{summary.mpm} as produced by the function of the same name}
  \item{filename}{path to the output file}
}
% TODO explain more extensively
\details{
  Polar (spherical) coordinates are added if the \code{summary.spm} object
  contains 2 (3) dimensions.
}
\value{
  the output is returned invisibly
}
\author{R. Verbeeck, T. Verbeke}
\seealso{\code{\link{summary.mpm}}}
\keyword{manip}