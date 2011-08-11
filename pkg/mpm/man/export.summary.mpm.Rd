\name{export.summary.mpm}
\alias{export.summary.mpm}
\title{Export the summary output for an mpm object to a text file...}
\usage{export.summary.mpm(x, filename="", ...)
}
\description{Export the summary output for an mpm object to a text file
Output the mpm summary to a tab-demimited file for processing by other
programs (Excel, Spotfire...)  If the filename is empty, return the data
instead of writing to file (useful for web services).}
\details{Polar (spherical) coordinates are added if the \code{summary.spm} object
contains 2 (3) dimensions.}
\value{the output is returned invisibly}
\author{Rudi Verbeeck, Tobias Verbeke}
\seealso{\code{\link{summary.mpm}}}
\keyword{manip}
\arguments{\item{x}{object of class \code{summary.mpm} as produced by the function of
the same name}
\item{filename}{prefix used to name the output file following <filename>_xyz.txt}
\item{...}{further arguments; currently none are used}
}

