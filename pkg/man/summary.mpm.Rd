\name{summary.spm}
\alias{summary.mpm}
\alias{print.summary.mpm}
\title{Summary Statistics for Spectral Map Analysis}
\description{
  Summary method for object of class \code{mpm}.
}
\usage{
  \method{summary}{mpm}(object, maxdim = 4, ...)
  \method{print}{summary.mpm}(x, digits = 2, what = c("columns", "rows", "all"), ...)
}
\arguments{
  \item{object}{an object of class \code{mpm} resulting from a call to \code{mpm}}
  \item{maxdim}{maximum number of principal factors to be
    reported. Defaults to \code{4}}
  \item{x}{object of class \code{summary.mpm}}
  \item{digits}{accuracy of printing}
  \item{what}{optional character string specifying whether
    \code{\dQuote{columns}}, \code{\dQuote{rows}}, or both
    (\code{\dQuote{all}}) are to appear in the printed report. Defaults
    to \code{\dQuote{columns}}.}
  \item{...}{further arguments to the (default) summary or print
    methods}
}
\details{
The function \code{summary.mpm} computes and returns a list of summary
statistics of the spectral map analysis given in \code{x}.
}
\value{
  An object of class \code{summary.mpm} with the following components:
  \item{call}{the call to \code{mpm}}
  \item{Vxy}{sum of eigenvalues}
  \item{VPF}{a matrix with on the first line the eigenvalues and on the
    second line the cumulative eigenvalues of each of the principal
    factors (\code{PRF1} to \code{PRFmaxdim}) followed by the residual % to ?
    eigenvalues and the total eigenvalue.}
  \item{Rows}{a data frame with summary statistics for the row-items, as
    described below.}
  \item{Columns}{a data frame with with summary statistics for the
    column-items, as described below.\cr\cr The \code{Rows} and
    \code{Columns} data frames contain the following columns:\cr\cr}
  \item{Posit}{binary indication of whether the row or column was
    positioned (\code{1}) or not (\code{0}).}
  \item{Weight}{weight applied to the row or column in the function
    \code{mpm}.}
  \item{PRF1-PRFmaxdim}{factor scores or loadings for the first
  \code{maxdim} factors using eigenvalue scaling.}
  \item{Resid}{residual score or loading not accounted for by the first
    \code{maxdim} factors.}
  \item{Norm}{length of the vector representing the row or column in factor space.}
  \item{Contrib}{contribution of row or column to the sum of eigenvalues.}
  \item{Accuracy}{accuracy of the representation of the row or column by
    means of the first \code{maxdim} principal factors.}}
\references{
  Wouters, L., Goehlmann, H., Bijnens, L., Kass, S.U., Molenberghs, G.,
  Lewi, P.J. (2003). Graphical exploration of gene expression data: a
  comparative study of three multivariate methods. \emph{Biometrics}
  \bold{59}, 1131-1140.
} 
\author{Luc Wouters}
\seealso{\code{\link{mpm}}, \code{\link{plot.mpm}}}

\examples{
  # Example 1 weighted spectral map analysis Golub data
  data(Golub)
  r.sma <- mpm(Golub[,1:39], row.weight = "mean", col.weight = "mean")
  # summary report
  summary(r.sma)
  # Example 2 using print function
  data(Famin81A)
  r.fam <- mpm(Famin81A, row.weight = "mean", col.weight = "mean")
  r.sum <- summary(r.fam)
  print(r.sum, what = "all")
}
\keyword{multivariate}

