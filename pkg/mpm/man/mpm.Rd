\name{mpm}
\alias{mpm}
\alias{print.mpm}
\title{Spectral Map Analysis}
\description{Produces an object of class \code{mpm} that allows for
  exploratory multivariate analysis of large data matrices, such as gene
  expression data from microarray experiments.
}
\usage{
mpm(data, logtrans = TRUE, logrepl = 1e-09, 
    center = c("double", "row", "column", "global", "none"),
    normal = c("global", "row", "column", "none"),
    closure = c("none", "row", "column", "global", "double"), 
    row.weight = c("constant", "mean", "median", "max", "logmean", "RW"),
    col.weight = c("constant", "mean", "median", "max", "logmean", "CW"),
    CW = rep(1, ncol(data) - 1),
    RW = rep(1, nrow(data)),
    pos.row = rep(FALSE, nrow(data)),
    pos.column = rep(FALSE, ncol(data) - 1))

\method{print}{mpm}(x, digits = 3, ...)
}
\arguments{
  \item{data}{a data frame with the row descriptors in the first
    column. For microarray data rows indicate genes and columns
    biological samples.}
  \item{logtrans}{an optional logical value. If \code{TRUE}, data are
    first transformed to logarithms (base e) before the other
    operations. Non-positive numbers are replaced by \code{logrepl}. If
    \code{FALSE}, data are left unchanged. Defaults to \code{TRUE}.}
  \item{logrepl}{an optional numeric value that replaces non-positive
    numbers in log-transformations. Defaults to \code{1e-9}.}
  \item{closure}{optional character string specifying the closure
    operation that is carried out on the optionally log-transformed data
    matrix. If \kbd{"double"}.data are divided by row- and
    column-totals. If \kbd{"row"} data are divided by
    row-totals. If \kbd{"column"} data are divided by
    column-totals. If \kbd{"none"} no closure is carried
    out. Defaults to \kbd{"none"}.}
  \item{center}{optional character string specifying the centering
    operation that is carried out on the optionally log-transformed,
    closed data matrix. If \kbd{"double"} both row- and
    column-means are subtracted. If \kbd{"row"} row-means are
    subtracted. If \kbd{"column"} column-means are subtracted. If
    \kbd{"none"} the data are left uncentered. Defaults to
    \kbd{"double"}.}
  \item{normal}{optional character string specifying the normalization
    operation that is carried out on the optionally log-transformed,
    closed, and centered data matrix. If \kbd{"global"} the data
    are normalized using the global standard deviation. If
    \kbd{"row"} data are divided by the standard deviations of
    the respective row. If \kbd{"column"} data are divided by
    their respective column standard deviation. If \kbd{"none"}
    no normalization is carried out. Defaults to \kbd{"global"}.}
  \item{row.weight}{optional character string specifying the weights of
    the different rows in the analysis. This can be
    \kbd{"constant"}, \kbd{"mean"}, \kbd{"median"},
    \kbd{"max"}, \kbd{"logmean"}, or \kbd{"RW"}. If
    \kbd{"RW"} is specified, weights must be supplied in the
    vector \kbd{RW}. In other cases weights are computed from the
    data. Defaults to \kbd{"constant"}, i.e. constant weighting.}
  \item{col.weight}{optional character string specifying the weights of
    the different columns in the analysis. This can be
    \kbd{"constant"}, \kbd{"mean"}, \kbd{"median"},
    \kbd{"max"}, \kbd{"logmean"}, or \kbd{"CW"}. If
    \kbd{"CW"} is specified, weights must be supplied in the
    vector \code{CW}. In other cases weights are computed from the
    data. Defaults to \kbd{"constant"}, i.e. constant weighting.}
  \item{CW}{optional numeric vector with external column
    weights. Defaults to 1 (constant weights).}
  \item{RW}{optional numeric vector with external row weights. Defaults
    to 1 (constant weights).}
  \item{pos.row}{logical vector indicating rows that are not to be
    included in the analysis but must be positioned on the projection
    obtained with the remaining rows. Defaults to \code{FALSE}.}
  \item{pos.column}{logical vector indicating columns that are not to be
    included in the analysis but must be positioned on the projection
    obtained with the remaining columns. Defaults to \code{FALSE}.}
  \item{x}{object of class \kbd{"mpm"} to be printed.}
  \item{digits}{number of digits to be printed. Defaults to \code{3}.}
  \item{...}{further arguments passed to the (default) print method.}
}
\details{The function \code{mpm} presents a unified approach to
  exploratory multivariate analysis encompassing principal component
  analysis, correspondence factor analysis, and spectral map
  analysis. The algorithm computes projections of high dimensional data
  in an orthogonal space. The resulting object can subsequently be used
  in the construction of biplots (i.e. \code{plot.mpm}).

The projection of the pre-processed data matrix in the orthogonal space
is calculated using the \code{La.svd} function.}
\value{An object of class \code{mpm} representing the projection of data
  after the different operations of transformation, closure, centering,
  and normalization in an orthogonal space. Generic functions
  \code{plot} and \code{summary} have methods to show the results of the
  analysis in more detail.
The object consists of the following components:
  \item{TData}{matrix with the data after optional log-transformation,
    closure, centering and normalization.}
  \item{row.names}{character vector with names of the row elements as
    supplied in the first column of the original data matrix}
  \item{col.names}{character vector with the names of columns obtained
    from the column names from the original data matrix}
  \item{closure}{closure operation as specified in the function call}
  \item{center}{centering operation as specified in the function call}
  \item{normal}{normalization operation as specified in the function call}
  \item{row.weight}{type of weighting used for rows as specified in the function call}
  \item{col.weight}{type of weighting used for columns as specified in
    the function call}
  \item{Wn}{vector with calculated weights for rows}
  \item{Wp}{vector with calculated weights for columns}
  \item{RM}{vector with row means of original data}
  \item{CM}{vector with column means of original data}
  \item{pos.row}{logical vector indicating positioned rows as specified
    in the function call}
  \item{pos.column}{logical vector indicating positioned columns as
    specified in the function call}
  \item{SVD}{list with components returned by \code{La.svd}}
  \item{eigen}{eigenvalues for each orthogonal factor from obtained from
    the weighted singular value decomposition}
  \item{contrib}{contributions of each factor to the total variance of
    the pre-processed data, i.e. the eigenvalues as a fraction of the
    total eigenvalue.}
  \item{call}{the matched call.}
}
\references{
Wouters, L., Goehlmann, H., Bijnens, L., Kass, S.U., Molenberghs, G.,
Lewi, P.J. (2003). Graphical exploration of gene expression data: a
comparative study of three multivariate methods. \emph{Biometrics}
\bold{59}, 1131-1140.
} 
\author{Luc Wouters}
\note{Principal component analysis is defined as the projection onto an
  orthogonal space of the column-centered and column-normalized data. In
  correspondence factor analysis the data are pre-processed by double
  closure, double centering, and global normalization. Orthogonal
  projection is carried out using the weighted singular value
  decomposition. Spectral map analysis is in essence a principal
  component analysis on the log-transformed, double centered and global
  normalized data. Weighted spectral map analysis has been proven to be
  successful in the detection of patterns in gene expression data
  (Wouters et al., 2003).}

\seealso{\code{\link{plot.mpm}}, \code{\link{summary.mpm}}}

\examples{
  data(Golub)
  # Principal component analysis
  r.pca <- mpm(Golub[,1:39], center = "column", normal = "column")
  # Correspondence factor analysis
  r.cfa <- mpm(Golub[,1:39],logtrans = FALSE, row.weight = "mean",
             col.weight = "mean", closure = "double")
  # Weighted spectral map analysis
  r.sma <- mpm(Golub[,1:39], row.weight = "mean", col.weight = "mean")
}
\keyword{multivariate}

