\name{plot.mpm}
\alias{plot.mpm}
\title{Spectral Map Plot of Multivariate Data...}
\usage{\method{plot}{mpm}(x, scale=c("singul", "eigen", "uvr", "uvc"), dim=c(1, 2), zoom=rep(1,
    2), show.row=c("all", "position"), show.col=c("all", "position"),
    col.group=rep(1, length(x$col.names)), colors=c("orange1", "red",
    rainbow(length(unique(col.group)), start = 2/6, end = 4/6)),
    col.areas=TRUE, col.symbols=c(1, rep(2,
    length(unique(col.group)))), sampleNames=TRUE, rot=rep(-1,
    length(dim)), labels, label.tol=1, label.col.tol=1, lab.size=0.725,
    col.size=10, row.size=10, do.smoothScatter=FALSE, do.plot=TRUE, ...)
}
\description{Spectral Map Plot of Multivariate Data
Produces a spectral map plot (biplot) of an object of class \code{mpm}}
\details{Spectral maps are special types of biplots with the area of the symbols
proportional to some measure, usually the row or column mean value and an
identification of row- and column-items.  For large matrices, such as gene
expression data, where there is an abundance of rows, this can obscure the
plot.  In this case, the argument \code{label.tol} can be used to select the
most informative rows, i.e. rows that are most distant from the center of
the plot.  Only these row-items are then labeled and represented as circles
with their areas proportional to the marginal mean value.  For the
column-items it can be useful to apply some grouping specified by
\code{col.group}. Examples of groupings are different pathologies, such as
specified in \code{Golub.grp}}
\value{An object of class \code{plot.mpm} that has the following
components: \item{Rows}{a data frame with the X and Y coordinates of the
rows and an indication \code{Select} of whether the row was selected
according to \code{label.tol}} \item{Columns}{a data frame with the X and
Y coordinates of the columns}}
\note{\code{value} is returned invisibly, but is available for further use
when an explicit assignment is made}
\author{Luc Wouters}
\seealso{\code{\link{mpm}}, \code{\link{summary.mpm}}}
\references{Wouters, L., Goehlmann, H., Bijnens, L., Kass, S.U.,
Molenberghs, G., Lewi, P.J. (2003). Graphical exploration of gene
expression data: a comparative study of three multivariate methods.
\emph{Biometrics} \bold{59}, 1131-1140.}
\keyword{multivariate}
\keyword{hplot}
\arguments{\item{x}{object of class \code{mpm} a result of a call to \code{mpm}.}
\item{scale}{optional character string specifying the type of factor scaling
of the biplot. This can be either \kbd{"singul"} (singular value scaling),
\kbd{"eigen"} (eigenvalue scaling), \kbd{"uvr"} (unit row-variance
scaling), \kbd{"uvc"} (unit column-variance scaling). The latter is of
particular value when analyzing large matrices, such as gene expression
data. Singular value scaling \kbd{"singul"} is customary in spectral map
analysis. Defaults to \kbd{"singul"}.}
\item{dim}{optional principal factors that are plotted along the horizontal
and vertical axis. Defaults to \code{c(1,2)}.}
\item{zoom}{optional zoom factor for row and column items. Defaults to
\code{c(1,1)}.}
\item{show.row}{optional character string indicating whether all rows
(\kbd{"all"}) are to be plotted or just the positioned rows
\kbd{"position"}.}
\item{show.col}{optional character string indicating whether all columns
(\kbd{"all"}) are to be plotted or just the positioned columns
\kbd{"position"}.}
\item{col.group}{optional vector (character or numeric) indicating the
different groupings of the columns, e.g. \code{Golub.grp}. Defaults to 1.}
\item{colors}{vector specifying the colors for the annotation of the plot;
the first two elements concern the rows; the third till the last element
concern the columns; the first element will be used to color the unlabeled
rows; the second element for the labeled rows and the remaining elements
to give different colors to different groups of columns.}
\item{col.areas}{logical value indicating whether columns should be plotted
as squares with areas proportional to their marginal mean and colors
representing the different groups (\code{TRUE}), or with symbols
representing the groupings and identical size (\code{FALSE}). Defaults to
\code{TRUE}.}
\item{col.symbols}{vector of symbols when \code{col.areas=FALSE} corresponds
to the \code{pch} argument of the function \code{plot}.}
\item{sampleNames}{Either a logical vector of length one or a character
vector of length equal to the number of samples in the dataset. If a
logical is provided, sample names will be displayed on the plot
(\code{TRUE}; default) or not (\code{FALSE}); if a character vector is
provided, the names provided will be used to label the samples instead of
the default column names.}
\item{rot}{rotation of plot. Defaults to \code{c(-1,-1)}.}
\item{labels}{character vector to be used for labeling points on the graph;
if \code{NULL}, the row names of \code{x} are used instead}
\item{label.tol}{numerical value specifying either the percentile
(\code{label.tol<=1}) of rows or the number of rows (\code{label.tol>1})
most distant from the plot-center (0,0) that are labeled and are plotted
as circles with area proportional to the marginal means of the original
data.}
\item{label.col.tol}{numerical value specifying either the percentile
(\code{label.col.tol<=1}) of columns or the number of columns (\code{label.col.tol>1})
most distant from the plot-center (0,0) that are labeled and are plotted
as circles with area proportional to the marginal means of the original
data.}
\item{lab.size}{size of identifying labels for row- and column-items as
\code{cex} parameter of the \code{text} function}
\item{col.size}{size in mm of the column symbols}
\item{row.size}{size in mm of the row symbols}
\item{do.smoothScatter}{use smoothScatter or not instead of plotting
individual points}
\item{do.plot}{produce a plot or not}
\item{...}{further arguments to \code{eqscaleplot} which draws the canvas
for the plot; useful for adding a \code{main} or a custom \code{sub}}
}
\examples{# Weighted spectral map analysis
data(Golub) # Gene expression data of leukemia patients
data(Golub.grp) # Pathological classes coded as 1, 2, 3
r.sma <- mpm(Golub[,1:39], row.weight = "mean", col.weight = "mean")
# Spectral map biplot with result
r <- plot(r.sma, label.tol = 20, scale = "uvc",
col.group = (Golub.grp)[1:38], zoom = c(1,1.2), col.size = 5)
Golub[r$Rows$Select, 1] # 20 most extreme genes}

