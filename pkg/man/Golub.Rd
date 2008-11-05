\name{Golub}
\alias{Golub}
\alias{Golub.grp}
\docType{data}
\title{ Golub (1999) Data}
\description{
  Golub et al. (1999) data on gene expression profiles of 38 
  patients suffering from acute leukemia and a validation sample
  of 34 patients. % TODO add details / MIT1 and MIT2 group etc. 
}
\usage{data(Golub)}
\format{
  The expression data are available in data frame \code{Golub} with 5327 observations on the following 73 variables.
  \describe{
    \item{\code{Gene}}{a character vector with gene identifiers} % TODO which ID ?
    \item{\code{1}}{gene expression data for sample 1}
    \item{\code{2}}{gene expression data for sample 2}
    \item{\code{3}}{gene expression data for sample 3}
    \item{\code{4}}{gene expression data for sample 4}
    \item{\code{5}}{gene expression data for sample 5}
    \item{\code{6}}{gene expression data for sample 6}
    \item{\code{7}}{gene expression data for sample 7}
    \item{\code{8}}{gene expression data for sample 8}
    \item{\code{9}}{gene expression data for sample 9}
    \item{\code{10}}{gene expression data for sample 10}
    \item{\code{11}}{gene expression data for sample 11}
    \item{\code{12}}{gene expression data for sample 12}
    \item{\code{13}}{gene expression data for sample 13}
    \item{\code{14}}{gene expression data for sample 14}
    \item{\code{15}}{gene expression data for sample 15}
    \item{\code{16}}{gene expression data for sample 16}
    \item{\code{17}}{gene expression data for sample 17}
    \item{\code{18}}{gene expression data for sample 18}
    \item{\code{19}}{gene expression data for sample 19}
    \item{\code{20}}{gene expression data for sample 20}
    \item{\code{21}}{gene expression data for sample 21}
    \item{\code{22}}{gene expression data for sample 22}
    \item{\code{23}}{gene expression data for sample 23}
    \item{\code{24}}{gene expression data for sample 24}
    \item{\code{25}}{gene expression data for sample 25}
    \item{\code{26}}{gene expression data for sample 26}
    \item{\code{27}}{gene expression data for sample 27}
    \item{\code{34}}{gene expression data for sample 34}
    \item{\code{35}}{gene expression data for sample 35}
    \item{\code{36}}{gene expression data for sample 36}
    \item{\code{37}}{gene expression data for sample 37}
    \item{\code{38}}{gene expression data for sample 38}
    \item{\code{28}}{gene expression data for sample 28}
    \item{\code{29}}{gene expression data for sample 29}
    \item{\code{30}}{gene expression data for sample 30}
    \item{\code{31}}{gene expression data for sample 31}
    \item{\code{32}}{gene expression data for sample 32}
    \item{\code{33}}{gene expression data for sample 33}
    \item{\code{39}}{gene expression data for sample 39}
    \item{\code{40}}{gene expression data for sample 40}
    \item{\code{42}}{gene expression data for sample 42}
    \item{\code{47}}{gene expression data for sample 47}
    \item{\code{48}}{gene expression data for sample 48}
    \item{\code{49}}{gene expression data for sample 49}
    \item{\code{41}}{gene expression data for sample 41}
    \item{\code{43}}{gene expression data for sample 43}
    \item{\code{44}}{gene expression data for sample 44}
    \item{\code{45}}{gene expression data for sample 45}
    \item{\code{46}}{gene expression data for sample 46}
    \item{\code{70}}{gene expression data for sample 70}
    \item{\code{71}}{gene expression data for sample 71}
    \item{\code{72}}{gene expression data for sample 72}
    \item{\code{68}}{gene expression data for sample 68}
    \item{\code{69}}{gene expression data for sample 69}
    \item{\code{67}}{gene expression data for sample 67}
    \item{\code{55}}{gene expression data for sample 55}
    \item{\code{56}}{gene expression data for sample 56}
    \item{\code{59}}{gene expression data for sample 59}
    \item{\code{52}}{gene expression data for sample 52}
    \item{\code{53}}{gene expression data for sample 53}
    \item{\code{51}}{gene expression data for sample 51}
    \item{\code{50}}{gene expression data for sample 50}
    \item{\code{54}}{gene expression data for sample 54}
    \item{\code{57}}{gene expression data for sample 57}
    \item{\code{58}}{gene expression data for sample 58}
    \item{\code{60}}{gene expression data for sample 60}
    \item{\code{61}}{gene expression data for sample 61}
    \item{\code{65}}{gene expression data for sample 65}
    \item{\code{66}}{gene expression data for sample 66}
    \item{\code{63}}{gene expression data for sample 63}
    \item{\code{64}}{gene expression data for sample 64}
    \item{\code{62}}{gene expression data for sample 62}
  }
  
  The classes are in a separate numeric vector \code{Golub.grp} with values
  \code{1} for the 38 ALL B-Cell samples, \code{2} for the 9 ALL T-Cell samples
  and \code{3} for the 25 AML samples.
}
\details{
  The original data of Golub et al. (1999) were preprocessed 
  as follows: genes that were called 'absent' in all samples 
  were removed from the data sets, since these measurements 
  are considered unreliable by the manufacturer of the technology. 
  Negative measurements in the data were set to 1. 
  
  The resulting data frame contains 5327 genes of the 6817 
  originally reported by Golub et al. (1999).
}
\note{Luc Wouters et al. (2003), p. 1134 contains a typo 
  concerning the sample sizes of AML- and ALL-type  and erroneously reported }
\source{
  Golub, T. R., Slonim, D. K., Tamayo, P., et al. (1999). Molecular
  classification of cancer: Class discovery and class prediction by 
  gene expression monitoring. Science 286, 531 -- 537.
}
\references{
  Luc Wouters et al. (2003). Graphical Exploration of Gene Expression Data: 
  A Comparative Study of Three Multivariate Methods, Biometrics, 59, 1131-1139.
}
\keyword{datasets}
