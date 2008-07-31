###
### Multivariate Projection Methods, mpm
###
### Copyright 2003-2008 Luc Wouters <wouters_luc@telenet.be>
###
### This file is part of the mpm package for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

### Web service enabled by Rudi Verbeeck <rverbeec@prdbe.jnj.com>
### Included smoothScatter functionality from geneplotter package (BioConductor).

mpm <- function(data, # input data; column 1 contains row names
	logtrans = TRUE,    # logtransform input data or not (`reexpression')
  logrepl = 1e-9,     # replace input data <= 0 with this value before taking the log
	center = c("double", "row", "column", "global", "none"),
	normal = c("global", "row", "column", "none"),
	closure = c("none", "row", "column", "global", "double"),
	row.weight = c("constant", "mean", "median", "max", "logmean", "RW"),
	col.weight = c("constant", "mean", "median", "max", "logmean", "CW"),
	CW = rep(1, ncol(data)-1),  # Column weight vector (default all 1) 
	RW = rep(1, nrow(data)),    # Row weight vector (default all 1) 
	pos.row = rep(FALSE, nrow(data)),            # Positioned rows and columns are not taken into account during calculation,
	pos.column = rep(FALSE, ncol(data) - 1)){    # but are still plotted at the correct location
  
  
  ### Error checking and argument matching
  if (missing(data))
      stop("Argument \"data\" is missing, with no default")
  NData <- as.matrix(data[, -1]) # drop column 1 with row names
  if (any(is.na(NData)) || !is.numeric(NData))
      stop("Data must be numeric without NA's")
  if (length(pos.row) != dim(NData)[1])
      stop("Length of pos.row argument not equal to number of rows in table")
  if (length(pos.column) != dim(NData)[2])
      stop("Length of pos.column argument not equal to number of columns in table")
  if (length(RW) != nrow(NData))
      stop("Length of RW not equal to number of rows in table")
  if (length(CW) != ncol(NData))
      stop ("Length of CW not equal to number of columns in table")
    
    
  # Parse other arguments, store in variables
  # If no value is specified, the first option in the list is taken as the default value  
  center <- match.arg(center)
  normal <- match.arg(normal)
  closure <- match.arg(closure)
  row.weight <- match.arg(row.weight)
  col.weight <- match.arg(col.weight)
  logrepl <- max(1e-9, logrepl)
  
  ### Get row-descriptors as character variable from column 1
  Row.Names <- as.character(data[, 1]) # TV: simplified
  
  ###################
  ### Positioning ###
  ################### 
  
  #### Determine positioned rows-columns, set to NA 
  RData <- NData
  RData[pos.row, ] <- NA
  RData[, pos.column] <- NA
  
  ####################
  ### Reexpression ###
  ####################
  
  ### Logarithmic transform with replacement of non-positive numbers
  
  ## Define logtransform function
  logtransf <- function(x, logrepl, comp)
  {
    if (any(x <= 0, na.rm = TRUE)) # Test if there is non-positive data
    # (RV: added na.rm=T: test fails when positioning is used)
    {
      warning(paste("Non-positive data replaced by",
              logrepl, "in computing", comp, "in: spectralmap.\n"),
          call. = FALSE)
      # Replace non-positive data by logrepl
      x[x <= 0] <- logrepl
    }
    return(log(x))
  }

  LData <- if (logtrans) logtransf(NData, logrepl, comp = "logarithms") else NData # just copy the data from the previous step
  
  ### means of original data matrix
  RM <- rowMeans(NData[, !pos.column]) ### transformation does not 
                                       ### have impact on this !!
  CM <- colMeans(NData[!pos.row, ])   
  
  ### define weights
  Wn <- pmax(0, switch(row.weight,
      constant = rep(1, length(RM)),
  	  mean = apply(RData, 1, mean, na.rm = TRUE),
  	  median = apply(RData, 1, median, na.rm = TRUE),
  	  max = apply(RData, 1, max,na.rm = TRUE),
  	  logmean = apply(logtransf(RData, logrepl, comp = "logmean weights"),
                    1, mean, na.rm = TRUE),
  	  RW = RW))

  Wp <- pmax(0,switch(col.weight,
  	  constant = rep(1, length(CM)),
  	  mean = apply(RData, 2, mean, na.rm = TRUE),
  	  median = apply(RData, 2, median, na.rm = TRUE),
  	  max = apply(RData,2,max,na.rm=TRUE),
  	  logmean = apply(logtransf(RData, logrepl, comp = "logmean weights"),
                      2, mean, na.rm = TRUE) ,
  	  CW = CW))

  Wn[pos.row] <- 0
  Wp[pos.column] <- 0
  
  Wn <- Wn / sum(Wn)   # normalize weights to unit sum
  Wp <- Wp / sum(Wp)   
  
  ###############
  ### Closure ###
  ###############
  
  if (closure != "none" && any(LData < 0)) 
    warning("Closure operation with non-positive data")
  
  Tn <- rowSums(LData[, !pos.column], na.rm = TRUE) # row totals
  Tp <- colSums(LData[!pos.row, ], na.rm = TRUE)    # column totals
  # Positioned rows (in Tn) / columns (in Tp) are not excluded to preserve the dimensions
  # Consequently, Tn <> Tp <> Tt
  Tt <- sum(LData[!pos.row, !pos.column], na.rm = TRUE)   # global total
  
  ClData <- switch(closure,
      none = LData,
      row = sweep(LData, 1, Tn, "/"),
      column = sweep(LData, 2, Tp, "/"),
      global = LData / Tt,
      double = Tt * sweep(sweep(LData, 1, Tn, "/"), 2, Tp, "/"))
  
  if (any(!is.finite(ClData))) # TV: a bit indirect ?, changed to is.finite 
    stop("Division by 0 in closure operation")
    
  #################
  ### Centering ### 
  #################
  
  ## using weighted means (note: positioned data have weight 0)
  
  Mp <- colSums(sweep(ClData, 1, Wn, "*"))  # weighted row means. Use "sum" as Wn is normalized.
  Mn <- rowSums(sweep(ClData, 2, Wp, "*"))  # weighted column means
  Mg <- sum(Mp * Wp) # weighted global mean (weigh each row by Wn, each column by Wp and add everything together)

  CData <- switch(center,
  	double = Mg + sweep(sweep(ClData, 2, Mp), 1, Mn),
  	row = sweep(ClData, 1, Mn),
  	column = sweep(ClData, 2, Mp), # RV: originally the code swept across LData instead of ClData
  	global = ClData - Mg,
  	none = ClData)

  #######################
  ### Standardization ###
  #######################
  
  Vp <- colSums(sweep(CData^2, 1, Wn, "*"))
  Vn <- rowSums(sweep(CData^2, 2, Wp, "*"))
  Vg <- sum(Vp * Wp)
  
  SData <- switch(normal,
  	  global = CData / sqrt(Vg),
  	  row = sweep(CData, 1, sqrt(Vn), "/"),
  	  column = sweep(CData, 2, sqrt(Vp), "/"), # RV: originally the code swept across LData instead of ClData
  	  none = CData)
  
  #####################
  ### Factorization ###
  #####################
  
  WData <- sweep(sweep(SData, 1, sqrt(Wn), "*"), 2, sqrt(Wp), "*") #  weighted data matrix
  svd.res <- La.svd(WData)       # Singular Value Decomposition
  eigen <- svd.res$d^2           # Eigenvalues
  contrib <- eigen / sum(eigen)  # Contributions (normalized singular values)
  
  ### return
  r <- list(TData = SData,
          	row.names = Row.Names,
          	col.names = names(data)[-1],
          	closure = closure,
          	center = center,
          	normal = normal,
          	row.weight = row.weight,
          	col.weight = col.weight,
          	eigen = eigen,
          	contrib = contrib,
          	Rm = RM,
          	Cm = CM,
          	Wn = Wn,
          	Wp = Wp,
          	SVD = svd.res,
          	pos.column = pos.column,
          	pos.row = pos.row,
          	call = match.call())
          class(r) <- "mpm"
  return(r)
}

#print.plot.spm <- function(x) {invisible(x)}

