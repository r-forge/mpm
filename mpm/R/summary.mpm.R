summary.mpm <- function(object, maxdim = 4, ...){
  # scale <- match.arg(scale) ## TV
  if (!inherits(object, "mpm")) 
    stop("Use only with 'mpm' objects")
  svd.res <- object$SVD

  ### extract data from object
  Wn <- object$Wn     # row weights
  Wp <- object$Wp     # column weights
  Z <- object$TData   # raw data after standardisation
  La <- object$eigen  # eigen values
  d <- svd.res$d      # singular values
  U <- svd.res$u      # left singular vectors  (columns of U)
  V <- t(svd.res$vt)  # right singular vectors (columns of V)

  nf <- length(d) # number of extracted factors
  
  # eigenvalues and contributions
  Vxy <- sum(La)   # sum of all eigenvalues
  Las <- La / Vxy  # eigenvalues scaled to unit sum (proportions)
  VLa <- c(Las[1:min(maxdim, nf)], rep(0, max(0, maxdim - nf))) # proportions for retained factors only
  VQ <- if (maxdim < nf) sum(Las[(maxdim + 1):nf]) else 0  # sum of proportions for all non-retained factors
  VPf <- rbind(c(VLa,         VQ, sum(VLa, VQ)), # summary presentation
               c(cumsum(VLa), VQ, sum(VLa, VQ))) # summary presentation (cumulative)
  
  # scores and loadings (projections)	
  S <- crossprod(t(as.matrix(sweep(Z, 2, sqrt(Wp), "*"))), V) / sqrt(Vxy)
  L <- crossprod(as.matrix(sweep(Z, 1, sqrt(Wn),"*")), U) / sqrt(Vxy)
  nS <- nrow(S)
  nL <- nrow(L)
  
  # residual scores and loadings not accounted for by first maxdim dimensions
  if (maxdim < nf) {
    SQr <- sqrt(rowSums(S[, (maxdim+1):nf, drop = FALSE]^2))
    SQc <- sqrt(rowSums(L[, (maxdim+1):nf, drop = FALSE]^2))
  } else {
    SQr <- rep(0, nS)
    SQc <- rep(0, nL)
  }
  
  # row and column norms
  SYr <- sqrt(rowSums(S^2))
  SYc <- sqrt(rowSums(L^2))
  
  # Contributions to global weighted sums of squares
  CPr <- (Wn * SYr^2) / sum(Wn * SYr^2)
  CPc <- (Wp * SYc^2) / sum(Wp * SYc^2)
  
  # first maxdim dominant factors
  S[, d < 1E-6] <- 0
  L[, d < 1E-6] <- 0
  Prf <- cbind(S[, 1:min(maxdim, nf)], 
      matrix(0, nrow = nS, ncol = max(0, maxdim - nf)))
  Pcf <- cbind(L[, 1:min(maxdim, nf)], 
      matrix(0, nrow = nL, ncol = max(0, maxdim - nf)))
  
  # Accuracies in maxdim space
  APr <- rowSums(Prf^2) / rowSums(S^2)
  APc <- rowSums(Pcf^2) / rowSums(L^2)
  
  # Collect row and column results
  Rows <- as.data.frame(cbind(object$pos.row, Wn, Prf, SQr, SYr, CPr, APr))
  dimnames(Rows) <- list(object$row.names, c("Posit", "Weight", paste("Prf", 1:maxdim, sep = ""),
                  "Resid", "Norm", "Contrib", "Accuracy"))
  
  Cols <- as.data.frame(cbind(object$pos.column, Wp, Pcf, SQc, SYc, CPc, APc))
  dimnames(VPf) <- list(c("Individual", "Cumulative"),
      c(paste("PRF", 1:maxdim, sep = ""), "Residual", "Total"))
  names(Cols) <- c("Posit", "Weight", paste("Pcf", 1:maxdim, sep = ""),
                   "Resid", "Norm", "Contrib", "Accuracy")
  r <- list(call = object$call, Vxy = Vxy, VPF = VPf, 
            Rows = Rows, Columns = Cols)
  class(r) <- "summary.mpm"
  return(r)
}

