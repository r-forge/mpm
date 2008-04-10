print.mpm <- function(x, digits = 3, ...){
  cat("Call:\n")
  dput(x$call)
  cat("\nContributions:\n")
  print(round(x$contrib, digits), ...)
  cat("\n", length(x$pos.column), " columns and ", length(x$pos.row), " rows.\n")
  cat("\n", sum(x$pos.column), " columns and ", sum(x$pos.row), " rows positioned.\n")
  invisible(x)
}
