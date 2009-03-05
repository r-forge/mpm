print.summary.mpm <- function(x, digits = 2, 
    what = c("columns", "rows", "all"), ...){
  if (missing(x)) stop("Argument \"x\" is missing, with no default")
  what <- match.arg(what)
  cat("\nCall:\n")
  dput(x$call)
  cat("\n", length(x$Columns$Posit), " columns and ", length(x$Row$Posit), " rows.\n")
  cat("\n", sum(x$Columns$Posit), " columns and ", sum(x$Rows$Posit), " rows positioned.\n")
  cat("\nContributions:\n")
  cx <- format(round(x$VPF,digits), digits = digits)
  print(cx, quote = FALSE,...)
  if (what == "columns" || what == "all"){
    cat("\nColumns:\n")
    cx <- format(round(x$Columns, digits), digits = digits)
    print(cx, quote = FALSE, ...)
  }
  if (what == "rows" || what == "all"){
    cat("\nRows:\n")
    cx <- format(round(x$Rows, digits), digits = digits)
    print(cx, quote = FALSE, ...)
  }
  invisible(x)
}
