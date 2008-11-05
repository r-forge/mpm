.onAttach <- function(libname, pkgname){
  message(paste("\nmpm version ", packageDescription("mpm")$Version, 
          "\n", sep = ""))
}

