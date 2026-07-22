.onAttach <- function(libname, pkgname) {
  
  version <- utils::packageVersion(pkgname)
  
  packageStartupMessage(paste0(
    "bkmrGS version ", version, "\n",
    "Bayesian Kernel Machine Regression with Group Seperable Kernel Functionality\n",
    "For guided examples, see vignettes('bkmrGS_Vignette'). This package is based on the bkmr package version 0.2.2."
  ))
}

release_questions <- function() {
  c(
    "Have you updated the vignette and posted to GitHub?"
  )
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables("Liu_data")
}