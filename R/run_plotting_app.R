#' Launch the group-separable BKMR Shiny App
#'
#' @export
plotting_app <- function() {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2'.")
  }
  
  appDir <- system.file("shiny", "PlottingApp.R", package = "bkmrGS")
  
  if (appDir == "") {
    stop("Could not find app directory.")
  }
  
  shiny::runApp(appDir)
}