#' Launch the BKMR Shiny App
#'
#' @export
plotting_app <- function() {
  
  appDir <- system.file("shiny", "PlottingApp.R", package = "bkmrGS")
  
  if (appDir == "") {
    stop("Could not find app directory.")
  }
  
  shiny::runApp(appDir)
}