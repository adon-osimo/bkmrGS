#' Launch the BKMR Shiny App
#'
#' @export
run_plotting_app <- function() {
  
  appDir <- system.file("shiny", "App4.R", package = "bkmrGS")
  
  if (appDir == "") {
    stop("Could not find app directory.")
  }
  
  shiny::runApp(appDir)
}