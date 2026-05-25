#' Trace plot
#' 
#' Fits a trace plot of the given parameter, along with the posterior mean estimate 
#'
#' @inheritParams ExtractEsts
#' @param par which parameter to plot
#' @param comp which component of the parameter vector to plot
#' @param main title
#' @param xlab x axis label
#' @param ylab y axis label
#' @param ... other arguments to pass onto the plotting function
#' @export
#' @import graphics
#' @details For guided examples, see vignette(VignettteDraft)
#' 
#' @return No return value, generates plot
#' 
#' @examples
#' set.seed(111) 
#' fitkm <- bkmrGS(y ~ h(Log_Lead, Log_Manganese, Log_Arsenic, mod = Sex) +
#'                           Age + Gestation + Delivery + Birth_order +
#'                           Education_parent1 + Education_parent2 + Smoking +
#'                           HOME_emotional + HOME_avoid + HOME_careg + 
#'                           HOME_env + HOME_play + HOME_stim + 
#'                           Energy,
#'                         data = Liu_data, 
#'                         iter = 10, 
#'                         verbose = FALSE) 
#'                         
#' TracePlot(fit = fitkm, par = "beta")
#' TracePlot(fit = fitkm, par = "lambda")
#' TracePlot(fit = fitkm, par = "sigsq.eps")
#' TracePlot(fit = fitkm, par = "r", comp = 1)
TracePlot <- function(fit, par, comp = 1, sel = NULL, main = "", xlab = "iteration", ylab = NULL, ...) {
    samps <- ExtractSamps(fit, sel = sel)[[par]]
    if (!is.null(ncol(samps))) {
        nm <- colnames(samps)[comp]
        samps <- samps[, comp]
    } else {
        nm <- par
    }
    
    if(is.null(ylab)){
      ylab <- paste0(as.character(nm), " value")
    }
    
    main <- paste0(main, "\nTrace Plot for ", nm)
    plot(samps, type = "l", main = main, xlab = xlab,  ylab = ylab, ...)
    abline(h = mean(samps), col = "#D9782D", lwd = 2) #originally blue but I made it green, then orange, go CSU
}
