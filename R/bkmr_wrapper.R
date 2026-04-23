#wrapper for kmbayes
#' Fit Bayesian kernel machine regression
#'
#' Fits the Bayesian kernel machine regression (BKMR) model with group-separable kernel using Markov chain Monte Carlo (MCMC) methods.
#'
#' @export
#'
#' @param formula a formula
#' @param modifier a vector categorical values of length \code{n} that may modify the exposure-response associations. Control level ordering by using class \code{factor}, or use default level orders with any vector class.
#' @param iter number of iterations to run the sampler
#' @param family a description of the error distribution and link function to be used in the model. Currently implemented for \code{gaussian} and \code{binomial} families.
#' @param id optional vector (of length \code{n}) of grouping factors for fitting a model with a random intercept. If NULL then no random intercept will be included.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate diagnostic information during the model fitting.
#' @param starting.values list of starting values for each parameter. If not specified default values will be chosen.
#' @param control.params list of parameters specifying the prior distributions and tuning parameters for the MCMC algorithm. If not specified default values will be chosen.
#' @param varsel TRUE or FALSE: indicator for whether to conduct variable selection on the Z variables in \code{h}
#' @param groups optional vector (of length \code{M}) of group indicators for fitting hierarchical variable selection if varsel=TRUE. If varsel=TRUE without group specification, component-wise variable selections will be performed.
#' @param knots optional matrix of knot locations for implementing the Gaussian predictive process of Banerjee et al. (2008). Currently only implemented for models without a random intercept.
#' @param ztest optional vector indicating on which variables in Z to conduct variable selection (the remaining variables will be forced into the model).
#' @param rmethod for those predictors being forced into the \code{h} function, the method for sampling the \code{r[m]} values. Takes the value of 'varying' to allow separate \code{r[m]} for each predictor; 'equal' to force the same \code{r[m]} for each predictor; or 'fixed' to fix the \code{r[m]} to their starting values
#' @param est.h TRUE or FALSE: indicator for whether to sample from the posterior distribution of the subject-specific effects h_i within the main sampler. This will slow down the model fitting.
#' @param kernel.method When \code{modifier = "TRUE"}, use \code{kernel.method = "one"} for the standard approach or \code{kernel.method = "two"} for a group-separable approach with a categorical modifier
#' @param gs.tau TRUE or FALSE: indicator for whether to use group-specific tau parameters, only available for the group-separable method (\code{kernel.method = "two"})
#' @param gs.sig TRUE or FALSE: indicator for whether to estimate separate error variance terms for each group. Only available for the group-separable method (\code{kernel.method = "two"})
#' @return an object of class "bkmrfit" (containing the posterior samples from the model fit), which has the associated methods:
#' \itemize{
#'   \item \code{\link{print}} (i.e., \code{\link{print.bkmrfit}}) 
#'   \item \code{\link{summary}} (i.e., \code{\link{summary.bkmrfit}})
#' }
#' 
#' @seealso For guided examples, see vignette(bkmrGSOverview)
#' @references Bobb, JF, Valeri L, Claus Henn B, Christiani DC, Wright RO, Mazumdar M, Godleski JJ, Coull BA (2015). Bayesian Kernel Machine Regression for Estimating the Health Effects of Multi-Pollutant Mixtures. Biostatistics 16, no. 3: 493-508.
#' @references Banerjee S, Gelfand AE, Finley AO, Sang H (2008). Gaussian predictive process models for large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(4), 825-848.
#' @import utils
#' 
#' @examples
#' ## First generate data set
#' y <- ex_data$y
#' Z <- ex_data$Z
#' modifier <- ex_data$X$Sex
#' X_full <- ex_data$X[,-2] #remove Sex from the covariate matrix because it is the modifier
#' #create design matrix to account for factor variables, remove the intercept column
#' X <- model.matrix(~., data=X_full)[,-1] 
#' 
#' ## Fit model 
#' ## Using only 10 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, modifier = modifier, X = X, iter = 10, verbose = FALSE) 
kmbayesWrapper <- function(formula, data,
                           iter = 1000, family = "gaussian", id = NULL, verbose = TRUE, 
                           starting.values = NULL, control.params = NULL, varsel = FALSE, 
                           groups = NULL, knots = NULL, ztest = NULL, rmethod = "varying",
                           est.h = FALSE, kernel.method = "two",
                           gs.tau = TRUE, gs.sig = FALSE) {
  
  #seperate out the formula? 
  #checks
  # So this all works, now my goal is to be able to from the data frame
  # make sure there is only 1 h() 
  # pull out the response OR reject if no response is found in the data frame
  # pull out the Z and modifier (mod) from either h or kernel OR reject if no 
  # Z or modifier (mod) is in there, OR if there are more than 1 h()
  # pull out the X OR reject if no X is found
  # ensure that the data frame is a data frame or matrix
  # ensure that the size checks pass
  
  #seperate formula
  formula <- as.formula(formula)
  terms_obj <- terms(formula)
  term_labels <- attr(terms_obj, "term.labels")
  
  #extract response
  response <- all.vars(formula)[attr(terms_obj, "response")]
  
  #check response
  if(length(response) != 1){
    stop("Response needs to be only one column")
  }
  
  if(!(response %in% names(data))){
    stop("Response " + as.character(mod) + " is not in the data frame")
  }
  
  # Extract h()
  h_term <- term_labels[grepl("^h\\(", term_labels)]
  
  #check h()
  if (length(h_term) > 1) {
    stop("Exactly one h() term required.")
  }
  
  if (length(h_term) < 1){
    stop("Need to declare exactly one h().")
  }
  
  #parse h()  
  h_call <- str2lang(h_term)
  args <- as.list(h_call)[-1]
  arg_names <- names(args)
  
  #seperate Z
  Z_vars <- args[arg_names == "" | is.null(arg_names)]
  Z_names <- sapply(Z_vars, deparse)
  
  #check Z
  if (length(Z_names) == 0) {
    stop("No exposure variables supplied inside h().")
  }
  
  missing_Z <- setdiff(Z_names, names(data))
  
  if (length(missing_Z) > 0) {
    stop(
      sprintf(
        "Variables in kernel() not found in data: %s",
        paste(missing_Z, collapse = ", ")
      )
    )
  }
  
  #seperate mod/modifier
  mod_name <- "" 
  if("mod" %in% arg_names){
    mod_name <- deparse(args[["mod"]])
  } 
  else if ("modifier" %in% arg_names){
    mod_name <- deparse(args[["modifier"]])
  }
  else NULL
  
  #check modifer 
  if(length(mod_name) != 1){
    stop("Modifier can only be one column")
  }
  
  if (!is.null(mod_name)) {
    if (!(mod_name %in% names(data))) {
      stop(paste0("Modifier '", mod_name, "' not found in data."))
    }
  }
  
  #get X 
  X_terms <- term_labels[!grepl("^h\\(", term_labels)]
  
  #check X
  if (length(X_terms) == 0) {
    stop("No covaraites supplied inside h().")
  }
  
  missing_X <- setdiff(X_terms, names(data))
  
  if (length(missing_X) > 0) {
    stop(
      sprintf(
        "Covariates not found in data: %s",
        paste(missing_X, collapse = ", ")
      )
    )
  }
  
  #passed all checks 
  
  y <- data[[response]]
  X <- data[c(X_terms, mod_name)]
  Z <- data[Z_names]
  modifier <- mod_name
  
  return(kmbayesBlocked(y, Z, X, modifier, 
                        iter, family, id, verbose, starting.values, control.params, varsel, 
                        groups, knots, ztest, rmethod,
                        est.h, kernel.method,
                        gs.tau, gs.sig))
    

    
}
