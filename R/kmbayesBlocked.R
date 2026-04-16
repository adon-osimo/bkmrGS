get_block_ids <- function(modifier){
  groups <- apply(modifier, 1, paste, collapse = "")
  split(seq_along(groups), groups)
}

makeVcompsBlocked <- function(r, lambda, Z, data.comps, modifier = NULL) {
  #Check if data.comps is null
  if (is.null(data.comps$knots)) {
    #Check if we have a modifier
    #if so, use the sped up version levering block structure
    if(!is.null(modifier)){
      block_ids <- get_block_ids(modifier)
      G <- length(block_ids)
      Z_list <- lapply(block_ids, function(id) Z[id, ,drop = FALSE])
    }
    
    else{
      block_ids <- list(1:nrow(Z))
      G <- 1
      Z_list <- list(Z)
    }
    
    #If modifier = null, run with standard
    scalar_lambda <- lambda_scalar(mod1 = modifier, mod2 = modifier, lambda = lambda, data.comps = data.comps)
    scalar_lambda_list <- lambda_scalar_blocked(modifier, lambda, data.comps, block_ids)
    
    #for model with random intercept
    #if (data.comps$randint) { #changed from data.comps$nlambda == 2 upon GS-tau model implementation
    #  V <- V + lambda[length(lambda)-1]*data.comps$crossTT
    #}
    #not sure what this does really 
    Vinv_list <- vector("list", G)
    logdetVinv_list <- vector("list", G)
    logdetVinv <- 0
    
    Vinv = matrix(0, nrow = nrow(Z), ncol = nrow(Z))
    
    for(g in 1:G){
      
      id <- block_ids[[g]]
      Z_g <- Z_list[[g]]
      
      Kpart_g <- makeKpart(r, Z_g)
      K_g <- exp(-Kpart_g)
      scalar_lambda_g <- scalar_lambda_list[[g]]
      
      lambda_val <- scalar_lambda_g[1,1]
      V_g <- diag(1, length(id)) + lambda_val * K_g
      
      cholV_g <- chol(V_g) #recode for speed in Rcpp?
      
      Vinv_list[[g]] <- chol2inv(cholV_g) #recode for speed in Rcpp?
      
      logdetVinv_list[[g]] <- - 2* sum(log(diag(cholV_g)))
      logdetVinv <- logdetVinv + logdetVinv_list[[g]]
      
      Vinv[id,id] <- Vinv_list[[g]]
    }
    
    Vcomps <- list(Vinv = Vinv, Vinv_list = Vinv_list, logdetVinv = logdetVinv, scalar_lambda = scalar_lambda, logdetVinv_list = logdetVinv_list)
  } else {## predictive process approach
    ## note: currently does not work with random intercept model
    nugget <- 0.001
    n0 <- nrow(Z)
    n1 <- nrow(data.comps$knots)
    nall <- n0 + n1
    # Kpartall <- makeKpart(r, rbind(Z, data.comps$knots))
    # Kall <- exp(-Kpartall)
    # K0 <- Kall[1:n0, 1:n0 ,drop=FALSE]
    # K1 <- Kall[(n0+1):nall, (n0+1):nall ,drop=FALSE]
    # K10 <- Kall[(n0+1):nall, 1:n0 ,drop=FALSE]
    K1 <- exp(-makeKpart(r, data.comps$knots))
    if(!is.null(modifier)){
      K1 <- block_kernel(mod_vec1 = modifier,
                         mod_vec2 = modifier,
                         K = K1)
    }
    K10 <- exp(-makeKpart(r, data.comps$knots, Z))
    if(!is.null(modifier)){
      K10 <- block_kernel(mod_vec1 = modifier, 
                          mod_vec2 = modifier,
                          K = K10)
    }
    Q <- K1 + diag(nugget, n1, n1)
    scalar_lambda <- lambda_scalar(mod1 = modifier, 
                                   mod2 = modifier, 
                                   lambda = lambda,
                                   data.comps = data.comps)
    R <- Q + scalar_lambda*tcrossprod(K10)
    cholQ <- chol(Q)
    cholR <- chol(R)
    Qinv <- chol2inv(cholQ)
    Rinv <- chol2inv(cholR)
    Vinv <- diag(1, n0, n0) - scalar_lambda*t(K10) %*% Rinv %*% K10
    logdetVinv <- 2*sum(log(diag(cholQ))) - 2*sum(log(diag(cholR)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv, cholR = cholR, Q = Q, K10 = K10, Qinv = Qinv, Rinv = Rinv)
  }
  Vcomps
}

lambda_scalar_blocked <- function(modifier, lambda, data.comps, block_ids) {
  scalar_list <- vector("list", length(block_ids))
  
  if(data.comps$gs.tau){
    if(data.comps$randint){
      lambda <- lambda[1:(length(lambda)-1)]
    }
    
    mod_levels <- apply(modifier, 1, paste, collapse = "")
    
    for(g in seq_along(block_ids)) {
      ids <- block_ids[[g]]
      lambda_block <- lambda[mod_levels[ids]]
      
      n <- length(ids)
      scalar_list[[g]] <- matrix(rep(lambda_block, each = n), nrow = n)
    }
    
  } else {
    for(g in seq_along(block_ids)) {
      n <- length(block_ids[[g]])
      scalar_list[[g]] <- matrix(lambda[1], n, n)
    }
  }
  
  return(scalar_list)
}

#' Fit Bayesian kernel machine regression
#'
#' Fits the Bayesian kernel machine regression (BKMR) model with group-separable kernel using Markov chain Monte Carlo (MCMC) methods.
#'
#' @export
#'
#' @param y a vector of outcome data of length \code{n}.
#' @param Z an \code{n}-by-\code{M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n}-by-\code{K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
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
#' ## First, seperate the data set
#' y <- ex_data$y
#' Z <- ex_data$Z
#' X <- ex_data$X
#' 
#' ## Fit model 
#' ## Using only 10 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayesBlocked(y = y, Z = Z, modifier = "Sex", X = X, iter = 10, verbose = FALSE) 
kmbayesBlocked <- function(y, Z, X = NULL, 
                    modifier = NULL, #added by DD, 
                    iter = 1000, family = "gaussian", id = NULL, verbose = TRUE, 
                    starting.values = NULL, control.params = NULL, varsel = FALSE, 
                    groups = NULL, knots = NULL, ztest = NULL, rmethod = "varying",
                    est.h = FALSE, kernel.method = "two",
                    gs.tau = TRUE, gs.sig = FALSE) {
  #browser()
  missingX <- is.null(X)
  if (missingX) X <- matrix(0, length(y), 1)
  hier_varsel <- !is.null(groups)
  
  #I don't think this line is necessary as if modifier is null
  #we already run standard bkmr, just slightly slower
  #if(is.null(modifier)){
  #  bkmr::kmbayes(y, Z, X, family, id, verbose, 
  #                starting.values, control.params, varsel, 
  #                groups, knots, ztest, rmethod,
  #                est.h , kernel.method,
  #                gs.tau, gs.sig)
  #}
  
  #convert X to model matrix and split out modifier
  if(!is.null(modifier)){
    
    if(is.character(modifier)){

      col <- names(X) %in% modifier
    
      if(!any(col)){
        stop("Modifier name not found in X")
      }
      
      if(sum(col) != 1){
        stop("Modifier must correspond to exactly one column")
      }
      orig_modifier <- X[, col]
      modifier <- as.factor(orig_modifier)
      
      X_full <- X[, !col, drop = FALSE]
    }
    
    else if(is.double(modifier)){
      
      if(length(modifier) != 1 || modifier < 1 || modifier > ncol(X)){
        warning("Modifier is not in matrix")
      }
      
      col = modifier
      orig_modifier <- X[, modifier]
      modifier <- as.factor(orig_modifier)
       
      X_full <- X[, -orig_modifier, drop = FALSE]
        
    }
    else{
      #throw an exception? 
      warning("Modifier is neither a name nor a value in the matrix")
    }
    
    modifier <- as.matrix(model.matrix(~modifier)[,-1])
    X <- as.matrix(model.matrix(~., data = X_full)[,-1])
  }
  
  #if(!is.null(modifier)){
  #  orig_modifier <- modifier
  #  modifier <- as.factor(modifier)
  #  modifier <- as.matrix(model.matrix(~modifier)[,-1])
  #}
  
  #convert modifier to factor and construct contrast matrix
  
  #make sure kernel.method is provided properly
  if(kernel.method == "one"){
    kern_modifier <- NULL
  }else if (kernel.method == "two"){
    if(is.null(modifier)){
      warning("Group-separable aproach is only for BKMR with effect modification. Setting kernel.method to 'one'.")
      kernel.method <- "one"
    }else{
      kern_modifier <- modifier
    }
  }else{
    stop("Invalid kernel.method. Only 'one' and 'two' are supported.")
  }
  
  #make sure that gs.tau is FALSE when kernel.method is one
  if(kernel.method == "one" & gs.tau){
    stop("Group-specific tau method only available for group-separable model. Change gs.tau to FALSE or change kernel.method to 'two'.")
  }
  
  #make sure that gs.sig is FALSE when kernel.method is one
  if(kernel.method == "one" & gs.sig){
    stop("Group-specific residual variance method only available for group-separable model. Change gs.sig to FALSE or change kernel.method to 'two'.")
  }
  
  ##Argument check 1, required arguments without defaults
  ##check vector/matrix sizes
  stopifnot (length(y) > 0, is.numeric(y), anyNA(y) == FALSE)
  if(!is.null(modifier)){
    stopifnot (length(modifier) > 0, is.numeric(modifier), anyNA(modifier) == FALSE)
  }
  if (!inherits(Z, "matrix"))  Z <- as.matrix(Z)
  stopifnot (is.numeric(Z), nrow(Z) == length(y), anyNA(Z) == FALSE)
  if (!inherits(X, "matrix"))  X <- as.matrix(X)
  stopifnot (is.numeric(X), nrow(X) == length(y), anyNA(X) == FALSE) 
  
  ##Argument check 2: for those with defaults, write message and reset to default if invalid
  if (iter < 1) {
    message ("invalid input for iter, resetting to default value 1000")
    nsamp <- 1000
  } else {
    nsamp <- iter
  }
  if (!family %in% c("gaussian", "binomial")) {
    stop("family", family, "not yet implemented; must specify either 'gaussian' or 'binomial'")
  }
  if (family == "binomial") {
    message("Fitting probit regression model")
    if (!all(y %in% c(0, 1))) {
      stop("When family == 'binomial', y must be a vector containing only zeros and ones")
    }
  }
  if (rmethod != "varying" & rmethod != "equal" & rmethod != "fixed") {
    message ("invalid value for rmethod, resetting to default varying")
    rmethod <- "varying"
  }
  if (verbose != FALSE & verbose != TRUE) {
    message ("invalid value for verbose, resetting to default FALSE")
    verbose <- FALSE
  }
  if (varsel != FALSE & varsel != TRUE) {
    message ("invalid value for varsel, resetting to default FALSE")
    varsel <- FALSE
  }
  
  ##Argument check 3: the rest id (below) znew, knots, groups, ztest
  if (!is.null(id)) { 
    stopifnot(length(id) == length(y), anyNA(id) == FALSE)
    if (!is.null(knots)) { 
      message ("knots cannot be specified with id, resetting knots to null")
      knots<-NA
    }
  }
  # if (!is.null(Znew)) { 
  #   if (!inherits(Znew, "matrix"))  Znew <- as.matrix(Znew)
  #   stopifnot(is.numeric(Znew), ncol(Znew) == ncol(Z), anyNA(Znew) == FALSE)
  # }
  if (!is.null(knots)) { 
    if (!inherits(knots, "matrix"))  knots <- as.matrix(knots)
    stopifnot(is.numeric(knots), ncol(knots )== ncol(Z), anyNA(knots) == FALSE)
  }
  if (!is.null(groups)) { 
    if (varsel == FALSE) {
      message ("groups should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(groups), length(groups) == ncol(Z), anyNA(groups) == FALSE)
    }
  }
  if (!is.null(ztest)) { 
    if (varsel == FALSE) {
      message ("ztest should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(ztest), length(ztest) <= ncol(Z), anyNA(ztest) == FALSE, max(ztest) <= ncol(Z) )
    }
  }
  
  if(is.null(colnames(Z))){
    colnames(Z) <- paste0("z",1:ncol(Z))
  }
  
  
  #add the modifier, if given, to the exposure matrix and covariate matrix
  if(!is.null(modifier)){
    
    if(kernel.method == "one"){ #add the modifier into the kernel for the one-kernel interaction approach
      Z <- cbind(Z, modifier) 
    }else{ #do not add the modifier in the kernel for the two-kernel approach
      Z <- Z
    }
    if(is.null(colnames(X))){#always add the modifier into the covaraites when modifier is provided
      colnames(X) <- paste0("X",1:ncol(X))
    }
    X <- cbind(X, modifier) 
  }
  
  ## start JB code
  if (!is.null(id)) { ## for random intercept model
    randint <- TRUE
    id <- as.numeric(as.factor(id))
    nid <- length(unique(id))
    if(gs.tau){
      nlambda <- length(levels(as.factor(orig_modifier))) + 1
    }else{
      nlambda <- 2
    }
    
    ## matrix that multiplies the random intercept vector
    TT <- matrix(0, length(id), nid)
    for (i in 1:nid) {
      TT[which(id == i), i] <- 1
    }
    crossTT <- tcrossprod(TT)
    rm(TT, nid)
  } else {
    randint <- FALSE
    if(gs.tau){
      nlambda <- length(levels(as.factor(orig_modifier)))
    }else{
      nlambda <- 1
    }
    crossTT <- 0
  }
  
  #set levels of the modifier
  if(is.null(modifier)){
    lvls <- NULL
  }else{
    lvls <- unique(apply(modifier, 1, paste, collapse = ""))
  }
  
  #set the number of residual variance terms
  if(gs.sig){
    n.sigsq.eps <- length(lvls)
  }else{
    n.sigsq.eps <- 1
  }
  
  data.comps <- list(randint = randint, nlambda = nlambda, crossTT = crossTT, gs.tau = gs.tau, levels = lvls, gs.sig = gs.sig)
  if (!is.null(knots)) data.comps$knots <- knots
  rm(randint, nlambda, crossTT)
  
  ## create empty matrices to store the posterior draws in
  chain <- list(h.hat = matrix(0, nsamp, nrow(Z)),
                beta = matrix(0, nsamp, ncol(X)),
                lambda = matrix(NA, nsamp, data.comps$nlambda),
                sigsq.eps = matrix(NA, nsamp, n.sigsq.eps),
                r = matrix(NA, nsamp, ncol(Z)),
                acc.r = matrix(0, nsamp, ncol(Z)),
                acc.lambda = matrix(0, nsamp, data.comps$nlambda),
                delta = matrix(1, nsamp, ncol(Z))
  )
  if(gs.tau){
    if(data.comps$randint){
      colnames(chain$lambda) <- c(lvls, "int")
    }else{
      colnames(chain$lambda) <- lvls
    }
  }
  if(gs.sig){
    colnames(chain$lambda) <- lvls
  }
  
  if (varsel) {
    chain$acc.rdelta <- rep(0, nsamp)
    chain$move.type <- rep(0, nsamp)
  }
  if (family == "binomial") {
    chain$ystar <- matrix(0, nsamp, length(y))
  }
  
  # ## components to predict h(Znew)
  # if (!is.null(Znew)) {
  #   if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
  #   if (inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)
  #   if (ncol(Z) != ncol(Znew)) {
  #     stop("Znew must have the same number of columns as Z")
  #   }
  #   ##Kpartall <- as.matrix(dist(rbind(Z,Znew)))^2
  #   chain$hnew <- matrix(0,nsamp,nrow(Znew))
  #   colnames(chain$hnew) <- rownames(Znew)
  # }
  
  modtest <- FALSE #modtest used to be a test on the modifier for the modifier-in-kernel model. removing since we do not recommend this method. easiest way to remove is to set to false
  
  if (varsel) {
    if (is.null(ztest) & !modtest) {
      ztest <- 1:(ncol(Z) - !is.null(modifier))
    }
    if(modtest){
      ztest <- c(ztest, (ncol(Z) - !is.null(modifier))+1)
    }
    rdelta.update <- rdelta.comp.update
  } else {
    ztest <- NULL
  }
  
  ## control parameters
  control.params.default <- list(lambda.jump = rep(10, data.comps$nlambda), mu.lambda = rep(10, data.comps$nlambda), sigma.lambda = rep(10, data.comps$nlambda), a.p0 = 1, b.p0 = 1, r.prior = "invunif", a.sigsq = 1e-3, b.sigsq = 1e-3, mu.r = 5, sigma.r = 5, r.muprop = 1, r.jump = 0.1, r.jump1 = 2, r.jump2 = 0.1, r.a = 0, r.b = 100)
  if (!is.null(control.params)){
    control.params <- modifyList(control.params.default, as.list(control.params))
    validateControlParams(varsel, family, id, control.params)
  } else {
    control.params <- control.params.default
  }
  
  control.params$r.params <- with(control.params, list(mu.r = mu.r, sigma.r = sigma.r, r.muprop = r.muprop, r.jump = r.jump, r.jump1 = r.jump1, r.jump2 = r.jump2, r.a = r.a, r.b = r.b))
  
  ## components if grouped model selection is being done
  if (!is.null(groups)) {
    if (!varsel) {
      stop("if doing grouped variable selection, must set varsel = TRUE")
    }
    rdelta.update <- rdelta.group.update
    control.params$group.params <- list(groups = groups, sel.groups = sapply(unique(groups), function(x) min(seq_along(groups)[groups == x])), neach.group = sapply(unique(groups), function(x) sum(groups %in% x)))
  }
  
  ## specify functions for doing the Metropolis-Hastings steps to update r
  e <- environment()
  rfn <- set.r.MH.functions(r.prior = control.params$r.prior)
  rprior.logdens <- rfn$rprior.logdens
  environment(rprior.logdens) <- e
  rprop.gen1 <- rfn$rprop.gen1
  environment(rprop.gen1) <- e
  rprop.logdens1 <- rfn$rprop.logdens1
  environment(rprop.logdens1) <- e
  rprop.gen2 <- rfn$rprop.gen2
  environment(rprop.gen2) <- e
  rprop.logdens2 <- rfn$rprop.logdens2
  environment(rprop.logdens2) <- e
  rprop.gen <- rfn$rprop.gen
  environment(rprop.gen) <- e
  rprop.logdens <- rfn$rprop.logdens
  environment(rprop.logdens) <- e
  rm(e, rfn)
  
  ## initial values
  starting.values0 <- list(h.hat = 1, beta = NULL, sigsq.eps = NULL, r = 1, lambda = 10, delta = 1)
  if (is.null(starting.values)) {
    starting.values <- starting.values0
  } else {
    starting.values <- modifyList(starting.values0, starting.values)
    validateStartingValues(varsel, y, X, Z, starting.values, rmethod)
  }
  if (family == "gaussian") {
    if (is.null(starting.values$beta) | is.null(starting.values$sigsq.eps)) {
      lmfit0 <- lm(y ~ Z + X)
      if (is.null(starting.values$beta)) {
        coefX <- coef(lmfit0)[grep("X", names(coef(lmfit0)))]
        starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
      }
      if (is.null(starting.values$sigsq.eps)) {
        starting.values$sigsq.eps <- summary(lmfit0)$sigma^2
      }
    } 
  } else if (family == "binomial") {
    starting.values$sigsq.eps <- 1 ## always equal to 1
    if (is.null(starting.values$beta) | is.null(starting.values$ystar)) {
      probitfit0 <- try(glm(y ~ Z + X, family = binomial(link = "probit")))
      if (!inherits(probitfit0, "try-error")) {
        if (is.null(starting.values$beta)) {
          coefX <- coef(probitfit0)[grep("X", names(coef(probitfit0)))]
          starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
        }
        if (is.null(starting.values$ystar)) {
          #prd <- predict(probitfit0)
          #starting.values$ystar <- ifelse(y == 1, abs(prd), -abs(prd))
          starting.values$ystar <- ifelse(y == 1, 1/2, -1/2)
        }
      } else {
        starting.values$beta <- 0
        starting.values$ystar <- ifelse(y == 1, 1/2, -1/2)
      }
    } 
  }
  
  ##print (starting.values)
  ##truncate vectors that are too long
  if (length(starting.values$h.hat) > length(y)) {
    starting.values$h.hat <- starting.values$h.hat[1:length(y)]
  }
  if (length(starting.values$beta) > ncol(X)) {
    starting.values$beta <- starting.values$beta[1:ncol(X)]
  }
  if (length(starting.values$delta) > ncol(Z)) {
    starting.values$delta <- starting.values$delta[1:ncol(Z)]
  }
  if (varsel==FALSE & rmethod == "equal" & length(starting.values$r) > 1) {
    starting.values$r <- starting.values$r[1] ## this should only happen if rmethod == "equal"
  } else if (length(starting.values$r) > ncol(Z)) {
    starting.values$r <- starting.values$r[1:ncol(Z)]
  }
  
  chain$h.hat[1, ] <- starting.values$h.hat
  chain$beta[1, ] <- starting.values$beta
  chain$lambda[1, ] <- starting.values$lambda
  chain$sigsq.eps[1, ] <- rep(starting.values$sigsq.eps, n.sigsq.eps)
  chain$r[1, ] <- starting.values$r
  if (varsel) {
    chain$delta[1,ztest] <- starting.values$delta
  }
  if (family == "binomial") {
    chain$ystar[1, ] <- starting.values$ystar
    chain$sigsq.eps[] <- starting.values$sigsq.eps ## does not get updated
  }
  if (!is.null(groups)) {
    ## make sure starting values are consistent with structure of model
    if (!all(sapply(unique(groups), function(x) sum(chain$delta[1, ztest][groups == x])) == 1)) {
      # warning("Specified starting values for delta not consistent with model; using default")
      starting.values$delta <- rep(0, length(groups))
      starting.values$delta[sapply(unique(groups), function(x) min(which(groups == x)))] <- 1
    }
    chain$delta[1,ztest] <- starting.values$delta
    chain$r[1,ztest] <- ifelse(chain$delta[1,ztest] == 1, chain$r[1,ztest], 0)
  }
  chain$est.h <- est.h
  
  ## components
  Vcomps <- makeVcompsBlocked(r = chain$r[1, ], lambda = chain$lambda[1, ], Z = Z, data.comps = data.comps, modifier = kern_modifier)
  
  block_ids <- get_block_ids(kern_modifier)
  
  X_list <- lapply(block_ids, function(id) X[id, , drop = FALSE])
  y_list <- lapply(block_ids, function(id) y[id])
  Z_list <- lapply(block_ids, function(id) Z[id, , drop = FALSE])
  
  
  ## start sampling ####
  chain$time1 <- Sys.time()
  for (s in 2:nsamp) {
    
    ## continuous version of outcome (latent outcome under binomial probit model)
    if (family == "gaussian") {
      ycont <- y
    } else if (family == "binomial") {
      if (est.h) {
        chain$ystar[s,] <- ystar.update(y = y, X = X, beta = chain$beta[s - 1,], h = chain$h[s - 1, ])
      } else {
        chain$ystar[s,] <- ystar.update.noh(y = y, X = X, beta = chain$beta[s - 1,], Vinv = Vcomps$Vinv, ystar = chain$ystar[s - 1, ])
      }
      ycont <- chain$ystar[s, ]
    }
    
    ## generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)
    
    ## beta
    if (!missingX) {
      #gs.sig method: make sigsq.eps a vector (not matrix) of length n
      if(gs.sig){
        sigs <- chain$sigsq.eps[s - 1, ]
        names(sigs) <- lvls
        sigsq.eps.betaup <- sigs[apply(modifier, 1, paste, collapse = "")]
      }else{
        sigsq.eps.betaup <- rep(chain$sigsq.eps[s - 1, ], nrow(X)) #if single sigsq.eps, make vector of the same value
      }
      chain$beta[s,] <- beta.update.block(X_list = X_list, Vinv_list = Vcomps$Vinv_list, y_list = y_list, sigsq.eps = sigsq.eps.betaup)
    }
    
    ## \sigma_\epsilon^2
    if (family == "gaussian") {
      if(n.sigsq.eps==1){ #gs.sig == FALSE
        chain$sigsq.eps[s,1] <- sigsq.eps.update.block(y_list = y_list, X_list = X_list, beta = chain$beta[s,], Vinv_list = Vcomps$Vinv_list, a.eps = control.params$a.sigsq, b.eps = control.params$b.sigsq)
      }else{ #gs.sig method: call this only on the data for group p
        for(p in 1:n.sigsq.eps){
          grp <- lvls[p]
          grp_idx <- which(apply(modifier, 1, paste, collapse = "")==grp)
          grp_coords <- expand.grid(grp_idx,grp_idx)
          Vinv_grp <- matrix(Vcomps$Vinv[as.matrix(grp_coords)], ncol = length(grp_idx))
          chain$sigsq.eps[s,p] <- sigsq.eps.update(y = ycont[grp_idx], X = X[grp_idx,], beta = chain$beta[s,], Vinv = Vinv_grp, a.eps = control.params$a.sigsq, b.eps = control.params$b.sigsq)
        }
      }
    }
    
    ## lambda
    #gs.sig method: only pass data and sigsq.eps for each group
    lambdaSim <- chain$lambda[s - 1,]
    for (comp in 1:data.comps$nlambda) {
      varcomps <- lambda.update.block(r = chain$r[s - 1,], delta = chain$delta[s - 1,], lambda = lambdaSim, whichcomp = comp, y = y, X = X, Z = Z, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s, ], Vcomps = Vcomps, data.comps = data.comps, control.params = control.params, modifier = kern_modifier)
      lambdaSim <- varcomps$lambda
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.lambda[s,comp] <- varcomps$acc
      }
    }
    chain$lambda[s,] <- lambdaSim
    
    ## r
    rSim <- chain$r[s - 1,]
    comp <- which(!1:ncol(Z) %in% ztest)
    if (length(comp) != 0) {
      if (rmethod == "equal") { ## common r for those variables not being selected
        varcomps <- r.update.block(r = rSim, whichcomp = comp, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = y, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s, ], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens, modifier = kern_modifier)
        rSim <- varcomps$r
        if (varcomps$acc) {
          Vcomps <- varcomps$Vcomps
          chain$acc.r[s, comp] <- varcomps$acc
        }
      } else if (rmethod == "varying") { ## allow a different r_m
        for (whichr in comp) {
          varcomps <- r.update.block(r = rSim, whichcomp = whichr, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = y, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s, ], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens, modifier = kern_modifier)
          rSim <- varcomps$r
          if (varcomps$acc) {
            Vcomps <- varcomps$Vcomps
            chain$acc.r[s, whichr] <- varcomps$acc
          }
        }
      }
    }
    ## for those variables being selected: joint posterior of (r,delta)
    if (varsel) {
      varcomps <- rdelta.update.block(r = rSim, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = y, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s, ], Vcomps = Vcomps, Z = Z, ztest = ztest, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
      chain$delta[s,] <- varcomps$delta
      rSim <- varcomps$r
      chain$move.type[s] <- varcomps$move.type
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.rdelta[s] <- varcomps$acc
      }
    }
    chain$r[s,] <- rSim
    
    ###################################################
    ## generate posterior sample of h(z) from its posterior P(h | beta, sigsq.eps, lambda, r, y)
    
    if (est.h) {
      hcomps <- h.update.block(lambda = chain$lambda[s,], Vcomps = Vcomps, sigsq.eps = chain$sigsq.eps[s, ], y = ycont, X = X, beta = chain$beta[s,], r = chain$r[s,], Z = Z, data.comps = data.comps, modifier = modifier, kernel.method = kernel.method)
      chain$h.hat[s,] <- hcomps$hsamp
      if (!is.null(hcomps$hsamp.star)) { ## GPP
        Vcomps$hsamp.star <- hcomps$hsamp.star
      }
      rm(hcomps)
    }
    
    # ###################################################
    # ## generate posterior samples of h(Znew) from its posterior P(hnew | beta, sigsq.eps, lambda, r, y)
    # 
    # if (!is.null(Znew)) {
    #   chain$hnew[s,] <- newh.update(Z = Z, Znew = Znew, mod_new = mod_new, Vcomps = Vcomps, lambda = chain$lambda[s,], sigsq.eps = chain$sigsq.eps[s, ], r = chain$r[s,], y = ycont, X = X, beta = chain$beta[s,], data.comps = data.comps, modifier = modifier, kernel.method = kernel.method)
    # }
    
    ###################################################
    ## print details of the model fit so far
    opts <- set_verbose_opts(
      verbose_freq = control.params$verbose_freq, 
      verbose_digits = control.params$verbose_digits,
      verbose_show_ests = control.params$verbose_show_ests
    )
    print_diagnostics(verbose = verbose, opts = opts, curr_iter = s, tot_iter = nsamp, chain = chain, varsel = varsel, hier_varsel = hier_varsel, ztest = ztest, Z = Z, groups = groups)
    
  }
  control.params$r.params <- NULL
  chain$time2 <- Sys.time()
  chain$iter <- nsamp
  chain$family <- family
  chain$starting.values <- starting.values
  chain$control.params <- control.params
  if(!is.null(modifier)){#added by DD
    chain$X <- X[,1:(ncol(X)-ncol(modifier))]#added by DD
    if(kernel.method == "one"){
      chain$Z <- Z[,1:(ncol(Z)-ncol(modifier))]#added by DD
    }else{
      chain$Z <- Z
    }
    chain$modifier <- as.factor(orig_modifier) #added by DD
  }else{#added by DD
    chain$X <- X#added by DD
    chain$Z <- Z#added by DD
  }#added by DD
  chain$y <- y
  chain$ztest <- ztest
  chain$data.comps <- data.comps
  #if (!is.null(Znew)) chain$Znew <- Znew
  if (!is.null(groups)) chain$groups <- groups
  chain$varsel <- varsel
  chain$kernel.method <- kernel.method
  chain$gs.tau <- gs.tau
  chain$gs.sig <- gs.sig
  class(chain) <- c("bkmrfit", class(chain))
  chain
}

#' Print basic summary of BKMR model fit
#'
#' \code{print} method for class "bkmrfit"
#'
#' @param x an object of class "bkmrfit"
#' @param digits the number of digits to show when printing
#' @param ...	further arguments passed to or from other methods.
#'  
#' @export
#' 
#' @return No return value, prints basic summary of fit to console
#' 
#' @examples
#' ## First generate data set
#' y <- ex_data$y
#' Z <- ex_data$Z
#' modifier <- ex_data$X$Sex
#' X_full <- ex_data$X[,-2] #remove Sex from the covariate matrix because it is the modifier
#' #create design matrix to account for factor variables, remove the intercept column
#' X <- model.matrix(~., data=X_full)[,-1] 
#' ## Fit model 
#' ## Using only 10 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, modifier = modifier, X = X, iter = 10, verbose = FALSE) 
#' fitkm
print.bkmrfit <- function(x, digits = 5, ...) {
  cat("Fitted object of class 'bkmrfit'\n")
  cat("Iterations:", x$iter, "\n")
  cat("Outcome family:", x$family, ifelse(x$family == "binomial", "(probit link)", ""), "\n")
  cat("Model fit on:", as.character(x$time2), "\n")
}

#' Summarizing BKMR model fits
#'
#' \code{summary} method for class "bkmrfit"
#'
#' @param object an object of class "bkmrfit"
#' @param q quantiles of posterior distribution to show
#' @param digits the number of digits to show when printing
#' @param show_ests logical; if \code{TRUE}, prints summary statistics of posterior distribution
#' @param show_MH logical; if \code{TRUE}, prints acceptance rates from the Metropolis-Hastings algorithm
#' @param ...	further arguments passed to or from other methods.
#'  
#' @export
#' 
#' @return No return value, prints more detailed summary of fit to console
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
#' summary(fitkm)
summary.bkmrfit <- function(object, q = c(0.025, 0.975), digits = 5, show_ests = TRUE, show_MH = TRUE, ...) {
  x <- object
  elapsed_time <- difftime(x$time2, x$time1)
  
  print(x, digits = digits)  
  cat("Running time: ", round(elapsed_time, digits), attr(elapsed_time, "units"), "\n")
  cat("\n")
  
  if (show_MH) {
    cat("Acceptance rates for Metropolis-Hastings algorithm:\n")
    accep_rates <- data.frame()
    ## lambda
    nm <- "lambda"
    rate <- colMeans(x$acc.lambda[2:x$iter, ,drop = FALSE])
    if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
    accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
    ## r_m
    if (!x$varsel) {
      nm <- "r"
      rate <- colMeans(x$acc.r[2:x$iter, , drop = FALSE])
      if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
    } else {
      nm <- "r/delta (overall)"
      rate <- mean(x$acc.rdelta[2:x$iter])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      ##
      nm <- "r/delta  (move 1)"
      rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 1])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      ##
      nm <- "r/delta  (move 2)"
      rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 2])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      if (!is.null(x$groups)) {
        nm <- "r/delta  (move 3)"
        rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 3])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      }
    }
    print(accep_rates)
  }
  if (show_ests) {
    sel <- with(x, seq(floor(iter/2) + 1, iter))
    cat("\nParameter estimates (based on iterations ", min(sel), "-", max(sel), "):\n", sep = "")
    ests <- ExtractEsts(x, q = q, sel = sel)
    if (!is.null(ests$h)) {
      ests$h <- ests$h[c(1,2,nrow(ests$h)), ]
    }
    if (!is.null(ests$ystar)) {
      ests$ystar <- ests$ystar[c(1,2,nrow(ests$ystar)), ]
    }
    summ <- with(ests, rbind(beta, sigsq.eps, r, lambda))
    if (!is.null(ests$h)) {
      summ <- rbind(summ, ests$h)
    }
    if (!is.null(ests$ystar)) {
      summ <- rbind(summ, ests$ystar)
    }
    summ <- data.frame(param = rownames(summ), round(summ, digits))
    rownames(summ) <- NULL
    print(summ)
    if (x$varsel) {
      cat("\nPosterior inclusion probabilities:\n")
      pips <- ExtractPIPs(x)
      pips[, -1] <- round(pips[, -1], digits)
      print(pips)
    }
  }
  return()
}
