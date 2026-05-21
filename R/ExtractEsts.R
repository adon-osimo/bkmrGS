#' Compute summary statistics
#'
#' @param q vector of quantiles
#' @param s vector of posterior samples
#' 
#' @noRd
SummarySamps <- function(s, q = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    qs <- quantile(s, q)
    names(qs) <- paste0("q_", 100*q)
    std <- sd(s)
    ess <- effective_sample_size(s)
    summ <- c(mean = mean(s), sd = std, mean_se = mean_standard_error(std, ess), qs, ess = ess)
    summ <- matrix(summ, nrow = 1, dimnames = list(NULL, names(summ)))
}

#' Compute the effective sample size
#'
#' @param x vector for specified parameter
#'
effective_sample_size <- function(x){
  n <- length(x)
  
  max_lag <- floor(n/2)
  
  acf_vals <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]
  
  positive_acf <- acf_vals[acf_vals > 0]
  
  ess <- n / (1 + 2 * sum(positive_acf))
  
  return(ess)
}

#' Compute the mean standard error, like in Rstan and Stan
#' Take the standard deviation estimate and divide it by the square root of the ess
#'
#' @param std value for standard deviation
#' @param ess value for effecitve sample size
#'
mean_standard_error <- function(std, ess){
  ret <- NULL
  if(!is.null(ess)){
    ret <- std / sqrt(ess)
  }
  return(ret)
}


#' Extract summary statistics
#'
#' Obtain summary statistics of each parameter from the BKMR fit
#'
#' @param fit An object containing the results returned by a the \code{kmbayes} function 
#' @param q vector of quantiles
#' @param sel logical expression indicating samples to keep; defaults to keeping the second half of all samples 
#'
#' @export
#' 
#' @return a list where each component is a data frame containing the summary statistics of the posterior distribution of one of the parameters (or vector of parameters) being estimated
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
#' 
#' ests <- ExtractEsts(fitkm)
#' names(ests)
#' ests$beta
ExtractEsts <- function(fit, q = c(0.025, 0.25, 0.5, 0.75, 0.975), sel = NULL) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(sel)) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter))
    }
    sigsq.eps <- SummarySamps(fit$sigsq.eps[sel], q = q)
    rownames(sigsq.eps) <- "sigsq.eps"
    
    r <- t(apply(fit$r[sel, , drop = FALSE], 2, SummarySamps, q = q))
    rownames(r) <- paste0("r", 1:nrow(r))
    
    beta <- t(apply(fit$beta[sel, , drop = FALSE], 2, SummarySamps, q = q))
    
    lambda <- t(apply(fit$lambda[sel, ,drop = FALSE], 2, SummarySamps, q = q))
    if (nrow(lambda) > 1) {
      rownames(lambda) <- paste0("lambda", 1:nrow(lambda))
    } else {
      rownames(lambda) <- "lambda"
    }
    
    if (fit$est.h) {
      h <- t(apply(fit$h.hat[sel, ], 2, SummarySamps, q = q))
      rownames(h) <- paste0("h", 1:nrow(h))
    }
    
    if (!is.null(fit$hnew)) {
      hnew <- t(apply(fit$hnew[sel, ], 2, SummarySamps, q = q))
      rownames(hnew) <- paste0("hnew", 1:nrow(hnew))
    }
    
    if (!is.null(fit$ystar)) {
      ystar <- t(apply(fit$ystar[sel, ], 2, SummarySamps, q = q))
      rownames(ystar) <- paste0("ystar", 1:nrow(ystar))
    }
  }
  
  if (nrow(beta) > 1) {
    rownames(beta) <- paste0("beta", 1:nrow(beta))
  } else {
    rownames(beta) <- "beta"
  }
  
  colnames(beta) <- colnames(sigsq.eps)
  colnames(r) <- colnames(sigsq.eps)
  colnames(lambda) <- colnames(sigsq.eps)
  if (fit$est.h) {
    colnames(h) <- colnames(sigsq.eps)
  }
  if (!is.null(fit$hnew)) {
    colnames(hnew) <- colnames(sigsq.eps)
  }
  if (!is.null(fit$ystar)) {
    colnames(ystar) <- colnames(sigsq.eps)
  }
  
  ret <- list(sigsq.eps = data.frame(sigsq.eps), beta = beta, lambda = lambda, r = r)
  if (fit$est.h) ret$h <- h
  if (!is.null(fit$hnew)) ret$hnew <- hnew
  if (!is.null(fit$ystar)) ret$ystar <- ystar
  
  ret
}

#' Extract samples
#'
#' Extract samples of each parameter from the BKMR fit
#'
#' @inheritParams ExtractEsts
#'
#' @export
#' @return a list where each component contains the posterior samples of one of the parameters (or vector of parameters) being estimated
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
#' 
#' samps <- ExtractSamps(fitkm)
ExtractSamps <- function(fit, sel = NULL) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(sel)) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter))
    }
    
    sigsq.eps <- fit$sigsq.eps[sel]
    sig.eps <- sqrt(sigsq.eps)
    r <- fit$r[sel, , drop = FALSE]
    beta <- fit$beta[sel, ]
    lambda <- fit$lambda[sel, ]
    tau <- lambda*sigsq.eps
    h <- fit$h.hat[sel, ]
    if (!is.null(fit$hnew)) hnew <- fit$hnew[sel, ]
    if (!is.null(fit$ystar)) ystar <- fit$ystar[sel, ]
  }
  
    if (!is.null(ncol(beta))) colnames(beta) <- paste0("beta", 1:ncol(beta))
    if (!is.null(ncol(lambda))) colnames(lambda) <- paste0("lambda", 1:ncol(lambda))
    colnames(r) <- paste0("r", 1:ncol(r))
    colnames(h) <- paste0("h", 1:ncol(h))
    if (!is.null(fit$hnew)) colnames(hnew) <- paste0("hnew", 1:ncol(hnew))
    if (!is.null(fit$ystar)) colnames(ystar) <- paste0("ystar", 1:ncol(ystar))
    
    res <- list(sigsq.eps = sigsq.eps, sig.eps = sig.eps, r = r, beta = beta, lambda = lambda, tau = tau, h = h)
    if (!is.null(fit$hnew)) res$hnew <- hnew
    if (!is.null(fit$ystar)) res$ystar <- ystar
    res
}


#' Extract summary statistics
#'
#' Obtain summary statistics of each parameter from the BKMR fit
#'
#' @param fit An object containing the results returned by a the \code{kmbayes} function
#' @param par The parameter of interest that you want returned
#' @param q vector of quantiles
#' @param sel logical expression indicating samples to keep; defaults to keeping the burin in removed samples
#' @param digits input for the amount of output digits required 
#'
#' @export
#' 
#' @return a list where each component is a data frame containing the summary statistics of the posterior distribution of one of the parameters (or vector of parameters) being estimated
#' 
#' @examples
#' ## Fit model 
#' ## Using only 10 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayesWrapper(y ~ h(Log_Lead, Log_Manganese, Log_Arsenic, mod = Sex) +
#'                           Age + Gestation + Delivery + Birth_order +
#'                           Education_parent1 + Education_parent2 + Smoking +
#'                           HOME_emotional + HOME_avoid + HOME_careg + 
#'                           HOME_env + HOME_play + HOME_stim + 
#'                           Energy, 
#'                         data = Liu_data, 
#'                         iter = 10, 
#'                         verbose = FALSE) 
#' 
#' ests <- GetEsts(fitkm)
#' names(ests)
#' ests$beta

GetEsts <- function(fit, par = NULL, q = c(0.025, 0.25, 0.5, 0.75, 0.975), sel = NULL, digits = 5) {
  if(length(par) > 1){
    stop("par must be exactly one parameter")
  }
  
  if(is.null(par) || par %in% c("beta", "sigsq.eps", "r", "lambda", "h", "hnew", "ystar") ){
    if (inherits(fit, "bkmrfit")) {
      if (is.null(sel)) {
        sel <- with(fit, seq(burnin + 1, iter))
      }
      sigsq.eps <- round(SummarySamps(fit$sigsq.eps[sel], q = q), digits)
      rownames(sigsq.eps) <- "sigsq.eps"
      
      r <- round(t(apply(fit$r[sel, , drop = FALSE], 2, SummarySamps, q = q)), digits)
      rownames(r) <- paste0("r", 1:nrow(r))
      
      beta <- round(t(apply(fit$beta[sel, , drop = FALSE], 2, SummarySamps, q = q)), digits)
      
      lambda <- round(t(apply(fit$lambda[sel, ,drop = FALSE], 2, SummarySamps, q = q)), digits)
      if (nrow(lambda) > 1) {
        rownames(lambda) <- paste0("lambda", 1:nrow(lambda))
      } else {
        rownames(lambda) <- "lambda"
      }
      
      if (fit$est.h) {
        h <- round(t(apply(fit$h.hat[sel, ], 2, SummarySamps, q = q)), digits)
        rownames(h) <- paste0("h", 1:nrow(h))
      }
      
      if (!is.null(fit$hnew)) {
        hnew <- round(t(apply(fit$hnew[sel, ], 2, SummarySamps, q = q)), digits)
        rownames(hnew) <- paste0("hnew", 1:nrow(hnew))
      }
      
      if (!is.null(fit$ystar)) {
        ystar <- round(t(apply(fit$ystar[sel, ], 2, SummarySamps, q = q)), digits)
        rownames(ystar) <- paste0("ystar", 1:nrow(ystar))
      }
    }
    
    if (nrow(beta) > 1) {
      rownames(beta) <- paste0("beta", 1:nrow(beta))
    } else {
      rownames(beta) <- "beta"
    }
    
    colnames(beta) <- colnames(sigsq.eps)
    colnames(r) <- colnames(sigsq.eps)
    colnames(lambda) <- colnames(sigsq.eps)
    if (fit$est.h) {
      colnames(h) <- colnames(sigsq.eps)
    }
    if (!is.null(fit$hnew)) {
      colnames(hnew) <- colnames(sigsq.eps)
    }
    if (!is.null(fit$ystar)) {
      colnames(ystar) <- colnames(sigsq.eps)
    }
    
    ret <- list(sigsq.eps = sigsq.eps, beta = beta, lambda = lambda, r = r)
    if (fit$est.h) ret$h <- h
    if (!is.null(fit$hnew)) ret$hnew <- hnew
    if (!is.null(fit$ystar)) ret$ystar <- ystar
    if(!is.null(par)){
      
      #return beta estimates if it is beta
      if(par == "beta"){
        ret <- ret$beta
      }
        
        #return sigsq.eps estimates if it is sigsq.eps
      if(par == "sigsq.eps"){
        ret <- ret$sigsq.eps
  
      }
        
        #return r estimates if it is r
      if(par == "r"){
        ret <- ret$r
      }
        
        #return lambda estimates if it is lambda
      if(par == "lambda"){
        ret <- ret$lambda
      }
        
        #return est.h estimates if it is est.h
      if(par == "est.h"){
        if (!is.null(fit$est.h)) {
          ret <- ret$est.h
        }
        else{
          stop('"est.h" is not a valid parameter')
        }
      }
        
        #return hnew estimates if it is hnew
      if(par == "hnew"){
        if (!is.null(fit$hnew)) {
          ret <- ret$hnew
        }
        else{
          stop('"hnew" is not a valid parameter')
        }
      }
        
        #return ystar estimates if it is ystar
      if(par == "ystar"){
        if (!is.null(fit$ystar)) {
          ret <- ret$ystar
        }
        else{
          stop('"ystar" is not a valid parameter')
        }
      }
    }
  }
  else{
    stop('Invalid parameter supllied, must be one of "beta", "sigsq.eps", "r", "lambda", "h", "hnew", or "ystar" ')
  }
  cat("\nParameter estimates (based on iterations ", min(sel), "-", max(sel), "):\n", sep = "")
  ret
}

#' Extract samples
#'
#' Extract samples of each parameter from the BKMR fit
#'
#' @inheritParams GetEsts
#'
#' @export
#' @return a list where each component contains the posterior samples of one of the parameters (or vector of parameters) being estimated
#' 
#' @examples
#' ## Fit model 
#' ## Using only 10 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayesWrapper(y ~ h(Log_Lead, Log_Manganese, Log_Arsenic, mod = Sex) +
#'                           Age + Gestation + Delivery + Birth_order +
#'                           Education_parent1 + Education_parent2 + Smoking +
#'                           HOME_emotional + HOME_avoid + HOME_careg + 
#'                           HOME_env + HOME_play + HOME_stim + 
#'                           Energy, 
#'                         data = Liu_data, 
#'                         iter = 10, 
#'                         verbose = FALSE) 
#' 
#' samps <- GetSamps(fitkm)
GetSamps <- function(fit, par = NULL, sel = NULL) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(sel)) {
      sel <- with(fit, seq(burnin + 1, iter))
    }
    
    sigsq.eps <- fit$sigsq.eps[sel]
    sig.eps <- sqrt(sigsq.eps)
    r <- fit$r[sel, , drop = FALSE]
    beta <- fit$beta[sel, ]
    lambda <- fit$lambda[sel, ]
    tau <- lambda*sigsq.eps
    h <- fit$h.hat[sel, ]
    if (!is.null(fit$hnew)) hnew <- fit$hnew[sel, ]
    if (!is.null(fit$ystar)) ystar <- fit$ystar[sel, ]
  }
  
  if (!is.null(ncol(beta))) colnames(beta) <- paste0("beta", 1:ncol(beta))
  colnames(r) <- paste0("r", 1:ncol(r))
  colnames(h) <- paste0("h", 1:ncol(h))
  if (!is.null(fit$hnew)) colnames(hnew) <- paste0("hnew", 1:ncol(hnew))
  if (!is.null(fit$ystar)) colnames(ystar) <- paste0("ystar", 1:ncol(ystar))
  
  res <- list(sigsq.eps = sigsq.eps, sig.eps = sig.eps, r = r, beta = beta, lambda = lambda, tau = tau, h = h)
  if (!is.null(fit$hnew)) res$hnew <- hnew
  if (!is.null(fit$ystar)) res$ystar <- ystar
  
  if(is.null(par)){
    return(res)
  }
  
  if(par %in% names(res)){
    return(res[[par]])
  }
  
  else{
    stop('Invalid parameter supllied, must be one of "beta", "sigsq.eps", "r", "lambda", "h", "hnew", or "ystar"')
  }
}
