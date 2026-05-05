#' Compute summary statistics
#'
#' @param q vector of quantiles
#' @param s vector of posterior samples
#' 
#' @noRd
SummarySamps <- function(s, q = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    qs <- quantile(s, q)
    names(qs) <- paste0("q_", 100*q)
    summ <- c(mean = mean(s), sd = sd(s), qs)
    summ <- matrix(summ, nrow = 1, dimnames = list(NULL, names(summ)))
}

#'Compute the effective sample size
#'
#'@param par vector for specified parameter
#'
effective_sample_size <- function(x){
  n <- length(x)
  
  max_lag <- floor(n/2)
  
  acf_vals <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]
  
  positive_acf <- acf_vals[acf_vals > 0]
  
  ess <- n / (1 + 2 * sum(positive_acf))
  
  return(ess)
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
      sel <- with(fit, seq(burnin + 1, iter))
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

ExtractEsts2 <- function(fit, par, q = c(0.025, 0.25, 0.5, 0.75, 0.975), sel = NULL) {
  if(is.null(par)){
    stop("Need to supply one parameter")
  }
  
  if(length(par) != 1){
    stop("par must be exactly one parameter")
  }
  
  if(par %in% c("beta", "sigsq.eps", "r", "lambda", "h", "hnew", "ystar")){
    if (inherits(fit, "bkmrfit")) {
      if (is.null(sel)) {
        sel <- with(fit, seq(burnin + 1, iter))
      }
      
      #return beta estimates if it is beta
      if(par == "beta"){
        beta <- t(apply(fit$beta[sel, , drop = FALSE], 2, SummarySamps, q = q))
        
        if (nrow(beta) > 1) {
          rownames(beta) <- paste0("beta", 1:nrow(beta))
        } else {
          rownames(beta) <- "beta"
        }
        
        names <- SummarySamps(fit$beta[sel, , drop = FALSE], q = q)
        
        colnames(beta) <- colnames(names)
        
        effective_sample_size <- apply(fit$beta[sel, , drop = FALSE], 2, effective_sample_size)
        return(cbind(beta, effective_sample_size))
      }
      
      #return sigsq.eps estimates if it is sigsq.eps
      if(par == "sigsq.eps"){
        sigsq.eps <- SummarySamps(fit$sigsq.eps[sel], q = q)
        rownames(sigsq.eps) <- "sigsq.eps"
        
        #effective_sample_size <- apply(fit$sigsq.eps[sel,], 2, effective_sample_size)
        return(cbind(sigsq.eps, effective_sample_size = effective_sample_size(fit$sigsq.eps[sel])))
      }
      
      #return r estimates if it is r
      if(par == "r"){
        r <- t(apply(fit$r[sel, , drop = FALSE], 2, SummarySamps, q = q))
        rownames(r) <- paste0("r", 1:nrow(r))
        
        colnames(r) <- colnames(SummarySamps(fit$r[sel], q= q))
        
        effective_sample_size <- apply(fit$r[sel, , drop = FALSE], 2, effective_sample_size)
        return(cbind(r, effective_sample_size))
        
      }
      
      #return lambda estimates if it is lambda
      if(par == "lambda"){
        lambda <- t(apply(fit$lambda[sel, ,drop = FALSE], 2, SummarySamps, q = q))
        if (nrow(lambda) > 1) {
          rownames(lambda) <- paste0("lambda", 1:nrow(lambda))
        } else {
          rownames(lambda) <- "lambda"
        }
        
        colnames(lambda) <- colnames(SummarySamps(fit$lambda[sel], q = q))
        
        effective_sample_size <- apply(fit$lambda[sel, , drop = FALSE], 2, effective_sample_size)
        return(cbind(lambda, effective_sample_size))
      }
      
      #return est.h estimates if it is est.h
      if(par == "est.h"){
        if (fit$est.h) {
          h <- t(apply(fit$h.hat[sel, ], 2, SummarySamps, q = q))
          rownames(h) <- paste0("h", 1:nrow(h))
          colnames(h) <- colnames(SummarySamps(fit$h.hat[sel,], q=q))
          effective_sample_size <- apply(fit$h.hat[sel, , drop = FALSE], 2, effective_sample_size)
          return(cbind(h, effective_sample_size))
        }
        else{
          stop('"est.h" is not a valid parameter')
        }
      }
      
      #return hnew estimates if it is hnew
      if(par == "hnew"){
        if (!is.null(fit$hnew)) {
          hnew <- t(apply(fit$hnew[sel, ], 2, SummarySamps, q = q))
          rownames(hnew) <- paste0("hnew", 1:nrow(hnew))
          colnames(hnew) <- SummarySamps(fit$hnew[sel, ], q = q)
          effective_sample_size <- apply(fit$hnew[sel, , drop = FALSE], 2, effective_sample_size)
          return(cbind(hnew, effective_sample_size))
        }
        else{
          stop('"hnew" is not a valid parameter')
        }
      }
      
      #return ystar estimates if it is ystar
      if(par == "ystar"){
        if (!is.null(fit$ystar)) {
          ystar <- t(apply(fit$ystar[sel, ], 2, SummarySamps, q = q))
          rownames(ystar) <- paste0("ystar", 1:nrow(ystar))
          colnames(ystar) <- colnames(fit$ystar[sel, ], q = q)
          effective_sample_size <- apply(fit$ystar[sel, , drop = FALSE], 2, effective_sample_size)
          return(cbind(ystar, effective_sample_size))
        }
      }
    }
  }
  
  else{
    stop('Invalid parameter supllied, must be one of "beta", "sigsq.eps", "r", "lambda", "h", "hnew", or "ystar" ')
  }
  
}

#' Extract samples
#'
#' Extract samples of each parameter from the BKMR fit
#'
#' @inheritParams ExtractEsts2
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
ExtractSamps2 <- function(fit, par, sel = NULL) {
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
  
  if(par %in% names(res)){
    return(res[[par]])
  }
  
  else{
    stop('Invalid parameter supllied, must be one of "beta", "sigsq.eps", "r", "lambda", "h", "hnew", or "ystar"')
  }
}
