makeVcomps <- function(r, lambda, Z, data.comps, modifier = NULL) {
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
    
    #for model with random intercept
    #if (data.comps$randint) { #changed from data.comps$nlambda == 2 upon GS-tau model implementation
    #  V <- V + lambda[length(lambda)-1]*data.comps$crossTT
    #}
    #not sure what this does really 
    Vinv_list <- vector("list", G)
    logdetVinv <- 0
    
    Vinv = matrix(0, nrow = nrow(Z), ncol = nrow(Z))
    
    for(g in 1:G){
      
      id <- block_ids[[g]]
      Z_g <- Z_list[[g]]
      
      Kpart_g <- makeKpart(r, Z_g)
      K_g <- exp(-Kpart_g)
      
      V_g <- diag(1, length(id)) + scalar_lambda * K_g
      
      cholV_g <- chol(V_g) #recode for speed in Rcpp?
      
      Vinv_list[[g]] <- chol2inv(cholV_g) #recode for speed in Rcpp?
      
      logdetVinv <- logdetVinv - 2* sum(log(diag(cholV_g)))
      
      Vinv[id,id] = Vinv_list[[g]]
    }
    
    Vcomps <- list(Vinv = Vinv, Vinv_list = Vinv_list, logdetVinv = logdetVinv, scalar_lambda = scalar_lambda)
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