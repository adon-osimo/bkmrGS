#' 
#' @export
generate_code <- function(
    fit_name = "fitkm",
    exposure = "Overall",
    analysis = "Group Specific",
    centered = FALSE,
    sel = "NULL",
    m.fixed = NULL,
    qs = NULL,
    q.fixed = NULL,
    qs.diff = NULL,
    mod.diff = NULL,
    qs.fixed = NULL){
  
  
  code <- c(
    build_header(exposure, analysis),
    build_data_extract(fit_name),
    build_user_changeable_data(exposure, analysis, m.fixed, qs, q.fixed,
                               qs.diff, mod.diff, qs.fixed),
    build_prediction_function(sel),
    build_summary_function(analysis),
    build_point_constructor(exposure, analysis),
    build_iteration_engine(exposure, analysis),
    build_output_formatter(exposure, analysis),
    build_plotting(exposure, analysis)
  )
  
  paste(code, collapse = "\n")
}

build_header <- function(exposure, analysis){
    
    c(paste0("# ", exposure, " ", analysis, " Analysis"), "")
}

build_data_extract <- function(fit_name){
  c("#Update this line with the name of your bkmrfit object!!",
    paste0("fit <- ", fit_name, ""))
}

build_user_changeable_data <- function(exposure, analysis, m.fixed, qs, q.fixed,
                                       qs.diff, mod.diff, qs.fixed){
  ret <- c("#ALL OF THESE ARE CHANGEABLE AT YOUR DISCRESION",
           "#BE SURE THAT WHEN YOU RUN THIS CODE YOU UPDATE THE FOLLOWING LINES",
           "#TO REFLECT WHAT VALUES YOU WANT TO SEE")
  
  if(exposure == "Overall"){
    ret <- c(ret, c("",
                    paste0("qs <- ", as.character(qs), " "),
                    "",
                    paste0("q.fixed <- ", as.character(q.fixed), " ")))
  }
  
  if (analysis == "Group Specific"){
    ret <- c(ret, c("",
                    paste0("m.fixed <- '", as.character(m.fixed), "' ")))
  }
  
  if(analysis == "Between Groups"){
    ret <- c(ret, c("",
                    paste0("mod.diff <- ", as.character(mod.diff), " ")))
  }
  
  if(exposure == "Single"){
    ret <- c(ret, c("",
                    paste0("qs.diff <-", as.character(qs.diff), " ")))
    
    if(analysis == "Group Specific"){
      ret <- c(ret, c("",
                      paste0("q.fixed <- ", as.character(q.fixed), " ")))
    }
    
    else {
      ret <- c(ret, c("",
                      paste0("qs.fixed <- ", as.character(qs.fixed), " ")))
    }
  }
  
  ret <- c(ret, c("",
                  "#THE REST OF THIS CODE IS JUST RUNNABLE, DO NOT EDIT",
                  "#ANY CODE BETWEEN NOW AND plot_obj"))
  
  ret

}

build_prediction_function <- function(sel){
  c("preds.fun <- function(znew, modnew) { #do not comment out this line",
    "ComputePostmeanHnew(fit = fit,",
    "y = fit$y,",
    "Z = fit$Z,",
    "X = fit$X,",
    "modifier = fit$modifier,",
    "Znew = znew,",
    "mod_new = modnew,",
    paste0("sel = ", as.character(sel), ","),
    "method = 'exact')}")
}

build_summary_function <- function(analysis){
  if(analysis == "Group Specific"){
    
    c(
      "# Overall summary function",
      "summary.fun <- function(point1, point2, modnew = NULL, preds.fun) {",
      "  cc <- c(-1, 1)",
      "  newz <- rbind(point1, point2)",
      "  preds <- preds.fun(newz, modnew)",
      "  diff <- drop(cc %*% preds$postmean)",
      "  diff.sd <- drop(sqrt(cc %*% preds$postvar %*% cc))",
      "  c(est = diff, sd = diff.sd)",
      "}"
    )
  
  } else if(analysis == "Between Groups") {
  
  c(
    "# Interaction summary function",
    "summary.fun <- function(newz.q1, newz.q2, modnew.1, modnew.2, preds.fun) {",
    "  newz <- rbind(newz.q1, newz.q2)",
    "  modnew <- c(modnew.1, modnew.2)",
    "  preds <- preds.fun(newz, modnew)",
    "  cc <- c(-1 * c(-1, 1), c(-1, 1))",
    "  int <- drop(cc %*% preds$postmean)",
    "  int.sd <- drop(sqrt(cc %*% preds$postvar %*% cc))",
    "  c(est = int, sd = int.sd)",
    "}"
    )
  
  }
}

#This is one that I am going to have to update
build_point_constructor <- function(exposure, analysis){
  
  ret <- c()
  
  #if(analysis == "Group Specific"){
  #  if(exposure == "Overall"){
  #    
  #    ret <- c(ret, c(
  #      "# Single-variable Overall settings",
  #      "which.z <- 1:ncol(Z)",
  #      "z.names <- colnames(Z)",
  #      ""
  #    ))
  #    
  #  }
  #  
  #  
  #  ret <- c(ret, c("if(!is.null(m.fixed)) {",
  #  "  modnew <- matrix(rep(m.fixed, 2), ncol = 1)",
  #  "  Z_for_quants <- Z[fit$modifier == m.fixed, , drop = FALSE]",
  #  "} else {",
  #  "  modnew <- NULL",
  #  "  Z_for_quants <- Z",
  #  "}",
  #  ""))
  #}
  
  if(analysis == "Between Groups"){
  
    ret <- c(ret, c("Z_support <- function(Z, mod.diff, modifier){",
    "  mins <- c()",
    "  maxs <- c()",
    "  for(idx in mod.diff){",
    "    mod_idx <- which(modifier==idx)",
    "    mins <- rbind(mins, apply(Z[mod_idx,], 2, min))",
    "    maxs <- rbind(maxs, apply(Z[mod_idx,], 2, max))",
    "  }",
    "  expo_min <- apply(mins, 2, max)",
    "  expo_max <- apply(maxs, 2, min)",
    "  Z_idx <- c()",
    "  for(z_col in 1:ncol(Z)){",
    "    sub_z <- Z[,z_col]",
    "    Z_idx <- cbind(Z_idx, sub_z > expo_min[z_col] & sub_z < expo_max[z_col])",
    "  }",
    "  Z_for_quants <- Z[which(rowSums(Z_idx)==ncol(Z)),]",
    "  return(Z_for_quants)",
    "}"))
    
    if(exposure == "Single"){
      
      ret <- c(ret, c(
        "which.z <- 1:ncol(Z)",
        "z.names <- colnames(Z)"
      ))
      
    }
    
    
    ret <- c(ret, c("# Restrict exposure support",
                    "Z_for_quants <- Z_support(",
                    "  Z = Z,",
                    "  mod.diff = mod.diff,",
                    "  modifier = fit$modifier",
                    ")",
                    "modnew.1 <- rep(mod.diff[1], 2)",
                    "modnew.2 <- rep(mod.diff[2], 2)"))
  
  }
    
  ret
}

build_iteration_engine <- function(exposure, analysis){
  if(exposure == "Overall" && analysis == "Group Specific"){
    
    c("if(m.fixed == 'All'){",
      "   iter_over <- c(unique(as.character(fit$modifier)))",
      "} else{",
      "   iter_over <- c(m.fixed)}",
      "results <- data.frame()",
      "for(m in iter_over){",
      "  m.fixed <- m",
      "  if(!is.null(m.fixed)) {",
      "    modnew <- matrix(rep(m.fixed, 2), ncol = 1)",
      "    Z_for_quants <- Z[fit$modifier == m.fixed, ]",
      "  } else {",
      "    modnew <- NULL",
      "    Z_for_quants <- Z",
      "  }",
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "  tmp_fn <- function(q) {",
      "    point2 <- apply(Z_for_quants, 2, quantile, q)",
      "    summary.fun(",
      "      point1 = point1,",
      "      point2 = point2,",
      "      modnew = modnew,",
      "      preds.fun = preds.fun",
      "    )",
      "  }",
      "  results_temp <- t(sapply(qs, tmp_fn))",
      "  ",
      "  results_temp <- data.frame(",
      "    quantile = qs,",
      "    results_temp",
      "  )",
      "  ",
      "  results_temp$modifier <- rep(m.fixed, nrow(results_temp))",
      "  results<- rbind(results, results_temp)",
      "}"
    )
    
  }
  
  else if(exposure == "Single" && analysis == "Group Specific"){
    
    c("if(m.fixed == 'All'){",
      "   iter_over <- c(unique(as.character(fit$modifier)))",
      "} else{",
      "   iter_over <- c(m.fixed)}",
      "results <- data.frame()",
      "for(m in iter_over){",
      "  m.fixed <- m",
      "  modnew <- matrix(rep(m.fixed, 2), ncol = 1)",
      "  Z_for_quants <- Z[fit$modifier == m.fixed, ]",
      "  results_tmp <-  lapply(seq_along(q.fixed), function(i) {",
      "    lapply(seq_along(which.z), function(j) {",
      "      point1 <- apply(Z_for_quants, 2, quantile, q.fixed[i])",
      "      point2 <- point1",
      "      point1[which.z[j]] <- quantile(",
      "        Z_for_quants[, which.z[j]],",
      "        qs.diff[1]",
      "      )",
      "      point2[which.z[j]] <- quantile(",
      "        Z_for_quants[, which.z[j]],",
      "        qs.diff[2]",
      "      )",
      "      out <- summary.fun(",
      "        point1 = point1,",
      "        point2 = point2,",
      "        modnew = modnew,",
      "        preds.fun = preds.fun",
      "      )",
      "      data.frame(",
      "        q.fixed = q.fixed[i],",
      "        variable = z.names[j],",
      "        t(out)",
      "      )",
      "    })",
      "  })",
      "  results_tmp <- do.call(",
      "    rbind,",
      "    unlist(results_tmp, recursive = FALSE)",
      "  )",
      "  results_tmp$variable <- factor(",
      "    results_tmp$variable,",
      "    levels = z.names",
      "  )",
      "  results_tmp$q.fixed <- as.factor(results_tmp$q.fixed)",
      "  results_tmp$modifier <- rep(m.fixed, nrow(results_tmp))",
      "  results <- rbind(results, results_tmp)",
      "}"
    )
    
  }
  
  else if(exposure == "Overall" && analysis == "Between Groups"){
    
    c("tmp_fn <- function(q) {",
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "  point2 <- apply(Z_for_quants, 2, quantile, q)",
      "  newz.q1 <- rbind(point1, point2)",
      "  newz.q2 <- rbind(point1, point2)",
      "  summary.fun(",
      "    newz.q1 = newz.q1,",
      "    newz.q2 = newz.q2,",
      "    modnew.1 = modnew.1,",
      "    modnew.2 = modnew.2,",
      "    preds.fun = preds.fun",
      "  )",
      "}",
      "results <- t(sapply(qs, tmp_fn))"
    )
    
  }
  
  else if(exposure == "Single" && analysis == "Between Groups"){
    
    c("results <- lapply(which.z, function(j) {",
      "  q.fixed <- qs.fixed[1]",
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "  point2 <- point1",
      "  point1[j] <- quantile(Z_for_quants[, j], qs.diff[1])",
      "  point2[j] <- quantile(Z_for_quants[, j], qs.diff[2])",
      "  newz.q1 <- rbind(point1, point2)",
      "  q.fixed <- qs.fixed[2]",
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "  point2 <- point1",
      "  point1[j] <- quantile(Z_for_quants[, j], qs.diff[1])",
      "  point2[j] <- quantile(Z_for_quants[, j], qs.diff[2])",
      "  newz.q2 <- rbind(point1, point2)",
      "  out <- summary.fun(",
      "    newz.q1 = newz.q1,",
      "    newz.q2 = newz.q2,",
      "    modnew.1 = modnew.1,",
      "    modnew.2 = modnew.2,",
      "    preds.fun = preds.fun",
      "  )"
    )
    
  }
  
}

#Need to change this to only have the results object returned
build_output_formatter <- function(exposure, analysis){

  if(analysis == "Group Specific"){
    
    c("results")
    
  }
  
  else if(exposure == "Overall" && analysis == "Between Groups"){
    
    c("results <- data.frame(",
      "  quantile = qs,",
      "  results",
      ")",
      "print(results)"
    )
    
  }
  
  else if(exposure == "Single" && analysis == "Between Groups"){
    c("  data.frame(",
    "    variable = z.names[j],",
    "    t(out)",
    "  )",
    "})",
    "results <- do.call(rbind, results)",
    "results")
  }
}

#Need to change this to only have the plot_obj returned
#I think this is good, just have to check that what I am calling on is correct
build_plotting <- function(exposure, analysis){
  if(exposure == "Overall" && analysis == "Group Specific"){
    c("plot_obj <- ggplot(results,",
      "       aes(x=quantile, y = est, ymin = est - 1.96*sd,",
      "           ymax = est + 1.96*sd, color = modifier, ",
      "           shape = modifier)) + ",
      "  geom_hline(yintercept = 0) +",
      "  geom_pointrange(position = position_dodge(width = 0.05), size = 0.5) +",
      "  theme_bw()+",
      "  theme(legend.position = 'bottom',",
      "        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ",
      "  xlab('Quantiles')+",
      "  ylab('Difference in Response \n per Change in Exposure') + ",
      "  labs(color='modifier_name', shape = 'modifier_name')" )
  }

  else if(exposure == "Single" && analysis == "Group Specific"){
    c("plot_obj <- ggplot(results,",
      "                   aes(x=q.fixed, y = est, ymin = est - 1.96*sd,",
      "                       ymax = est + 1.96*sd, color = modifier, ",
      "                       shape = modifier))+",
      "  geom_hline(yintercept=0)+",
      "  geom_pointrange(position = position_dodge(width = 0.05), size=0.5)+",
      "  facet_wrap(vars(variable))+",
      "  theme_bw()+",
      "  theme(legend.position = 'bottom',",
      "        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ",
      "        xlab('Quantile')+",
      "        ylab('Difference in Response \n per Change in Exposure') +",
      "        labs(color='modifier_name', shape = 'modifier_name')")
  }
  
  else if(exposure == "Overall" && analysis == "Between Groups"){
      c("plot_obj <- ggplot(results, aes(x=quantile, y = est,", 
      "                    ymin = est - 1.96*sd, ymax = est + 1.96*sd))+",
      "   geom_hline(yintercept=0)+",
      "   geom_pointrange(size=0.5)+",
      "   theme_bw()+",
      "   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+",
      "   xlab('CHANGE THIS X-AXIS')+",
      "   ylab('CHANGE THIS Y-AXIS')")
  }
  
  else if(exposure == "Single" && analysis == "Between Groups"){
    c("plot_obj <- ggplot(results, aes(x='', y = est,", 
      "                    ymin = est - 1.96*sd, ymax = est + 1.96*sd))+",
      "   geom_hline(yintercept=0)+",
      "   geom_pointrange(size=0.5)+",
      "   facet_wrap(vars(variable))+",
      "   theme_bw()+",
      "   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+",
      "   xlab('CHANGE THIS X-AXIS')+",
      "   ylab('CHANGE THIS Y-AXIS')")
  }
}









