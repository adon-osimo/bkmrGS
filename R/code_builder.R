#' 
#' @export
generate_code <- function(
    fit_name = "fitkm",
    comparison = "Overall",
    movement = "Group Specific",
    centered = FALSE,
    sel = "NULL",
    m.fixed = NULL,
    qs = NULL,
    q.fixed = NULL,
    qs.diff = NULL,
    mod.diff = NULL,
    qs.fixed = NULL){
  
  
  code <- c(
    build_header(comparison, movement),
    build_data_extract(fit_name),
    build_user_changeable_data(comparison, movement, m.fixed, qs, q.fixed,
                               qs.diff, mod.diff, qs.fixed),
    build_prediction_function(sel),
    build_summary_function(movement),
    build_point_constructor(comparison, movement, m.fixed, qs, q.fixed),
    build_iteration_engine(comparison, movement),
    build_output_formatter(comparison, movement),
    build_plotting(comparison, movement)
  )
  
  paste(code, collapse = "\n")
}

build_header <- function(
    comparison,
    movement){
    
    c(paste0("# ", movement, " ", comparison, " Analysis"), "")
}

build_data_extract <- function(fit_name){
  c("#Update this line with the name of your bkmrfit object!!",
    paste0("fit <- ", fit_name, ""),
    "#Change this if you want a different Z",
    "Z <- fit$Z",
    "#Change this if you want a different X",
    "X <- fit$X",
    "#y <- fit$y",
    "#modifier <- fit$modifier")
}

build_user_changeable_data <- function(comparison, movement, m.fixed, qs, q.fixed,
                                       qs.diff, mod.diff, qs.fixed){
  ret <- c("#ALL OF THESE ARE CHANGEABLE AT YOUR DISCRESION",
           "#BE SURE THAT WHEN YOU RUN THIS CODE YOU UPDATE THE FOLLOWING LINES",
           "#TO REFLECT WHAT VALUES YOU WANT TO SEE")
  
  if(comparison == "Overall"){
    ret <- c(ret, c(paste0("qs <- ", as.character(qs), " "),
                    paste0("q.fixed <- ", as.character(q.fixed), " ")))
  }
  
  if (movement == "Group Specific"){
    ret <- c(ret, c(paste0("m.fixed <- '", as.character(m.fixed), "' ")))
  }
  
  if(movement == "Between Groups"){
    ret <- c(ret, c(paste0("mod.diff <- ", as.character(mod.diff), " ")))
  }
  
  if(comparison == "Single"){
    ret <- c(ret, c(paste0("qs.diff <-", as.character(qs.diff), " ")))
    
    if(movement == "Group Specific"){
      ret <- c(ret, c(paste0("q.fixed <- ", as.character(q.fixed), " ")))
    }
    
    else {
      ret <- c(ret, c(paste0("qs.fixed <- ", as.character(qs.fixed), " ")))
    }
  }
  
  ret <- c(ret, c("#THE REST OF THIS CODE IS JUST RUNNABLE, DO NOT EDIT",
                  "#ANY CODE BETWEEN NOW AND plot_obj"))
  
  ret

}

build_prediction_function <- function(sel){
  c("preds.fun <- function(znew, modnew) { #do not comment out this line",
    "ComputePostmeanHnew(fit = fit,",
    "y = fit$y,",
    "Z = Z,",
    "X = X,",
    "modifier = fit$modifier,",
    "Znew = znew,",
    "mod_new = modnew,",
    paste0("sel = ", as.character(sel), ","),
    "method = 'exact')}")
}

build_summary_function <- function(movement){
  if(movement == "Group Specific"){
    
    c(
      "# Overall summary function",
      "summary.fun <- function(point1, point2, modnew = NULL, preds.fun) {",
      "",
      "  cc <- c(-1, 1)",
      "  newz <- rbind(point1, point2)",
      "  preds <- preds.fun(newz, modnew)",
      "",
      "  diff <- drop(cc %*% preds$postmean)",
      "  diff.sd <- drop(sqrt(cc %*% preds$postvar %*% cc))",
      "",
      "  c(est = diff, sd = diff.sd)",
      "}",
      ""
    )
  
  } else if(movement == "Between Groups") {
  
  c(
    "# Interaction summary function",
    "summary.fun <- function(newz.q1, newz.q2, modnew.1, modnew.2, preds.fun) {",
    "",
    "  newz <- rbind(newz.q1, newz.q2)",
    "  modnew <- c(modnew.1, modnew.2)",
    "",
    "  preds <- preds.fun(newz, modnew)",
    "",
    "  cc <- c(-1 * c(-1, 1), c(-1, 1))",
    "",
    "  int <- drop(cc %*% preds$postmean)",
    "  int.sd <- drop(sqrt(cc %*% preds$postvar %*% cc))",
    "",
    "  c(est = int, sd = int.sd)",
    "}",
    ""
    )
  
  }
}

build_point_constructor <- function(comparison, movement, m.fixed, qs, q.fixed){
  
  ret <- c()
  
  if(movement == "Group Specific"){
    
    if(comparison == "Overall"){
      
      ret <- c(ret, c(
        "# Single-variable Overall settings",
        "which.z <- 1:ncol(Z)",
        "z.names <- colnames(Z)",
        ""
      ))
      
    }
    
    
    ret <- c(ret, c("",
    paste0("m.fixed <- '", m.fixed, "'"),
    "",
    "if(!is.null(m.fixed)) {",
    "  modnew <- matrix(rep(m.fixed, 2), ncol = 1)",
    "  Z_for_quants <- Z[fit$modifier == m.fixed, , drop = FALSE]",
    "} else {",
    "  modnew <- NULL",
    "  Z_for_quants <- Z",
    "}",
    ""))
  }
  
  if(movement == "Between Groups"){
  
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
    
    if(comparison == "Single"){
      
      ret <- c(ret, c(
        "# Interaction settings",
        "",
        "which.z <- 1:ncol(Z)",
        "z.names <- colnames(Z)",
        ""
      ))
      
    }
    
    
    ret <- c(ret, c("# Restrict exposure support",
                    "Z_for_quants <- Z_support(",
                    "  Z = Z,",
                    "  mod.diff = mod.diff,",
                    "  modifier = fit$modifier",
                    ")",
                    "",
                    "# Modifier values",
                    "modnew.1 <- rep(mod.diff[1], 2)",
                    "modnew.2 <- rep(mod.diff[2], 2)",
                    ""))
  
  }
    
  ret
}

build_iteration_engine <- function(comparison, movement){
  if(comparison == "Overall" && movement == "Group Specific"){
    
    c(
      "# Run Overall Overall analysis",
      "point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "tmp_fn <- function(q) {",
      "",
      "  point2 <- apply(Z_for_quants, 2, quantile, q)",
      "",
      "  summary.fun(",
      "    point1 = point1,",
      "    point2 = point2,",
      "    modnew = modnew,",
      "    preds.fun = preds.fun",
      "  )",
      "}",
      "",
      "results <- t(sapply(qs, tmp_fn))",
      ""
    )
    
  }
  
  else if(comparison == "Single" && movement == "Group Specific"){
    
    c(
      "# Run Single-variable Overall analysis",
      "",
      "results <- lapply(seq_along(q.fixed), function(i) {",
      "",
      "  lapply(seq_along(which.z), function(j) {",
      "",
      "    point1 <- apply(Z_for_quants, 2, quantile, q.fixed[i])",
      "    point2 <- point1",
      "",
      "    point1[which.z[j]] <- quantile(",
      "      Z_for_quants[, which.z[j]],",
      "      qs.diff[1]",
      "    )",
      "",
      "    point2[which.z[j]] <- quantile(",
      "      Z_for_quants[, which.z[j]],",
      "      qs.diff[2]",
      "    )",
      "",
      "    out <- summary.fun(",
      "      point1 = point1,",
      "      point2 = point2,",
      "      modnew = modnew,",
      "      preds.fun = preds.fun",
      "    )",
      "",
      "    data.frame(",
      "      q.fixed = q.fixed[i],",
      "      variable = z.names[j],",
      "      t(out)",
      "    )",
      "",
      "  })",
      "",
      "})",
      "",
      "results <- do.call(",
      "  rbind,",
      "  unlist(results, recursive = FALSE)",
      ")",
      ""
    )
    
  }
  
  else if(comparison == "Overall" && movement == "Between Groups"){
    
    c(
      "# Run Overall interaction analysis",
      "",
      "tmp_fn <- function(q) {",
      "",
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "  point2 <- apply(Z_for_quants, 2, quantile, q)",
      "",
      "  newz.q1 <- rbind(point1, point2)",
      "  newz.q2 <- rbind(point1, point2)",
      "",
      "  summary.fun(",
      "    newz.q1 = newz.q1,",
      "    newz.q2 = newz.q2,",
      "    modnew.1 = modnew.1,",
      "    modnew.2 = modnew.2,",
      "    preds.fun = preds.fun",
      "  )",
      "}",
      "",
      "results <- t(sapply(qs, tmp_fn))",
      ""
    )
    
  }
  
  else if(comparison == "Single" && movement == "Between Groups"){
    
    c(
      "# Run Single-variable interaction analysis",
      "results <- lapply(which.z, function(j) {",
      "",
      "  # First comparison",
      "  q.fixed <- qs.fixed[1]",
      "",
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "  point2 <- point1",
      "",
      "  point1[j] <- quantile(Z_for_quants[, j], qs.diff[1])",
      "  point2[j] <- quantile(Z_for_quants[, j], qs.diff[2])",
      "",
      "  newz.q1 <- rbind(point1, point2)",
      "",
      "  # Second comparison",
      "  q.fixed <- qs.fixed[2]",
      "",
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "  point2 <- point1",
      "",
      "  point1[j] <- quantile(Z_for_quants[, j], qs.diff[1])",
      "  point2[j] <- quantile(Z_for_quants[, j], qs.diff[2])",
      "",
      "  newz.q2 <- rbind(point1, point2)",
      "",
      "  out <- summary.fun(",
      "    newz.q1 = newz.q1,",
      "    newz.q2 = newz.q2,",
      "    modnew.1 = modnew.1,",
      "    modnew.2 = modnew.2,",
      "    preds.fun = preds.fun",
      "  )",
      ""
    )
    
  }
  
}

build_output_formatter <- function(comparison, movement){

  if(comparison == "Overall" && movement == "Group Specific"){
    
    c(
      "# Format results",
      "results <- data.frame(",
      "  quantile = qs,",
      "  results",
      ")",
      "",
      "results"
    )
    
  }
  
  else if(comparison == "Single" && movement == "Group Specific"){
    
    c(
      "# Format results",
      "results$variable <- factor(",
      "  results$variable,",
      "  levels = z.names",
      ")",
      "",
      "results$q.fixed <- as.factor(results$q.fixed)",
      "",
      "print(results)"
    )
    
  }
  
  else if(comparison == "Overall" && movement == "Between Groups"){
    
    c(
      "# Format results",
      "results <- data.frame(",
      "  quantile = qs,",
      "  results",
      ")",
      "",
      "print(results)"
    )
    
  }
  
  else if(comparison == "Single" && movement == "Between Groups"){
    c("  data.frame(",
    "    variable = z.names[j],",
    "    t(out)",
    "  )",
    "",
    "})",
    "",
    "results <- do.call(rbind, results)",
    "")
  }
}

#needs a work beyone overall and group specific
build_plotting <- function(comparison, movement){
  if(comparison == "Overall" && movement == "Group Specific"){
    c("results$modifier <- rep(m.fixed, nrow(results))",
      "resultsfinal <- data.frame()",
      "for(m in unique(as.character(fit$modifier))){",
      "  m.fixed <- m",
      "  ",
      "  if(!is.null(m.fixed)) {",
      "    modnew <- matrix(rep(m.fixed, 2), ncol = 1)",
      "    Z_for_quants <- Z[fit$modifier == m.fixed, ]",
      "  } else {",
      "    modnew <- NULL",
      "    Z_for_quants <- Z",
      "  }",
      "  ",
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "",
      "  tmp_fn <- function(q) {",
      "    ",
      "    point2 <- apply(Z_for_quants, 2, quantile, q)",
      "    ",
      "    summary.fun(",
      "      point1 = point1,",
      "      point2 = point2,",
      "      modnew = modnew,",
      "      preds.fun = preds.fun",
      "    )",
      "  }",
      "  ",
      "  results <- t(sapply(qs, tmp_fn))",
      "  ",
      "  # Format results",
      "  results <- data.frame(",
      "    quantile = qs,",
      "    results",
      "  )",
      "  ",
      "  results$modifier <- rep(m.fixed, nrow(results))",
      "  resultsfinal <- rbind(resultsfinal, results)",
      "}",
      "",
      "plot_obj <- ggplot(resultsfinal,",
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
  
  else if(comparison == "Overall" && movement == "Between Groups"){
    c("results$modifier <- rep(m.fixed, nrow(results))",
      "#You must rerun your code, with the varying levels of m.fixed", 
      "#for example if your modifier is split by group ('group_1', 'group_2', 'group_3') and originally `m.fixed = 'group_1'`",
      "#you must rerun the code setting `m.fixed = 'group_2'` and create a new results (results_tmp_1)",
      "#repeat this process with however many modifier levels you have, changing X, m.fixed, and results_tmp_X accordingly",
      "results_tmp_X <- #put new code results here",
      "results_tmp_X$modifier <- rep(m.fixed, nrow(results_tmp_X))",
      "results <- rbind(results, results_tmp_X)",
      "",
      "plot_obj <- ggplot(results,", 
      "       aes(x=q.fixed, y = est, ymin = est - 1.96*sd,", 
      "           ymax = est + 1.96*sd, color = modifier, ",
      "           shape = modifier))+",
      "   geom_hline(yintercept=0)+",
      "   geom_pointrange(position = position_dodge(width = 0.05), size=0.5)+",
      "   facet_wrap(vars(variable))+",
      "   theme_bw()+",
      "   theme(legend.position = 'bottom',",
      "         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+",
      "   xlab('CHANGE THIS X-AXIS')+",
      "   ylab('CHANGE THIS Y-AXIS') +",
      "   labs(color='modifier_name', shape = 'modifier_name') +",
      "   scale_color_manual(values=c('color_1', ...))")
  }
  
  else if(comparison == "Single" && movement == "Group Specific"){
      c("plot_obj <- ggplot(results, aes(x=quantile, y = est,", 
      "                    ymin = est - 1.96*sd, ymax = est + 1.96*sd))+",
      "   geom_hline(yintercept=0)+",
      "   geom_pointrange(size=0.5)+",
      "   theme_bw()+",
      "   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+",
      "   xlab('CHANGE THIS X-AXIS')+",
      "   ylab('CHANGE THIS Y-AXIS')")
  }
  
  else if(comparison == "Single" && movement == "Between Groups"){
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









