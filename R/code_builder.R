#' Creates a text file of plotting code 
#' 
#' Takes a specified exposure type, and analysis and returns copy and pasteable code text, which when executed will return results and plot objects
#' @param fit_name The name of a current fitted object in the environment
#' @param exposure exposure type for what change in exposure is being measured. Options are 'Overall', 'Single' and 'Bivariate'; defaults to 'Overall'
#' @param analysis analysis type, describes what difference we are looking at. Options are 'Group Specific', 'Between Groups', and 'Response Function'; defaults to 'Group Specific'
#' @param centered defaults to NULL
#' @param sel defaults to 'NULL'
#' @param m.fixed defaults to NULL
#' @param qs defaults to NULL
#' @param q.fixed defaults to NULL
#' @param qs.diff defaults to NULL
#' @param mod.diff defaults to NULL
#' @param ngrid defaults to NULL
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
    ngrid = NULL){
  
  
  code <- c(
    build_header(exposure, analysis),
    build_data_extract(fit_name),
    build_user_changeable_data(exposure, analysis, m.fixed, qs, q.fixed,
                               qs.diff, mod.diff, ngrid),
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
                                       qs.diff, mod.diff, ngrid){
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
                    paste0("q.fixed <- ", as.character(q.fixed), " ")))
    
    if(analysis == "Response Function"){
      ret <- c(ret, c("",
                      paste0("ngrid <- ", as.character(ngrid), ""),
                      "",
                      "centered <- TRUE"))
    }
    
    else{
      ret <- c(ret, c("",
                    paste0("qs.diff <-", as.character(qs.diff), " ")))
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
  
  if(analysis == "Group Specific"){
    
    if(exposure == "Single"){
      ret <- c(ret, c(
              "# Single-variable Overall settings",
              "Z <- fit$Z",
              "which.z <- 1:ncol(Z)",
              "z.names <- colnames(Z)",
              ""))
    }
    
  }
  
  if(analysis == "Response Function"){
    if(exposure == "Single"){
      ret <- c(ret,c("Z <- fit$Z",
                     "modifier <- fit$modifier",
                     "z.names <- paste0('z', 1:ncol(Z))",
                     "which.mod <- levels(modifier)",
                     "which.z = 1:ncol(Z)")
      )
    }
  }
  
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
        "Z <- fit$Z",
        "which.z <- 1:ncol(Z)",
        "z.names <- colnames(Z)"
      ))
      
    }
    
    
    ret <- c(ret, c("# Restrict exposure support",
                    "Z_for_quants <- Z_support(",
                    "  Z = fit$Z,",
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
    
    c("Z <- fit$Z",
      "if(m.fixed == 'All'){",
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
      "      point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
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
      "        q.fixed = q.fixed,",
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
      "  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)",
      "  point2 <- point1",
      "  point1[j] <- quantile(Z_for_quants[, j], qs.diff[1])",
      "  point2[j] <- quantile(Z_for_quants[, j], qs.diff[2])",
      "  newz.q1 <- rbind(point1, point2)",
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
  
  else if(exposure == "Single" && analysis == "Response Function"){
    c("results <- dplyr::tibble()",
      "for(i in which.z) {",
      "  whichz <- i",
      "  colnames(Z) <- paste0('z', 1:ncol(Z))",
      "  Z_for_quants <- Z_support(Z = Z, mod.diff = levels(modifier), modifier = modifier)",
      "  ord <- c(whichz, setdiff(1:ncol(Z_for_quants), whichz))",
      "  z1 <- seq(min(Z_for_quants[,ord[1]]), max(Z_for_quants[,ord[1]]), length = ngrid)",
      "  z.others <- lapply(2:ncol(Z_for_quants), function(x) quantile(Z[,ord[x]], q.fixed))",
      "  z.all <- c(list(z1), z.others)",
      "  newz.grid <- expand.grid(z.all)",
      "  colnames(newz.grid) <- colnames(Z)[ord]",
      "  newz.grid <- newz.grid[,colnames(Z)]",
      "  if(!is.null(which.mod)){",
      "    mod_new <- rep(which.mod, each = nrow(newz.grid))",
      "    if(length(which.mod) > 1){",
      "      newz.grid_tmp <- newz.grid",
      "      z1_tmp <- z1",
      "      for(k in 2:length(which.mod)){",
      "        newz.grid <- rbind(newz.grid, newz.grid_tmp)",
      "        z1 <- c(z1, z1_tmp)",
      "      }",
      "    }",
      "  }",
      "  mindists <- rep(NA,nrow(newz.grid))",
      "  for (j in seq_along(mindists)) {",
      "    pt <- as.numeric(newz.grid[j, colnames(Z)[ord[1]]])",
      "    dists <- fields::rdist(matrix(pt, nrow = 1), Z_for_quants[, colnames(Z)[ord[1]]])",
      "    mindists[j] <- min(dists)",
      "  }",
      "  preds <- preds.fun(znew = newz.grid, modnew = mod_new)",
      "  preds.plot <- preds$postmean",
      "  se.plot <- sqrt(diag(preds$postvar))",
      "  if(centered) {",
      "    preds.plot <- preds.plot - mean(preds.plot)", 
      "  }",
      "  res <- dplyr::tibble(z = z1, modifier = mod_new, est = preds.plot, se = se.plot)",
      "  df0 <- dplyr::mutate(res, variable = z.names[i]) %>% dplyr::select_at(c('variable', 'z', 'modifier', 'est', 'se'))",
      "  results <- dplyr::bind_rows(results, df0)",
      "}",
      "results$variable <- factor(results$variable, levels = z.names[which.z])")
  }
  
}

#Need to change this to only have the results object returned
build_output_formatter <- function(exposure, analysis){

  if(exposure == "Overall" && analysis == "Between Groups"){
    
    c("results <- data.frame(",
      "  quantile = qs,",
      "  results",
      ")",
      "results"
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
  
  else{
    
    c("results")
    
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
  
  else if(exposure == "Single" && analysis == "Response Function"){
    c(
      "plot_obj <- ggplot(results, aes(z, est, ymin = est - 1.96*se, ",
      "                             ymax = est + 1.96*se, color = modifier, fill = modifier, linetype = modifier)) + ",
      "  geom_smooth(stat = 'identity') + ",
      "  geom_hline(yintercept=0)+",
      "  facet_wrap(vars(variable)) +",
      "  xlab('Xlab') +",
      "  ylab('Ylab') +",
      "  theme_classic() +",
      "  theme(legend.position = 'bottom')+",
      "  labs(color='modifier_name', fill = 'modifier_name', linetype = 'modifier_name')"
    )
  }
}



#' Creates a text file of interpretation code 
#' 
#' Takes a specified exposure type, and analysis and returns copy and pasteable code text, which when executed will return results and plot objects
#' @param exposure exposure type for what change in exposure is being measured. Options are 'Overall', 'Single' and 'Bivariate'; defaults to 'Overall'
#' @param analysis analysis type, describes what difference we are looking at. Options are 'Group Specific', 'Between Groups', and 'Response Function'; defaults to 'Group Specific'
#' 
#' @export
generate_interpretation <- function(exposure, analysis){
  
  ret <- c()

  if(exposure == "Overall" && analysis == "Group Specific"){
    ret <- c(ret, c("For a fixed modifier group \\(w\\), it evaluates the change in the outcome when the entire exposure mixture is set to a specified quantile level, relative to a reference quantile level. Mathematically this corresponds to:",
                    ""  ,
                    "$$h_{w}(Z^{qs}) - h_{w}(Z^{q.fixed})$$",
                    ""  ,
                    "where \\(h_w (Z)\\) denotes the group-specific exposure-response function,\\(Z^{qs}\\) represents all exposures jointly set to the specified quantiles, and\\(Z^{q.fixed}\\) represents the exposures set to the reference quantiles. ",
                    "",
                    "The estimated value (`est`) can be interpreted as the expected change in the outcome for modifier group \\(w\\) when the exposure mixture shifts jointly from the reference quantile to the specified quantile. Graphically, we see these estimates as a point, surrounded by a Wald confidence interval. "))
  }
  
  if(exposure == "Single" && analysis == "Group Specific"){
    ret <- c(ret, c("For a fixed modifier group \\(w\\), it evaluates the change in the outcome when the a single exposure is set to a specified quantile level, relative to all other exposures being set to a reference quantile level. Mathematically this corresponds to:",
                    ""  ,
                    "$$h_{w}(Z_1^{q.fixed}, . . . , Z_k^{qs.diff[1]},. . .,Z_n^{q.fixed}) - h_{w}(Z_1^{q.fixed}, . . ., Z_k^{qs.diff[2]},. . . ,Z_n^{q.fixed})$$",
                    ""  ,
                    "where \\(h_w (Z)\\) denotes the group-specific exposure-response function, \\(Z_k^{qs.diff}\\) represents a specific exposure (\\k\\) set to the specified quantiles, and \\(Z_k^{q.fixed}\\) represents a specific exposure (\\k\\) set to the reference quantiles. ",
                    "",
                    "The estimated value (`est`) can be interpreted as the expected changein the outcome for modifier group \\(w\\) when the isolated exposure shiftsfrom the first difference quantile to the second difference quantile, while holding every other exposure at the reference quantile. Graphically,we see these estimates as a point, surrounded by a Wald confidence interval."))
  }
  
  if(exposure == "Overall" && analysis == "Between Groups"){
    ret <- c(ret, c("This function evaluates the overall mixture effect differs between modifier groups when the entire exposure mixture is shifted from a reference quantile level to a specified quantile level. Mathematically this corresponds to:",
                    ""  ,
                    "$$(h_{w_1}(Z^{qs}) - h_{w_1}(Z^{q.fixed})) - (h_{w_2}(Z^{qs}) - h_{w_2}(Z^{q.fixed}))$$",
                    ""  ,
                    "where \\(h_w (Z)\\) denotes the group-specific exposure-response function jointly set to the specified quantile level, and \\(Z^{qs}\\) represents all exposures jointly set to the specified quantiles, and \\(Z^{q.fixed}\\) represents the exposures set to the reference quantiles.",
                    "",
                    "The estimated effect (`est`) measures the difference in overall mixture effects between modifier groups  \\(w_1\\) and \\(w_2\\). Specifically, it represents how much larger the expected change in the outcome associated with shifting the entire exposure mixture from \\(q.fixed\\) to \\(qs\\) is in group \\(w_1\\) compared with group \\(w_2\\)"
    ))
  }

  
  if(exposure == "Single" && analysis == "Between Groups"){
    ret <- c(ret, c("This function evaluates how the effect of a single exposure differs between modifier groups while holding all other exposures fixed at a reference quantile level. Mathematically, this corresponds to",
                    "",
                    "$$(h_{w_1}(Z_1^{q.fixed}, . . . , Z_k^{qs.diff[1]},...,Z_n^{q.fixed}) - h_{w_1}(Z_1^{q.fixed}, . . ., Z_k^{qs.diff[1]},...,Z_n^{q.fixed}) - $$",
                    "$$(h_{w_2}(Z_1^{q.fixed}, . . ., Z_k^{qs.diff[1]},...,Z_n^{q.fixed}) - h_{w_2}(Z_1^{q.fixed}, . . ., Z_k^{qs.diff[1]},...,Z_n^{q.fixed})).$$",
                    "" ,
                    "Here, \\(h_w(Z)\\) denotes the group-specific exposure-response function, \\(z_m^{qs_1}\\) and \\(z_m^{qs_2}\\) represent the selected exposure set to the lower and upper quantiles specified by `qs.diff`, and \\(Z_{-m}^{q.fixed}\\) represents all remaining exposures held fixed at the reference quantile level.",
                    "",
                    "The estimated effect (`est`) measures the difference in the single-exposure effect between modifier groups \\(w_1\\) and \\(w_2\\). Specifically, it represents how much larger (or smaller) the expected change in the outcome associated with increasing the selected exposure from \\(qs_1\\) to \\(qs_2\\) is in group \\(w_1\\) compared with group \\(w_2\\), while all other exposures remain fixed.",
                    ""))
  }

  if(exposure == "Single" && analysis == "Response Function"){
    ret <- c(ret, c(
      "This function estimates the exposure-response relationship for a single exposure within each modifier group while holding all remaining exposures fixed at a reference quantile level. Mathematically, this corresponds to:",
      "",
      "$$h_w(Z_1^{q.fixed},...,Z_k,..., Z_m^{q.fixed})$$",
      "",
      "where \\(h_w(Z)\\) denotes the group-specific exposure-response function, \\(Z_k\\) is the exposure being varied across its observed range, and \\(Z_{m}^{q.fixed}\\) represents all remaining exposures held fixed at the reference quantile level.",
      "",
      "The resulting curves describe the expected outcome as the selected exposure changes while all other exposures remain fixed. Separate curves are estimated for each modifier group, allowing direct visualization of differences in exposure-response relationships across groups.",
      "",
      "Differences in the shape, slope, or magnitude of the curves suggest that the association between the selected exposure and the outcome varies across modifier groups."
    ))
    
  }
  
  paste(ret, collapse = "<br/>")
  
}

