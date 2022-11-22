
### Quantile Function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.995), names = TRUE)
}


### State space model (common model)
ssm_common <- function(cell_list, mvtime, warmup=1000, sampling=1000, thin=3,
                       res_name, ex_name, unit1, unit2,
                       shade = TRUE, ps = 7, theme_plot = "bw",
                       out = out){
  
  ## Set up
  
  # Create output directories
  if(file.exists(paste0(out, "/csv"))==F){
    dir.create(paste0(out, "/csv"), recursive=T)
  }
  if(file.exists(paste0(out, "/pdf"))==F){
    dir.create(paste0(out, "/pdf"), recursive=T)
  }
  if(file.exists(paste0(out, "/diagnosis"))==F){
    dir.create(paste0(out, "/diagnosis"), recursive=T)
  }
  
  # Adjust data.frame
  mvtime <- mvtime[,1:4]
  names(mvtime)[1:4] <- c("cell", "each", "visual", "predicted")
  distance <- NULL
  for(i in 1:length(cell_list)){
    distance <- c(distance, as.numeric(cell_list[[i]][max(which(cell_list[[i]]$ex == 0)),-(1:2)]))
  }
  mvtime$distance <- distance
  
  # Remove infinity from cell_list and mvtime
  null_cell <- mvtime[is.infinite(rowSums(mvtime)),]$cell
  null_each <- mvtime[is.infinite(rowSums(mvtime)),]$each
  if(length(null_cell) > 0){
    for(i in 1:length(null_cell)){
      cell_list[[null_cell[i]]] <- cell_list[[null_cell[i]]][,-(null_each[i]+2)]
      mvtime <- mvtime[!is.infinite(rowSums(mvtime)),]
    }
  }
  
  # Load stan model
  model <- cmdstan_model("stan_model/common_model_221101.stan")
  
  
  ## Execution of the Bayesian inference
  
  # for loop of cells
  for(i in 1:length(cell_list)){
    
    ## Bayesian estimation of the states
    data_list <- list(
      N = nrow(cell_list[[i]])-1,
      N_ex = length(which(cell_list[[i]]$ex==1)),
      N_each = ncol(cell_list[[i]])-2,
      ex = cell_list[[i]]$ex[-1],
      Y = apply(cell_list[[i]][,-(1:2)], 2, diff),
      start = mvtime$predicted[mvtime$cell==i]
    )
    
    # Execute MCMC
    fit <- model$sample(
      data = data_list,
      seed = 1,
      iter_warmup = warmup*thin,
      iter_sampling = sampling*thin,
      chains = 4,
      parallel_chains = 4,
      refresh = floor(warmup/2.5*thin),
      #show_messages = F,
      #sig_figs = 4,
      output_dir = "/tmp",
      output_basename = paste0("cell", i),
      adapt_delta = 0.95,
      thin = thin
    )
    
    # 99% Bayesian credible intervals
    outcsv_name <- list.files("/tmp")
    outcsv_name <- outcsv_name[grep(paste0("cell", i), outcsv_name)]
    tmp_csv_w <- NULL
    tmp_csv_b_ex <- NULL
    tmp_csv_alpha <- NULL
    tmp_csv_b_ex_each <- NULL
    tmp_csv_alpha_each <- NULL
    tmp_csv_dist <- NULL
    tmp_csv_s_w <- NULL
    tmp_csv_s_b_ex <- NULL
    tmp_csv_s_Y <- NULL
    
    for(k in 1:length(outcsv_name)){
      tmp_csv <- as.data.frame(fread(cmd = paste0("grep -v '^#' ", "/tmp/", outcsv_name[k])))
      tmp_csv_w <- rbind(tmp_csv_w, tmp_csv[,str_starts(names(tmp_csv), "w\\.")])
      tmp_csv_b_ex <- rbind(tmp_csv_b_ex, tmp_csv[,str_starts(names(tmp_csv), "b_ex\\.")])
      tmp_csv_alpha <- rbind(tmp_csv_alpha, tmp_csv[,str_starts(names(tmp_csv), "alpha\\.")])
      tmp_csv_b_ex_each <- rbind(tmp_csv_b_ex_each, tmp_csv[,str_starts(names(tmp_csv), "b_ex_each")])
      tmp_csv_alpha_each <- rbind(tmp_csv_alpha_each, tmp_csv[,str_starts(names(tmp_csv), "alpha_each")])
      tmp_csv_dist <- rbind(tmp_csv_dist, tmp_csv[,str_starts(names(tmp_csv), "dist")])
      tmp_csv_s_w <- c(tmp_csv_s_w, tmp_csv[,str_starts(names(tmp_csv), "s_w")])
      tmp_csv_s_b_ex <- c(tmp_csv_s_b_ex, tmp_csv[,str_starts(names(tmp_csv), "s_b_ex")])
      tmp_csv_s_Y <- c(tmp_csv_s_Y, tmp_csv[,str_starts(names(tmp_csv), "s_Y")])
    }
    
    df_w <- as.data.frame(t(apply(tmp_csv_w, 2, quantile99)))
    df_b_ex <- as.data.frame(t(apply(tmp_csv_b_ex, 2, quantile99)))
    df_alpha <- as.data.frame(t(apply(tmp_csv_alpha, 2, quantile99)))
    df_b_ex_each <- as.data.frame(t(apply(tmp_csv_b_ex_each, 2, quantile99)))
    df_alpha_each <- as.data.frame(t(apply(tmp_csv_alpha_each, 2, quantile99)))
    df_dist <- as.data.frame(t(apply(tmp_csv_dist, 2, quantile99)))
    df_s <- t(data.frame(s_w = quantile99(tmp_csv_s_w),
                         s_b_ex = quantile99(tmp_csv_s_b_ex),
                         s_Y = quantile99(tmp_csv_s_Y)))
    df_s <- cbind(data.frame(s_name = row.names(df_s)), df_s)
    
    # Save output
    colnames(df_w) <- paste0("w_", colnames(df_w))
    colnames(df_b_ex) <- paste0("b_ex_", colnames(df_b_ex))
    colnames(df_alpha) <- paste0("alpha_", colnames(df_alpha))
    colnames(df_b_ex_each) <- paste0("b_ex_each_", colnames(df_b_ex_each))
    colnames(df_alpha_each) <- paste0("alpha_each_", colnames(df_alpha_each))
    colnames(df_dist) <- paste0("dist_", colnames(df_dist))
    df_b_ex <- as.data.frame(rbind(matrix(NA, nrow = data_list$N - data_list$N_ex, ncol = 21), as.matrix(df_b_ex)))
    
    for(l in 1: data_list$N_each){
      eval(parse(text = paste0(
        "df_alpha_each", l, " <- df_alpha_each[str_ends(row.names(df_alpha_each), '\\\\.", l, "'),]
           df_b_ex_each", l, " <- df_b_ex_each[str_ends(row.names(df_b_ex_each), '\\\\.", l, "'),]
           colnames(df_alpha_each", l, ") <- paste0(colnames(df_alpha_each", l, "), '_", l, "')
           colnames(df_b_ex_each", l, ") <- paste0(colnames(df_b_ex_each", l, "), '_", l, "')
           df_b_ex_each", l, " <- as.data.frame(rbind(matrix(NA, nrow = data_list$N - data_list$N_ex, ncol = 21), as.matrix(df_b_ex_each", l, ")))"
      )))
    }
    
    dfs <- cbind(data.frame(time = cell_list[[i]]$time[-1]),
                 data.frame(Y = data_list$Y),
                 df_w, df_alpha, df_dist, df_b_ex)
    
    for(l in 1: data_list$N_each){
      eval(parse(text = paste0(
        "dfs <- cbind(dfs, df_alpha_each", l, ")
           dfs <- cbind(dfs, df_b_ex_each", l, ")"
      )))
    }
    
    fwrite(dfs, file = paste0(out, "/csv/ssm_common_cell", i, ".csv"))
    fwrite(df_s, file = paste0(out, "/csv/ssm_common_cell", i, "_sd.csv"))
    
    
    ## Diagnosis of MCMC
    
    # Check of Rhat
    bayesplot::color_scheme_set("viridisC")
    bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_size = ps+2, base_family = "sans"))
    suppressWarnings(g <- bayesplot::mcmc_rhat(bayesplot::rhat(fit)))
    suppressWarnings(ggsave(paste0(out, "/diagnosis/ssm_common_rhat_cell", i, ".pdf"),
                            g, height = ps*20, width = ps*20, units = "mm"))
    max_rhat <- names(which.max(bayesplot::rhat(fit)))
    
    # Confirmation of convergence
    g <- bayesplot::mcmc_combo(
      fit$draws(),
      combo = c("dens_overlay", "trace"),
      widths = c(1, 1),
      pars = c("b_ex[1]", paste0("b_ex[", data_list$N_ex, "]"), "alpha[1]", paste0("alpha[", data_list$N, "]"), max_rhat),
      gg_theme = theme_classic(base_size = ps+2)
    )
    ggsave(paste0(out, "/diagnosis/ssm_common_combo_cell", i, ".pdf"),
           g, height = ps*20, width = ps*20, units = "mm")
    
    
    # Remove temporary files
    file.remove(paste0("/tmp/", outcsv_name))
    rm(fit, tmp_csv, tmp_csv_w, tmp_csv_b_ex, tmp_csv_alpha, tmp_csv_dist, 
       tmp_csv_b_ex_each, tmp_csv_alpha_each,
       tmp_csv_s_w, tmp_csv_s_b_ex, tmp_csv_s_Y,
       dfs)
  }
  
  
  ## Y range of plots
  ymax_dist_all <- NULL
  ymin_dist_all <- NULL
  ymax_alpha_all <- NULL
  ymin_alpha_all <- NULL
  ymax_b_ex_all <- NULL
  ymin_b_ex_all <- NULL
  ymax_w_all <- NULL
  ymin_w_all <- NULL
  
  # for loop of cells
  for(i in 1:length(cell_list)){
    
    # Load data
    dfs <- fread(file = paste0(out, "/csv/ssm_common_cell", i, ".csv"))
    
    ymax_dist <- max(dfs$`dist_97.5%`)
    ymax_dist_all <- cbind(ymax_dist_all, ymax_dist)
    ymin_dist <- min(dfs$`dist_2.5%`)
    ymin_dist_all <- cbind(ymin_dist_all, ymin_dist)
    
    ymax_alpha <- max(dfs$`alpha_97.5%`)
    ymax_alpha_all <- cbind(ymax_alpha_all, ymax_alpha)
    ymin_alpha <- min(dfs$`alpha_2.5%`)
    ymin_alpha_all <- cbind(ymin_alpha_all, ymin_alpha)
    
    ymax_b_ex <- max(dfs$`b_ex_97.5%`, na.rm = T)
    ymax_b_ex_all <- cbind(ymax_b_ex_all, ymax_b_ex)
    ymin_b_ex <- min(dfs$`b_ex_2.5%`, na.rm = T)
    ymin_b_ex_all <- cbind(ymin_b_ex_all, ymin_b_ex)
    
    ymax_w <- max(dfs$`w_97.5%`)
    ymax_w_all <- cbind(ymax_w_all, ymax_w)
    ymin_w <- min(dfs$`w_2.5%`)
    ymin_w_all <- cbind(ymin_w_all, ymin_w)
  }
  
  
  ## Plotting
  
  # for loop of cells
  for(i in 1:length(cell_list)){
    
    # Load data
    dfs <- fread(file = paste0(out, "/csv/ssm_common_cell", i, ".csv"))
    
    # Shade
    if(shade == T){
      alpha = 0.3
    }else{
      alpha = 0
    }
    
    # Theme
    if(theme_plot == "bw"){
      theme_plot2 <- theme_bw(base_size = ps)
    }else if(theme_plot == "light"){
      theme_plot2 <- theme_light(base_size = ps)
    }else if(theme_plot == "classic"){
      theme_plot2 <- theme_classic(base_size = ps)
    }else if(theme_plot == "gray"){
      theme_plot2 <- theme_gray(base_size = ps)
    }else if(theme_plot == "dark"){
      theme_plot2 <- theme_dark(base_size = ps)
    }else if(theme_plot == "test"){
      theme_plot2 <- theme_test(base_size = ps)
    }else if(theme_plot == "minimal"){
      theme_plot2 <- theme_minimal(base_size = ps)
    }else if(theme_plot == "void"){
      theme_plot2 <- theme_void(base_size = ps)
    }
    
    # label
    if(unit1=="meter"){
      label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (m))))
      label_alpha <- bquote(atop("Velocity of movement", (m/.(unit2))))
      label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (m/.(unit2))))
      label_random <- bquote(atop("Random fluctuation", (m/.(unit2))))
    }else if(unit1=="centimeter"){
      label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (cm))))
      label_alpha <- bquote(atop("Velocity of movement", (cm/.(unit2))))
      label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (cm/.(unit2))))
      label_random <- bquote(atop("Random fluctuation", (cm/.(unit2))))
    }else if(unit1=="millimeter"){
      label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (mm))))
      label_alpha <- bquote(atop("Velocity of movement", (mm/.(unit2))))
      label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (mm/.(unit2))))
      label_random <- bquote(atop("Random fluctuation", (mm/.(unit2))))
    }else if(unit1=="micrometer"){
      label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (mu*m))))
      label_alpha <- bquote(atop("Velocity of movement", (mu*m/.(unit2))))
      label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (mu*m/.(unit2))))
      label_random <- bquote(atop("Random fluctuation", (mu*m/.(unit2))))
    }else if(unit1=="nanometer"){
      label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (nm))))
      label_alpha <- bquote(atop("Velocity of movement", (nm/.(unit2))))
      label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (nm/.(unit2))))
      label_random <- bquote(atop("Random fluctuation", (nm/.(unit2))))
    }else{
      label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " (",  .(unit1), ")")))
      label_alpha <- bquote(atop("Velocity of movement", (.(unit1)/.(unit2))))
      label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), paste("(", .(unit1), " / ", .(unit2), ")")))
      label_random <- bquote(atop("Random fluctuation", paste("(", .(unit1), " / ", .(unit2), ")")))
    }
    
    # Shade of the x-axis min and max
    shade_xmin <- min(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
    shade_xmax <- max(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
    
    # distance
    ymax_dist <- max(ymax_dist_all)
    ymin_dist <- min(ymin_dist_all)
    yrange_dist <- (ymax_dist - ymin_dist)
    yceiling_dist <-  ymax_dist + yrange_dist * 0.05
    yfloor_dist <- ymin_dist - yrange_dist * 0.05
    
    g_dist <- ggplot(data = dfs, aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = yfloor_dist, ymax = yceiling_dist, alpha = alpha, fill = "gray50") +
      geom_ribbon(aes(ymin = `dist_2.5%`, ymax = `dist_97.5%`), alpha = 0.5) +
      geom_line(aes(y = `dist_50%`), size=0.5) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_plot2 +
      theme(legend.position = "none",
            axis.title=element_text(size = ps),
            axis.title.x=element_blank(),
            axis.text = element_text(size = ps),
            plot.title = element_text(size = ps, face = "bold")) +
      labs(title = paste("Cell ", i, sep=""),
           y = label_y)
    
    # Velocity
    ymax_alpha <- max(ymax_alpha_all)
    ymin_alpha <- min(ymin_alpha_all)
    yrange_alpha <- (ymax_alpha - ymin_alpha)
    yceiling_alpha <-  ymax_alpha + yrange_alpha * 0.05
    yfloor_alpha <- ymin_alpha - yrange_alpha * 0.05
    
    g_velocity <- ggplot(data = dfs, aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = yfloor_alpha, ymax = yceiling_alpha, alpha = alpha, fill = "gray50") +
      geom_ribbon(aes(ymin = `alpha_2.5%`, ymax = `alpha_97.5%`), alpha = 0.5) +
      geom_line(aes(y = `alpha_50%`), size = 0.5) +
      geom_hline(yintercept = 0, linetype="dashed") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_plot2 +
      theme(legend.position = "none",
            axis.title=element_text(size = ps),
            axis.title.x=element_blank(),
            axis.text = element_text(size = ps),
            plot.title = element_blank()) +
      labs(y = label_alpha)
    
    # beta_ex
    ymax_b <- max(ymax_b_ex_all)
    ymin_b <- min(ymin_b_ex_all)
    yrange_b <- (ymax_b - ymin_b)
    yceiling_b <-  ymax_b + yrange_b * 0.05
    yfloor_b <- ymin_b - yrange_b * 0.05
    
    g_b_ex <- ggplot(data = dfs, aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = yfloor_b, ymax = yceiling_b, alpha = alpha, fill = "gray50") +
      geom_ribbon(aes(ymin = `b_ex_2.5%`, ymax = `b_ex_97.5%`), alpha = 0.5) +
      geom_line(aes(y = `b_ex_50%`), size = 0.5) +
      geom_hline(yintercept = 0, linetype="dashed") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_plot2 +
      theme(legend.position = "none",
            axis.title=element_text(size = ps),
            axis.title.x=element_blank(),
            axis.text = element_text(size = ps),
            plot.title = element_blank()) +
      labs(y = label_beta)
    
    # Random movement
    ymax_w <- max(ymax_w_all)
    ymin_w <- min(ymin_w_all)
    yrange_w <- (ymax_w - ymin_w)
    yceiling_w <-  ymax_w + yrange_w * 0.05
    yfloor_w <- ymin_w - yrange_w * 0.05
    
    g_w <- ggplot(data = dfs, aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = yfloor_w, ymax = yceiling_w, alpha = alpha, fill = "gray50") +
      geom_ribbon(aes(ymin = `w_2.5%`, ymax = `w_97.5%`), alpha = 0.5) +
      geom_line(aes(y = `w_50%`), size = 0.5) +
      geom_hline(yintercept = 0, linetype="dashed") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_plot2 +
      theme(legend.position = "none",
            axis.title=element_text(size = ps),
            axis.text = element_text(size = ps),
            plot.title = element_blank()) +
      labs(x = "Time (min)",
           y = label_random)
    
    # Integrate plots
    g <- g_dist + g_velocity + g_b_ex + g_w +
      plot_layout(ncol = 1, heights = c(1, 1, 1, 1))
    suppressWarnings(
      ggsave(paste0(out, "/pdf/ssm_common_cell", i, ".pdf"),
             g, height = ps*20*4/4, width = ps*10*1.2, units = "mm")
    )
    
    
    ## Release memory
    gc(reset = T);gc(reset = T)
  }
  
}
