
### Quantile Function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.995), names = TRUE)
}


### State space model (individual model)
ssm_individual <- function(cell_list, visual, warmup=1000, sampling=1000, thin=3,
                           ex_sign = "negative", res_name, ex_name, unit1, unit2,
                           shade = TRUE, start_line = TRUE, ps = 7, theme_plot = "bw",
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
  
  # Prepare a container for movement time
  df_mv <- data.frame(NULL)
  
  # Load stan model
  model <- cmdstan_model("stan_model/individual_model_221101.stan")
  
  
  ## Execution of the Bayesian inference
  
  # for loop of cells
  for(i in 1:length(cell_list)){
    
    # for loop of response variable
    for (j in 1:(ncol(cell_list[[i]])-2)){
      
      
      ## Bayesian estimation of the states
      data_list <- list(
        N = nrow(cell_list[[i]])-1,
        N_ex = length(which(cell_list[[i]]$ex==1)),
        ex = cell_list[[i]]$ex[-1],
        Y = diff(cell_list[[i]][,j+2])
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
        output_basename = paste0("cell", i, "_", res_name, j),
        adapt_delta = 0.95,
        thin = thin
      )
      
      # 99% Bayesian credible intervals
      outcsv_name <- list.files("/tmp")
      outcsv_name <- outcsv_name[grep(paste0("cell", i, "_", res_name, j, "-"), outcsv_name)]
      tmp_csv_w <- NULL
      tmp_csv_b_ex <- NULL
      tmp_csv_alpha <- NULL
      tmp_csv_s_w <- NULL
      tmp_csv_s_b_ex <- NULL
      tmp_csv_s_Y <- NULL
      
      for(k in 1:length(outcsv_name)){
        tmp_csv <- as.data.frame(fread(cmd = paste0("grep -v '^#' ", "/tmp/", outcsv_name[k])))
        tmp_csv_w <- rbind(tmp_csv_w, tmp_csv[,str_starts(names(tmp_csv), "w")])
        tmp_csv_b_ex <- rbind(tmp_csv_b_ex, tmp_csv[,str_starts(names(tmp_csv), "b_ex")])
        tmp_csv_alpha <- rbind(tmp_csv_alpha, tmp_csv[,str_starts(names(tmp_csv), "alpha")])
        tmp_csv_s_w <- c(tmp_csv_s_w, tmp_csv[,str_starts(names(tmp_csv), "s_w")])
        tmp_csv_s_b_ex <- c(tmp_csv_s_b_ex, tmp_csv[,str_starts(names(tmp_csv), "s_b_ex")])
        tmp_csv_s_Y <- c(tmp_csv_s_Y, tmp_csv[,str_starts(names(tmp_csv), "s_Y")])
      }
      
      df_w <- as.data.frame(t(apply(tmp_csv_w, 2, quantile99)))
      df_b_ex <- as.data.frame(t(apply(tmp_csv_b_ex, 2, quantile99)))
      df_alpha <- as.data.frame(t(apply(tmp_csv_alpha, 2, quantile99)))
      df_s <- t(data.frame(s_w = quantile99(tmp_csv_s_w),
                           s_b_ex = quantile99(tmp_csv_s_b_ex),
                           s_Y = quantile99(tmp_csv_s_Y)))
      df_s <- cbind(data.frame(s_name = row.names(df_s)), df_s)
      
      # Save output
      colnames(df_w) <- paste0("w_", colnames(df_w))
      colnames(df_b_ex) <- paste0("b_ex_", colnames(df_b_ex))
      colnames(df_alpha) <- paste0("alpha_", colnames(df_alpha))
      df <- cbind(data.frame(time = cell_list[[i]]$time[-1]),
                  data.frame(Y = diff(cell_list[[i]][,j+2])),
                  df_w, df_alpha,
                  as.data.frame(rbind(matrix(NA, nrow = data_list$N - data_list$N_ex, ncol = 21), as.matrix(df_b_ex))))
      fwrite(df, file = paste0(out, "/csv/ssm_individual_cell", i, "_", res_name, j, ".csv"))
      fwrite(df_s, file = paste0(out, "/csv/ssm_individual_cell", i, "_", res_name, j, "_sd.csv"))
      
      
      ## Movement time
      ex_period <- max(cell_list[[i]]$time[cell_list[[i]]$ex == 1], na.rm = T) - min(cell_list[[i]]$time[cell_list[[i]]$ex == 1], na.rm = T) + 1
      
      if(ex_sign == "positive"){  # ex_sign == "positive"
        
        if(length(which(df$`b_ex_2.5%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_2.5%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_2.5%` > 0)][length(which(df$`b_ex_2.5%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_5%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_5%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_5%` > 0)][length(which(df$`b_ex_5%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_10%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_10%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_10%` > 0)][length(which(df$`b_ex_10%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_15%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_15%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_15%` > 0)][length(which(df$`b_ex_15%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_20%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_5%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_5%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_5%` > 0)][length(which(df$`b_ex_5%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_10%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_10%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_10%` > 0)][length(which(df$`b_ex_10%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_15%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_15%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_15%` > 0)][length(which(df$`b_ex_15%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_20%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_10%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_10%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_10%` > 0)][length(which(df$`b_ex_10%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_15%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_15%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_15%` > 0)][length(which(df$`b_ex_15%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_20%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_15%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_15%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_15%` > 0)][length(which(df$`b_ex_15%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_20%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_20%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_25%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
          move_time <- end_time - start_time
          
        }else{
          start_time <- Inf
          end_time <- Inf
          move_time <- Inf
        }
        
      }else{  # ex_sign == "negative"
        
        if(length(which(df$`b_ex_97.5%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_97.5%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_97.5%` < 0)][length(which(df$`b_ex_97.5%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_95%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_95%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_95%` < 0)][length(which(df$`b_ex_95%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_90%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_90%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_90%` < 0)][length(which(df$`b_ex_90%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_85%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_85%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_85%` < 0)][length(which(df$`b_ex_85%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_80%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_95%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_95%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_95%` < 0)][length(which(df$`b_ex_95%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_90%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_90%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_90%` < 0)][length(which(df$`b_ex_90%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_85%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_85%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_85%` < 0)][length(which(df$`b_ex_85%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_80%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_90%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_90%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_90%` < 0)][length(which(df$`b_ex_90%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_85%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_85%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_85%` < 0)][length(which(df$`b_ex_85%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_80%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_85%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_85%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_85%` < 0)][length(which(df$`b_ex_85%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_80%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_80%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/5){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }
          
        }else if(length(which(df$`b_ex_75%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
          move_time <- end_time - start_time
          
        }else{
          start_time <- Inf
          end_time <- Inf
          move_time <- Inf
        }
        
      }
      
      mv_time <- data.frame(start_time = start_time, end_time = end_time, move_time = move_time)
      df_mv <- rbind(df_mv, mv_time)
      
      
      ## Plotting
      
      # Visual
      if(is.data.frame(visual)){
        vis <- filter(visual, cell == i & index == j)$time
      }else{vis <- NULL}
      
      # Shade
      if(shade == T){
        alpha = 0.3
      }else{
        alpha = 0
      }
      
      # Start time
      if(start_line == T){
        col1 = "orange"
        col2 = "aquamarine3"
      }else{
        col1 = "transparent"
        col2 = "transparent"
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
      
      # Location of text
      text_x1 <- max(cell_list[[i]]$time) - (max(cell_list[[i]]$time) - min(cell_list[[i]]$time)) * 0.13
      text_x2 <- max(cell_list[[i]]$time) - (max(cell_list[[i]]$time) - min(cell_list[[i]]$time)) * 0.16
      
      # Distance
      ymax <- max(cell_list[[i]][,j+2])
      ymin <- min(cell_list[[i]][,j+2])
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05
      
      g_dist <- ggplot(data = cell_list[[i]]) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_line(aes(x = time, y = cell_list[[i]][,j+2]), size=0.5) +
        geom_vline(xintercept = mv_time$start_time, linetype="solid", col = col1) +
        geom_vline(xintercept = vis, linetype="dashed", col = col2) +
        annotate("text", x=text_x1, y=yceiling-yrange*0.08, label="Statistical", col=col1, size = ps/ggplot2::.pt) +
        annotate("text", x=text_x2, y=yceiling-yrange*0.2, label="Visual", col=col2, size = ps/ggplot2::.pt) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_plot2 +
        theme(legend.position = "none",
              axis.title=element_text(size = ps),
              axis.title.x=element_blank(),
              axis.text = element_text(size = ps),
              plot.title = element_text(size = ps, face = "bold")) +
        labs(title = paste("Cell ", i, ", ", res_name, " ", j, sep=""),
             y = label_y)
      
      # alpha
      ymax <- max(df$Y)
      ymin <- min(df$Y)
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05
      
      g_alpha <- ggplot(data = df, aes(x = time)) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_ribbon(aes(ymin = `alpha_2.5%`, ymax = `alpha_97.5%`), alpha = 0.5) +
        geom_line(aes(y = `alpha_50%`), size = 0.5) +
        geom_point(aes(y = Y), alpha = 0.5, size=0.5) +
        geom_vline(xintercept = mv_time$start_time, linetype="solid", col = col1) +
        geom_vline(xintercept = vis, linetype="dashed", col = col2) +
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
      ymax <- max(df$`b_ex_97.5%`, na.rm = T)
      ymin <- min(df$`b_ex_2.5%`, na.rm = T)
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05
      
      g_b_ex <- ggplot(data = df, aes(x = time)) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_ribbon(aes(ymin = `b_ex_2.5%`, ymax = `b_ex_97.5%`), alpha = 0.5) +
        geom_line(aes(y = `b_ex_50%`), size = 0.5) +
        geom_vline(xintercept = mv_time$start_time, linetype="solid", col = col1) +
        geom_vline(xintercept = vis, linetype="dashed", col = col2) +
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
      ymax <- max(df$`w_97.5%`)
      ymin <- min(df$`w_2.5%`)
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05
      
      g_w <- ggplot(data = df, aes(x = time)) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_ribbon(aes(ymin = `w_2.5%`, ymax = `w_97.5%`), alpha = 0.5) +
        geom_line(aes(y = `w_50%`), size = 0.5) +
        geom_vline(xintercept = mv_time$start_time, linetype="solid", col = col1) +
        geom_vline(xintercept = vis, linetype="dashed", col = col2) +
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
      g <- g_dist + g_alpha + g_b_ex + g_w +
        plot_layout(ncol = 1, heights = c(1, 1, 1, 1))
      suppressWarnings(
        ggsave(paste0(out, "/pdf/ssm_individual_cell", i, "_", res_name, j, ".pdf"),
               g, height = ps*20*4/4, width = ps*10*1.2, units = "mm")
      )
      
      
      ## Diagnosis of MCMC
      
      # Check of Rhat
      bayesplot::color_scheme_set("viridisC")
      bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_size = ps+2, base_family = "sans"))
      g <- bayesplot::mcmc_rhat(bayesplot::rhat(fit))
      ggsave(paste0(out, "/diagnosis/ssm_individual_rhat_cell", i, "_", res_name, j, ".pdf"),
             g, height = ps*20, width = ps*20, units = "mm")
      max_rhat <- names(which.max(bayesplot::rhat(fit)))
      
      # Confirmation of convergence
      g <- bayesplot::mcmc_combo(
        fit$draws(),
        combo = c("dens_overlay", "trace"),
        widths = c(1, 1),
        pars = c("b_ex[1]", paste0("b_ex[", data_list$N_ex, "]"), "alpha[1]", paste0("alpha[", data_list$N, "]"), max_rhat),
        gg_theme = theme_classic(base_size = ps+2)
      )
      ggsave(paste0(out, "/diagnosis/ssm_individual_combo_cell", i, "_", res_name, j, ".pdf"),
             g, height = ps*20, width = ps*20, units = "mm")
      
      
      ## Remove temporary files
      file.remove(paste0("/tmp/", outcsv_name))
      
      ## Release memory
      rm(fit, tmp_csv, tmp_csv_w, tmp_csv_b_ex, tmp_csv_alpha,
         tmp_csv_s_w, tmp_csv_s_b_ex, tmp_csv_s_Y,
         df)
      gc(reset = T);gc(reset = T)
    }
  }
  
  
  ## Save movement time
  if(is.data.frame(visual)){
    df_visual_mv <- cbind(visual, df_mv)
  }else{
    df_visual_mv <- df_mv
  }
  
  fwrite(df_visual_mv, file = paste0(out, "/csv/ssm_individual_mvtime.csv"))
}
