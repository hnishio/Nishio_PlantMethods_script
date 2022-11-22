
##### State space model (individual model) #####

### Quantile Function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.995), names = TRUE)
}



### State space model with Bayesian inference
calctime_ssm_individual <- function(cell_list, visual, warmup=1000, sampling=1000, thin=3,
                                    ex_sign = "negative", res_name, ex_name, unit1, unit2,
                                    shade = TRUE, start_line = TRUE, ps = 7, theme_plot = "bw",
                                    out = out){
  
  ## Set up
  
  # Prepare a container for movement time
  df_mv <- data.frame(NULL)
  
  # Load stan model
  model <- cmdstan_model("stan_model/Individual_model_221101.stan")
  
  
  ## Execution of the model
  
  # for loop of cells
  for(i in 1:length(cell_list)){
    
    # for loop of chloroplasts
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
        output_basename = paste0("cell", i, "_", "chl", j),
        adapt_delta = 0.95,
        thin = thin
      )
      
      # 99% Bayesian credible intervals
      outcsv_name <- list.files("/tmp")
      outcsv_name <- outcsv_name[grep(paste0("cell", i, "_", "chl", j, "-"), outcsv_name)]
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
      
      # Save output
      colnames(df_w) <- paste0("w_", colnames(df_w))
      colnames(df_b_ex) <- paste0("b_ex_", colnames(df_b_ex))
      colnames(df_alpha) <- paste0("alpha_", colnames(df_alpha))
      df <- cbind(data.frame(time = cell_list[[i]]$time[-1]),
                  data.frame(Y = diff(cell_list[[i]][,j+2])),
                  df_w, df_alpha,
                  as.data.frame(rbind(matrix(NA, nrow = data_list$N - data_list$N_ex, ncol = 21), as.matrix(df_b_ex))))
      
      
      ## Movement time
      ex_period <- max(cell_list[[i]]$time[cell_list[[i]]$ex == 1], na.rm = T) - min(cell_list[[i]]$time[cell_list[[i]]$ex == 1], na.rm = T) + 1
      
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
      
      mv_time <- data.frame(start_time = start_time, end_time = end_time, move_time = move_time)
      df_mv <- rbind(df_mv, mv_time)
      
      
      ## Release memory
      file.remove(paste0("/tmp/", outcsv_name))
      rm(fit, tmp_csv, tmp_csv_w, tmp_csv_b_ex, tmp_csv_alpha,
         tmp_csv_s_w, tmp_csv_s_b_ex, tmp_csv_s_Y,
         df)
      gc(reset = T);gc(reset = T)
    }
  }
  
  
  ## Save movement time
  fwrite(df_mv, file = paste0(out, "/calctime_ssm_individual_mvtime.csv"))
}
