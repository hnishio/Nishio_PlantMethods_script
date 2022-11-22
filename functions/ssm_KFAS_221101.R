
### Quantile Function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.995), names = TRUE)
}


### State space model (Kalman filter)
ssm_KFAS <- function(cell_list, visual,
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
  
  # Prepare a container for movement time
  df_mv <- data.frame(NULL)
  df_s_all <- data.frame(NULL)
  
  
  ## Execution of the Kalman filter
  
  # for loop of cells
  for(i in 1:length(cell_list)){
    
    # for loop of the response variable
    for (j in 1:(ncol(cell_list[[i]])-2)){
      
      
      ## Estimation of the states by the Kalman filter
      
      # Definition of model
      ex = cell_list[[i]]$ex[-1]
      Y = diff(cell_list[[i]][,j+2])
      modReg <- SSModel(
        H = NA,
        as.numeric(Y) ~ 
          SSMregression(~ ex, Q = NA, a1 = 0)
      )
      modReg["P1inf", 2,2] <- 0
      modReg["P1inf", 1,1] <- 0
      
      # Estimation of parameters
      fitReg <- fitSSM(modReg, inits = c(0,0))
      
      # System noise and observation error
      df_s <- data.frame(cell = i,
                         idx = j,
                         s_b_ex = sqrt(exp(fitReg$optim.out$par)[1]),
                         s_Y = sqrt(exp(fitReg$optim.out$par)[2]))
      df_s_all <- rbind(df_s_all, df_s)
      
      # Confidence interval of state, coefficient
      interval_state_50 <- as.data.frame(predict(fitReg$model, interval = "confidence", level = 0.50))
      interval_beta_50 <- as.data.frame(predict(fitReg$model, states = "regression", interval = "confidence", level = 0.50))
      interval_state_60 <- as.data.frame(predict(fitReg$model, interval = "confidence", level = 0.60))
      interval_beta_60 <- as.data.frame(predict(fitReg$model, states = "regression", interval = "confidence", level = 0.60))
      interval_state_70 <- as.data.frame(predict(fitReg$model, interval = "confidence", level = 0.70))
      interval_beta_70 <- as.data.frame(predict(fitReg$model, states = "regression", interval = "confidence", level = 0.70))
      interval_state_80 <- as.data.frame(predict(fitReg$model, interval = "confidence", level = 0.80))
      interval_beta_80 <- as.data.frame(predict(fitReg$model, states = "regression", interval = "confidence", level = 0.80))
      interval_state_90 <- as.data.frame(predict(fitReg$model, interval = "confidence", level = 0.90))
      interval_beta_90 <- as.data.frame(predict(fitReg$model, states = "regression", interval = "confidence", level = 0.90))
      interval_state_95 <- as.data.frame(predict(fitReg$model, interval = "confidence", level = 0.95))
      interval_beta_95 <- as.data.frame(predict(fitReg$model, states = "regression", interval = "confidence", level = 0.95))
      interval_state_99 <- as.data.frame(predict(fitReg$model, interval = "confidence", level = 0.99))
      interval_beta_99 <- as.data.frame(predict(fitReg$model, states = "regression", interval = "confidence", level = 0.99))
      
      df <- as.data.frame(cbind(interval_state_99$lwr, interval_state_95$lwr, interval_state_90$lwr, interval_state_80$lwr, interval_state_70$lwr, interval_state_60$lwr, interval_state_50$lwr, interval_state_90$fit, interval_state_50$upr, interval_state_60$upr, interval_state_70$upr, interval_state_80$upr, interval_state_90$upr, interval_state_95$upr, interval_state_99$upr,
                                interval_beta_99$lwr, interval_beta_95$lwr, interval_beta_90$lwr, interval_beta_80$lwr, interval_beta_70$lwr, interval_beta_60$lwr, interval_beta_50$lwr, interval_beta_90$fit, interval_beta_50$upr, interval_beta_60$upr, interval_beta_70$upr, interval_beta_80$upr, interval_beta_90$upr, interval_beta_95$upr, interval_beta_99$upr))
      names(df) <- c("alpha_0.5%",	"alpha_2.5%",	"alpha_5%",	"alpha_10%",	"alpha_15%",	"alpha_20%",	"alpha_25%",	"alpha_50%", "alpha_75%",	"alpha_80%",	"alpha_85%",	"alpha_90%",	"alpha_95%",	"alpha_97.5%", "alpha_99.5%",
                     "b_ex_0.5%",	"b_ex_2.5%",	"b_ex_5%",	"b_ex_10%",	"b_ex_15%",	"b_ex_20%",	"b_ex_25%",	"b_ex_50%", "b_ex_75%",	"b_ex_80%",	"b_ex_85%",	"b_ex_90%",	"b_ex_95%",	"b_ex_97.5%", "b_ex_99.5%")
      
      # Save output
      df <- cbind(data.frame(time = cell_list[[i]]$time[-1],
                             Y = Y),
                  df)
      
      fwrite(df, file = paste0(out, "/csv/ssm_KFAS_cell", i, "_", res_name, j, ".csv"))
      
      
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
              axis.text = element_text(size = ps),
              plot.title = element_blank()) +
        labs(x = "Time (min)", y = label_beta)
      
      # Integrate plots
      g <- g_dist + g_alpha + g_b_ex +
        plot_layout(ncol = 1, heights = c(1, 1, 1))
      suppressWarnings(
        ggsave(paste0(out, "/pdf/ssm_KFAS_cell", i, "_", res_name, j, ".pdf"),
               g, height = ps*20*3/4, width = ps*10*1.2, units = "mm")
      )
      
    }
  }
  
  
  ## Save movement time
  if(is.data.frame(visual)){
    df_visual_mv <- cbind(visual, df_mv)
  }else{
    df_visual_mv <- df_mv
  }
  
  fwrite(df_visual_mv, file = paste0(out, "/csv/ssm_KFAS_mvtime.csv"))
  
  
  ## Save sd of system noise and observation error
  fwrite(df_s_all, file = paste0(out, "/csv/ssm_KFAS_sd.csv"))
  
}
