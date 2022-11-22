
### Quantile Function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.995), names = TRUE)
}


### State space model (individual model)
calctime_ssm_KFAS <- function(cell_list, visual,
                              ex_sign = "negative", res_name, ex_name, unit1, unit2,
                              shade = TRUE, start_line = TRUE, ps = 7, theme_plot = "bw",
                              out = out){
  
  ## Set up
  
  # Prepare a container for movement time
  df_mv <- data.frame(NULL)
  
  
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
      
    }
  }
  
  
  ## Save movement time
  fwrite(df_mv, file = paste0(out, "/calctime_ssm_KFAS_mvtime.csv"))
  
}
