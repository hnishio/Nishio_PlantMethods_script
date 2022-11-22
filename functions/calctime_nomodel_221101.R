
### Estimation of movement without model
nomodel <- function(cell_list, visual,
                    ex_sign = "negative", fold = 2, res_name, ex_name, unit1, unit2,
                    shade = TRUE, start_line = TRUE, ps = 7, theme_plot = "bw",
                    out = out){
  
  # Prepare a container for movement time
  df_mv <- data.frame(NULL)
  
  # for loop of cells
  for(i in 1:length(cell_list)){
    
    # for loop of the response variable
    for(j in 1:(ncol(cell_list[[i]])-2)){
      
      Y = cell_list[[i]][,j+2]
      
      ## Estimate the start of movement 
      ## definition1: approaching to light
      st_index <- min(which(cell_list[[i]]$ex == 1)) : max(which(cell_list[[i]]$ex == 1))
      st_index1 <- st_index[Y[st_index] > Y[st_index+1]] # decrease at the next point
      st_index1 <- st_index1[-length(st_index1)]
      
      ## definition2: moving average of differences (10 points) are negative
      diff <- diff(Y[st_index]) #difference between adjacent data
      sma <- NULL
      for(k in 1:(length(diff)-9)){
        sma[k] <- mean(diff[k:(k+9)]) # sma of next 10 differences
      }
      st_index2 <- st_index[which(sma < 0)]
      
      ## definition3: difference from the data after 5 tp is larger than fold*5/N_ex times of max - min
      N_ex <- sum(cell_list[[i]]$ex == 1)
      st_index3 <- st_index[(Y[st_index] - Y[st_index+5]) > 
                              (fold*5/N_ex)*(max(Y[st_index]) - min(Y[st_index]))]
      st_index3 <- st_index3[-((length(st_index3)-4):length(st_index3))]
      
      ## Overlap of def1, def2 and def3
      st_index1_2 <- intersect(st_index1, st_index2)
      st_index1_2_3 <- intersect(st_index1_2, st_index3)
      st_index_first <- min(st_index1_2_3, na.rm = T)
      st_index_last <- max(st_index1_2_3, na.rm = T)
      
      start_time <- cell_list[[i]]$time[st_index_first] #time at st_index_first
      end_time <- cell_list[[i]]$time[st_index_last] #time at st_index_last
      move_time <- end_time - start_time
      
      mv_time <- data.frame(start_time = start_time, end_time = end_time, move_time = move_time)
      df_mv <- rbind(df_mv, mv_time)
      
    }
  }
  
  ## Save movement time
  fwrite(df_mv, file = paste0(out, "/calctime_nomodel_mvtime.csv"))
  
}

