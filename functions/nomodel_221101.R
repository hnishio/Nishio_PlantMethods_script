
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
      
      if(ex_sign == "positive"){  # positive
        
        ## Estimate the start of movement 
        ## definition1: approaching to light
        st_index <- min(which(cell_list[[i]]$ex == 1)) : max(which(cell_list[[i]]$ex == 1))
        st_index1 <- st_index[Y[st_index] < Y[st_index+1]] # increase at the next point
        st_index1 <- st_index1[-length(st_index1)]
        
        ## definition2: moving average of differences (10 points) are positive
        diff <- diff(Y[st_index]) #difference between adjacent data
        sma <- NULL
        for(k in 1:(length(diff)-9)){
          sma[k] <- mean(diff[k:(k+9)]) # sma of next 10 differences
        }
        st_index2 <- st_index[which(sma > 0)]
        
        ## definition3: difference from the data after 5 tp is larger than fold*5/N_ex times of max - min
        N_ex <- sum(cell_list[[i]]$ex == 1)
        st_index3 <- st_index[(Y[st_index+5] - Y[st_index]) > 
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
        
      }else{  # negative
        
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
      }else if(unit1=="centimeter"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (cm))))
      }else if(unit1=="millimeter"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (mm))))
      }else if(unit1=="micrometer"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (mu*m))))
      }else if(unit1=="nanometer"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (nm))))
      }else{
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " (",  .(unit1), ")")))
      }
      
      # Shade of the x-axis min and max
      shade_xmin <- min(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
      shade_xmax <- max(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
      
      # Location of text
      text_x1 <- max(cell_list[[i]]$time) - (max(cell_list[[i]]$time) - min(cell_list[[i]]$time)) * 0.13
      text_x2 <- max(cell_list[[i]]$time) - (max(cell_list[[i]]$time) - min(cell_list[[i]]$time)) * 0.16
      
      # Distance
      df_g <- data.frame(x = cell_list[[i]]$time, y = Y)
      ymax <- max(df_g$y)
      ymin <- min(df_g$y)
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05
      
      g_dist <- ggplot(data = df_g) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_line(aes(x = x, y = y), size=0.5) +
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
             x = "Time (min)",
             y = label_y)
      
      # Integrate plots
      suppressWarnings(
        ggsave(paste0(out, "/nomodel_cell", i, "_", res_name, j, ".pdf"),
               g_dist, height = ps*20*1/4, width = ps*10*1.2, units = "mm")
      )
      
    }
  }
      
      ## Save movement time
      if(is.data.frame(visual)){
        df_visual_mv <- cbind(visual, df_mv)
      }else{
        df_visual_mv <- df_mv
      }
      fwrite(df_visual_mv, file = paste0(out, "/nomodel_mvtime.csv"))
      
}

