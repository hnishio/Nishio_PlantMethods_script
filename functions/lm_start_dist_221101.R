
##### Linear regression (x: distance, y: start time) #####
lm_start_dist <- function(cell_list, mvtime,
                          ex_name, unit1, unit2, title = NULL, ps = 7,
                          theme_plot = "bw"){
  
  ## Error message
  if(length(cell_list) != length(unique(mvtime$cell))){
    stop(paste("The length of 'cell_list' should be the same as the unique identifier at 'mvtime$cell'"))
  }
  
  
  ## Adjust data.frame
  df <- mvtime[,1:4]
  names(df)[1:4] <- c("cell", "index", "visual", "predicted")
  distance <- NULL
  for(i in 1:length(cell_list)){
    distance <- c(distance, as.numeric(cell_list[[i]][max(which(cell_list[[i]]$ex == 0)),-(1:2)]))
  }
  df$distance <- distance
  df <- df[!is.infinite(rowSums(df)),]
  
  
  ## label
  if(unit1=="meter"){
    label_y <- bquote(paste("Distance from ", .(ex_name), " ", (m)))
    label_x <- bquote(paste("Start time ", (.(unit2))))
  }else if(unit1=="centimeter"){
    label_y <- bquote(paste("Distance from ", .(ex_name), " ", (cm)))
    label_x <- bquote(paste("Start time ", (.(unit2))))
  }else if(unit1=="millimeter"){
    label_y <- bquote(paste("Distance from ", .(ex_name), " ", (mm)))
    label_x <- bquote(paste("Start time ", (.(unit2))))
  }else if(unit1=="micrometer"){
    label_y <- bquote(paste("Distance from ", .(ex_name), " ", (mu*m)))
    label_x <- bquote(paste("Start time ", (.(unit2))))
  }else if(unit1=="nanometer"){
    label_y <- bquote(paste("Distance from ", .(ex_name), " ", (nm)))
    label_x <- bquote(paste("Start time ", (.(unit2))))
  }else{
    label_y <- bquote(paste("Distance from ", .(ex_name), " (", .(unit1), ")"))
    label_x <- bquote(paste("Start time ", (.(unit2))))
  }
  
  
  ## Theme
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
  
  
  ## Title of the plots
  if(is.null(title) == T){
    title <- paste0("Cell ", as.character(unique(df$cell)))
    title <- c(title, "All cells")
  }
  
  ## Plotting for each cell
  glist <- list()
  for(i in 1:length(cell_list)){
    data <- df[df$cell==i,]
    
    range_x <- diff(range(df$predicted))
    min_axis_x <- min(df$predicted)
    max_axis_x <- max(df$predicted)
    range_y <- diff(range(df$distance))
    min_axis_y <- min(df$distance) - range_y*0.1
    max_axis_y <- max(df$distance) + range_y*0.1
    
    glist[[i]] <- ggplot(data, aes(x=predicted, y=distance)) +
      geom_smooth(method="lm", color = "steelblue", fill = "steelblue") +
      geom_point(size=0.8, alpha=0.5) +
      stat_regline_equation(label.x=min_axis_x+range_x*0.05, label.y=max_axis_y-range_y*0.05, size=ps/ggplot2::.pt) +
      stat_cor(aes(label=..rr.label..), digits = 2, label.x=min_axis_x+range_x*0.05, label.y=max_axis_y-range_y*0.17, size=ps/ggplot2::.pt) +
      coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
      theme_plot2 +
      theme(plot.title = element_text(size=ps, face = "bold"),
            axis.title=element_text(size=ps), 
            axis.text=element_text(size=ps),
            plot.tag = element_text(size = 12, face = "bold")) + 
      labs(title = title[i], x=label_x, 
           y = label_y)
  }
  
  
  ## Plotting for all cells
  data <- df
  
  range_x <- diff(range(df$predicted))
  min_axis_x <- min(df$predicted)
  max_axis_x <- max(df$predicted)
  range_y <- diff(range(df$distance))
  min_axis_y <- min(df$distance) - range_y*0.1
  max_axis_y <- max(df$distance) + range_y*0.1
  
  glist[[length(cell_list) + 1]] <- ggplot(data, aes(x=predicted, y=distance)) +
    geom_smooth(method="lm", color = "steelblue", fill = "steelblue") +
    geom_point(size=0.8, alpha=0.5) +
    stat_regline_equation(label.x=min_axis_x+range_x*0.05, label.y=max_axis_y-range_y*0.05, size=ps/ggplot2::.pt) +
    stat_cor(aes(label=..rr.label..), digits = 2, label.x=min_axis_x+range_x*0.05, label.y=max_axis_y-range_y*0.17, size=ps/ggplot2::.pt) +
    coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
    theme_plot2 +
    theme(plot.title = element_text(size=ps, face = "bold"),
          axis.title=element_text(size=ps), 
          axis.text=element_text(size=ps),
          plot.tag = element_text(size = 12, face = "bold")) + 
    labs(title = title[i+1], x=label_x, 
         y = label_y)
  
  
  ## Signaling speed
  
  # Speed of all cells
  lm_allcells <- lm(formula = distance~predicted, data = df)
  
  # Speed of each cell
  ydiff <- 1 / (length(cell_list) + 1)
  g_speed <- ggplot() + theme_void()
  
  for(i in 1:length(cell_list)){
    data <- df[df$cell==i,]
    lm_cell <- lm(formula = distance~predicted, data = df[df$cell==i,])
    
    # label
    if(unit1=="meter"){
      label_all <- bquote(paste("All cells:  ", .(round(lm_allcells$coefficients[2], 2)), " ", (m/.(unit2)), sep = ""))
      label_each <- bquote(paste("Cell ", .(i), ":  ", .(round(lm_cell$coefficients[2], 2)), " ", (m/.(unit2)), sep = ""))
    }else if(unit1=="centimeter"){
      label_all <- bquote(paste("All cells:  ", .(round(lm_allcells$coefficients[2], 2)), " ", (cm/.(unit2)), sep = ""))
      label_each <- bquote(paste("Cell ", .(i), ":  ", .(round(lm_cell$coefficients[2], 2)), " ", (cm/.(unit2)), sep = ""))
    }else if(unit1=="millimeter"){
      label_all <- bquote(paste("All cells:  ", .(round(lm_allcells$coefficients[2], 2)), " ", (mm/.(unit2)), sep = ""))
      label_each <- bquote(paste("Cell ", .(i), ":  ", .(round(lm_cell$coefficients[2], 2)), " ", (mm/.(unit2)), sep = ""))
    }else if(unit1=="micrometer"){
      label_all <- bquote(paste("All cells:  ", .(round(lm_allcells$coefficients[2], 2)), " ", (mu*m/.(unit2)), sep = ""))
      label_each <- bquote(paste("Cell ", .(i), ":  ", .(round(lm_cell$coefficients[2], 2)), " ", (mu*m/.(unit2)), sep = ""))
    }else if(unit1=="nanometer"){
      label_all <- bquote(paste("All cells:  ", .(round(lm_allcells$coefficients[2], 2)), " ", (nm/.(unit2)), sep = ""))
      label_each <- bquote(paste("Cell ", .(i), ":  ", .(round(lm_cell$coefficients[2], 2)), " ", (nm/.(unit2)), sep = ""))
    }else{
      
    }
    g_speed <- g_speed + 
      annotate("text", x = 0.5, y = 1 - ydiff*i,
               label = label_each, 
               size = ps/ggplot2::.pt)
  }
  
  # Speed of all cells
  g_speed <- g_speed + 
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip='off') +
    annotate("text", x=0.5, y=1.0,
             label="Signal transfer speed", 
             size=ps/ggplot2::.pt, fontface = 2) +
    annotate("text", x=0.5, y=0.0,
             label = label_all, 
             size=ps/ggplot2::.pt)
  
  glist[[length(cell_list) + 2]] <- g_speed
  
  return(glist)
}
