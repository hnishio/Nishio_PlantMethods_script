
##### Visualization of distance from ex #####
distance_vis <- function(cell_list, res_name, ex_name, unit1, unit2,
                         shade = TRUE, ps = 7, theme_plot = "bw",
                         out = out){
  
  # Set axis limits
  xmaxs <- NULL
  xmins <- NULL
  ymaxs <- NULL
  ymins <- NULL
  for(i in 1:length(cell_list)){
    xmax <- max(cell_list[[i]]$time)
    xmin <- min(cell_list[[i]]$time)
    ymax <- max(cell_list[[i]][,3:ncol(cell_list[[i]])])
    ymin <- min(cell_list[[i]][,3:ncol(cell_list[[i]])])
    xmaxs <- cbind(xmaxs, xmax)
    xmins <- cbind(xmins, xmin)
    ymaxs <- cbind(ymaxs, ymax)
    ymins <- cbind(ymins, ymin)
  }
  
  xmax <- max(xmaxs)
  xmin <- min(xmins)
  xrange <- xmax - xmin  
  ymax <- max(ymaxs)
  ymin <- min(ymins)
  yrange <- ymax - ymin
  
  max_yaxis <- ymax + yrange*0.05
  min_yaxis <- ymin - yrange*0.05
  max_xaxis <- xmax
  min_xaxis <- xmin
  
  
  ## Plotting
  
  # List of results
  glist <- list()
  
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
  
  # for loop of cells
  for(i in 1:length(cell_list)){
    
    # Shade of the x-axis min and max
    shade_xmin <- min(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
    shade_xmax <- max(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
    
    g <- ggplot(data = cell_list[[i]], aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = min_yaxis, ymax = max_yaxis, alpha = alpha, fill = "gray50") +
      coord_cartesian(xlim=c(min_xaxis, max_xaxis), ylim=c(min_yaxis, max_yaxis), clip='off') +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_plot2 +
      theme(legend.position = "none",
            axis.title = element_text(size = ps),
            axis.text = element_text(size = ps),
            plot.title = element_text(size = ps, face = "bold")) +
      labs(title = paste("Cell ", i, sep=""),
           x = "Time (min)", 
           y = label_y)
    
    # for loop of response variable
    for(j in 1:(ncol(cell_list[[i]])-2)){
      g <- g + eval(parse(text = paste("geom_line(aes(y = distance", j, "), col = as.vector(viridis(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
    }
    
    glist[[i]] <- g
  }
  
  return(glist)
  
}

