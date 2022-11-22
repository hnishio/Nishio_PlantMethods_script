
### Linear regression between the distance at time 0 and the coefficient of an explanatory variable
lm_dist_beta <- function(cell_list, mvtime, ssm_path,
                         ex_sign = "negative", ssm_method, res_name, ex_name, 
                         unit1, unit2, ps = 7,
                         theme_plot = "bw"){
  
  mvtime <- mvtime[,1:4]
  names(mvtime)[1:4] <- c("cell", "chl", "visual", "predicted")
  #mvtime <- mvtime[!is.infinite(rowSums(mvtime)),]
  distance <- NULL
  for(i in 1:length(cell_list)){
    distance <- c(distance, as.numeric(cell_list[[i]][max(which(cell_list[[i]]$ex == 0)),-(1:2)]))
  }
  mvtime$distance <- distance
  
  
  
  if(ssm_method == "Bayes"){
    
    ## Prepare a container for movement time
    vec_most_b <- numeric()
    vec_mean_b <- numeric()
    vec_most_alpha <- numeric()
    vec_mean_alpha <- numeric()
    vec_mean_w <- numeric()
    vec_s_w <- numeric()
    vec_s_b_ex <- numeric()
    vec_s_Y <- numeric()
    
    
    ## Get values
    if(ex_sign == "positive"){  # positive
      
      # for loop of cells
      for(i in 1:length(cell_list)){
        
        # for loop of explanatory variable
        for (j in 1:(ncol(cell_list[[i]])-2)){
          outcsv_name <- list.files(paste0(ssm_path, "/csv"))
          outcsv_name2 <- outcsv_name[grep(paste0("cell", i, "_", res_name, j, "\\."), outcsv_name)]
          outcsv_name3 <- outcsv_name[grep(paste0("cell", i, "_", res_name, j, "_sd\\."), outcsv_name)]
          
          eval(parse(text = paste0(
            "df_", i, "_", j, " <- as.data.frame(fread('", ssm_path,"/csv/", outcsv_name2, "'))
      
      most_b <- max(df_", i, "_", j, "$`b_ex_50%`, na.rm = T)
      vec_most_b <- c(vec_most_b, most_b)

      mean_b <- mean(df_", i, "_", j, "$`b_ex_50%`, na.rm = T)
      vec_mean_b <- c(vec_mean_b, mean_b)

      most_alpha <- max(df_", i, "_", j, "$`alpha_50%`[df_", i, "_", j, "$time>=1])
      vec_most_alpha <- c(vec_most_alpha, most_alpha)

      mean_alpha <- mean(df_", i, "_", j, "$`alpha_50%`[df_", i, "_", j, "$time>=1])
      vec_mean_alpha <- c(vec_mean_alpha, mean_alpha)

      mean_w <- mean(abs(df_", i, "_", j, "$`w_50%`[df_", i, "_", j, "$time>=1]))
      vec_mean_w <- c(vec_mean_w, mean_w)
      
      df_", i, "_", j, "_sd <- as.data.frame(fread('", ssm_path,"/csv/", outcsv_name3, "'))
      
      s_w <- df_", i, "_", j, "_sd$`50%`[df_", i, "_", j, "_sd$s_name=='s_w']
      s_b_ex <- df_", i, "_", j, "_sd$`50%`[df_", i, "_", j, "_sd$s_name=='s_b_ex']
      s_Y <- df_", i, "_", j, "_sd$`50%`[df_", i, "_", j, "_sd$s_name=='s_Y']
      vec_s_w <- c(vec_s_w, s_w)
      vec_s_b_ex <- c(vec_s_b_ex, s_b_ex)
      vec_s_Y <- c(vec_s_Y, s_Y)"
          )))
          
        }
      }
      
    }else{  # negative
      
      # for loop of cells
      for(i in 1:length(cell_list)){
        
        # for loop of explanatory variable
        for (j in 1:(ncol(cell_list[[i]])-2)){
          outcsv_name <- list.files(paste0(ssm_path, "/csv"))
          outcsv_name2 <- outcsv_name[grep(paste0("cell", i, "_", res_name, j, "\\."), outcsv_name)]
          outcsv_name3 <- outcsv_name[grep(paste0("cell", i, "_", res_name, j, "_sd\\."), outcsv_name)]
          
          eval(parse(text = paste0(
            "df_", i, "_", j, " <- as.data.frame(fread('", ssm_path,"/csv/", outcsv_name2, "'))
      
      most_b <- min(df_", i, "_", j, "$`b_ex_50%`, na.rm = T)
      vec_most_b <- c(vec_most_b, most_b)

      mean_b <- mean(df_", i, "_", j, "$`b_ex_50%`, na.rm = T)
      vec_mean_b <- c(vec_mean_b, mean_b)

      most_alpha <- min(df_", i, "_", j, "$`alpha_50%`[df_", i, "_", j, "$time>=1])
      vec_most_alpha <- c(vec_most_alpha, most_alpha)

      mean_alpha <- mean(df_", i, "_", j, "$`alpha_50%`[df_", i, "_", j, "$time>=1])
      vec_mean_alpha <- c(vec_mean_alpha, mean_alpha)

      mean_w <- mean(abs(df_", i, "_", j, "$`w_50%`[df_", i, "_", j, "$time>=1]))
      vec_mean_w <- c(vec_mean_w, mean_w)
      
      df_", i, "_", j, "_sd <- as.data.frame(fread('", ssm_path,"/csv/", outcsv_name3, "'))
      
      s_w <- df_", i, "_", j, "_sd$`50%`[df_", i, "_", j, "_sd$s_name=='s_w']
      s_b_ex <- df_", i, "_", j, "_sd$`50%`[df_", i, "_", j, "_sd$s_name=='s_b_ex']
      s_Y <- df_", i, "_", j, "_sd$`50%`[df_", i, "_", j, "_sd$s_name=='s_Y']
      vec_s_w <- c(vec_s_w, s_w)
      vec_s_b_ex <- c(vec_s_b_ex, s_b_ex)
      vec_s_Y <- c(vec_s_Y, s_Y)"
          )))
          
        }
      }
      
    }
    
    
    df_new <- cbind(mvtime, 
                    data.frame(most_b = vec_most_b,
                               mean_b = vec_mean_b,
                               most_alpha = vec_most_alpha,
                               mean_alpha = vec_mean_alpha,
                               mean_w = vec_mean_w,
                               s_w = vec_s_w,
                               s_b_ex = vec_s_b_ex,
                               s_Y = vec_s_Y)
    )
    
    
    ## label
    if(unit1=="meter"){
      label_x <- bquote(paste("Distance from ", .(ex_name)," (m)", sep = ""))
      label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (m/.(unit2)))))
      label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (m/.(unit2)))))
      label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (m/.(unit2)))))
      label_y_meanvelo <- bquote(atop("Mean velocity of", paste("movement ", (m/.(unit2)))))
      label_y_mostvelo <- bquote(atop(paste("Most ", .(ex_sign), " velocity of"), paste("movement ", (m/.(unit2)))))
    }else if(unit1=="centimeter"){
      label_x <- bquote(paste("Distance from ", .(ex_name)," (cm)", sep = ""))
      label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (cm/.(unit2)))))
      label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (cm/.(unit2)))))
      label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (cm/.(unit2)))))
      label_y_meanvelo <- bquote(atop("Mean velocity of", paste("movement ", (cm/.(unit2)))))
      label_y_mostvelo <- bquote(atop(paste("Most ", .(ex_sign), " velocity of"), paste("movement ", (cm/.(unit2)))))
    }else if(unit1=="millimeter"){
      label_x <- bquote(paste("Distance from ", .(ex_name)," (mm)", sep = ""))
      label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (mm/.(unit2)))))
      label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (mm/.(unit2)))))
      label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (mm/.(unit2)))))
      label_y_meanvelo <- bquote(atop("Mean velocity of", paste("movement ", (mm/.(unit2)))))
      label_y_mostvelo <- bquote(atop(paste("Most ", .(ex_sign), " velocity of"), paste("movement ", (mm/.(unit2)))))
    }else if(unit1=="micrometer"){
      label_x <- bquote(paste("Distance from ", .(ex_name)," (", mu, "m)", sep = ""))
      label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (mu*m/.(unit2)))))
      label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (mu*m/.(unit2)))))
      label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (mu*m/.(unit2)))))
      label_y_meanvelo <- bquote(atop("Mean velocity of", paste("movement ", (mu*m/.(unit2)))))
      label_y_mostvelo <- bquote(atop(paste("Most ", .(ex_sign), " velocity of"), paste("movement ", (mu*m/.(unit2)))))
    }else if(unit1=="nanometer"){
      label_x <- bquote(paste("Distance from ", .(ex_name)," (nm)", sep = ""))
      label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (nm/.(unit2)))))
      label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (nm/.(unit2)))))
      label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (nm/.(unit2)))))
      label_y_meanvelo <- bquote(atop("Mean velocity of", paste("movement ", (nm/.(unit2)))))
      label_y_mostvelo <- bquote(atop(paste("Most ", .(ex_sign), " velocity of"), paste("movement ", (nm/.(unit2)))))
    }else{
      label_x <- bquote(paste("Distance from ", .(ex_name)," (", .(unit1), ")", sep = ""))
      label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (.(unit1)/.(unit2)))))
      label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (.(unit1)/.(unit2)))))
      label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (.(unit1)/.(unit2)))))
      label_y_meanvelo <- bquote(atop("Mean velocity of", paste("movement ", (.(unit1)/.(unit2)))))
      label_y_mostvelo <- bquote(atop(paste("Most ", .(ex_sign), " velocity of"), paste("movement ", (.(unit1)/.(unit2)))))
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
    
    
    ##### Linear regression (x: distance, y: mean_b) #####
    model <- siegel_regression(mean_b ~ distance, data = df_new)
    #summary(model)
    
    newx = seq(min(df_new$distance), max(df_new$distance), by = 1)
    suppressWarnings(conf_interval <- predict(model, newdata=data.frame(x=newx), interval="confidence", level = 0.95))
    conf_interval2 <- as.data.frame(cbind(df_new$distance, conf_interval)[order(df_new$distance, decreasing = F),])
    names(conf_interval2)[1] <- "distance"
    
    
    ## Plotting
    r2 <- formatC(summary(model)$r.squared, digits=2)
    p <- formatC(summary(model)$coefficients[2,4], digits = 2, format = "e")
    r2lab <- bquote(paste(R^2, " = ", .(r2), sep=""))
    plab <- bquote(paste(italic(P), " = ", .(p), sep=""))
    
    min_axis_y <- min(df_new$mean_b) - diff(range(df_new$mean_b))*0.1
    max_axis_y <- max(df_new$mean_b) + diff(range(df_new$mean_b))*0.1
    min_axis_x <- min(df_new$distance)
    max_axis_x <- max(df_new$distance)
    
    g1 <- ggplot() +
      geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue") +
      geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
      geom_point(data = df_new, aes(x=distance, y=mean_b), size=0.8, alpha=0.5) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.03, 
               label=r2lab, size=ps/ggplot2::.pt) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.13, 
               label=plab, size=ps/ggplot2::.pt) +
      coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
      theme_plot2 +
      theme(plot.title = element_text(size=ps, face = "bold"),
            axis.title=element_text(size=ps), 
            axis.text=element_text(size=ps),
            plot.tag = element_text(size = 12, face = "bold")) + 
      labs(x = label_x, y = label_y_meancoef)
    
    
    
    
    
    ##### Linear regression (x: distance, y: most_b) #####
    model <- siegel_regression(most_b ~ distance, data = df_new)
    #summary(model)
    
    newx = seq(min(df_new$distance), max(df_new$distance), by = 1)
    suppressWarnings(conf_interval <- predict(model, newdata=data.frame(x=newx), interval="confidence", level = 0.95))
    conf_interval2 <- as.data.frame(cbind(df_new$distance, conf_interval)[order(df_new$distance, decreasing = F),])
    names(conf_interval2)[1] <- "distance"
    
    
    ## Plotting
    r2 <- formatC(summary(model)$r.squared, digits=2)
    p <- formatC(summary(model)$coefficients[2,4], digits = 2, format = "e")
    r2lab <- bquote(paste(R^2, " = ", .(r2), sep=""))
    plab <- bquote(paste(italic(P), " = ", .(p), sep=""))
    
    min_axis_y <- min(df_new$most_b) - diff(range(df_new$most_b))*0.1
    max_axis_y <- max(df_new$most_b) + diff(range(df_new$most_b))*0.1
    min_axis_x <- min(df_new$distance)
    max_axis_x <- max(df_new$distance)
    
    g2 <- ggplot() +
      geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue") +
      geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
      geom_point(data = df_new, aes(x=distance, y=most_b), size=0.8, alpha=0.5) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.03, 
               label=r2lab, size=ps/ggplot2::.pt) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.13, 
               label=plab, size=ps/ggplot2::.pt) +
      coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
      theme_plot2 +
      theme(plot.title = element_text(size=ps, face = "bold"),
            axis.title=element_text(size=ps), 
            axis.text=element_text(size=ps),
            plot.tag = element_text(size = 12, face = "bold")) + 
      labs(x = label_x, y = label_y_mostcoef)
    
    
    
    
    ##### Linear regression (x: distance, y: s_b_ex) #####
    model <- siegel_regression(s_b_ex ~ distance, data = df_new)
    #summary(model)
    
    newx = seq(min(df_new$distance), max(df_new$distance), by = 1)
    suppressWarnings(conf_interval <- predict(model, newdata=data.frame(x=newx), interval="confidence", level = 0.95))
    conf_interval2 <- as.data.frame(cbind(df_new$distance, conf_interval)[order(df_new$distance, decreasing = F),])
    names(conf_interval2)[1] <- "distance"
    
    
    ## Plotting
    r2 <- formatC(summary(model)$r.squared, digits=2)
    p <- formatC(summary(model)$coefficients[2,4], digits = 2, format = "e")
    r2lab <- bquote(paste(R^2, " = ", .(r2), sep=""))
    plab <- bquote(paste(italic(P), " = ", .(p), sep=""))
    
    min_axis_y <- min(df_new$s_b_ex) - diff(range(df_new$s_b_ex))*0.1
    max_axis_y <- max(df_new$s_b_ex) + diff(range(df_new$s_b_ex))*0.1
    min_axis_x <- min(df_new$distance)
    max_axis_x <- max(df_new$distance)
    
    g3 <- ggplot() +
      geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue") +
      geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
      geom_point(data = df_new, aes(x=distance, y=s_b_ex), size=0.8, alpha=0.5) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.03, 
               label=r2lab, size=ps/ggplot2::.pt) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.13, 
               label=plab, size=ps/ggplot2::.pt) +
      coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
      theme_plot2 +
      theme(plot.title = element_text(size=ps, face = "bold"),
            axis.title=element_text(size=ps), 
            axis.text=element_text(size=ps),
            plot.tag = element_text(size = 12, face = "bold")) + 
      labs(x = label_x, y = label_y_sdcoef)
    
    
    
    
    
    ##### Linear regression (x: distance, y: mean_alpha) #####
    model <- siegel_regression(mean_alpha ~ distance, data = df_new)
    summary(model)
    
    newx = seq(min(df_new$distance), max(df_new$distance), by = 1)
    suppressWarnings(conf_interval <- predict(model, newdata=data.frame(x=newx), interval="confidence", level = 0.95))
    conf_interval2 <- as.data.frame(cbind(df_new$distance, conf_interval)[order(df_new$distance, decreasing = F),])
    names(conf_interval2)[1] <- "distance"
    
    
    ## Plotting
    r2 <- formatC(summary(model)$r.squared, digits=2)
    p <- formatC(summary(model)$coefficients[2,4], digits = 2, format = "e")
    r2lab <- bquote(paste(R^2, " = ", .(r2), sep=""))
    plab <- bquote(paste(italic(P), " = ", .(p), sep=""))
    
    min_axis_y <- min(df_new$mean_alpha) - diff(range(df_new$mean_alpha))*0.1
    max_axis_y <- max(df_new$mean_alpha) + diff(range(df_new$mean_alpha))*0.1
    min_axis_x <- min(df_new$distance)
    max_axis_x <- max(df_new$distance)
    
    g4 <- ggplot() +
      geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue") +
      geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
      geom_point(data = df_new, aes(x=distance, y=mean_alpha), size=0.8, alpha=0.5) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.03, 
               label=r2lab, size=ps/ggplot2::.pt) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.13, 
               label=plab, size=ps/ggplot2::.pt) +
      coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
      theme_plot2 +
      theme(plot.title = element_text(size=ps, face = "bold"),
            axis.title=element_text(size=ps), 
            axis.text=element_text(size=ps),
            plot.tag = element_text(size = 12, face = "bold")) + 
      labs(x = label_x, y = label_y_meanvelo)
    
    
    
    
    
    ##### Linear regression (x: distance, y: most_alpha) #####
    model <- siegel_regression(most_alpha ~ distance, data = df_new)
    summary(model)
    
    newx = seq(min(df_new$distance), max(df_new$distance), by = 1)
    suppressWarnings(conf_interval <- predict(model, newdata=data.frame(x=newx), interval="confidence", level = 0.95))
    conf_interval2 <- as.data.frame(cbind(df_new$distance, conf_interval)[order(df_new$distance, decreasing = F),])
    names(conf_interval2)[1] <- "distance"
    
    
    ## Plotting
    r2 <- formatC(summary(model)$r.squared, digits=2)
    p <- formatC(summary(model)$coefficients[2,4], digits = 2, format = "e")
    r2lab <- bquote(paste(R^2, " = ", .(r2), sep=""))
    plab <- bquote(paste(italic(P), " = ", .(p), sep=""))
    
    min_axis_y <- min(df_new$most_alpha) - diff(range(df_new$most_alpha))*0.1
    max_axis_y <- max(df_new$most_alpha) + diff(range(df_new$most_alpha))*0.1
    min_axis_x <- min(df_new$distance)
    max_axis_x <- max(df_new$distance)
    
    g5 <- ggplot() +
      geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue") +
      geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
      geom_point(data = df_new, aes(x=distance, y=most_alpha), size=0.8, alpha=0.5) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.03, 
               label=r2lab, size=ps/ggplot2::.pt) +
      annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.13, 
               label=plab, size=ps/ggplot2::.pt) +
      coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
      theme_plot2 +
      theme(plot.title = element_text(size=ps, face = "bold"),
            axis.title=element_text(size=ps), 
            axis.text=element_text(size=ps),
            plot.tag = element_text(size = 12, face = "bold")) + 
      labs(x = label_x, y = label_y_mostvelo)
    
    return(list(g1, g2, g3, g4, g5))
    
  }else if(ssm_method == "KFAS"){
    
    
    ## Prepare a container for movement time
    vec_most_b <- numeric()
    vec_mean_b <- numeric()
    
    
    ## Get values
    if(ex_sign == "positive"){  # positive
      
      # for loop of cells
      for(i in 1:length(cell_list)){
        
        # for loop of explanatory variable
        for (j in 1:(ncol(cell_list[[i]])-2)){
          outcsv_name <- list.files(paste0(ssm_path, "/csv"))
          outcsv_name2 <- outcsv_name[grep(paste0("cell", i, "_", res_name, j, "\\."), outcsv_name)]

          eval(parse(text = paste0(
            "df_", i, "_", j, " <- as.data.frame(fread('", ssm_path,"/csv/", outcsv_name2, "'))
      
           most_b <- max(df_", i, "_", j, "$`b_ex_50%`, na.rm = T)
           vec_most_b <- c(vec_most_b, most_b)

           mean_b <- mean(df_", i, "_", j, "$`b_ex_50%`, na.rm = T)
           vec_mean_b <- c(vec_mean_b, mean_b)"
          )))
          
        }
      }
      
    }else{  # negative
      
      # for loop of cells
      for(i in 1:length(cell_list)){
        
      # for loop of explanatory variable
      for (j in 1:(ncol(cell_list[[i]])-2)){
        outcsv_name <- list.files(paste0(ssm_path, "/csv"))
        outcsv_name2 <- outcsv_name[grep(paste0("cell", i, "_", res_name, j, "\\."), outcsv_name)]

        eval(parse(text = paste0(
          "df_", i, "_", j, " <- as.data.frame(fread('", ssm_path,"/csv/", outcsv_name2, "'))
      
           most_b <- min(df_", i, "_", j, "$`b_ex_50%`, na.rm = T)
           vec_most_b <- c(vec_most_b, most_b)

           mean_b <- mean(df_", i, "_", j, "$`b_ex_50%`, na.rm = T)
           vec_mean_b <- c(vec_mean_b, mean_b)"
        )))
        
      }
      }
    }
  

  ## Get SD
  df_sd <- as.data.frame(fread(paste0(ssm_path,"/csv/ssm_KFAS_sd.csv")))
  
  
  ## Prepare data.frame
  df_new <- cbind(mvtime, 
                  data.frame(most_b = vec_most_b,
                             mean_b = vec_mean_b,
                             s_b_ex = df_sd$s_b_ex,
                             s_Y = df_sd$s_Y)
  )
  
  
  ## label
  if(unit1=="meter"){
    label_x <- bquote(paste("Distance from ", .(ex_name)," (m)", sep = ""))
    label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (m/.(unit2)))))
    label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (m/.(unit2)))))
    label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (m/.(unit2)))))
  }else if(unit1=="centimeter"){
    label_x <- bquote(paste("Distance from ", .(ex_name)," (cm)", sep = ""))
    label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (cm/.(unit2)))))
    label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (cm/.(unit2)))))
    label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (cm/.(unit2)))))
  }else if(unit1=="millimeter"){
    label_x <- bquote(paste("Distance from ", .(ex_name)," (mm)", sep = ""))
    label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (mm/.(unit2)))))
    label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (mm/.(unit2)))))
    label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (mm/.(unit2)))))
  }else if(unit1=="micrometer"){
    label_x <- bquote(paste("Distance from ", .(ex_name)," (", mu, "m)", sep = ""))
    label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (mu*m/.(unit2)))))
    label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (mu*m/.(unit2)))))
    label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (mu*m/.(unit2)))))
  }else if(unit1=="nanometer"){
    label_x <- bquote(paste("Distance from ", .(ex_name)," (nm)", sep = ""))
    label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (nm/.(unit2)))))
    label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (nm/.(unit2)))))
    label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (nm/.(unit2)))))
  }else{
    label_x <- bquote(paste("Distance from ", .(ex_name)," (", .(unit1), ")", sep = ""))
    label_y_meancoef <- bquote(atop("Mean coefficient of", paste(.(ex_name), " ", (.(unit1)/.(unit2)))))
    label_y_mostcoef <- bquote(atop(paste("Most ", .(ex_sign)," coefficient of"), paste(.(ex_name), " ", (.(unit1)/.(unit2)))))
    label_y_sdcoef <- bquote(atop("SD of coefficient of", paste(.(ex_name), " ", (.(unit1)/.(unit2)))))
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
  
  
  ##### Linear regression (x: distance, y: mean_b) #####
  model <- siegel_regression(mean_b ~ distance, data = df_new)
  #summary(model)
  
  newx = seq(min(df_new$distance), max(df_new$distance), by = 1)
  suppressWarnings(conf_interval <- predict(model, newdata=data.frame(x=newx), interval="confidence", level = 0.95))
  conf_interval2 <- as.data.frame(cbind(df_new$distance, conf_interval)[order(df_new$distance, decreasing = F),])
  names(conf_interval2)[1] <- "distance"
  
  
  ## Plotting
  r2 <- formatC(summary(model)$r.squared, digits=2)
  p <- formatC(summary(model)$coefficients[2,4], digits = 2, format = "e")
  r2lab <- bquote(paste(R^2, " = ", .(r2), sep=""))
  plab <- bquote(paste(italic(P), " = ", .(p), sep=""))
  
  min_axis_y <- min(df_new$mean_b) - diff(range(df_new$mean_b))*0.1
  max_axis_y <- max(df_new$mean_b) + diff(range(df_new$mean_b))*0.1
  min_axis_x <- min(df_new$distance)
  max_axis_x <- max(df_new$distance)
  
  g1 <- ggplot() +
    geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue") +
    geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
    geom_point(data = df_new, aes(x=distance, y=mean_b), size=0.8, alpha=0.5) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.03, 
             label=r2lab, size=ps/ggplot2::.pt) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.13, 
             label=plab, size=ps/ggplot2::.pt) +
    coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
    theme_plot2 +
    theme(plot.title = element_text(size=ps, face = "bold"),
          axis.title=element_text(size=ps), 
          axis.text=element_text(size=ps),
          plot.tag = element_text(size = 12, face = "bold")) + 
    labs(x = label_x, y = label_y_meancoef)
  
  
  
  
  
  ##### Linear regression (x: distance, y: most_b) #####
  model <- siegel_regression(most_b ~ distance, data = df_new)
  #summary(model)
  
  newx = seq(min(df_new$distance), max(df_new$distance), by = 1)
  suppressWarnings(conf_interval <- predict(model, newdata=data.frame(x=newx), interval="confidence", level = 0.95))
  conf_interval2 <- as.data.frame(cbind(df_new$distance, conf_interval)[order(df_new$distance, decreasing = F),])
  names(conf_interval2)[1] <- "distance"
  
  
  ## Plotting
  r2 <- formatC(summary(model)$r.squared, digits=2)
  p <- formatC(summary(model)$coefficients[2,4], digits = 2, format = "e")
  r2lab <- bquote(paste(R^2, " = ", .(r2), sep=""))
  plab <- bquote(paste(italic(P), " = ", .(p), sep=""))
  
  min_axis_y <- min(df_new$most_b) - diff(range(df_new$most_b))*0.1
  max_axis_y <- max(df_new$most_b) + diff(range(df_new$most_b))*0.1
  min_axis_x <- min(df_new$distance)
  max_axis_x <- max(df_new$distance)
  
  g2 <- ggplot() +
    geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue") +
    geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
    geom_point(data = df_new, aes(x=distance, y=most_b), size=0.8, alpha=0.5) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.03, 
             label=r2lab, size=ps/ggplot2::.pt) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.13, 
             label=plab, size=ps/ggplot2::.pt) +
    coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
    theme_plot2 +
    theme(plot.title = element_text(size=ps, face = "bold"),
          axis.title=element_text(size=ps), 
          axis.text=element_text(size=ps),
          plot.tag = element_text(size = 12, face = "bold")) + 
    labs(x = label_x, y = label_y_mostcoef)
  
  
  
  
  ##### Linear regression (x: distance, y: s_b_ex) #####
  model <- siegel_regression(s_b_ex ~ distance, data = df_new)
  #summary(model)
  
  newx = seq(min(df_new$distance), max(df_new$distance), by = 1)
  suppressWarnings(conf_interval <- predict(model, newdata=data.frame(x=newx), interval="confidence", level = 0.95))
  conf_interval2 <- as.data.frame(cbind(df_new$distance, conf_interval)[order(df_new$distance, decreasing = F),])
  names(conf_interval2)[1] <- "distance"
  
  
  ## Plotting
  r2 <- formatC(summary(model)$r.squared, digits=2)
  p <- formatC(summary(model)$coefficients[2,4], digits = 2, format = "e")
  r2lab <- bquote(paste(R^2, " = ", .(r2), sep=""))
  plab <- bquote(paste(italic(P), " = ", .(p), sep=""))
  
  min_axis_y <- min(df_new$s_b_ex) - diff(range(df_new$s_b_ex))*0.1
  max_axis_y <- max(df_new$s_b_ex) + diff(range(df_new$s_b_ex))*0.1
  min_axis_x <- min(df_new$distance)
  max_axis_x <- max(df_new$distance)
  
  g3 <- ggplot() +
    geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue") +
    geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
    geom_point(data = df_new, aes(x=distance, y=s_b_ex), size=0.8, alpha=0.5) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.03, 
             label=r2lab, size=ps/ggplot2::.pt) +
    annotate("text", x=max_axis_x - (max_axis_x - min_axis_x)*0.25, y=max_axis_y - (max_axis_y - min_axis_y)*0.13, 
             label=plab, size=ps/ggplot2::.pt) +
    coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
    theme_plot2 +
    theme(plot.title = element_text(size=ps, face = "bold"),
          axis.title=element_text(size=ps), 
          axis.text=element_text(size=ps),
          plot.tag = element_text(size = 12, face = "bold")) + 
    labs(x = label_x, y = label_y_sdcoef)
  
  return(list(g1, g2, g3))
  
}

}

