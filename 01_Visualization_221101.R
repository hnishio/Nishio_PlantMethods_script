
# Set working directory
setwd("~/Nishio_PlantMethods_script")

# Load packages
library(ggplot2)
library(patchwork)
library(viridis)

# Load plot functions
source("functions/distance_vis_221101.R")

# Create an output directory
out <- "01_Visualization_221101"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}


# Load data of chloroplasts
cell1 <- read.csv("data/cell1.csv")
cell2 <- read.csv("data/cell2.csv")
cell3 <- read.csv("data/cell3.csv")
cell4 <- read.csv("data/cell4.csv")
cell_list <- list(cell1, cell2, cell3, cell4)

# Plotting
glist <- distance_vis(cell_list = cell_list,
                      res_name = "chloroplast", ex_name = "microbeam", 
                      unit1 = "micrometer", unit2 = "min", 
                      out = out)

g <- glist[[1]] + glist[[2]] + glist[[3]] + glist[[4]] + 
  plot_layout(ncol = 2, heights = c(1, 1))
ggsave(paste0(out, "/distance_ex_221101.pdf"),
       g, height = 120, width = 180, units = "mm")

