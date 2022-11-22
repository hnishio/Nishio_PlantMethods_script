
# Set working directory
setwd("~/Nishio_PlantMethods_script")

# Load packages
library(ggplot2)
library(patchwork)
library(data.table)

# Load plot functions
source("functions/Plot_functions_221101.R")

# Create an output directory
out <- "09_KFAS_visual_vs_predicted_221101"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data
df <- as.data.frame(fread("08_KFAS_chloroplast_velocity_221101/csv/ssm_KFAS_mvtime.csv"))
df <- df[,1:4]
names(df)[1:4] <- c("cell", "index", "visual", "predicted")
df <- df[!is.infinite(rowSums(df)),]


## Plotting
g <- cor_vis_pred(df, "Kalman filter")
ggsave(paste0(out, "/KFAS_vis_vs_pred.pdf"),
       g, height = 50, width = 45, units = "mm")

