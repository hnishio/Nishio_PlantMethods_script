
# Set working directory
setwd("~/Nishio_PlantMethods_script")

# Load packages
library(ggplot2)
library(patchwork)
library(data.table)
library(RobustLinearReg)

# Load plot functions
source("functions/lm_dist_beta_221101.R")

# Create output directories
out <- "05_CmdStan_beta_vs_distance_221101"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data
cell1 <- read.csv("data/cell1.csv")
cell2 <- read.csv("data/cell2.csv")
cell3 <- read.csv("data/cell3.csv")
cell4 <- read.csv("data/cell4.csv")
visual <- read.csv("data/visual_judgement.csv")
cell_list <- list(cell1, cell2, cell3, cell4)

mvtime <- as.data.frame(fread("02_CmdStan_chloroplast_velocity_221101/csv/ssm_individual_mvtime.csv"))

glist <- lm_dist_beta(cell_list = cell_list, mvtime = mvtime, 
                      ssm_path = "02_CmdStan_chloroplast_velocity_221101",
                      ssm_method = "Bayes", res_name = "chloroplast",
                      ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

g <- (glist[[1]] + labs(tag = "A")) + (glist[[2]] + labs(tag = "B")) + (glist[[3]] + labs(tag = "C")) + 
  (glist[[4]] + labs(tag = "D")) + (glist[[5]] + labs(tag = "E")) + 
  plot_layout(nrow = 2)
suppressWarnings(ggsave(paste0(out, "/individual_lm_dist_beta.pdf"),
                        g, height = 104, width = 168, units = "mm"))
