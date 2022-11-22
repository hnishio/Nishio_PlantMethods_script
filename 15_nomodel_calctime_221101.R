
# Set working directory
setwd("~/Nishio_PlantMethods_script")

# Load packages
library(patchwork)
library(data.table)
library(tictoc)

# Load SSM functions
source("functions/calctime_nomodel_221101.R")

# Create output directories
out <- "15_nomodel_calctime_221101"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Load data
cell1 <- read.csv("data/cell1.csv")
cell2 <- read.csv("data/cell2.csv")
cell3 <- read.csv("data/cell3.csv")
cell4 <- read.csv("data/cell4.csv")
visual <- read.csv("data/visual_judgement.csv")
cell_list <- list(cell2)

# Predict movement
tic()
nomodel(cell_list = cell_list, visual = visual, 
        res_name = "chloroplast", ex_name = "microbeam", 
        unit1 = "micrometer", unit2 = "min", out = out)
tictoc_time <- toc()
calc_time <- data.frame(
  calc_time = as.numeric(tictoc_time$toc - tictoc_time$tic))
fwrite(calc_time, file = paste0(out, "/nomodel_calctime.csv"))

