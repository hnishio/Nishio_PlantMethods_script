
# Set working directory
setwd("~/Nishio_PlantMethods_script")

# Load packages
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(data.table)
library(cmdstanr)
library(bayesplot)
library(tictoc)
set_cmdstan_path("~/cmdstan/")

# Load SSM functions
source("functions/calctime_ssm_individual_221101.R")

# Create output directories
out <- "04_CmdStan_calctime_221101"
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

# Execution of state-space modeling
tic()
calctime_ssm_individual(cell_list = cell_list, visual = visual, 
                        res_name = "chloroplasts", ex_name = "microbeam", 
                        unit1 = "micrometer", unit2 = "min")
tictoc_time <- toc()
calc_time <- data.frame(
  calc_time = as.numeric(tictoc_time$toc - tictoc_time$tic))
fwrite(calc_time, file = paste0(out, "/CmdStan_calctime.csv"))

