
# Set working directory
setwd("~/Nishio_PlantMethods_script")

# For the first time usage, install cmdstanr by
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# See the details here: https://mc-stan.org/cmdstanr/articles/cmdstanr.html

# Load packages
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(data.table)
library(cmdstanr)
library(bayesplot)
set_cmdstan_path("~/cmdstan/")

# Load SSM functions
source("functions/ssm_common_221101.R")

# Create output directories
out <- "07_CmdStan_chloroplast_velocity_common_221101"
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

# Execution of state-space modeling
ssm_common(cell_list = cell_list, mvtime = mvtime,
           res_name = "chloroplast", ex_name = "microbeam", 
           unit1 = "micrometer", unit2 = "min",
           out = out)

