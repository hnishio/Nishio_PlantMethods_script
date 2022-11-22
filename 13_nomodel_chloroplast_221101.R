
# Set working directory
setwd("~/Nishio_PlantMethods_script")

# Load packages
library(ggplot2)
library(dplyr)
library(patchwork)
library(data.table)

# Load plot functions
source("functions/nomodel_221101.R")

# Create an output directory
out <- "13_nomodel_chloroplast_221101/"
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

# Predict movement
nomodel(cell_list = cell_list, visual = visual, 
         res_name = "chloroplast", ex_name = "microbeam", 
         unit1 = "micrometer", unit2 = "min", 
         out = out)

