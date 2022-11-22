
# Set working directory
setwd("~/Nishio_PlantMethods_script")

# Load packages
library(ggplot2)
library(patchwork)
library(data.table)
library(ggpubr)
library(exactRankTests)

# Load plot functions
source("functions/lm_dist_start_221101.R")
source("functions/lm_start_dist_221101.R")

# Create an output directory
out <- "12_KFAS_signal_speed_221101"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}



##### Linear regression (x: distance, y: start time) #####

# Load data
cell1 <- read.csv("data/cell1.csv")
cell2 <- read.csv("data/cell2.csv")
cell3 <- read.csv("data/cell3.csv")
cell4 <- read.csv("data/cell4.csv")
visual <- read.csv("data/visual_judgement.csv")
cell_list <- list(cell1, cell2, cell3, cell4)
mvtime <- as.data.frame(fread("08_KFAS_chloroplast_velocity_221101/csv/ssm_KFAS_mvtime.csv"))

# dist vs. start
glist <- lm_dist_start(cell_list = cell_list, mvtime = mvtime,
                       ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

g <- glist[[1]] + glist[[2]] + glist[[3]] + glist[[4]] + glist[[5]] + 
  plot_layout(ncol = 3)

ggsave(paste0(out, "/KFAS_lm_dist_start.pdf"),
       g, height = 110, width = 50*3, units = "mm")





##### Linear regression (x: start time, y: distance) #####

# Load data
cell1 <- read.csv("data/cell1.csv")
cell2 <- read.csv("data/cell2.csv")
cell3 <- read.csv("data/cell3.csv")
cell4 <- read.csv("data/cell4.csv")
visual <- read.csv("data/visual_judgement.csv")
cell_list <- list(cell1, cell2, cell3, cell4)
mvtime <- as.data.frame(fread("08_KFAS_chloroplast_velocity_221101/csv/ssm_KFAS_mvtime.csv"))

# dist vs. start
glist <- lm_start_dist(cell_list = cell_list, mvtime = mvtime,
                       ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

g <- glist[[1]] + glist[[2]] + glist[[3]] + glist[[4]] + glist[[5]] + glist[[6]] + 
  plot_layout(ncol = 3)

suppressWarnings(ggsave(paste0(out, "/KFAS_lm_start_dist.pdf"),
                        g, height = 110, width = 50*3, units = "mm"))





##### Mean and SD of distance #####

# Load data
cell1 <- read.csv("data/cell1.csv")
cell2 <- read.csv("data/cell2.csv")
cell3 <- read.csv("data/cell3.csv")
cell4 <- read.csv("data/cell4.csv")
visual <- read.csv("data/visual_judgement.csv")
cell_list <- list(cell1, cell2, cell3, cell4)
mvtime <- as.data.frame(fread("08_KFAS_chloroplast_velocity_221101/csv/ssm_KFAS_mvtime.csv"))

# Adjust data.frame
df <- mvtime[,1:4]
names(df)[1:4] <- c("cell", "chl", "visual", "predicted")
distance <- NULL
for(i in 1:length(cell_list)){
  distance <- c(distance, as.numeric(cell_list[[i]][max(which(cell_list[[i]]$ex == 0)),-(1:2)]))
}
df$distance <- distance
df <- df[!is.infinite(rowSums(df)),]

# Boxplot of distance
df_mean_sd <- data.frame(category = c("Cell 1", "Cell 2", "Cell 3", "Cell 4"),
                         mean =c(mean(df[df$cell==1,"distance"], na.rm = T), mean(df[df$cell==2,"distance"], na.rm = T),
                                 mean(df[df$cell==3,"distance"], na.rm = T), mean(df[df$cell==4,"distance"], na.rm = T)),
                         sd =c(sd(df[df$cell==1,"distance"], na.rm = T), sd(df[df$cell==2,"distance"], na.rm = T),
                               sd(df[df$cell==3,"distance"], na.rm = T), sd(df[df$cell==4,"distance"], na.rm = T)))

g <- ggplot(df_mean_sd, aes(x = category)) +
  geom_linerange(aes(ymin = mean-sd, ymax = mean+sd)) +
  geom_point(aes(y = mean)) +
  theme_bw() +
  theme(plot.title=element_text(size=7),
        axis.title=element_text(size=7), 
        axis.title.x=element_blank(),
        axis.text=element_text(size=7),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y=expression(paste("Distance (", mu, "m)")))

ggsave(paste0(out, "/distance.pdf"),
       g, height = 35, width = 40, units = "mm")

# Boxplot of predicted start time
df_mean_sd <- data.frame(category = c("Cell 1", "Cell 2", "Cell 3", "Cell 4"),
                         mean =c(mean(df[df$cell==1,"predicted"], na.rm = T), mean(df[df$cell==2,"predicted"], na.rm = T),
                                 mean(df[df$cell==3,"predicted"], na.rm = T), mean(df[df$cell==4,"predicted"], na.rm = T)),
                         sd =c(sd(df[df$cell==1,"predicted"], na.rm = T), sd(df[df$cell==2,"predicted"], na.rm = T),
                               sd(df[df$cell==3,"predicted"], na.rm = T), sd(df[df$cell==4,"predicted"], na.rm = T)))

g <- ggplot(df_mean_sd, aes(x = category)) +
  geom_linerange(aes(ymin = mean-sd, ymax = mean+sd)) +
  geom_point(aes(y = mean)) +
  theme_bw() +
  theme(plot.title=element_text(size=7),
        axis.title=element_text(size=7), 
        axis.title.x=element_blank(),
        axis.text=element_text(size=7),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y="Start time (min)")

ggsave(paste0(out, "/KFAS_predicted_starttime.pdf"),
       g, height = 35, width = 40, units = "mm")





##### Comparison between signaling time and start time #####

# Load data
cell1 <- read.csv("data/cell1.csv")
cell2 <- read.csv("data/cell2.csv")
cell3 <- read.csv("data/cell3.csv")
cell4 <- read.csv("data/cell4.csv")
visual <- read.csv("data/visual_judgement.csv")
cell_list <- list(cell1, cell2, cell3, cell4)
mvtime <- as.data.frame(fread("08_KFAS_chloroplast_velocity_221101/csv/ssm_KFAS_mvtime.csv"))

## Adjust data.frame
df <- mvtime[,1:4]
names(df)[1:4] <- c("cell", "chl", "visual", "predicted")
distance <- NULL
for(i in 1:length(cell_list)){
  distance <- c(distance, as.numeric(cell_list[[i]][max(which(cell_list[[i]]$ex == 0)),-(1:2)]))
}
df$distance <- distance
df <- df[!is.infinite(rowSums(df)),]

# Signaling time
lm_cell1 <- lm(formula = distance~predicted, data = df[df$cell==1,])
lm_cell2 <- lm(formula = distance~predicted, data = df[df$cell==2,])
lm_cell3 <- lm(formula = distance~predicted, data = df[df$cell==3,])
lm_cell4 <- lm(formula = distance~predicted, data = df[df$cell==4,])
lm_allcells <- lm(formula = distance~predicted, data = df)
summary(lm_cell4)

signaling_time_cell1 <- df[df$cell==1,]$distance / lm_cell1$coefficients[2]
start_time_cell1 <- df[df$cell==1,]$predicted
signaling_time_cell2 <- df[df$cell==2,]$distance / lm_cell2$coefficients[2]
start_time_cell2 <- df[df$cell==2,]$predicted
signaling_time_cell3 <- df[df$cell==3,]$distance / lm_cell3$coefficients[2]
start_time_cell3 <- df[df$cell==3,]$predicted
signaling_time_cell4 <- df[df$cell==4,]$distance / lm_cell4$coefficients[2]
start_time_cell4 <- df[df$cell==4,]$predicted
signaling_time_allcells <- df$distance / lm_allcells$coefficients[2]
start_time_allcells <- df$predicted


## Boxplot
df_cat_signaling <- data.frame(
  cell=c(rep("Cell 2", length(start_time_cell2)),
         rep("Cell 3", length(start_time_cell3)), 
         rep("Cell 4", length(start_time_cell4)),
         rep("Cell 2", length(signaling_time_cell2)),
         rep("Cell 3", length(signaling_time_cell3)), 
         rep("Cell 4", length(signaling_time_cell4))),
  STorSIG = c(rep("Total reaction time", length(start_time_cell2)+length(start_time_cell3)+length(start_time_cell4)), 
              rep("Signaling time", length(signaling_time_cell2)+length(signaling_time_cell3)+length(signaling_time_cell4))),
  time=c(start_time_cell2, start_time_cell3, start_time_cell4,
         signaling_time_cell2, signaling_time_cell3, signaling_time_cell4))
df_cat_signaling$STorSIG <- factor(df_cat_signaling$STorSIG, levels=c("Signaling time", "Total reaction time"))

# Wilcoxon test
t2 <- wilcox.exact(signaling_time_cell2, start_time_cell2, "less", paired = T)
t3 <- wilcox.exact(signaling_time_cell3, start_time_cell3, "less", paired = T)
t4 <- wilcox.exact(signaling_time_cell4, start_time_cell4, "less", paired = T)
pvalues <- c(t2$p.value, t3$p.value, t4$p.value)
p.adjust(pvalues, method = "BH")

# Plotting
set.seed(1)
g_sig_st <- ggplot(df_cat_signaling, aes(x = cell, y = time, fill = STorSIG)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(size = 0.02, alpha = 0.3) +
  geom_point(size = 0.02, alpha = 0.3, position=position_jitterdodge()) +
  scale_fill_manual(values=c("darkturquoise", "gray70"))+
  coord_cartesian(ylim = c(0, 57), clip = "off") +
  annotate("text", x = 2, y = 56, label="Hypothesis: signaling time < total reaction time", size = 7/ggplot2::.pt) +
  annotate("text", x = 1:3, y = rep(51.5, 3), label=c("NS","NS","NS"), size = 7/ggplot2::.pt) +
  annotate("segment", x = 0.7, xend = 1.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 1.7, xend = 2.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  annotate("segment", x = 2.7, xend = 3.3, y = 49, yend = 49, size = 0.5,
           arrow = arrow(ends = "both", length = unit(0, "mm"))) +
  #annotate("text", x = 2.5, y = 40, label="Signaling time", size = 7/ggplot2::.pt, col = "darkturquoise") +
  #annotate("text", x = 2.5, y = 30, label="Total reaction time", size = 7/ggplot2::.pt, col = "gray70") +
  theme_bw() +
  theme(plot.title=element_text(size=7),
        axis.title.y=element_text(size=7), 
        axis.title.x=element_blank(),
        axis.text=element_text(size=7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.75, 0.65),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        #legend.margin=margin(0, 0, 0, -5),
        legend.spacing.x = unit(0.7, 'mm'),
        legend.key.width = unit(2, "mm"),
        plot.tag = element_text(size = 12, face = "bold")) + 
  labs(y="Time (min)")


ggsave(paste0(out, "/KFAS_signalingtime_starttime.pdf"),
       g_sig_st, height = 70, width = 68, units = "mm")

