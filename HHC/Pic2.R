# This script was used for plot subtype comparisons
library(xlsx)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggradar)        # Radar plot
library(RColorBrewer)
library(patchwork)


# Set working directory
setwd("E:/WorkingSpace/R/Paper/Depression/Subtypes/HHC")

# Load data --------------------------------------------------------------------
home <- "E:/WorkingSpace/R/Paper/Depression/Subtypes/HHC/HHC.xlsx"
patient <- read.xlsx2(home, 1,  stringsAsFactors = FALSE,
  colClasses = rep(c("character", "numeric"), times = c(8, 29)))
# Modify variables
patient$Gender <- factor(patient$Gender, labels = c("Male", "Female"))
patient$Sites <- factor(patient$Sites, levels = 1:15)
patient$Level1 <- factor(patient$Level1, levels = c(2, 1))
patient$Level2 <- factor(patient$Level2, levels = c(3, 2, 1, 4))

hamd_vars <- c("Depressed\nMood", "Guilt", "Suicide", "Early\nInsomnia", 
  "Middle\nInsomnia", "Late\nInsomnia", "Work\nInterests", "Retardation", 
  "Agitation", "Psychic\nAnxiety", "Somatic\nAnxiety", "Gastrointestinal", 
  "General\nSomatic", "Loss of\nLibido", "Hypochondriasis", "Weight\nLoss", 
  "Loss of\nInsight")
rm(home)

# Plot subtype comparisons -----------------------------------------------------
# Set the parameters of radar plot
radar_plot1 <- function(data, colors, labels) {
  ggradar(data, grid.max = 1.5, grid.mid = 0,  grid.min = -1.5, 
    values.radar = c(-1.5, 0, 1.5), axis.labels = labels, grid.line.width = 1.5, 
    gridline.min.linetype = 3, gridline.mid.linetype = 1, gridline.max.linetype = 1,  
    gridline.mid.colour = "grey", grid.label.size = 45, axis.label.offset = 1.1, 
    gridline.label.offset = 1.1, group.line.width = 5, group.point.size  = 5, 
    group.colours = colors, axis.label.size = 36,  legend.position = "bottom", 
    legend.text.size = 1, background.circle.colour = "white")
}

pca_levels <- subset(patient, select = c(Level1:Level2, Guilt:Interests.Loss))
pca_level1 <- pca_levels %>% select(-Level2) %>% group_by(Level1) %>% summarise_all(mean)
pca_level2 <- pca_levels %>% group_by(Level1, Level2) %>% summarise_all(mean)

# 1 PCA components--------------------------------------------------------------
# level 1
radar_plot1(pca_level1, c("#B30000","#2166AC"),
  labels = c(
    "Guilt", 
    "          ***\nInsomnia", 
    "                             ***\nSomatic and Anxiety", 
    "              ***\nInterests\nLoss"
  )
)
ggsave("Radar1_12_pca.pdf", width = 120, height = 80, unit = "cm")

# level 2
# SAI+ vs SAI-
radar_plot1(pca_level2[c(1,2), -1], c("#4393C3", "#92C5DE"),
  labels = c(
    "      ***\nGuilt",
    "          ***\nInsomnia",
    "                             ***\nSomatic and Anxiety",
    "              ***\nInterests\nLoss"
  )
)
ggsave("Radar2_32_pca.pdf", width = 120, height = 80, unit = "cm")
# LOI+ vs LOI-
radar_plot1(pca_level2[c(3,4), -1], c("#EF3B2C", "#FC9272"),
  labels = c("Guilt***", "Insomnia***", "Somatic and Anxiety***", "***\nInterests Loss"))
ggsave("Radar2_14_pca.pdf", width = 120, height = 80, unit = "cm")

# Set the parameters of radar plot ---------------------------------------------
radar_plot2 <- function(data, colors, labels) {
  ggradar(data, grid.max = 7, grid.mid = 3.5,  grid.min = 0, values.radar = c(0, 3.5, 7),
    axis.labels = labels, grid.line.width = 1.5, gridline.min.linetype = 3, 
    gridline.mid.linetype = 1, gridline.max.linetype = 1, gridline.mid.colour = "grey", 
    grid.label.size = 15, axis.label.offset = 1.1, gridline.label.offset = 1, 
    group.line.width = 5, group.point.size  = 5, group.colours = colors,
    axis.label.size = 15, plot.legend = TRUE, legend.position = "right", 
    legend.text.size = 15, background.circle.colour = "white")
}

# 2 Sub-dimensions -------------------------------------------------------------
subd_levels <- subset(patient, select = c(Level1:Level2, CD:NVSM))
subd_level1 <- subd_levels %>% select(-Level2) %>% group_by(Level1) %>% summarise_all(mean)
subd_level2 <- subd_levels %>% group_by(Level1, Level2) %>% summarise_all(mean)

# level 1
# LOI vs SAI
radar_plot2(subd_level1, c("#B30000","#2166AC"),
  labels = c("CD**", "ANX*", "NVSM***"))
ggsave("Radar1_12_subd.pdf", width = 60, height = 60, units = "cm")

# level 2
# SAI+ vs SAI-
radar_plot2(subd_level2[c(1,2), -1], c("#4393C3", "#92C5DE"),
  labels = c("CD***", "ANX***", "NVSM***"))
ggsave("Radar2_32_subd.pdf",  width = 60, height = 60, units = "cm")
# LOI+ vs LOI-
radar_plot2(subd_level2[c(3,4), -1], c("#EF3B2C", "#FC9272"),
  labels = c("CD***", "ANX***", "NVSM***"))
ggsave("Radar2_14_subd.pdf",  width = 60, height = 60, units = "cm")

# 3 Symptoms -------------------------------------------------------------------
with(patient, ftable(Level1, Level2))
hamd_levels <- patient[, c(7, 8, 13:29)]
hamd_levels$Level1 <- factor(hamd_levels$Level1, level = 1:2, 
  labels = c("LOI", "SAI"))
hamd_levels$Level2 <- factor(hamd_levels$Level2, level = c(2, 3, 1, 4), 
  labels = c("SAI+", "SAI-", "LOI+", "LOI-"))

hamd_level1 <- hamd_levels %>% select(-Level2) %>% group_by(Level1) %>% summarise_all(mean) %>% melt() 
hamd_level1$variable <- as.numeric(hamd_level1$variable)

hamd_level2_LOI <- hamd_levels %>% group_by(Level1, Level2) %>% summarise_all(mean) %>% 
  subset(Level1 == "LOI") %>% melt()
hamd_level2_LOI$variable <- as.numeric(hamd_level2_LOI$variable)

hamd_level2_SAI <- hamd_levels %>% group_by(Level1, Level2) %>% summarise_all(mean) %>% 
  subset(Level1 == "SAI") %>% melt()
hamd_level2_SAI$variable <- as.numeric(hamd_level2_SAI$variable)

# Plot theme -------------------------------------------------------------------
themes <-  theme(text = element_text(size = 30),
  axis.text.y = element_text(size = 7, face = "bold"),
  legend.position = c(0.9, 0.9),
  legend.title = element_blank(),
  axis.text.x = element_text(size = 17, face = "bold", color = "black"),
  axis.text.y.left = element_text(size = 20),
  axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
  panel.background = element_rect(fill = "white"))

p1 <- ggplot(hamd_level1, aes(variable, value, color = Level1)) + 
  geom_line(size = 2) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = 1:17, labels = hamd_vars) + 
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  labs(x = NULL, y = NULL) +
  scale_colour_manual(values = c("#B30000", "#2166AC")) +
  themes

p2 <- ggplot(hamd_level2_LOI, aes(variable, value, color = Level2)) + 
  geom_line(size = 2) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = 1:17, labels = hamd_vars) + 
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  labs(x = NULL, y = NULL) + 
  scale_colour_manual(values = c("#EF3B2C", "#FC9272")) +
  themes

p3 <- ggplot(hamd_level2_SAI, aes(variable, value, color = Level2)) + 
  geom_line(size = 2) +
  geom_point(size = 4) +
  scale_x_continuous(breaks = 1:17, labels = hamd_vars) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  labs(x = NULL, y = NULL) +
  scale_colour_manual(values = c("#4393C3", "#92C5DE")) +
  themes

p1 + p2 + p3 + plot_layout(ncol = 1)
ggsave("Figure S5.pdf", width = 70, height = 70, units = "cm")