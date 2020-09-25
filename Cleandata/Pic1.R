# This script is used to plot site information
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(xlsx)

# Set working directory
setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/Cleandata")

# Loading Data -----------------------------------------------------------------
home <- "Data.xlsx"
patient <- read.xlsx2(home, sheetIndex = 1, stringsAsFactors = FALSE,
  colClasses = c(rep("character", 6), rep("numeric", 29)))

patient$Gender <- factor(patient$Gender, labels = c("Male", "Female"))
patient$Sites <- factor(patient$Sites, levels = 1:15)

rm(home)

## Site information -------------------------------------------------------------
# plot theme -------------------------------------------------------------------
themes <- theme(
  legend.position  = c(0.9, 0.9),
  legend.title     = element_blank(),
  text             = element_text(size = 25),
  axis.text.x      = element_text(size = 20, face = "bold", color = "black"),
  axis.text.y.left = element_text(size = 15, face = "bold"),
  axis.line        = element_line(color = "black", size = 1, linetype = "solid"),
  panel.background = element_rect(fill = "white")
  )

# Sample distribution per site and Gender distribution per site; Figure S2
table1 <- patient %>% select(Sites) %>% group_by(Sites) %>% summarise(Count = n())

p1 <- ggplot(table1, aes(Sites, Count, fill = Sites)) +
  geom_bar(stat = "identity", show.legend = FALSE, width = 0.55) +
  labs(x = "Sites", y = "Count") +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  geom_text(label = table1$Count, vjust = 00, alpha = 0.5, size = 8) +
  scale_fill_manual(values = c(brewer.pal(11, "Set3"), brewer.pal(4, "Set2"))) +
  themes

p2 <- ggplot(patient, aes(Sites, fill = Gender)) +
  geom_bar(stat = "count", width = 0.7, alpha = 0.8) +
  labs(y = "") +
  scale_y_continuous(breaks = seq(0, 150, 5)) +
  scale_fill_brewer(palette = "Set1") +
  themes

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S2.pdf", width = 50, height = 60, units = "cm")