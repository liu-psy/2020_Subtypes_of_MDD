# The script was an integrated version

################################################################################
################################# Data clean ###################################
################################################################################
library(xlsx)
library(dplyr)
library(mice)                 # Imputation
library(psych)                # PCA
library(reshape2)
library(ggplot2)
library(corrplot)


# Set working directory
setwd("E:/WorkingSpace/R/Paper/Depression/Subtypes/All")

# Loading Data -----------------------------------------------------------------
# Data of REST-meta-MDD
home1 <- "E:/WorkingSpace/R/Paper/Depression/Subtypes/All/Depression.xlsx"
# Data of SWU
home2 <- "E:/WorkingSpace/R/Paper/Depression/Subtypes/All/SWU_Data.xlsx"

patient_rest_project <- read.xlsx2(home1, sheetIndex = 1, stringsAsFactors = FALSE,
  colClasses = c(rep("character", 2), rep("numeric", 2), rep("character", 2), rep("numeric", 27)))
# Exclude irrelevant variables
patient_rest_project <- subset(patient_rest_project, select = -c(HAMD, HAMA, X18:X24))
patient_swu <- read.xlsx2(home2, sheetIndex = 1, stringsAsFactors = FALSE,
  colClasses = c(rep("character", 2), rep("numeric", 2), "character", "numeric",
    "character", rep("numeric", 17)))

# Combine Datasets -------------------------------------------------------------
# Patients without HAMD-17 were excluded
max <- apply(patient_rest_project[, 8:24], 1, max, na.rm = TRUE)
patient_rest_project <- patient_rest_project[which(max != -Inf), ]

# Replace the original data of SWU(S20)
patient <- rbind(patient_rest_project[1:grep("S21", patient_rest_project$ID)[1] - 1, ], patient_swu,
  patient_rest_project[grep("S21", patient_rest_project$ID)[1] : nrow(patient_rest_project), ])

rm(home1, home2, patient_swu, patient_rest_project, max)

# Data clean -------------------------------------------------------------------
# Remove duplicated patients
patient <- patient[!duplicated(patient %>% select(X1:X17)), ]

# Control age range (18-65)
patient <- subset(patient, Age >= 18 & Age <= 65)
range(patient$Age)

# Set impossible values as NA
hamd <- subset(patient, select = c(X1:X17))
sum(is.na(hamd))

make_na <- function(data) {
  # Those items socres range 0-4
  group1 <- c(1, 2, 3, 7, 8, 9, 10, 11, 15)
  # Those items scores range 0-2
  group2 <- c(1:17)[-group1]
  
  data[, group1] <- apply(data[, group1], 2, function(x) {x[!x %in% 0:4] <- NA; x})
  data[, group2] <- apply(data[, group2], 2, function(x) {x[!x %in% 0:2] <- NA; x})
  cat("Total NAs:", sum(is.na(data)))
  return(data)
}
hamd <- make_na(hamd)

# Check by eyes
lapply(hamd, unique)
# Check the missing pattern
md.pattern(hamd)

# Interpolation using random forest
imputed_hamd <- complete(mice(hamd, m = 10, method = "rf", seed = 7))
# Check again
lapply(imputed_hamd, unique)

patient[, 8:24] <- imputed_hamd

rm(hamd, imputed_hamd, make_na)

# Control HAMD-17 total score, HAMD-17 total socre > 7
patient$HAMD <- rowSums(subset(patient, select = c(X1:X17)))
patient <- subset(patient, HAMD > 7)
range(patient$HAMD)

# The patients whose HAMD-17 total score was greater than 3 standard deviations were excluded
patient$Scale <- scale(patient$HAMD)
patient <- subset(patient, select = -c(Scale), subset = abs(patient$Scale) <= 3)

# Modify the site variable
get_sites <- function(data) {
  name_list <- strsplit(data$ID, "-", fixed = FALSE)
  sites <- vector(length = length(name_list))
  for (i in seq(name_list)) {
    sites[i] <- name_list[[i]][1]
  }
  return(sites)
}
sites <- get_sites(patient)
# 15 sites
site_list <- unique(sites)
patient$Sites <- as.numeric(factor(sites, levels = site_list))

rm(site_list, sites, get_sites)

# Add Extra variables ----------------------------------------------------------
# Add head motion of fMRI data path of every patient
headmotion <- function(type) {
  home <- "F:/MDD/RealignParameter/"
  id <- patient$ID
  path <- paste0(home, id)
  list <- vector("list", length = nrow(patient))
  for (i in seq(list)) {
    list[i] <- read.table(paste0(path[i], "/", type, id[i], ".txt"))
  }
  patient$HeadMotion <- sapply(list, mean)
  return(patient)
}
# Jenkinson
patient <- headmotion("FD_Jenkinson_")

# Add fMRI file dir, Weighted DC
add_dir <- function(type) {
  home <- "F:/MDD/DegreeCentrality_FunImgARCWF/"
  filename <- paste0("DegreeCentrality_Positive", type, "SumBrainMap_")
  dir <- paste0(home, "f", filename, patient$ID, ".nii")
  patient$Dir <- dir
  return(patient)
}
patient <- add_dir("Weighted")

rm(add_dir, headmotion)

# Correlation between symptoms -------------------------------------------------
hamd <- subset(patient, select = c(X1:X17))
cor_hamd <- round(cor(hamd), 3)

# Correlation range
range(cor_hamd[upper.tri(cor_hamd)])

# Plot Correlation matrix; Figure S11
pdf("Figure S11.pdf", width = 10, height = 8)
corrplot(cor_hamd, method = "square")
dev.off()

# PCA for reducing demensions --------------------------------------------------
# Parallel analysis
PA <- fa.parallel(hamd, fa = "pc", sim = FALSE, n.iter = 10000, quant = 0.95,
  main = "Parallel Analysis")
# 4 components come from parallel analysis
PA$ncomp

threshold <- apply(PA$values, 2, quantile, prob = 0.95)
PA$pc.values > threshold

# Plot for parallel analysis; Figure S4
pa_table <- data.frame(t(rbind(PA$pc.values, PA$pc.simr)))
names(pa_table) <- c("Acutual Data", "Resample Data")
pa_table <- melt(pa_table)
pa_table$Components <- rep(1:17, times = 2)

ggplot(pa_table, aes(Components, value, color = variable)) +
  geom_line(size = 2, linetype = rep(1:2, each = 17)) +
  geom_point(size = 5) +
  scale_x_continuous(breaks = 1:17) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  labs(y = "Eigen value") +
  theme(title = element_text(size = 30),
    text = element_text(size = 30),
    legend.position = c(0.8, 0.8),
    legend.title = element_blank(),
    axis.text.y = element_text(size = 30),
    axis.text.x = element_text(size = 25),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    panel.background = element_rect(fill = "white"))
ggsave("Figure S4.pdf", width = 15, height = 8)

# promax rotation
PCA <- principal(hamd, nfactors = 4, rotate = "promax", scores = TRUE,
  oblique.scores = FALSE)
# Check the correlation matrix of components,
# and all the correlation coefficients are smaller than 0.32
round(PCA$r.scores, 2)

# varimax rotation
PCA <- principal(hamd, nfactors = 4, rotate = "varimax", scores = TRUE,
  oblique.scores = FALSE)
PCA$loadings
PCA$r.scores
describe(PCA$communality)

hamd_vars <- c("Depressed Mood", "Guilt", "Suicide", "Early Insomnia",
  "Middle Insomnia",  "Late Insomnia", "Work Interests", "Retardation",
  "Agitation", "Psychic Anxiety", "Somatic Anxiety", "Gastrointestinal",
  "General Somatic", "Loss of Libido", "Hypochondriasis", "Weight Loss",
  "Loss of Insight")
components <- c("Guilt", "Insomnia", "Somatic and Anxiety", "Interests Loss")

# Modify components position
loading_matrix <- PCA$loadings[, c(4, 1, 2, 3)]
colnames(loading_matrix) <- components

# Plot loading matrix, Figure 1
# This code is inspired by https://rpubs.com/danmirman/plotting_factor_analysis
loading <- as.data.frame(loading_matrix)
loading$vars <- hamd_vars %>% factor(levels = rev(hamd_vars))
loading <- melt(loading, id = "vars")
loading$variable <- factor(loading$variable)

ggplot(loading, aes(vars, value, fill = value)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~ variable, nrow = 1) +
  labs(x = NULL, y = "Loading Strength") +
  coord_flip() +
  scale_fill_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0,
    guide = FALSE) +
  theme(
    title = element_text(size = 20),
    text = element_text(size = 20),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    panel.background = element_rect(fill = "white")
  )
ggsave("Figure 1.pdf", width = 15, height = 8)

patient <- cbind(patient, PCA$scores[, c(4, 1, 2, 3)])
names(patient)[29:32] <- components

# Add Sub-dimensions -----------------------------------------------------------
patient <- patient %>% mutate("CD" = (X1 + X7), "ANX" = (X9 + X10 + X11 + X15),
  "NVSM" = (X6 + X12))

# Modify Data ------------------------------------------------------------------
patient <- patient %>% select(ID, Gender, Sites, Episode, Medication, Dir, Age,
  Education, Month, HeadMotion, everything())

# Remove redundant variables
rm(PA, PCA, loading, loading_matrix, components, threshold, hamd, cor_hamd, pa_table)

##################################### Test1 ####################################
library(effsize)     # Cohen's d
library(sjstats)     # Crame's V

patient$Gender <- factor(patient$Gender, levels = c(1, 2), labels = c("Male", "Female"))
patient$Sites <- factor(patient$Sites, levels = 1:15)
names(patient)[11:27] <- paste(1:17, hamd_vars)

table(patient$Gender)
table(patient$Sites)
table(patient$Sites, patient$Gender)
table(patient$Sites, patient$Episode)
table(patient$Sites, patient$Episode)
table(patient$Sites, patient$Medication)

patient %>% select(Age, Education, HAMD, Sites) %>% group_by(Sites) %>%
  summarise_all(mean)
patient %>% select(Age, Education, HAMD, Sites) %>% group_by(Sites) %>%
  summarise_all(sd)

lapply(c("Age", "Education", "HAMD"), function(x) describe(patient[[x]]))

# Gender test ------------------------------------------------------------------
Male <- patient[patient$Gender == "Male", ]
Female <- patient[patient$Gender == "Female", ]

# Gender Distribution
chisq.test(table(patient$Gender, patient$Sites)) 

# Multiply t-tests
m_ttest <- function(data1, data2, x) {
  results <- matrix(nrow = length(x), ncol = 7)
  colnames(results) <- c("Stat", "df", "P", "Lower", "Upper", "d", "P(corrected)")
  rownames(results) <- x
  for (i in seq(x)) {
    t <- t.test(data1[[x[i]]], data2[[x[i]]])
    d <- cohen.d(data1[[x[i]]], data2[[x[i]]])
    results[i, 1:6] <- c(abs(round(t$statistic, 2)), abs(round(t$parameter)),
      abs(round(t$p.value, 3)), round(t$conf.int[1], 2),
      round(t$conf.int[2], 2), abs(round(d$estimate, 3)))
  }
  # fdr correction for p-vaules
  results[, 7] <- round(p.adjust(results[, 3], method = "fdr"), 3)
  results <- results[, c(1:3,7,4,5,6)]
  return(results)
}

# Gender differences
test_vars <- names(patient)[c(7,8, 11:28)]
output <- m_ttest(Male, Female, test_vars)

write.xlsx2(output, "Gender Differences.xlsx")

# modify HAMD-17
names(patient)[11:27] <- paste0("X", 1:17)

# Remove redundant variables
rm(output, Female, Male, test_vars, m_ttest)

##################################### Pic1 #####################################
library(RColorBrewer)
library(patchwork)

## Site information -------------------------------------------------------------
# plot theme -------------------------------------------------------------------
themes <-  theme(text = element_text(size = 25),
  legend.position = c(0.9, 0.9),
  legend.title = element_blank(),
  axis.text.x = element_text(size = 20, face = "bold", color = "black"),
  axis.text.y.left = element_text(size = 15, face = "bold"),
  axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
  panel.background = element_rect(fill = "white"))

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
  scale_y_continuous(breaks = seq(0, 150, 5)) +
  labs(y = "") +
  scale_fill_brewer(palette = "Set1") +
  themes

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S2.pdf", width = 50, height = 60, units = "cm")

# Remove redundant variables
rm(p1, p2, table1, themes)

################################################################################
####################### Hybrid Hierarchical Clustering #########################
################################################################################
library(hybridHclust)          # Hybrid Hierarchical Clustering
library(dynamicTreeCut)        # Dynamic Tree Cut algorithm
library(gplots)                # Heatmap
library(dendextend)            # Dendrogram
library(clusteval)             # Robustness of clusering

# Euclidean Distant Matrix
hamdDist <- patient %>% select(X1:X17) %>% scale() %>% dist() %>% as.matrix()

set.seed(1234)
HHC <- hybridHclust(hamdDist, trace = T)

# DynamicTreeCut ---------------------------------------------------------------
cut <- cutreeHybrid(HHC, distM = hamdDist)
# k = 4 is the best cut desicions
unique(cut$labels)

# Visual inspection ------------------------------------------------------------
# k = 2 looks good
plot(HHC, labels = FALSE, hang = -1)

# cut trees --------------------------------------------------------------------
# k = 2, 4
labels <- cutree(HHC, k = c(2, 4))
table(labels[, 1])
table(labels[, 2])

# Robustness of HHC ------------------------------------------------------------
labels2 <- cutree(hclust(as.dist(hamdDist), method = "average"), k = 2)
labels4 <- cutree(hclust(as.dist(hamdDist), method = "average"), k = 4)

# jaccard index
round(cluster_similarity(labels1 = labels[, 1], labels2 = labels2, similarity = "jaccard"), 2)
round(cluster_similarity(labels1 = labels[, 2], labels2 = labels4, similarity = "jaccard"), 2)
# rand index
round(cluster_similarity(labels1 = labels[, 1], labels2 = labels2, similarity = "rand"), 2)
round(cluster_similarity(labels1 = labels[, 2], labels2 = labels4, similarity = "rand"), 2)

# Plot Dendrogram; Figure 2A ---------------------------------------------------
pdf("Dendrogram.pdf", width = 10, height = 10)
HHC %>% as.dendrogram() %>% set("branches_k_color", c("#2166AC", "#B30000"), k = 2) %>%
  set("branches_k_color", c("#92C5DE", "#4393C3", "#EF3B2C", "#FC9272"), k = 4) %>%
  set("branches_k_lty", rep("twodash", 2), k = 2) %>% set("branches_lwd", 10) %>%
  set("labels", rep("", nrow(hamdDist))) %>% plot()
abline(h = 1.5e5, lwd = 4, lty = 4)
text(20, 1.7e5, "k = 4", cex = 2, font = 2)
abline(h = 3.5e5, lwd = 4, lty = 4)
text(20, 3.8e5, "k = 2", cex = 2, font = 2)
dev.off()

# Heatmap of Eculidean distant matrix; Figure 2B -------------------------------
color <- rep(c("#92C5DE", "#4393C3", "#EF3B2C", "#FC9272"), times = c(227, 178, 116, 129))
color <- color[order(HHC$order)]

pdf("Heatmap.pdf", width = 15, height = 18)
heatmap.2(hamdDist, Rowv = as.dendrogram(HHC), Colv = as.dendrogram(HHC),
  distfun = function(x) x, hclustfun = hybridHclust, dendrogram = "none",
  symm = TRUE, revC = TRUE, col = bluered(16), trace = "none", tracecol = NULL,
  labRow = NA, labCol = NA, keysize = 1, density.info = "none", key.title = NA, key.xlab = NA,
  key.xtickfun = function() {list(side = 1, padj = 0.8, labels = seq(0, 12, 2),
    at = seq(0, 1, length.out = 7), cex.axis = 6, tck = 0.5, tcl = 1, font = 2,
    lwd.ticks = 4, mtext("Euclidean Distance", side = 1, padj = 2, cex = 5, font = 2))},
  key.par = list(mai = c(2, 0, 0, 0.65)), ColSideColors = color,
  lmat = rbind(c(4, 1), c(3, 2), c(0, 5)), lhei = c(0.2, 4, 0.8), lwid = c(0.2, 4))
dev.off()

# Modify data ------------------------------------------------------------------
patient <- data.frame(patient, Level1 = labels[, 1], Level2 = labels[, 2], check.names = FALSE) %>%
  select(names(patient)[1:6], Level1, Level2, everything())

# Remove redundant variables
rm(cut, hamdDist, HHC, labels, labels2, labels4, color)

################################# Test 2 #######################################
# Set character to factor
patient[, c(7:8)][] <- lapply(patient[, c(7:8)], factor)

# Modify variable names --------------------------------------------------------
# HAMD-17
names(patient)[13:29] <- paste(1:17, hamd_vars)
# PCA
pca_vars <- names(patient)[31:34]
# Sub-dimensions
sub_vars <- names(patient)[35:37]
# symptoms
sym_vars <- names(patient)[13:29]

# Check subgroups --------------------------------------------------------------
# 1vs2 in level1
with(patient, ftable(Level1))
# 1vs4 2vs3 in level2
with(patient, ftable(Level1, Level2))

# Gender distribution ----------------------------------------------------------
table1 <- table(patient$Level1, patient$Gender)
table2 <- table(patient$Level2, patient$Gender)

lapply(list(table1, table2[c(1,4), ], table2[c(2,3), ]), chisq.test)
# Effect size, Cramer's V
lapply(list(table1, table2[c(1,4), ], table2[c(2,3), ]), cramer)

# Multiply t-tests -------------------------------------------------------------
mtest <- function(Level, cluster1, cluster2, x) {
  group1 <- patient[patient[[Level]] %in% cluster1, ]
  group2 <- patient[patient[[Level]] %in% cluster2, ]
  results <- matrix(nrow = length(x), ncol = 7)
  colnames(results) <- c("Stat", "df", "P", "Lower", "Upper", "d", "P(corrected)")
  rownames(results) <- x
  
  for (i in seq(x)) {
    t <- t.test(group1[[x[i]]], group2[[x[i]]])
    d <- cohen.d(group1[[x[i]]], group2[[x[i]]])
    results[i, 1:6] <- c(round(t$statistic, 1), abs(round(t$parameter)),
      abs(round(t$p.value, 3)), round(t$conf.int[1], 2),
      round(t$conf.int[2], 2), abs(round(d$estimate, 2)))
  }
  # fdr correction for p-vaules
  results[, 7] <- round(p.adjust(results[, 3], method = "fdr"), 3)
  results <- results[, c(1:3,7,4,5,6)]
  return(results)
}

# Two main subtypes ------------------------------------------------------------
# LOI vs SAI
mtest("Level1", 1, 2, pca_vars)
mtest("Level1", 1, 2, sub_vars)
mtest("Level1", 1, 2, sym_vars)
mtest("Level1", 1, 2, c("Age", "Education", "HAMD"))

# Four minor subtpyes ----------------------------------------------------------
# LOI+ vs LOI-
mtest("Level2", 1, 4, pca_vars)
mtest("Level2", 1, 4, sub_vars)
mtest("Level2", 1, 4, sym_vars)
mtest("Level2", 1, 4, c("Age", "Education", "HAMD"))
# SAI+ vs SAI-
mtest("Level2", 2, 3, pca_vars)
mtest("Level2", 2, 3, sub_vars)
mtest("Level2", 2, 3, sym_vars)
mtest("Level2", 2, 3, c("Age", "Education", "HAMD"))

# modify HAMD-17
names(patient)[13:29] <- paste0("X", 1:17)

############################# Pic 2 ############################################
library(ggradar)        # Radar plot

patient$Level1 <- factor(patient$Level1, levels = c(2, 1))
patient$Level2 <- factor(patient$Level2, levels = c(3, 2, 1, 4))

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

pca_levels <- subset(patient, select = c(Level1:Level2, 31:34))
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
vars <- c("Depressed\nMood", "Guilt", "Suicide", "Early\nInsomnia", 
  "Middle\nInsomnia", "Late\nInsomnia", "Work\nInterests", "Retardation", 
  "Agitation", "Psychic\nAnxiety", "Somatic\nAnxiety", "Gastrointestinal", 
  "General\nSomatic", "Loss of\nLibido", "Hypochondriasis", "Weight\nLoss", 
  "Loss of\nInsight")

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
  scale_x_continuous(breaks = 1:17, labels = vars) + 
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  labs(x = NULL, y = NULL) +
  scale_colour_manual(values = c("#B30000", "#2166AC")) +
  themes

p2 <- ggplot(hamd_level2_LOI, aes(variable, value, color = Level2)) + 
  geom_line(size = 2) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = 1:17, labels = vars) + 
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  labs(x = NULL, y = NULL) + 
  scale_colour_manual(values = c("#EF3B2C", "#FC9272")) +
  themes

p3 <- ggplot(hamd_level2_SAI, aes(variable, value, color = Level2)) + 
  geom_line(size = 2) +
  geom_point(size = 4) +
  scale_x_continuous(breaks = 1:17, labels = vars) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  labs(x = NULL, y = NULL) +
  scale_colour_manual(values = c("#4393C3", "#92C5DE")) +
  themes

p1 + p2 + p3 + plot_layout(ncol = 1)
ggsave("Figure S5.pdf", width = 70, height = 70, units = "cm")

# Remove redundant variables
rm(p1, p2, p3, pca_level1, pca_level2, pca_levels, subd_level1, subd_level2,
  subd_levels, hamd_level1, hamd_level2_LOI, hamd_level2_SAI, hamd_levels, 
  themes, table1, table2, radar_plot1, radar_plot2, sym_vars, vars, pca_vars,
  sub_vars, mtest)

################################################################################
########################### Network Analysis ###################################
################################################################################
    
# This part was mainly derived from https://github.com/Melovainio/Network_depression       
library(qgraph)
library(igraph)
library(EstimateGroupNetwork)
library(mgm)
library(NetworkComparisonTest)
library(bootnet)
library(NetworkToolbox)
library(networktools)

# Extract symptom of subtypes
LOI <- subset(patient, Level1 == 1, X1:X17)
SAI <- subset(patient, Level1 == 2, X1:X17)

### Network Analysis -----------------------------------------------------------
# Estimated polychoric correlations among symptoms -----------------------------
LOI_cor <- cor_auto(LOI)
SAI_cor <- cor_auto(SAI)

# Estimate Joint Graphical Lasso model
net <- EstimateGroupNetwork(list(LOI, SAI), inputType = "list.of.dataframes",
  method = "crossvalidation", criterion = "ebic", gamma = 0.5,
  simplifyOutput = FALSE, seed = 100, ncores = 16)
net$network

# 1 for LOI; 2 for SAI
nw1 <- getWmat(qgraph(net$network[[1]], sampleSize = nrow(LOI), DoNotPlot = TRUE))
nw2 <- getWmat(qgraph(net$network[[2]], sampleSize = nrow(SAI), DoNotPlot = TRUE))

Max <- max(c(nw1, nw2))
Max
# layout of nodes
L <- averageLayout(nw1, nw2)

# Preidictability Networks -----------------------------------------------------
# Mixed Graphical Models -------------------------------------------------------
# 1 for LOI; 2 for SAI
set.seed(1)
fit1 <- mgm(LOI, type = rep('g', 17), level = rep(1, 17), lambdaSel = 'CV',
  ruleReg = 'OR')
fit2 <- mgm(SAI, type = rep('g', 17), level = rep(1, 17), lambdaSel = 'CV',
  ruleReg = 'OR')

pred1 <- predict(object = fit1, data = LOI, errorCon = 'R2')
pred2 <- predict(object = fit2, data = SAI, errorCon = 'R2')

pred1$error$R2
hamd_vars[order(pred1$error$R2, decreasing = TRUE)]
pred2$error$R2
hamd_vars[order(pred2$error$R2, decreasing = TRUE)]

# Average node predictability
round(mean(pred1$error$R2), 2)
round(mean(pred2$error$R2), 2)

# Plot networks; Figure 3A
pdf("Figure 3a.pdf", width = 14, height = 8)
par(mfrow = c(1, 2))
gr1 <- qgraph(net$network[[1]], layout = L, title = "LOI", title.cex = 2.5,
  maximum = Max,  theme = "Hollywood", pieColor = "#FC8D62", pie = pred1$error$R2,
  border.width = 2, vsize = 9, label.cex = 1, tuning = 0.25)
gr2 <- qgraph(net$network[[2]], layout = L, title = "SAI", title.cex = 2.5,
  maximum = Max, theme = "Hollywood", pieColor = "#FC8D62", pie = pred2$error$R2,
  border.width = 2, vsize = 9, label.cex = 1, tuning = 0.25)
dev.off()

# Correlations joint networks with each other
cor(getWmat(gr1)[lower.tri(getWmat(gr1))], getWmat(gr2)[lower.tri(getWmat(gr2))],
  method = "spearman")
mean(getWmat(gr1))
mean(getWmat(gr2))

# Network Comparison -----------------------------------------------------------
q1 <- estimateNetwork(LOI, "EBICglasso", corMethod = "cor_auto")
q2 <- estimateNetwork(SAI, "EBICglasso", corMethod = "cor_auto")
plot(q1)
plot(q2)

# Time consumption! 10000 permutations
nct_12 <- NCT(q1, q2, it = 10, progressbar = TRUE, test.edges = TRUE,
  edges = 'all')
nct_12

res_nct <- matrix(c(1, 2, 0, 0, 0, 0, 0, 0), ncol = 8, byrow = TRUE)
colnames(res_nct) <- c("LOI subtype", "SAI subtype", "Difference in global strength",
  "p-value of global strength", "Max.diff invariance",  "p-value of invar.",
  "Global strength  Measure1", "Global strength  Measure2")
res_nct[1, 3:8] <- c(nct_12$glstrinv.real, nct_12$glstrinv.pval,
  nct_12$nwinv.real,  nct_12$nwinv.pval, nct_12$glstrinv.sep)
res_nct

plot(nct_12, what = "network")
# Plot results of global strength invariance test (not reliable with only 10 permutations!)
plot(nct_12, what = "strength")

# Community Struture -----------------------------------------------------------
# the walktrap-algorithm -------------------------------------------------------
graph.g1 <- as.igraph(gr1)
wc1 <- walktrap.community(graph.g1)
n.dim1 <- max(wc1$membership)

graph.g2 <- as.igraph(gr2)
wc2 <- walktrap.community(graph.g2)
n.dim2 <- max(wc2$membership)

# Figure S6
pdf("Figure S6.pdf", width = 14, height = 8)
par(mfrow = c(1, 2))
plot.ega1 <- qgraph(gr1, layout = "spring", vsize = 10,
  groups = as.factor(wc1$membership), overlay = TRUE, label.cex = 1,
  labels = colnames(LOI), maximum = Max)

plot.ega2 <- qgraph(gr2, layout = "spring", vsize = 10,
  groups = as.factor(wc2$membership), overlay = TRUE, label.cex = 1,
  labels = colnames(SAI), maximum = Max)
dev.off()

# check walktrap robustness. 10 times.
rand <- .Random.seed
walkch <- matrix(0, nrow = 17, ncol = 10)
for (i in 1:10) {
  set.seed(rand[i])
  walkch[, i] <- walktrap.community(graph.g1)$membership
}
walkch

for (i in 1:10) {
  set.seed(rand[i])
  walkch[, i] <- walktrap.community(graph.g2)$membership
}
walkch

# Sub-network structures -------------------------------------------------------
# Minimum spanning trees -------------------------------------------------------
md1 <- MaST(LOI, normal = TRUE, na.data = "none", depend = FALSE)
a <- qgraph(md1, layout = "spring", DoNotPlot = TRUE)
gmst1 <- as.igraph(a)
V(gmst1)$name <- hamd_vars
scale01 <- function(x) {(x - min(x))/(max(x) - min(x))}
vSizes <- (scale01(apply(LOI, 1, mean)) + 1.0) * 10
edgeweights <- gmst1$weight * 2.0

md2 <- MaST(SAI, normal = TRUE, na.data = "none", depend = FALSE)
b <- qgraph(md2, labels = colnames(SAI), DoNotPlot = TRUE)
gmst2 <- as.igraph(b)
vSizes <- (scale01(apply(SAI, 1, mean)) + 1.0) * 10
edgeweights <- gmst2$weight * 2.0

# Figure S7
pdf("Figure S7.pdf", width = 14, height = 14)
par(mfrow = c(2, 1))
dd <- qgraph(md1, layout = "spring", labels = TRUE, edge.labels = TRUE, label.cex = 2,
  edge.label.cex = 1.5, vsize = 7, esize = 9, label.color = "black",
  theme = "gimme" , borders = TRUE, title = "LOI", title.cex = 3)
qgraph(md2, layout = "spring", labels = TRUE, edge.labels = TRUE, label.cex = 2,
  edge.label.cex = 1.5, vsize = 7, esize = 9, label.color = "black",
  theme = "gimme" , borders = TRUE, title = "SAI", title.cex = 3)
dev.off()

# Estimate and plot centrality -------------------------------------------------
# 1 for LOI; 2 for SAI
# plot themes
themes <- theme(
  legend.title = element_blank(),
  legend.position = c(0.95, 0.95),
  legend.text = element_text(size = 15, face = "bold"),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.text.y = element_text(face = "bold", size = 12),
  axis.title.y = element_blank(),
  strip.text = element_text(size = 15),
  strip.background = element_blank(),
  strip.placement = "outside",
  axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = "black", size = 0.5)
)

# Standardized
centra1 <- centralityTable(gr1)
centra2 <- centralityTable(gr2)
centra2$graph <- "graph 2"

centra <- rbind(centra1, centra2)
centra$node <- factor(centra$node, levels = paste0("X", 1:17))
centra$graph <- factor(centra$graph, labels = c("LOI", "SAI"))

# Replace expected influence (step1) with expected influence (step2)
ef1 <- expectedInf(gr1, step = "both", directed = FALSE)$step2
ef2 <- expectedInf(gr2, step = "both", directed = FALSE)$step2
ef <- c(scale(ef1), scale(ef2))
centra[centra$measure == "ExpectedInfluence", ]$value <- ef

# Figure 3B
ggplot(centra, aes(node, value, group = graph, color = graph)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  facet_wrap(~measure, strip.position = "bottom", scales = "free") +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  themes
ggsave("Figure 3B.pdf", width = 14, height = 8)

# Unstandardized
us_centra1 <- centralityTable(gr1, standardized = FALSE)
us_centra2 <- centralityTable(gr2, standardized = FALSE)
us_centra2$graph <- "graph 2"

us_centra <- rbind(us_centra1, us_centra2)
us_centra$node <- factor(us_centra$node, levels = paste0("X", 1:17))
us_centra$graph <- factor(us_centra$graph, labels = c("LOI", "SAI"))

# Replace expected influence (step1) with expected influence (step2)
us_ef <- c(ef1, ef2)
us_centra[us_centra$measure == "ExpectedInfluence", ]$value <- us_ef

# plot Figure S8
ggplot(us_centra, aes(node, value, group = graph, color = graph)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  facet_wrap(~measure, strip.position = "bottom", scales = "free") +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  themes
ggsave("Figure S8.pdf", width = 14, height = 8)

# Correlation of the centrality of subtypes (Standardized) ---------------------
centra$node <- rep(hamd_vars, times = 8)
centra_list <- centra %>% group_by(graph, measure) %>% group_split()
# Add labels of each list
names(centra_list) <- paste0(
  rep(c("LOI", "SAI"), each = 4), "_",
  rep(c("Betweenness", "Closeness", "Strength", "ExpectedInfluence"), times = 2)
)

# Betweenness
cor(centra_list$LOI_Betweenness$value, centra_list$SAI_Betweenness$value)
# Clossness
cor(centra_list$LOI_Closeness$value, centra_list$SAI_Closeness$value)
# Strength
cor(centra_list$LOI_Strength$value, centra_list$SAI_Strength$value)
# ExpectedInfluence
cor(centra_list$LOI_ExpectedInfluence$value, centra_list$SAI_ExpectedInfluence$value)

# Describe nodes ---------------------------------------------------------------
decribe_node <- function(group) {
  lists <- centra_list[[group]][order(centra_list[[group]]$value, decreasing = TRUE), ]
  return(lists)
}
rank_nodes <- sapply(names(centra_list), decribe_node, simplify = FALSE)

# Clossness
rank_nodes$LOI_Closeness; rank_nodes$SAI_Closeness

# Strength
rank_nodes$LOI_Strength; rank_nodes$SAI_Strength

# Expected influence
rank_nodes$LOI_ExpectedInfluence; rank_nodes$SAI_ExpectedInfluence

# Betweenness
centra_list$LOI_Betweenness[
  order(centra_list$LOI_Betweenness$value - centra_list$SAI_Betweenness$value,
    decreasing = TRUE), 3]

# Stability estimates ----------------------------------------------------------
# making the test again in individual networks
# 1 for LOI; 2 for SAI
set.seed(1)
network1b <- estimateNetwork(LOI, default = "EBICglasso")
network2b <- estimateNetwork(SAI, default = "EBICglasso")
plot(network1b, layout = "spring", labels = TRUE)
plot(network2b, layout = "spring", labels = TRUE)

kboot1a <- bootnet(network1b, nBoots = 1000, nCores = 16)
kboot1b <- bootnet(network1b, nBoots = 1000, type = "case", nCores = 16)

kboot2a <- bootnet(network2b, nBoots = 1000, nCores = 16)
kboot2b <- bootnet(network2b, nBoots = 1000, type = "case", nCores = 16)

# Plot edge weight CI; Figure S9
p1 <- plot(kboot1a, labels = FALSE, order = "sample") +
  labs(x = "LOI") +
  theme(axis.title.x = element_text(face = "bold", size = 15))
p2 <- plot(kboot2a, labels = FALSE, order = "sample", legend = FALSE) +
  labs(x = "SAI") +
  theme(axis.title.x = element_text(face = "bold", size = 15))

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S9.pdf", width = 10, height = 12)

# Plot centrality stability (Strength); Figure S10
# Plot themes
themes <- theme(
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(face = "bold", size = 15),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(face = "bold", size = 15),
  axis.text.y = element_text(face = "bold", size = 12),
  axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = "black", size = 0.5)
)

p1 <- plot(kboot1b, legend = FALSE) +
  labs(x = "Sampled People (LOI)") +
  geom_line(size = 1.5, color = "#E41A1C") +
  geom_point(size = 2, color = "#E41A1C") +
  themes

p2 <- plot(kboot2b, legend = FALSE) +
  labs(x = "Sampled People (SAI)") +
  geom_line(size = 1.5, color = "#377EB8") +
  geom_point(size = 2, color = "#377EB8") +
  themes

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S10.pdf", width = 12, height = 10)

# Centrality stability coefficient
cs1 <- corStability(kboot1b)
cs2 <- corStability(kboot2b)
cs <- matrix(nrow = 2, ncol = 2)
cs[1, 1:2] <- round(cs1, digits = 3)
cs[2, 1:2] <- round(cs2, digits = 3)
colnames(cs) <- c("Edge", "Stregnth")
rownames(cs) <- c("LOI", "SAI")
cs

# Edge weights diff test; Figure S11
p1 <- plot(kboot1a, "edge", plot = "difference", onlyNonZero = TRUE,
  order = "sample", panels = FALSE) +
  labs(x = "LOI") +
  theme(axis.title.x = element_text(face = "bold", size = 15))

p2 <- plot(kboot2a, "edge", plot = "difference", onlyNonZero = TRUE,
  order = "sample", panels = FALSE) +
  labs(x = "SAI") +
  theme(axis.title.x = element_text(face = "bold", size = 15))

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S11.pdf", width = 8, height = 12)

# Centrality diff test (Strength); Figure S12
p1 <- plot(kboot1a, "strength", order = "sample", labels = TRUE, panels = FALSE) +
  labs(x = "LOI") +
  theme(axis.title.x = element_text(face = "bold", size = 15))
p2 <- plot(kboot2a, "strength", order = "sample", labels = TRUE, panels = FALSE) +
  labs(x = "SAI") +
  theme(axis.title.x = element_text(face = "bold", size = 15))

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S12.pdf", width = 8, height = 12)

################################################################################
############################# SPM file and ROI #################################
################################################################################
# Check subgroups --------------------------------------------------------------
# 1 for LOI;2 for SAI
with(patient, ftable(Level1))
# 1vs4  2vs3. 1 for LOI+; 4 for LOI-; 3 for SAI+; 4 for SAI-
with(patient, ftable(Level1, Level2))   

# Function factory -------------------------------------------------------------
extract <- function(vars) {
  extract_vars <- function(level, cluster) {
    dirs <- patient[patient[[level]] %in% cluster, ][vars]
    dirs
  }
}

# Extract SPM directory for every patient --------------------------------------
extract_dir <- extract("Dir")
Dirs <- list("1-1" = extract_dir("Level1", 1), "1-2" = extract_dir("Level1", 2),
  "2-1" = extract_dir("Level2", 1), "2-4" = extract_dir("Level2", 4),
  "2-2" = extract_dir("Level2", 2), "2-3" = extract_dir("Level2", 3))

# Extract covariates -----------------------------------------------------------
# covariates
covariates <- c("Age", "Gender", "Education", "Sites", "HeadMotion")
extract_covars <- extract(covariates)

Covariates <- vector("list", length = 3)
names(Covariates) <- c("1_12", "2_14", "2_23")
# Level 1
Covariates[[1]] <- rbind(extract_covars("Level1", 1), extract_covars("Level1", 2))                            # AONVA
# Level 2
Covariates[[2]] <- rbind(extract_covars("Level2", 1), extract_covars("Level2", 4))
Covariates[[3]] <- rbind(extract_covars("Level2", 2), extract_covars("Level2", 3))

# Fix the Site variable --------------------------------------------------------
with(patient, ftable(Level1, Sites))
with(patient, ftable(Level1, Level2, Sites))          # 1vs4 needs to be fixed

extract_sites <- extract("Sites")
Sites <- rbind(extract_sites("Level2", 1), extract_sites("Level2", 4))
fix_sites <- function(data) {
  data <- data[[1]]
  raw <- unique(data)
  Replace <- seq(raw)
  for (i in seq(raw)) {
    data <- replace(data, which(data == raw[i]), Replace[i])
  }
  # Abandon useless level
  data <- droplevels(data)
  data
}
Sites <- fix_sites(Sites)
Covariates$`2_14`$Sites <- Sites

# Output -----------------------------------------------------------------------
# Covariates for t-test
for (i in seq(Covariates)) {
  write.table(Covariates[[i]], paste0("Covariates", names(Covariates)[[i]], ".txt"), 
    quote = FALSE, row.names = FALSE, col.names = FALSE)
}                                                           
# SPM file directory
for (i in seq(Dirs)) {
  write.table(Dirs[[i]], paste0("Dirs", names(Dirs)[[i]], ".txt"), 
    quote = FALSE, row.names = FALSE, col.names = FALSE)
}

######################## Regression analysis ###################################
# ROIs
home2 <- "E:/WorkingSpace/R/Paper/Depression/Subtypes/All/"

ROIs <- matrix(nrow = nrow(patient), ncol = 17)
colnames(ROIs) <- paste0("ROI", seq(ncol(ROIs)))

# Add ROIs mean signals
for (i in seq(colnames(ROIs))) {
  data <- read.table(paste0(home2, "ROI", i, ".txt"))
  data <- data[, , drop = TRUE]
  ROIs[, i] <- data
}
# No NAs
sum(is.na(ROIs))

# Covariates
covariates <- c("Age", "Gender", "Education", "Sites", "HeadMotion")
# PCA components
pca_vars <- names(patient)[31:34]
# ROI
ROIs_vars <- cbind(ROIs, patient[, c(pca_vars, covariates, "Level1")])
# Convert factor to numeric
ROIs_vars[] <- lapply(ROIs_vars, as.numeric)
# Convert somatic and anxiety to SA
names(ROIs_vars)[20] <- "SA"
names(ROIs_vars)[21] <- "IL"

# Regression analysis --------------------------------------------------------
reg <- function(data) {
  reg_results <- data.frame()
  for (i in 1:17) {
    reg <- lm(data[, i] ~ Guilt + Insomnia + SA + IL + Age + Gender + Education + Sites + HeadMotion, 
      data = data)
    reg_result <- summary(reg)$coefficients[2:5, 3:4]
    reg_results  <- rbind(reg_results, reg_result)
  }
  
  reg_tb <- data.frame(rep(names(ROIs_vars)[1:17], each = 4), 
    rep(names(patient)[31:34], times = 17))
  reg_df <- cbind(reg_tb, reg_results)
  names(reg_df) <- c("ROIs", "Components", "Statistic", "P")
  rownames(reg_df) <- NULL
  
  # FDR correction for p-values
  reg_df$P_corrected <- p.adjust(reg_df$P, method = "fdr")
  # modify output
  reg_df$P <- round(reg_df$P, 3)
  reg_df$P_corrected <- round(reg_df$P_corrected, 3)
  reg_df$Statistic <- round(reg_df$Statistic, 2)
  
  return(reg_df)
}
result_all <- reg(ROIs_vars)

# Output -----------------------------------------------------------------------
write.xlsx2(result_all, "ROI_Regression.xlsx", row.names = FALSE)
