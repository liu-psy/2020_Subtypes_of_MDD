# This script was used to data clean
library(corrplot)
library(dplyr)
library(ggplot2)
library(mice)                 # Imputation
library(psych)                # PCA
library(reshape2)
library(xlsx)

# Set working directory
setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/Cleandata")

# Loading Data -----------------------------------------------------------------
# Data of REST-meta-MDD
home1 <- "Depression.xlsx"
# Data of SWU
home2 <- "SWU_Data.xlsx"

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
patient <- rbind(
  patient_rest_project[1:grep("S21", patient_rest_project$ID)[1] - 1, ], 
  patient_swu,
  patient_rest_project[grep("S21", patient_rest_project$ID)[1] : nrow(patient_rest_project), ]
)

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
names(pa_table) <- c("Acutual Data", "Resampled Data")
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

vars <- c("Depressed Mood", "Guilt", "Suicide", "Early Insomnia",
  "Middle Insomnia",  "Late Insomnia", "Work Interests", "Retardation",
  "Agitation", "Psychic Anxiety", "Somatic Anxiety", "Gastrointestinal",
  "General Somatic", "Loss of Libido", "Hypochondriasis", "Weight Loss",
  "Loss of Insight")
components <- c("Guilt", "Insomnia", "Somatic and Anxiety", "Interests Loss")

# Modify components position
loading_matrix <- PCA$loadings[, c(4, 1, 2, 3)]
colnames(loading_matrix) <- components
write.csv(loading_matrix, "Loading.csv", row.names = FALSE)

# Plot loading matrix, Figure 1
# This code is inspired by https://rpubs.com/danmirman/plotting_factor_analysis
loading <- as.data.frame(loading_matrix)
loading$vars <- vars %>% factor(levels = rev(vars))
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

# Output Data ------------------------------------------------------------------
patient <- patient %>% select(ID, Gender, Sites, Episode, Medication, Dir, Age,
  Education, Month, HeadMotion, everything())

write.xlsx2(patient, "DATA.xlsx", sheetName = "patient", row.names = FALSE)