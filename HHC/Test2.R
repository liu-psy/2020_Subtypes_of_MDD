# This script is used to compare the behavior of subtypes
library(xlsx)
library(dplyr)
library(sjstats)     # Cramer's V
library(effsize)     # Cohen's d


# Set working directory
setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/HHC")

# Load data --------------------------------------------------------------------
home <- "HHC.xlsx"
patient <- read.xlsx2(home, 1,  stringsAsFactors = FALSE,
  colClasses = rep(c("character", "numeric"), times = c(8, 29)))

# Set character to factor
patient[, c(2, 3, 7:8)][] <- lapply(patient[, c(2, 3, 7:8)], factor)

# Modify variable names --------------------------------------------------------
# HAMD-17
vars <- c("Depressed Mood", "Guilt", "Suicide", "Early Insomnia",  "Middle Insomnia",
  "Late Insomnia", "Work Interests", "Retardation", "Agitation", "Psychic Anxiety",
  "Somatic Anxiety", "Gastrointestinal", "General Somatic", "Loss of Libido",
  "Hypochondriasis", "Weight Loss", "Loss of Insight")
names(patient)[13:29] <- paste(1:17, vars)
# PCA
pca_vars <- names(patient)[31:34]
# Sub-dimensions
sub_vars <- names(patient)[35:37]
# symptoms
sym_vars <- names(patient)[13:29]

rm(home, vars)

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