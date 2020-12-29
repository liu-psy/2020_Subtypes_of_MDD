# This script was used to extract file path and covariates for SPM analysis
library(xlsx)
library(dplyr)


setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/ROI")
# Load data --------------------------------------------------------------------
home <- "E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/HHC/HHC.xlsx"
patient <- read.xlsx2(
  home,
  sheetIndex       = 1,
  stringsAsFactors = FALSE,
  check.names      = FALSE,
  colClasses       = rep(c("character", "numeric"), times = c(8, 29)),
)
patient$Sites <- factor(patient$Sites, levels = 1:15)

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
Dirs <- list(
  "1-1" = extract_dir("Level1", 1),
  "1-2" = extract_dir("Level1", 2),
  "2-1" = extract_dir("Level2", 1),
  "2-4" = extract_dir("Level2", 4),
  "2-2" = extract_dir("Level2", 2),
  "2-3" = extract_dir("Level2", 3)
)

# Extract covariates -----------------------------------------------------------
# covariates
covariates     <- c("Age", "Gender", "Education", "Sites", "HeadMotion")
extract_covars <- extract(covariates)

Covariates        <- vector("list", length = 3)
names(Covariates) <- c("1_12", "2_14", "2_23")
# Level 1
Covariates[[1]] <- rbind(
  extract_covars("Level1", 1), 
  extract_covars("Level1", 2)
)
# Level 2
Covariates[[2]] <- rbind(
  extract_covars("Level2", 1), 
  extract_covars("Level2", 4)
)
Covariates[[3]] <- rbind(
  extract_covars("Level2", 2), 
  extract_covars("Level2", 3)
)

# Fix the Site variable --------------------------------------------------------
with(patient, ftable(Level1, Sites))
with(patient, ftable(Level1, Level2, Sites))          # 1vs4 needs to be fixed

extract_sites <- extract("Sites")
Sites <- rbind(extract_sites("Level2", 1), extract_sites("Level2", 4))
fix_sites <- function(data) {
  data    <- data[[1]]
  raw     <- unique(data)
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
  write.table(
    Covariates[[i]],
    paste0("Covariates", names(Covariates)[[i]], ".txt"),
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}
# SPM file directory
for (i in seq(Dirs)) {
  write.table(
    Dirs[[i]],
    paste0("Dirs", names(Dirs)[[i]], ".txt"),
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}