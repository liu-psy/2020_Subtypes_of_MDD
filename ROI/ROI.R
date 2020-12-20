# This script was used to load ROIs signals and regression analysis
library(xlsx)
library(dplyr)


setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/ROI")
# Load data --------------------------------------------------------------------
home1 <- "E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/HHC/HHC.xlsx"
# ROIs
home2 <- "E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/ROI/"

patient <- read.xlsx2(
  home1,
  sheetIndex       = 1,
  stringsAsFactors = FALSE,
  check.names      = FALSE,
  colClasses       =rep(c("character", "numeric"), times = c(8, 29))
)
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
    reg_result   <- summary(reg)$coefficients[2:5, 3:4]
    reg_results  <- rbind(reg_results, reg_result)
  }

  reg_tb <- data.frame(
    rep(names(ROIs_vars)[1:17], each = 4),
    rep(names(patient)[31:34], times = 17)
  )
  reg_df           <- cbind(reg_tb, reg_results)
  names(reg_df)    <- c("ROIs", "Components", "Statistic", "P")
  rownames(reg_df) <- NULL

  # FDR correction for p-values
  reg_df$P_corrected <- p.adjust(reg_df$P, method = "fdr")
  # modify output
  reg_df$P           <- round(reg_df$P, 3)
  reg_df$P_corrected <- round(reg_df$P_corrected, 3)
  reg_df$Statistic   <- round(reg_df$Statistic, 2)

  return(reg_df)
}
result_all <- reg(ROIs_vars)

# Output -----------------------------------------------------------------------
write.xlsx2(result_all, "ROI_Regression.xlsx", row.names = FALSE)