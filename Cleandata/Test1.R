# This script was used to gender analysis
library(xlsx)
library(dplyr)
library(psych)
library(effsize)     # Cohen's d
library(sjstats)     # Crame's V


# Set working directory
setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/Cleandata")

# Loading Data -----------------------------------------------------------------
home <- "Data.xlsx"
patient <- read.xlsx2(home, 1,  stringsAsFactors = FALSE,
  colClasses = c(rep("character", 6), rep("numeric", 29)))

patient$Gender <- factor(patient$Gender, levels = c(1, 2), labels = c("Male", "Female"))
patient$Sites <- factor(patient$Sites, levels = 1:15)
hamd_vars <- c("Depressed Mood", "Guilt", "Suicide", "Early Insomnia", "Middle Insomnia",
  "Late Insomnia", "Work Interests", "Retardation", "Agitation", "Psychic Anxiety",
  "Somatic Anxiety", "Gastrointestinal", "General Somatic", "Loss of Libido",
  "Hypochondriasis", "Weight Loss", "Loss of Insight")
names(patient)[11:27] <- paste(1:17, hamd_vars)

rm(home)

# Describe Data ----------------------------------------------------------------
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