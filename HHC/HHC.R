# This script was used to Hybrid Hierarchical Clustering and plot
library(xlsx)
library(dplyr)
library(hybridHclust)          # Hybrid Hierarchical Clustering
library(dynamicTreeCut)        # Dynamic Tree Cut algorithm
library(gplots)                # Heatmap
library(dendextend)            # Dendrogram
library(clusteval)             # Robustness of clusering


# Set working directory and load data ------------------------------------------
setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/HHC")

home <- "E:/WorkingSpace/R/Paper/Depression/Subtypes/Cleandata/Data.xlsx"
patient <- read.xlsx2(home, 1, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = rep(c("character", "numeric"), times = c(6, 29)))

# Hybrid Hierarchnical Clustering ----------------------------------------------
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

# Output -----------------------------------------------------------------------
patient <- data.frame(patient, Level1 = labels[, 1], Level2 = labels[, 2], check.names = FALSE) %>%
  select(names(patient)[1:6], Level1, Level2, everything())
write.xlsx2(patient, "HHC.xlsx", sheetName = "patient", row.names = FALSE)