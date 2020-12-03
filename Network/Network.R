# This script was used to network analysis
library(bootnet)
library(dplyr)
library(ggplot2)
library(mgm)
library(NetworkComparisonTest)
library(patchwork)
library(qgraph)
library(RColorBrewer)
library(reshape2)
library(xlsx)


setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/Network")
# Load Data --------------------------------------------------------------------
home1    <- "E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/HHC/HHC.xlsx"
patient  <- read.xlsx2(home1, 1, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = rep(c("character", "numeric"), times = c(8, 29)))

# Load Loading matrix
home2    <- "E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/Cleandata/loading.csv"
loading  <- read.csv(home2)

# Get communities from loading matrix
community <- apply(loading, 1, which.max)
group     <- list(c(1:3, 16), c(4:6, 12), c(9:11, 13:15), c(7:8, 17))
# Color of communities
colors    <- brewer.pal(length(group), "Pastel1")

# Extract HAMD-17 of subtypes
LOI <- subset(patient, Level1 == 1, X1:X17)
SAI <- subset(patient, Level1 == 2, X1:X17)

# Items of HAMD-17
vars <- c("Depressed Mood", "Guilt", "Suicide", "Early Insomnia",
  "Middle Insomnia", "Late Insomnia", "Work Interests", "Retardation",
  "Agitation", "Psychic Anxiety", "Somatic Anxiety", "Gastrointestinal",
  "General Somatic", "Loss of Libido", "Hypochondriasis", "Weight Loss",
  "Loss of Insight")
labs     <- paste0("X", seq(LOI))
symptoms <- paste(seq(LOI), vars)

rm(home1, home2, loading)

## Network Analysis ------------------------------------------------------------
# Network estimation -----------------------------------------------------------
# glasso Models with EBIC ------------------------------------------------------
# 1 for LOI; 2 for SAI
gr1 <- LOI %>% cor_auto() %>% EBICglasso(n = nrow(LOI))
gr2 <- SAI %>% cor_auto() %>% EBICglasso(n = nrow(SAI))

# Layout of networks
Max <- max(c(gr1, gr2))
L   <- averageLayout(gr1, gr2)

# Predicability Networks -------------------------------------------------------
set.seed(1)
f1 <- function(subtype) {
  # Mixed Graphical Models
  fit    <- mgm(
    subtype,
    type = rep("g", ncol(subtype)),
    level = rep(1, ncol(subtype)),
    lambdaSel = "CV",
    ruleReg = "OR"
  )
  pred   <- predict(fit, data = subtype, errorCon = "R2")
  result <- pred$error$R2
  names(result) <- symptoms
  
  # Output in descending order according to centrality
  cat("\n\n\n", substitute(subtype), ":", "\n")
  print(result[order(result, decreasing = TRUE)])
  # Average node predictability
  cat("\n", substitute(subtype), "average node predictability:", round(mean(result), 2), "\n")
  
  return(result)
}

# 1 for LOI; 2 for SAI
predicability1 <- f1(LOI)
predicability2 <- f1(SAI)

# Plot networks; Figure 3A
f2 <- function(net, title, pred) {
  network <- qgraph(net, title = title, pie = pred, layout = L,
    title.cex = 2.5, maximum = Max, theme = "Hollywood",
    pieColor = "#FC8D62", border.width = 2, vsize = 9, label.cex = 1,
    groups = group, color = colors, labels = labs
  )
  return(network)
}

pdf("Figure 3A.pdf", width = 14, height = 8)
par(mfrow = c(1, 2))
gr1 <- f2(gr1, "LOI", predicability1)
gr2 <- f2(gr2, "SAI", predicability2)
dev.off()

# Estimate and plot node centrality --------------------------------------------
# plot themes
themes <- theme(
  legend.title       = element_blank(),
  legend.position    = c(0.95, 0.95),
  legend.text        = element_text(size = 15, face = "bold"),
  panel.grid.minor   = element_blank(),
  axis.text.x        = element_text(size = 12),
  axis.title.x       = element_blank(),
  axis.text.y        = element_text(face = "bold", size = 12),
  axis.title.y       = element_blank(),
  strip.text         = element_text(size = 15),
  strip.background   = element_blank(),
  strip.placement    = "outside",
  axis.line          = element_line(color = "black", size = 1, linetype = "solid"),
  panel.background   = element_rect(fill = "white"),
  panel.grid.major.y = element_line(color = "black", size = 0.5)
)

f3 <- function(Standardized) {
  # Get node centrality (Betweenness, Closeness, ExpectedInfluence, Strength)
  c1        <- centralityTable(gr1, standardized = Standardized)
  c2        <- centralityTable(gr2, standardized = Standardized)
  c2$graph  <- "graph 2"
  # Combine node centrality of two networks
  cc        <- rbind(c1, c2)
  cc$node   <- factor(cc$node, levels = labs)
  cc$graph  <- factor(cc$graph, labels = c("LOI", "SAI"))
  # plot node centrality
  p <- ggplot(cc, aes(node, value, group = graph, color = graph)) +
    geom_line(size = 1.5) +
    geom_point(size = 3) +
    facet_wrap(~measure, strip.position = "bottom", scales = "free") +
    scale_color_manual(values = c("#E41A1C", "#377EB8")) +
    themes
  print(p)
  return(cc)
}

# Standardized node centrality; Figure 3B
centra <- f3(Standardized = TRUE)
ggsave("Figure 3B.pdf", width = 14, height = 8)

# Unstandardized node centrality; Figure S6
un_centra <- f3(Standardized = FALSE)
ggsave("Figure S6.pdf", width = 14, height = 8)

# Correlation of the centrality between subtypes -------------------------------
centralities <- c("Betweenness", "Closeness", "ExpectedInfluence", "Strength")
centra$node  <- rep(vars, times = 8)

centra_list  <- centra %>%
  select(-type) %>%
  group_by(graph, measure) %>%
  select(graph, measure, node, value) %>%
  mutate(rank = rank(value, na.last = FALSE)) %>%
  group_split()
# Add labels of each list
names(centra_list) <- paste0(
  rep(c("LOI", "SAI"), each = 4),
  "_",
  rep(centralities, times = 2)
)
centra_list

f4 <- function(centrality) {
  var1 <- paste0("LOI_", centrality)
  x1   <- centra_list[[var1]]$value
  var2 <- paste0("SAI_", centrality)
  x2   <- centra_list[[var2]]$value
  
  # Compute perason correlation between two networks (NA is ignored)
  cors <- cor(x1, x2, use = "pairwise.complete.obs") %>%
    round(2)
  cat("Correlation coefficient -", centrality, ":", cors, "\n")
}
f4("Betweenness")
f4("Closeness")
f4("ExpectedInfluence")
f4("Strength")

# Mean rank of centralities
f5 <- function(subtype) {
  subtype_centralities <- paste0(subtype, "_", centralities)
  result <- centra_list[[subtype_centralities[1]]]$rank +
    centra_list[[subtype_centralities[2]]]$rank +
    centra_list[[subtype_centralities[3]]]$rank +
    centra_list[[subtype_centralities[4]]]$rank
  
  result         <- round(result / length(centralities), 1)
  names(result)  <- vars
  # Sorting results
  ordered_result <- result[order(result, decreasing = TRUE)]
  result_list    <- list("Origin" = result, "Ordered" = ordered_result)
  
  return(result_list)
}
f5("LOI")
f5("SAI")

# Network Comparison -----------------------------------------------------------
# 1 for LOI, 2 for SAI
node_centrality <- c("closeness", "betweenness", "strength", "expectedInfluence")

q1 <- estimateNetwork(LOI, "EBICglasso", corMethod = "cor_auto")
q2 <- estimateNetwork(SAI, "EBICglasso", corMethod = "cor_auto")

# NCT with 10000 iterations
nct <- NCT(q1, q2, it = 1000, progressbar = TRUE, test.edges = TRUE, 
  test.centrality = TRUE, centrality = node_centrality)
summary(nct)

# p-value of edge differences
nct$einv.pvals[nct$einv.pvals$`p-value` <=0.05, ]
# p-values of local nodes
nct$diffcen.pval
# Differences of local nodes; LOI - SAI
round(nct$diffcen.real, 3)

# Plot results of the network structure invariance test
plot(nct, what = "network")
# Plot results of global strength invariance test
plot(nct, what = "strength")

# Stability estimates ----------------------------------------------------------
# 1 for LOI; 2 for SAI
kboot1a <- bootnet(q1, nBoots = 1000, nCores = 16)
kboot1b <- bootnet(q1, nBoots = 1000, nCores = 16, type = "case",
  statistics = node_centrality)

kboot2a <- bootnet(q2, nBoots = 1000, nCores = 16)
kboot2b <- bootnet(q2, nBoots = 1000, nCores = 16, type = "case",
  statistics = node_centrality)

# Plot edge weight CI; Figure S7
f6 <- function(data, xlab) {
  plot(data, labels = FALSE, order = "sample") +
    labs(x = xlab) +
    theme(axis.title.x = element_text(face = "bold", size = 15))
}
p1 <- f6(kboot1a, "LOI")
p2 <- f6(kboot2a, "SAI")

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S7.pdf", width = 10, height = 12)

# Plot themes
themes <- theme(
  panel.grid.minor   = element_blank(),
  axis.text.x        = element_text(size = 12),
  axis.title.x       = element_text(face = "bold", size = 15),
  axis.title.y       = element_text(face = "bold", size = 15),
  axis.text.y        = element_text(face = "bold", size = 12),
  axis.line          = element_line(color = "black", size = 1, linetype = "solid"),
  panel.background   = element_rect(fill = "white"),
  panel.grid.major.y = element_line(color = "black", size = 0.5)
)

# Plot centrality stability; Figure S8
f7 <- function(data, xlab) {
  plot(data, statistics = "all") +
    labs(x = xlab) +
    geom_line(size = 1.5) +
    geom_point(size = 2) +
    themes
}
p1 <- f7(kboot1b, "Sampled People (LOI)")
p2 <- f7(kboot2b, "Sampled People (SAI)")

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S8.pdf", width = 12, height = 10)

# Centrality stability coefficient
cs_LOI <- corStability(kboot1b)
cs_SAI <- corStability(kboot2b)

# Edge weights diff test; Figure S9
f8 <- function(data, xlab) {
  plot(data, "edge", plot = "difference", onlyNonZero = TRUE,
    order = "sample", panels = FALSE) +
    labs(x = xlab) +
    theme(axis.title.x = element_text(face = "bold", size = 15))
}
p1 <- f8(kboot1a, "LOI")
p2 <- f8(kboot2a, "SAI")

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S9.pdf", width = 8, height = 12)

# Centrality diff test (Strength); Figure S10
f9 <- function(data, xlab) {
  plot(data, "strength", order = "sample", labels = TRUE, panels = FALSE) +
    labs(x = xlab) +
    theme(axis.title.x = element_text(face = "bold", size = 15))
}
p1 <- f9(kboot1a, "LOI")
p2 <- f9(kboot2a, "SAI")

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S10.pdf", width = 8, height = 12)