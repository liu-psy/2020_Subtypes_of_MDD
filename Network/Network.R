# This script was used to network analysis
# And this script was mainly derived from https://github.com/Melovainio/Network_depression
library(xlsx)
library(dplyr)
library(reshape2)
library(qgraph)
library(igraph)
library(EstimateGroupNetwork)
library(mgm)
library(NetworkComparisonTest)
library(bootnet)
library(NetworkToolbox)
library(networktools)
library(ggplot2)
library(patchwork)


setwd("E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/Network")
# Load Data --------------------------------------------------------------------
home    <- "E:/WorkingSpace/Project/2020_Symptom_Subtyping_MDD/HHC/HHC.xlsx"
patient <- read.xlsx2(home, 1, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = rep(c("character", "numeric"), times = c(8, 29)))

# Items of HAMD-17
vars <- c("Depressed Mood", "Guilt", "Suicide", "Early Insomnia", "Middle Insomnia",
  "Late Insomnia", "Work Interests", "Retardation", "Agitation", "Psychic Anxiety",
  "Somatic Anxiety", "Gastrointestinal", "General Somatic", "Loss of Libido",
  "Hypochondriasis", "Weight Loss", "Loss of Insight")

# Extract symptom of subtypes
LOI <- subset(patient, Level1 == 1, X1:X17)
SAI <- subset(patient, Level1 == 2, X1:X17)

### Network Analysis -----------------------------------------------------------
# Estimated polychoric correlations among symptoms -----------------------------
LOI_cor <- cor_auto(LOI)
SAI_cor <- cor_auto(SAI)

# Estimate Joint Graphical Lasso model
net <- EstimateGroupNetwork(
  list(LOI, SAI), 
  inputType = "list.of.dataframes",
  method = "crossvalidation", 
  criterion = "ebic", 
  gamma = 0.5,
  simplifyOutput = FALSE,
  seed = 100, 
  ncores = 16
)
LOI_net <- net$network[[1]]
SAI_net <- net$network[[2]]

# 1 for LOI; 2 for SAI
nw1 <- getWmat(qgraph(LOI_net, sampleSize = nrow(LOI), DoNotPlot = TRUE))
nw2 <- getWmat(qgraph(SAI_net, sampleSize = nrow(SAI), DoNotPlot = TRUE))

Max <- max(c(nw1, nw2))
Max
# layout of nodes
L <- averageLayout(nw1, nw2)

# Preidictability Networks -----------------------------------------------------
# Mixed Graphical Models -------------------------------------------------------
# 1 for LOI; 2 for SAI
set.seed(1)
fit1 <- mgm(LOI, type = rep('g', ncol(LOI)), level = rep(1, ncol(LOI)), 
  lambdaSel = 'CV', ruleReg = 'OR')
fit2 <- mgm(SAI, type = rep('g', ncol(SAI)), level = rep(1, ncol(SAI)), 
  lambdaSel = 'CV', ruleReg = 'OR')

pred1 <- predict(fit1, data = LOI, errorCon = 'R2')
pred2 <- predict(fit2, data = SAI, errorCon = 'R2')

pred1$error$R2
vars[order(pred1$error$R2, decreasing = TRUE)]
pred2$error$R2
vars[order(pred2$error$R2, decreasing = TRUE)]

# Average node predictability
round(mean(pred1$error$R2), 2)
round(mean(pred2$error$R2), 2)

# Plot networks; Figure 3A
f1 <- function(data, title) {
  network <- qgraph(data, layout = L, title = title, title.cex = 2.5,
    maximum = Max,  theme = "Hollywood", pieColor = "#FC8D62", pie = pred1$error$R2,
    border.width = 2, vsize = 9, label.cex = 1, tuning = 0.25)
  return(network)
}

pdf("Figure 3a.pdf", width = 14, height = 8)
par(mfrow = c(1, 2))
gr1 <- f1(LOI_net, "LOI")
gr2 <- f1(SAI_net, "SAI")
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
nct_12 <- NCT(q1, q2, it = 10000, progressbar = TRUE, test.edges = TRUE,
  edges = 'all')
nct_12

res_nct <- matrix(c(1, 2, 0, 0, 0, 0, 0, 0), ncol = 8, byrow = TRUE)
colnames(res_nct) <- c("LOI subtype", "SAI subtype", "Difference in global strength",
  "p-value of global strength", "Max.diff invariance",  "p-value of invar.",
  "Global strength  Measure1", "Global strength  Measure2")
res_nct[1, 3:8] <- c(nct_12$glstrinv.real, nct_12$glstrinv.pval,
  nct_12$nwinv.real,  nct_12$nwinv.pval, nct_12$glstrinv.sep)
res_nct

# Plot results of the network structure invariance test
plot(nct_12, what = "network")
# Plot results of global strength invariance test
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
# Figure S7
f2 <- function(data, title) {
  # Maximum Spanning Tree
  md <- MaST(data, normal = TRUE, na.data = "none", depend = FALSE)
  result <- qgraph(md, title = title, layout = "spring", labels = TRUE, 
    edge.labels = TRUE, label.cex = 2, edge.label.cex = 1.5, vsize = 7, 
    esize = 9, label.color = "black", theme = "gimme" , borders = TRUE, 
    title.cex = 3)
  return(result)
}

pdf("Figure S7.pdf", width = 14, height = 14)
par(mfrow = c(2, 1))
f2(LOI, "LOI")
f2(SAI, "SAI")
dev.off()

# Estimate and plot centrality -------------------------------------------------
# 1 for LOI; 2 for SAI
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

# Standardized
centra1       <- centralityTable(gr1)
centra2       <- centralityTable(gr2)
centra2$graph <- "graph 2"

centra        <- rbind(centra1, centra2)
centra$node   <- factor(centra$node, levels = paste0("X", 1:17))
centra$graph  <- factor(centra$graph, labels = c("LOI", "SAI"))

# Replace expected influence (step1) with expected influence (step2)
ef1 <- expectedInf(gr1, step = "both", directed = FALSE)$step2
ef2 <- expectedInf(gr2, step = "both", directed = FALSE)$step2
ef  <- c(scale(ef1), scale(ef2))
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
us_centra1       <- centralityTable(gr1, standardized = FALSE)
us_centra2       <- centralityTable(gr2, standardized = FALSE)
us_centra2$graph <- "graph 2"

us_centra        <- rbind(us_centra1, us_centra2)
us_centra$node   <- factor(us_centra$node, levels = paste0("X", 1:17))
us_centra$graph  <- factor(us_centra$graph, labels = c("LOI", "SAI"))

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
centra$node <- rep(vars, times = 8)
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
kboot1b <- bootnet(network1b, nBoots = 1000, nCores = 16, type = "case")

kboot2a <- bootnet(network2b, nBoots = 1000, nCores = 16)
kboot2b <- bootnet(network2b, nBoots = 1000, nCores = 16, type = "case")

# Plot edge weight CI; Figure S9
f3 <- function(data, xlab) {
  plot(data, labels = FALSE, order = "sample") +
    labs(x = xlab) +
    theme(axis.title.x = element_text(face = "bold", size = 15))
}
p1 <- f3(kboot1a, "LOI")
p2 <- f3(kboot2a, "SAI")

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S9.pdf", width = 10, height = 12)

# Plot centrality stability (Strength); Figure S10
# Plot themes
themes <- theme(
  panel.grid.minor   = element_blank(),
  axis.title.x       = element_text(face = "bold", size = 15),
  axis.text.x        = element_text(size = 12),
  axis.title.y       = element_text(face = "bold", size = 15),
  axis.text.y        = element_text(face = "bold", size = 12),
  axis.line          = element_line(color = "black", size = 1, linetype = "solid"),
  panel.background   = element_rect(fill = "white"),
  panel.grid.major.y = element_line(color = "black", size = 0.5)
)

# Figure S10
f4 <- function(data, xlab) {
  plot(data, legend = FALSE) +
    labs(x = xlab) +
    geom_line(size = 1.5, color = "#E41A1C") +
    geom_point(size = 2, color = "#E41A1C") +
    themes
}
p1 <- f4(kboot1b, "Sampled People (LOI)")
p2 <- f4(kboot2b, "Sampled People (SAI)")

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S10.pdf", width = 12, height = 10)

# Centrality stability coefficient
cs1 <- corStability(kboot1b)
cs2 <- corStability(kboot2b)
cs  <- matrix(nrow = 2, ncol = 2)
cs[1, 1:2] <- round(cs1, digits = 3)
cs[2, 1:2] <- round(cs2, digits = 3)
colnames(cs) <- c("Edge", "Stregnth")
rownames(cs) <- c("LOI", "SAI")
cs

# Edge weights diff test; Figure S11
f5 <- function(data, xlab) {
  plot(data, "edge", plot = "difference", onlyNonZero = TRUE,
    order = "sample", panels = FALSE) +
    labs(x = xlab) +
    theme(axis.title.x = element_text(face = "bold", size = 15))
}
p1 <- f5(kboot1a, "LOI")
p2 <- f5(kboot2a, "SAI")

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S11.pdf", width = 8, height = 12)

# Centrality diff test (Strength); Figure S12
f6 <- function(data, xlab) {
  plot(data, "strength", order = "sample", labels = TRUE, panels = FALSE) +
    labs(x = xlab) +
    theme(axis.title.x = element_text(face = "bold", size = 15))
}
p1 <- f6(kboot1a, "LOI")
p2 <- f6(kboot2a, "SAI")

p1 + p2 + plot_layout(ncol = 1)
ggsave("Figure S12.pdf", width = 8, height = 12)