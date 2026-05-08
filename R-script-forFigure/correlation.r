library(ggplot2)
library(vegan)
library(patchwork)
library(dplyr)
library(reshape2)
library(tidyverse)
library(randomForest)
library(pheatmap)
library(RColorBrewer)
library(forcats)
library(indicspecies)
library(VennDiagram)
library(linkET)
library(igraph)
library(ggcor)

otu <- read.delim("correlation.ampID.xls", row.names = 1, check.names = FALSE)
meta_raw <- read.delim("correlation.meta.xls", row.names = 1, check.names = FALSE)
group <- meta_raw$Group
meta <- meta_raw[, !(colnames(meta_raw) %in% c("Group"))]
common_samples <- intersect(colnames(otu), rownames(meta))
otu <- otu[, common_samples]
meta <- meta[common_samples, ]
otu <- otu[rowSums(otu) > 0, ]
otu_top <- otu[order(rowSums(otu), decreasing = TRUE)[1:50], ]
otu_top <- as.data.frame(otu_top)
otu.wilcox.biomarker <- data.frame(V1 = rownames(otu_top))
metabolites.wilcox.biomarker <- data.frame(V1 = colnames(meta)) 
meta_t <- as.data.frame(t(meta))  
meta_for_func <- as.data.frame(t(meta_t))  
meta_for_func <- as.data.frame(t(meta_for_func)) 
source("correlation.dotplot.R")
res <- correlation.dotplot(otu.wilcox.biomarker,
                           metabolites.wilcox.biomarker,
                           otu_top, meta_for_func,"rumen")
ggsave("top50_correlation.pdf", res[[3]], width = 12, height = 4)
write.table(res[[1]], "correlation_r_top50.txt", sep = "\t", quote = FALSE)
write.table(res[[2]], "correlation_p_top50.txt", sep = "\t", quote = FALSE)

