# Fig 2b Blood vs SVZ T cell DEG Heatmap SMART-Seq
# Matthew Buckley

rm(list = ls())
library(tidyverse)
library(reshape2)
library(viridis)
library(Seurat)
library(pvclust)
library(gplots)
library(corrplot)
library(heatmap.2x) #install.packages("devtools");devtools::install_github("TomKellyGenetics/heatmap.2x", ref="test")

setwd("~/Dropbox/DBN2019/smartseq/visuals")
load("../seurat/data/tcl_combo.rda")

#=================================================================================================
data <- as.matrix(tcl@data)
data <- data[, grepl("^Old", colnames(data))]
meta <- tcl@meta.data
meta <- meta[grepl("^Old", rownames(meta)), ]

loc_colors <- c("#C70039", "#188ad0")
rep_colors <- c("lightgrey", "darkgrey")
color <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = -1, option = "magma")
w = 8; h = 12
col_matrix <- rbind(loc_colors[factor(meta$Location)], rep_colors[factor(meta$Replicate)])
rownames(col_matrix) <- c("Location", "Replicate")

# Manual Genesets
# DEGs
load("../seurat/data/top50_DEGs.rda")
input <- data[rownames(data) %in% top50,] # SUBSET
input <- input[top50,] # REORDER

pdf(paste0("plots/Fig2b_DEG_magma_unscaled_", Sys.Date(), "Salah.pdf"), width = 5, height = 6)
p <- heatmap.2x(x = input,
     Rowv = TRUE,
     Colv = TRUE,
     revC = TRUE,
     scale = "none",
     hclustfun = hclust,
     dendrogram = "column",
     labCol = FALSE,
     key = TRUE,
     ColSideColors = col_matrix,
     cexCol = .9,
     offsetRow = 0,
     trace = "none",
     density.info = "none",
     col=color)
dev.off()

















