# Plotting and analysis of seurat cell cycle calls.
# See 10x/seurat/Seurat_All6.R
# Matthew Buckley

rm(list = ls())
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggthemes)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/tsne")

# Load Data
load(file = "../seurat/data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata

CELLTYPES <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts",  "Neurons",
 "OPC", "Oligodendrocytes","Endothelial", "Mural_cells","Microglia", "Macrophages", "T_cells")
d$Celltypes <- factor(d$Celltype, levels=CELLTYPES, ordered = T)
d$Age <- factor(d$Age, levels=c("y","o"), ordered=T)

###### Cell Cycle tSNEs #######

q <- ggplot(data = d, aes(x=tSNE_1, y=tSNE_2, color=Phase))
q <- q + geom_point(alpha = 0.4, size=.7)
q <- q + ggtitle("Cell Cycle")
q <- q + scale_color_tableau(palette="Tableau 10", type="regular")
ggsave(paste0("plots/SupFig4h_CellCycle_", Sys.Date(), ".pdf"), q, height = 7, width = 7)

# tSNE
q <- ggplot(data = d, aes(x=tSNE_1, y=tSNE_2, color=Phase))
q <- q + geom_point(alpha = 0.5, size=.8)
q <- q + facet_wrap(~Age)
q <- q + ggtitle("Cell Cycle")
q <- q + scale_color_tableau(palette="Tableau 10", type="regular")
ggsave(paste0("plots/CellCycle_Age_", Sys.Date(), ".pdf"), q, height = 7, width = 12)

