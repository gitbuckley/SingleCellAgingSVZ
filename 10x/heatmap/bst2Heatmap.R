# Bst2 Heatmap

library(Seurat)
library(tidyverse)
library(cowplot)
library(heatmap.2x)
library(viridis)

# Modify this path. All following are relative. 
setwd("~/Desktop/Dropbox/DBN2019/10x/heatmap")

# Load Data
load("../seurat/data/svz_celltypes_2019-01-31.rda")
load("../data/mouse.gene.hallmark.kegg.reactome.rda")

#============================================================================
svz__data <- as.matrix(svz@data)
svz__data[1:11,1:11]
rownames(svz__data) <- tolower(rownames(svz__data))
genelist<-tolower(hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE)
length(genelist); head(genelist)
genelist_name<-"HALLMARK_IFNG_RESPONSE"


# Get cells of interest column locations
cells_loc <- svz@meta.data$Celltype %in% c("aNSCs_NPCs", "Astrocytes_qNSCs") # 2463 TRUE

# Subset expression data to ifn genes and cells of interest
ifn_data <- svz__data[rownames(svz__data) %in% genelist, cells_loc ] # dim: 310 2463

# Ensure no rows have zero sum.
ifn_data <- ifn_data[rowSums(ifn_data)>1, ] # 246 2463

specific_genes <- c("b2m", "ifi27", "ifitm3", "ifit3", "ifit1", "xaf1", "bst2")
ifn_data_sp <- ifn_data[rownames(ifn_data) %in% specific_genes, ]

# Subset metadata to cells of interest
meta_data <- svz@meta.data[cells_loc, ]

# Heatmap
color <- magma(n = 20, direction = -1)
age_colors <- c("firebrick", "deepskyblue")

pdf(file = paste0("plots/SupFig5e_ifn_aNSC_qNSC_", Sys.Date(), "_.pdf"), width = 8, height = 6)
heatmap.2(x = ifn_data_sp,
          Rowv = T,
          Colv = T,
          key = T,
          hclustfun = hclust,
          revC=TRUE,
          dendrogram = "both",
          labCol = F,
          breaks=seq(0,3.5,length.out=21),
          ColSideColors = age_colors[factor(meta_data$Age)],
          trace = "none",
          col=color,
          main = "Ifng Response Gene Expression")
dev.off()

