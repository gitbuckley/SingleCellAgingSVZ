# Read in 6 gene matrices from 10X output using Seurat package.
# Combine, QC, assign clusters, save object, and make plotting dataframe.
# Perform Old/Young Differential expression using MAST on each cell type.
# Note: Takes >1 hour to run on a laptop.
# Matthew Buckley

rm(list=ls())
library(Seurat) # v2.3.4
library(dplyr) # v0.7.6
library(MAST) # v1.4.1
library(tidyverse) # v1.2.1
# See sessionInfo() at end of script for version information.

# Modify line below for code checking. All subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/seurat/")

# Read in 10X samples
y1.data <- Read10X(data.dir = "../data/Y1_Ben/filtered_gene_bc_matrices/mm10/")
o1.data <- Read10X(data.dir = "../data/O1_Ben/filtered_gene_bc_matrices/mm10/")
y2.data <- Read10X(data.dir = "../data/Y2_Ben/filtered_gene_bc_matrices/mm10/")
o2.data <- Read10X(data.dir = "../data/O2_Ben/filtered_gene_bc_matrices/mm10/")
y3.data <- Read10X(data.dir = "../data/Y3_MB/filtered_gene_bc_matrices/mm10/")
o3.data <- Read10X(data.dir = "../data/O3_MB/filtered_gene_bc_matrices/mm10/")

# Create data object
y1.obj <- CreateSeuratObject(raw.data = y1.data, project = "y1", min.cells = 0, min.genes = 0)

# Add samples to object
svz <- AddSamples(object = y1.obj, new.data = o1.data, add.cell.id = "o1")
svz <- AddSamples(object = svz, new.data = y2.data, add.cell.id = "y2")
svz <- AddSamples(object = svz, new.data = o2.data, add.cell.id = "o2")
svz <- AddSamples(object = svz, new.data = y3.data, add.cell.id = "y3")
svz <- AddSamples(object = svz, new.data = o3.data, add.cell.id = "o3")

# Inspect
svz; head(svz@meta.data); tail(svz@meta.data)
# 27998 genes across 15039 samples.

# Add labels
svz@meta.data$Age <- substr(svz@meta.data$orig.ident, 1, 1)
svz@meta.data$Replicate <- substr(svz@meta.data$orig.ident, 2, 2)

# Mitochondria Genes QC
mito.genes <- grep(pattern = "^mt-", x = rownames(x = svz@data), value = TRUE)
percent.mito <- Matrix::colSums(svz@raw.data[mito.genes, ])/Matrix::colSums(svz@raw.data)
svz@meta.data$percent.mito <- percent.mito
pdf("plots/ViolinQC_unfiltered.pdf", width = 16, height = 8)
VlnPlot(object = svz, features.plot = c("nGene", "nUMI", "percent.mito"), group.by = "orig.ident", point.size.use = 0.005)
dev.off()

# Filter
svz <- FilterCells(object = svz, subset.names = c("nGene", "percent.mito"), low.thresholds = c(400, -Inf), high.thresholds = c(4500, 0.15))
pdf("plots/ViolinQC_Filtered.pdf", width = 16, height = 8)
VlnPlot(object = svz, features.plot = c("nGene", "nUMI", "percent.mito"), group.by = "orig.ident", point.size.use = 0.005)
dev.off()

# Normalize
svz <- NormalizeData(object = svz, normalization.method = "LogNormalize", scale.factor = 10000)

# Find Variable Genes
svz <- FindVariableGenes(object = svz, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.75)
length(svz@var.genes) # 4125 genes

# Scale, Regress out nothing
svz <- ScaleData(object = svz, vars.to.regress = c())

# PCA
svz <- RunPCA(object = svz, pc.genes = svz@var.genes, do.print = FALSE)

# Cell cycle assignments based on seurat tutorial
cc.genes <- readLines(con = "../data/regev_lab_cell_cycle_genes.txt")
cc.genes <- tolower(cc.genes); cc.genes <- stringr::str_to_title(cc.genes)
s.genes <- cc.genes[1:43]; g2m.genes <- cc.genes[44:97]
svz <- CellCycleScoring(object = svz, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)

# Final cluster & tSNE choice.
svz <- FindClusters(object = svz, reduction.type = "pca", dims.use = 1:15, resolution = .25, print.output = 0, force.recalc = TRUE)
svz <- RunTSNE(object = svz, dims.use = 1:15, perplexity = 30, seed.use = 2)

pdf("plots/draft_tsne_age.pdf", width = 10, height = 10)
TSNEPlot(object = svz, group.by = "Age")
dev.off()

# Save expression data for Supplementary Tables
# m <- as.matrix(svz@raw.data)
# colnames(m)[!grepl("_", colnames(m))] <- paste0("y1_", colnames(m)[!grepl("_", colnames(m))])
# write.table(m, file="data/10x_RawExpressionMatrix.tsv", sep = "\t")

# m <- as.matrix(svz@data)
# colnames(m)[!grepl("_", colnames(m))] <- paste0("y1_", colnames(m)[!grepl("_", colnames(m))])
# write.table(m, file="data/10x_NormalizedExpressionMatrix.tsv", sep = "\t")

#=======================================================================================
# Find positive cluster marker genes
#=======================================================================================

svz.markers <- FindAllMarkers(object = svz, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25, type = "MAST")
save(svz.markers, file = "data/svz.markers.mast_All6.rda")

top5 <- svz.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
pdf("plots/SupFig1b_HeatmapCluster.pdf", width = 6.7, height = 3.9)
DoHeatmap(object = svz, genes.use = top5$gene, slim.col.label = TRUE, remove.key = TRUE, cex.row = 5)
dev.off()

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
new.cluster.ids <- c("Oligodendrocytes", "Endothelial","Microglia","Astrocytes_qNSCs",
                     "Neuroblasts","aNSCs_NPCs","Mural_cells","T_cells","OPC", "Macrophages", "Neurons")
svz@ident <- plyr::mapvalues(x = svz@ident, from = current.cluster.ids, to = new.cluster.ids)
svz <- AddMetaData(object = svz, metadata = svz@ident, col.name = "Celltype")

# Save Seurat object
save(svz, file = "data/svz_celltypes_2019-01-31.rda")

#=======================================================================================
# Convert Seurat object into metadata + data plotting dataframe
#=======================================================================================

embeds <- svz@dr$tsne@cell.embeddings
tSNE_1 <- embeds[,1] * -1 # Flip x-axis for aesthetics; does not affect any clustering/inferences.
tSNE_2 <- embeds[,2]

svz_alldata <- cbind(svz@meta.data, tSNE_1, tSNE_2, t(as.matrix(svz@data)))
svz_alldata <- svz_alldata[, !grepl("res.0*", colnames(svz_alldata))]
rownames(svz_alldata) <- NULL

# Save plotting dataframe
save(svz_alldata, file="data/svz_All6_Filtered_2019-01-31.rda")
write.csv(svz_alldata, file = "data/svz_All6_Filtered_2019-01-31.txt", quote = F)


#=======================================================================================
# Age-associated Differential Expression for each cell type
#=======================================================================================

## Append age to celltype label and add to metadata
cellAge <- paste0(svz@meta.data$Celltype, "_", svz@meta.data$Age)
m <- data.frame("Celltype_Age" = cellAge)
rownames(m) <- rownames(svz@meta.data)
svz <- AddMetaData(object = svz, metadata = m)

# Assign as main identity
svz <- SetAllIdent(svz, id="Celltype_Age")

CELLTYPES <- unique(svz@meta.data$Celltype)

mast_list <- vector(mode="list", length = 11)
names(mast_list) = CELLTYPES

# Compare young/old for each celltype, save each matrix in a list. 
for (CELLTYPE in CELLTYPES) {
	print(CELLTYPE)
	# Find  cluster marker genes
	svz.mast.de <- FindMarkers(object = svz,
								  ident.1 = paste0(CELLTYPE, "_o"),
								  ident.2 = paste0(CELLTYPE, "_y"),
								  only.pos = FALSE, 
								  min.pct = 0,
								  logfc.threshold = 0,
								  test.use = "MAST")
	svz.mast.de <- as.data.frame(svz.mast.de)
	svz.mast.de <- rownames_to_column(svz.mast.de, var = "gene")
	svz.mast.de$celltype <- rep(CELLTYPE, dim(svz.mast.de)[1])
	mast_list[[CELLTYPE]] <- svz.mast.de
}

# Combine all matrices into one dataframe
mast_df <- data.frame()
for (CELLTYPE in CELLTYPES) {
	print(CELLTYPE)
	mast_df  <- rbind(mast_df, mast_list[[CELLTYPE]])
}
dim(mast_df); head(mast_df) # 159516 by 7

# Save differential expression results
write.csv(mast_df, "../mast/data/mast_df.csv")

#=======================================================================================
# sessionInfo()

# R version 3.4.3 (2017-11-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
#  [1] forcats_0.3.0              stringr_1.3.1             
#  [3] purrr_0.2.5                readr_1.1.1               
#  [5] tidyr_0.8.1                tibble_2.0.1              
#  [7] tidyverse_1.2.1            MAST_1.4.1                
#  [9] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
# [11] matrixStats_0.54.0         Biobase_2.38.0            
# [13] GenomicRanges_1.30.3       GenomeInfoDb_1.14.0       
# [15] IRanges_2.12.0             S4Vectors_0.16.0          
# [17] BiocGenerics_0.24.0        dplyr_0.7.6               
# [19] Seurat_2.3.4               Matrix_1.2-14             
# [21] cowplot_0.9.3              ggplot2_3.1.0 

# End of Script
