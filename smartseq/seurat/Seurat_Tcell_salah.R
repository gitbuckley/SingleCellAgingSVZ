# Process and Visualize Raw Count T-Cell Gene Expression Matrix with Seurat
# Matthew Buckley 8/20/18

remove(list=ls())
library(Seurat)
library(magrittr)
library(MAST)
library(ggplot2)

# Modify this path. All subsequent paths are relative. 
setwd("~/Dropbox/DBN2019/smartseq/seurat/")

#==================================================================================================
# Basic seurat processing

data_mat_1 <- read.table(file = "../data/RawStarExprsMatrix_Exp1.txt") #45169 193
data_mat_2 <- read.table(file = "../data/RawStarExprsMatrix_Exp2.txt") #45118 70
data_mat_1[181] <- NULL # Remove undeteremined sample

# Create object and basic quality checks
tcl_1 <- CreateSeuratObject(raw.data = data_mat_1, min.cells = 2, min.genes = 500, project = "tcell1") # 15409 genes across 189 samples.
tcl_2 <- CreateSeuratObject(raw.data = data_mat_2, min.cells = 2, min.genes = 500, project = "tcell2") # 14619 genes across 70 samples.
tcl <- MergeSeurat(object1 = tcl_1, object2 = tcl_2, project = "ComboSmartSeq")  # 16567 genes across 259 samples.

# Generate metadata
names <- tcl@cell.names
n2 <- gsub("^[A-z0-9]*[\\.]", "", names); n2 # Trim everything before tissue.
tissue <- gsub("(^[A-z^CD]*)(.*)", "\\1", n2); tissue # Select tissue and ignore everything after.
tissue[tissue == "BloodCD"] <- "Blood" # Correct mislabel.
age <- gsub("([A-z]*)(.*)", "\\1", names); age # Just select initial letters that describe age.
mouse <- gsub("(^[A-z1-4]*)(.*)", "\\1", names); mouse 
replicate <- mouse
replicate[replicate != "Old4"] <- 1
replicate[replicate == "Old4"] <- 2

# Add metadata
meta_df <- data.frame("Age" = age, "Location" = tissue, "Mouse" = mouse, "Tissue" = tissue, "Replicate" = replicate)
row.names(meta_df) <- names
tcl <- AddMetaData(object = tcl, metadata = meta_df)
head(tcl@meta.data)
tcl <- SetAllIdent(object = tcl, id = "Replicate")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = tcl@data), value = TRUE)
percent.mito <- Matrix::colSums(tcl@raw.data[mito.genes, ])/Matrix::colSums(tcl@raw.data)
tcl <- AddMetaData(object = tcl, metadata = percent.mito, col.name = "percent.mito")
pdf("plots/qc_3plots_salah.pdf", height = 5, width = 10)
VlnPlot(object = tcl, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# Normalize
tcl <- NormalizeData(object = tcl, normalization.method = "LogNormalize", scale.factor = 10000)

# Scale data and regress out unwanted variation
tcl <- ScaleData(object = tcl, vars.to.regress = c("nUMI", "nGene"))

# Variable Genes
tcl <- FindVariableGenes(object = tcl, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 1)
length(x = tcl@var.genes) # 2066

# Linear Dimensional Reduction - PCA
tcl <- RunPCA(tcl, pc.genes = tcl@var.genes, do.print = T, pcs.print = 1:3, genes.print = 10)

# Cluster
tcl <- FindClusters(object = tcl, reduction.type = "pca", dims.use = 1:10, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc =T)
PrintFindClustersParams(object = tcl)

# Non linear dimensional reduction - tSNE
tcl <- RunTSNE(object = tcl, dims.use = 1:10, do.fast = TRUE, perplexity = 25, set.seed = 1)
TSNEPlot(object = tcl)

# How do I create a tSNE plot where cells are colored by replicate?  First,
# store the current identities in a new column of meta.data called CellType
tcl <- StashIdent(object = tcl, save.name = "ClusterID")
head(tcl@meta.data); tail(tcl@meta.data)
# Next, switch the identity class of all cells to reflect replicate ID
tcl <- SetAllIdent(object = tcl, id = "Tissue")
TSNEPlot(object = tcl)

# Differential Expression with MAST
tcl <- SetAllIdent(object = tcl, id = "Tissue")
markers.MAST <- FindMarkers(object = tcl, only.pos = FALSE, ident.1 = "SVZ", test.use = "MAST") # trhows an error: Error in inherits(data, "SingleCellAssay") :argument "data" is missing, with no default
markers.MAST$fdr <- p.adjust(markers.MAST$p_val, method = "fdr")
markers.MAST[c("Ifng", "Pdcd1"),c("fdr")] # 2.890102e-12 1.293227e-40
write.csv(markers.MAST, "data/Tcell_MAST_DE.csv")

SVZ.markers.MAST <- FindMarkers(object = tcl, only.pos = TRUE, ident.1 = "SVZ", test.use = "MAST")
SVZ.markers.MAST$fdr <- p.adjust(SVZ.markers.MAST$p_val, method = "fdr")
Blood.markers.MAST <- FindMarkers(object = tcl, only.pos = TRUE, ident.1 = "Blood", test.use = "MAST")
Blood.markers.MAST$fdr <- p.adjust(Blood.markers.MAST$p_val, method = "fdr")

# Get the top 50 gene names
top50 <- c(rownames(SVZ.markers.MAST)[1:25], rownames(Blood.markers.MAST)[1:25])
save(top50, file = "data/top50_DEGs.rda")

# Make draft heatmap
pdf("plots/Draft_Heatmap_top50Salah.pdf", height = 9, width = 9)
DoHeatmap(object = tcl, genes.use = top50, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

#----------------------------------
# SAVE
save(tcl, file = "data/tcl_combo.rda")
load("data/Tcell_MAST_DE.csv")
sessionInfo()
