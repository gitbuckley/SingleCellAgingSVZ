# Run basic GSEA analysis on MAST DEG data for each celltype
# Matthew Buckley 2/15/19

# Credit to Stephen Turner for a very nice tutorial.
# https://stephenturner.github.io/deseq-to-fgsea/

# See end of script for package versions.
library(tidyverse)
library(biomaRt)
library(fgsea)
library(DT)
library(org.Mm.eg.db)
library(msigdf)
library(NCmisc)
library(clusterProfiler)


setwd("~/Desktop/Dropbox/DBN2019/10x/gsea/")

# Load cell type differential expression results
mast <- read.csv("../mast/data/mast_df.csv")
mast$z <- p.to.Z(mast$p_val) * sign(mast$avg_logFC)
mast$z.adj <- p.to.Z(mast$p_val_adj) * sign(mast$avg_logFC)

# Load pathways. Genes are as human gene symbols.
pathways.hallmark <- gmtPathways("../data/h.all.v6.2.symbols.gmt")

# Get extra gene information, combine with DE results
gene_bitr <- bitr(mast$gene, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL", "ENSEMBL"), OrgDb = org.Mm.eg.db)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>%
  na_if("") %>% 
  na.omit()
gene_info <- inner_join(bm, gene_bitr, by = c("ensembl_gene_id" = "ENSEMBL"))
df <- inner_join(mast, gene_info, by = c("gene" = "SYMBOL"))

CELLTYPES <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")

gsea_all <- NULL
for (cell in CELLTYPES) {

	print(cell)

	# Filter for celltype, reduce to just gene and DE p value statistic.
	df2 <- df %>% dplyr::filter(celltype == cell)
	print(dim(df2))

	df2 <- df2 %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
	  na.omit() %>% 
	  distinct() %>% 
	  group_by(hsapiens_homolog_associated_gene_name) %>% 
	  summarize(stat=mean(z))

	# Genes ranked by Z score (not log fold change)
	ranks <- deframe(df2)

	# Run the Gene Set Enrichment Analysis
	fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

	fgseaResTidy <- fgseaRes %>%
	  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
	  as_tibble() %>%
	  arrange(desc(NES))
	print(head(fgseaResTidy))

	fgseaResTidy <- as.data.frame(fgseaResTidy)
	fgseaResTidy$`celltype` <- cell

	fname <- paste0("data/", cell, "_fgsea_", Sys.Date(),".txt")

	# Save individual tables
	write.table(fgseaResTidy, file=fname, sep = "\t")
	gsea_all <- rbind(gsea_all, fgseaResTidy)
}

# Save complete dataframe
fname <- paste0("data/ALL_fgsea_", Sys.Date(),".txt")
write.table(gsea_all, file=fname, sep = "\t")

# Optional: Inspect results in browser.
gsea_all %>% 
  arrange(-NES) %>% 
  DT::datatable()

sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
#  [1] clusterProfiler_3.6.0 DOSE_3.4.0            NCmisc_1.1.5         
#  [4] msigdf_5.2            org.Mm.eg.db_3.5.0    AnnotationDbi_1.40.0 
#  [7] IRanges_2.12.0        S4Vectors_0.16.0      Biobase_2.38.0       
# [10] BiocGenerics_0.24.0   DT_0.5                fgsea_1.4.1          
# [13] Rcpp_1.0.0            biomaRt_2.34.2        bindrcpp_0.2.2       
# [16] eulerr_5.1.0          reshape2_1.4.3        forcats_0.3.0        
# [19] stringr_1.3.1         dplyr_0.7.6           purrr_0.2.5          
# [22] readr_1.1.1           tidyr_0.8.1           tibble_2.0.1         
# [25] ggplot2_3.1.0         tidyverse_1.2.1   

