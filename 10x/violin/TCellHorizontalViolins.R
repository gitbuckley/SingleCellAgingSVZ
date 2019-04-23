# 10X Single-Cell T-cell Naive/Memory/Effect Characterization Analysis - Violin Plots
# Matthew Buckley

library(tidyverse)
library(reshape2)
library(viridis)
library(Seurat)

#==================================================================================================
# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/violin")

# Load Data
load(file = "../seurat/data/svz_All6_Filtered_2019-01-31.rda")
data <- svz_alldata

# Filter out non-t-cells
tcells <- data[data$Celltype == "T_cells",]
dim(tcells) # 255 28008.

#==================================================================================================
# Specify Genes of Interest
genes <- c("Cd3e", "Cd8a", "Cd4",  "Cd62L", "Cd44", "Cd69", "Xcl1","Itgal", "Itga4", "Ifng")

# Rename Sell as Cd62L
colnames(tcells)[colnames(tcells)=="Sell"] <- "Cd62L"

# Subset to genes of interest
tcells_sub <- tcells[, colnames(tcells) %in% c("Age", genes)]
dim(tcells_sub) # 255  11

# Tidy genes
tidyT <- gather(tcells_sub, key = "gene", value = "expression", Cd62L:Cd3e)
tidyT$gene <- factor(tidyT$gene, levels = rev(genes))

#==================================================================================================
# Plot vertical violins.

p <- ggplot(data = tidyT, aes(x = gene, y = expression)) +
			geom_violin(fill = "#e8000d", alpha = .75, scale = "width") +
			geom_point(alpha = 0.4, size = .3, position = position_jitter(width=.3, height=0)) +
			theme_classic() +
			labs(y = "Expression") +
			ggtitle("T-cell Phenotype") +
  			theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  			theme(axis.text.x=element_text(size=8)) +
  			theme(axis.text.y=element_text(size=8, margin = margin(r = 10))) +
  			theme(axis.ticks.y = element_blank()) +
  			theme(axis.title.y = element_text(margin = margin(r = 40))) +
  			theme(axis.title.x = element_text(size = 8)) +
  			coord_flip()
ggsave(paste0("plots/Fig1I_PreAnnotation_TCellViolin_", Sys.Date(), ".pdf"),
       p, height = 4, width = 2.7, bg = "transparent")
