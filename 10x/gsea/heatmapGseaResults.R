# Visualize Cell type GSEA Results

library(tidyverse)
library(scales) 

setwd("~/Desktop/Dropbox/DBN2019/10x/gsea/")

# HALLMARKS
d <- read.table("data/ALL_fgsea_2019-03-19.txt")
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
       "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
       "Microglia", "Macrophages", "T_cells")
d$celltype <- factor(d$celltype,  levels=CELLS, ordered=T)
d$pathway <- gsub("HALLMARK_", "", d$pathway)
d$pathway <- gsub("_", " ", d$pathway)
d$NES[d$padj > 0.05] <- NA # grey out non-significant
d <- d %>% filter(!celltype %in% c("Neurons", "Mural_cells"))

p <- ggplot(d, aes(x=celltype, y=pathway, fill=NES)) + geom_tile(colour="white",size=0.25) +
 scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), na.value="grey80") +
 theme_minimal() +
 labs(x="", y="") +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("plots/HM_heatmap_", Sys.Date(), ".pdf"), height=8, width=6)
