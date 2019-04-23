# Boxplots of gene detection by cell type for Extended Data Figure 1 A
# Matthew Buckley

library(Seurat)
library(dplyr)
library(tidyverse)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/violin")

# Load Data
d <- load(file = "../seurat/data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata


d$sample <- factor(d$orig.ident, levels=c("y1", "y2", "y3", "o1", "o2", "o3"), ordered=T)
sampleColors <- c("deepskyblue", "#ad42f4","slateblue", "darkorange", "#f41f1f", "firebrick")
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
d$Celltype <- factor(d$Celltype,  levels=CELLS, ordered=T)

p <- ggplot(d, aes(x=Celltype, y=nGene, fill=sample))
p <- p + geom_boxplot(alpha=0.5)
p<-p+ theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p+theme(axis.text.y=element_text(size=11))
p<-p+theme(axis.text.x=element_text(size=11))
p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
p<-p+theme(axis.title=element_text(size=13))
p<-p+theme(plot.title=element_text(size=13))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+labs(y="Genes Detected")
p<-p+theme(plot.title = element_text(hjust = 0.5,size=24))
p<-p+scale_fill_manual(values=sampleColors)
p<-p+theme(axis.title.x=element_blank())

ggsave("plots/SupFig1A_GenesDetectedbyCelltype.pdf", p, width=16,height=3)
