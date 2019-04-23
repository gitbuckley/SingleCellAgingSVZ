# Gene Expression Vertical Violin Plots by celltype
# Matthew Buckley

rm(list=ls())
library(Seurat)
library(tidyverse)
library(cowplot)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/violin")

# Load Data
load(file = "../seurat/data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata

# Randomize rows to mix up plotting order
d <- d[sample(nrow(d)),]

# Adjust factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
d$celltype_factor <- factor(d$Celltype,  levels=CELLS, ordered=T)

# Adjust colors
new_colors <- c("#03c03c", "#0054b4", "#966fd6", "#A52A2A",
                "#ffdf00", "#aec6cf", "#ffb347", "#e5aa70",
                "#db7093", "#2F4F4F", "#e8000d")
names(new_colors) <- levels(d$celltype_factor)

#######################
# Individual PLOTTING #
#######################

# Figure 1H

gene <- "Cd8a"
p_cd8 <- ggplot(d, aes(x=d$celltype_factor, y= d$Cd8a))
p_cd8 <- p_cd8 + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3), alpha=0.5)
p_cd8 <- p_cd8 + geom_violin(trim=T,scale="width",alpha=0.5,width=0.5)
p_cd8 <- p_cd8 + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p_cd8 <- p_cd8 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_cd8 <- p_cd8 + theme(axis.title.y=element_blank())
p_cd8 <- p_cd8 + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = -10)))
p_cd8 <- p_cd8 +labs(title=gene)
p_cd8 <- p_cd8 + theme(legend.position="none")
p_cd8 <- p_cd8 + scale_color_manual(values = c("darkred", "skyblue"))
#ggsave(paste0("All6_Indiv_", gene, "_nolabel_", Sys.Date(),".pdf"), plot=p_cd8, height = 1.5, width = 5)

gene <- "Cd3e"
p_cd3 <- ggplot(d, aes(x=d$celltype_factor, y= d$Cd3e))
p_cd3 <- p_cd3 + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3), alpha=0.5)
p_cd3 <- p_cd3 + geom_violin(trim=T,scale="width",alpha=0.5,width=0.5)
p_cd3 <- p_cd3 + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p_cd3 <- p_cd3 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_cd3 <- p_cd3 + theme(axis.title.y=element_blank())
p_cd3 <- p_cd3 + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = -10)))
p_cd3 <- p_cd3 +labs(title=gene)
p_cd3 <- p_cd3 + theme(legend.position="none")
p_cd3 <- p_cd3 + scale_color_manual(values = c("darkred", "skyblue"))
#ggsave(paste0("All6_Indiv_", gene, "_nolabel_", Sys.Date(),".pdf"), plot=p_cd3, height = 1.5, width = 5)

gene <- "Cd4"
p_cd4 <- ggplot(d, aes(x=d$celltype_factor, y= d$Cd4))
p_cd4 <- p_cd4 + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3), alpha=0.5)
p_cd4 <- p_cd4 + geom_violin(trim=T,scale="width",alpha=0.5,width=0.5)
p_cd4 <- p_cd4 + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p_cd4 <- p_cd4 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_cd4 <- p_cd4 + theme(axis.title.y=element_blank())
p_cd4 <- p_cd4 + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = -10)))
p_cd4 <- p_cd4 + labs(title=gene)
p_cd4 <- p_cd4 + theme(legend.position="none")
p_cd4 <- p_cd4 + scale_color_manual(values = c("darkred", "skyblue"))
#ggsave(paste0("All6_Indiv_", gene, "_nolabel_", Sys.Date(),".pdf"), plot=p_cd4, height = 1.5, width = 5)

gene <- "Ifng"
p_ifng <- ggplot(d, aes(x=d$celltype_factor, y= d$Ifng))
p_ifng <- p_ifng + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3), alpha=0.5)
p_ifng <- p_ifng + geom_violin(trim=T,scale="width",alpha=0.5,width=0.5)
p_ifng <- p_ifng+ theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p_ifng <- p_ifng + theme(axis.text.x=element_text(size=15))
p_ifng <- p_ifng + theme(axis.text.x=element_text(angle=45, hjust=1))
p_ifng <- p_ifng + theme(axis.title.y=element_blank())
p_ifng <- p_ifng + theme(axis.title.x=element_blank())
p_ifng <- p_ifng + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = -10)))
p_ifng <- p_ifng + labs(title=gene)
p_ifng <- p_ifng + theme(legend.position="none")
p_ifng <- p_ifng + scale_color_manual(values = c("darkred", "skyblue"))
#ggsave(paste0("All6_Indiv_", gene, "_nolabel_", Sys.Date(),".pdf"), plot=p_ifng, height = 3.05, width = 5)

pq <- plot_grid(p_cd3, p_cd8, p_cd4, p_ifng, ncol = 1, rel_heights = c(1,1,1,1.7), align = "v")
ggsave(paste0("plots/Fig1H2_4plots_", Sys.Date(), ".pdf"), pq, height = 7, width = 8.4)


#===============================================
# SupFig 2B

gene <- "Ifna1"
p_ifna1 <- ggplot(d, aes(x=d$celltype_factor, y= d$Ifna1))
p_ifna1 <- p_ifna1 + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3, height=0), alpha=0.5)
p_ifna1 <- p_ifna1 + geom_violin(trim=T,scale="width", alpha=0.5,width=0.5)
p_ifna1 <- p_ifna1 + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p_ifna1 <- p_ifna1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_ifna1 <- p_ifna1 + theme(axis.title.y=element_blank())
p_ifna1 <- p_ifna1 + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = -10)))
p_ifna1 <- p_ifna1 +labs(title=gene)
p_ifna1 <- p_ifna1 + theme(legend.position="none")
p_ifna1 <- p_ifna1 + scale_color_manual(values = c("darkred", "skyblue")) + expand_limits(y = c(0, 2))

gene <- "Ifnb1"
p_ifnb1 <- ggplot(d, aes(x=d$celltype_factor, y= d$Ifnb1))
p_ifnb1 <- p_ifnb1 + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3, height=0), alpha=0.5)
p_ifnb1 <- p_ifnb1 + geom_violin(trim=T,scale="width",alpha=0.5,width=0.5)
p_ifnb1 <- p_ifnb1 + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p_ifnb1 <- p_ifnb1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_ifnb1 <- p_ifnb1 + theme(axis.title.y=element_blank())
p_ifnb1 <- p_ifnb1 + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = -10)))
p_ifnb1 <- p_ifnb1 + labs(title=gene)
p_ifnb1 <- p_ifnb1 + theme(legend.position="none")
p_ifnb1 <- p_ifnb1 + scale_color_manual(values = c("darkred", "skyblue")) + expand_limits(y = c(0, 2))

gene <- "Ifng"
p_ifng <- ggplot(d, aes(x=d$celltype_factor, y= d$Ifng))
p_ifng <- p_ifng + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3), alpha=0.5)
p_ifng <- p_ifng + geom_violin(trim=T,scale="width",alpha=0.5,width=0.5)
p_ifng <- p_ifng+ theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p_ifng <- p_ifng + theme(axis.text.x=element_text(size=15))
p_ifng <- p_ifng + theme(axis.text.x=element_text(angle=45, hjust=1))
p_ifng <- p_ifng + theme(axis.title.y=element_blank())
p_ifng <- p_ifng + theme(axis.title.x=element_blank())
p_ifng <- p_ifng + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = -10)))
p_ifng <- p_ifng + labs(title=gene)
p_ifng <- p_ifng + theme(legend.position="none")
p_ifng <- p_ifng + scale_color_manual(values = c("darkred", "skyblue"))

pp <- plot_grid(p_ifna1, p_ifnb1, p_ifng, ncol = 1, rel_heights = c(1,1,1.7), align = "v")
ggsave(paste0("plots/SupFig2b_3plots_", Sys.Date(), ".pdf"), pp, height = 7.5, width = 8.4)


######################
# AUTOMATIC PLOTTING #
######################

# int_genes<-c("Cd8a","Ifng","Cd3e","Cd4", "Ifna1", "Ifnb1")

# for (gene in int_genes) {
# 	print(gene)
# 	p<-ggplot(d, aes(x=d$celltype_factor, y= d[,gene]))
# 	p <- p + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3), alpha=0.5)
# 	p <- p + geom_violin(trim=T,scale="width",alpha=0.5,width=0.5)
# 	p<-p+ theme_classic() + theme(
# 		axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
# 		axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
# 	p<-p+theme(axis.text.y=element_text(size=15))
# 	p<-p+theme(axis.text.x=element_text(size=15))
# 	p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
# 	p<-p+theme(axis.title=element_text(size=20))
# 	p<-p+theme(plot.title=element_text(size=20))
# 	p<-p+theme(axis.title.y=element_text(vjust=1))
# 	p<-p+theme(axis.title.x=element_text(vjust=-0.10))
# 	p<-p+theme(plot.title = element_text(hjust = 0.5,size=15))
# 	p<-p+theme(axis.title.x=element_blank())
# 	p<-p+labs(title=gene)
# 	p <- p + theme()
# 	p <- p + scale_color_manual(values = c("darkred", "skyblue"))
# 	print(p)
# 	ggsave(paste0("All6_Fig1i_", gene, "_", Sys.Date(),".pdf"), height = 3, width = 8)
# }

# # No labels or ticks
# for (gene in int_genes) {
# 	print(gene)
# 	p<-ggplot(d, aes(x=d$celltype_factor, y= d[,gene]))
# 	p <- p + geom_jitter(aes(color=Age), size = 0.7, position = position_jitter(width=.3), alpha=0.5)
# 	p <- p + geom_violin(trim=T,scale="width",alpha=0.5,width=0.5)
# 	p<-p+ theme_classic() + theme(
# 		axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
# 		axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
# 	p<-p+theme(axis.text.y=element_blank())
# 	p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
# 	p<-p+theme(axis.title=element_text(size=20))
# 	p<-p+theme(plot.title=element_text(size=20))
# 	p<-p+theme(axis.title.y=element_text(vjust=1))
# 	p<-p+theme(plot.title = element_text(hjust = 0.5,size=15))
# 	p<-p+labs(title=gene)
# 	p <- p + theme(legend.position="none")
# 	p <- p + scale_color_manual(values = c("darkred", "skyblue"))
# 	print(p)
# 	ggsave(paste0("All6_Fig1i_", gene, "_nolabel_", Sys.Date(),".pdf"), height = 2, width = 8)
# }

#================================================

sessionInfo()

# R version 3.4.3 (2017-11-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] bindrcpp_0.2.2  forcats_0.3.0   stringr_1.3.1   dplyr_0.7.6    
#  [5] purrr_0.2.5     readr_1.1.1     tidyr_0.8.1     tibble_2.0.1   
#  [9] tidyverse_1.2.1 Seurat_2.3.4    Matrix_1.2-14   cowplot_0.9.3  
# [13] ggplot2_3.1.0 

# End of Script

























