# Make pretty tSNE's for all 6 combined 10x replicates.
# Figures 1ab and SupFig 1c
# Matthew Buckley

library(Seurat)
library(tidyverse)
library(cowplot)

# Load Data
setwd("~/Desktop/Dropbox/BenNSCProject/10X/Seurat/")
d <- read.csv(file = "data/svz_All6_Filtered_2019-01-31.rda", header=T)
d <- d[,-1]
d[1:3,1:11]

# Randomize rows to mix up plotting order
d <- d[sample(nrow(d)),]

# Plot parameters
a <-  0.7
s <- .8
sampleColors <- c("#f41f1f", "darkorange", "firebrick", "deepskyblue", "#ad42f4","slateblue")
ageColors <- c("firebrick", "deepskyblue")
ben_colors <- c("#03c03c","#0054b4","#966fd6","#aec6cf","#ffdf00",
	"#ffb347","#e5aa70","#db7093","#e8000d", "#555f6d", "brown")
batchColors <- c("#03c03c","#0054b4","#966fd6")


### Sample
q <- ggplot(data = d, aes(x = tSNE_1, y = tSNE_2, color = orig.ident))
q <- q + geom_point(alpha=a, size=s)
q <- q + scale_colour_manual(values=sampleColors)
q <- q + guides(colour = guide_legend(override.aes = list(size=3)))
#q <- q + theme(legend.position="none")
q_sample <- q + ggtitle("")# + 
#		 theme(axis.line = element_blank(), axis.text = element_blank(),
#		 axis.ticks = element_blank(), axis.title = element_blank())
q_sample
ggsave(paste0("plots/All6/customeTSNE_LEG_sample_size.8_", Sys.Date(), ".pdf"), height = 7, width = 7)


### Replicate/Batch
q <- ggplot(data = d, aes(x = tSNE_1, y = tSNE_2, color = as.factor(Replicate)))
q <- q + geom_point(alpha=a, size=s)
q <- q + scale_colour_manual(values=batchColors)
q <- q + guides(colour = guide_legend(override.aes = list(size=3)))
#q <- q + theme(legend.position="none")
q_batch <- q + ggtitle("")# + 
#		 theme(axis.line = element_blank(), axis.text = element_blank(),
#		 axis.ticks = element_blank(), axis.title = element_blank())
q_batch
ggsave(paste0("plots/All6/customeTSNE_LEG_batch", Sys.Date(), ".pdf"), height = 7, width = 8)


### Age
q <- ggplot(data = d, aes(x = tSNE_1, y = tSNE_2, color = Age))
q <- q + geom_point(alpha=a, size=s)
q <- q + scale_colour_manual(values=ageColors)
q <- q + guides(colour = guide_legend(override.aes = list(size=3)))
q2 <- q + ggtitle("")# + 
#		 theme(axis.line = element_blank(), axis.text = element_blank(),
#		 axis.ticks = element_blank(), axis.title = element_blank())
q2
ggsave(paste0("plots/All6/customeTSNE_age_axis_legend", Sys.Date(), ".pdf"), height = 7, width = 8)

q <- ggplot(data = d, aes(x = tSNE_1, y = tSNE_2, color = Age))
q <- q + geom_point(alpha=a, size=s)
q <- q + scale_colour_manual(values=ageColors)
q <- q + theme(legend.position="none")
q2 <- q + ggtitle("")# + 
#		 theme(axis.line = element_blank(), axis.text = element_blank(),
#		 axis.ticks = element_blank(), axis.title = element_blank())
q2
ggsave(paste0("plots/All6/customeTSNE_age_axis_NOLegend", Sys.Date(), ".pdf"), height = 7, width = 6)

### Cell type
d$Celltypes <- factor(d$Celltype, levels=c("Astrocytes_qNSCs","aNSCs_NPCs","Neuroblasts",
	"Oligodendrocytes","OPC","Endothelial", "Mural_cells","Microglia","T_cells", 
	"Macrophages", "Neurons"), ordered = T)

q <- ggplot(data = d, aes(x = tSNE_1, y = tSNE_2, color = Celltypes))
q <- q + geom_point(alpha=a, size=s)
q <- q + scale_colour_manual(values=ben_colors)
#q <- q + guides(colour = guide_legend(override.aes = list(size=3)))
q <- q + theme(legend.position="none")
q3 <- q + ggtitle("")# + 
#		 theme(axis.line = element_blank(), axis.text = element_blank(),
#		 axis.ticks = element_blank(), axis.title = element_blank())
q3
ggsave(paste0("plots/All6/customeTSNE_celltype_axis_nolegend", Sys.Date(), ".pdf"), height = 7, width = 6)

plot_grid(q3, q2, ncol = 2, rel_widths = c(1, 1))
ggsave(paste0("plots/All6/customeTSNE_celltype_axis.pdf", height = 7, width = 12))


# Young or Old ONLY
q_o <- ggplot(data = dplyr::filter(d, Age == "o"), aes(x = tSNE_1, y = tSNE_2))
q_o <- q_o + geom_point(alpha=a, size=s, color = "firebrick")
q_o <- q_o + ggtitle("")
q_o
q_y <- ggplot(data = dplyr::filter(d, Age == "y"), aes(x = tSNE_1, y = tSNE_2))
q_y <- q_y + geom_point(alpha=a, size=s, color = "deepskyblue")
q_y <- q_y + ggtitle("")

pq <- plot_grid(q_y, q_o, ncol = 2)
pq
ggsave(paste0("plots/All6/Separated_Age", Sys.Date(), ".pdf"), pq, height = 9, width = 17)

fourplots <- plot_grid(q_batch, q_sample, q_y, q_o, ncol = 4)
fourplots
ggsave(paste0("plots/All6/Fourplots", Sys.Date(), ".pdf"), fourplots, height = 6, width = 17)
