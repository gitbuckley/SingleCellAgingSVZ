# Fig3 Interferon Gamma Response Signature Violins
# Matthew Buckley

rm(list=ls())
library(Seurat)
library(ggplot2)
library(ggpubr)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/violin")

# Load Data
load(file = "../seurat/data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata

#===================================================================================================
# Ifng response genes
# Figure 3B

load(file = "../data/mouse.gene.hallmark.kegg.reactome.rda")
genelist <- tolower(hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE)

#d <- as.data.frame(hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE)
#write.csv(d, "data/Hallmark_IfngResponseGenes_Mouse.csv")

# Subset data to signature genes
colnames(d) <- tolower(colnames(d))
ifn_data <- d[, colnames(d) %in% genelist] # 310 genes
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype")]
ifn_data$ifng_response <- rowSums(ifn_data)
ifn_data <- cbind(meta, ifn_data)

# Reorder factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
ifn_data$celltype <- factor(ifn_data$celltype,  levels=CELLS, ordered=T)
ifn_data$age <- factor(ifn_data$age, levels=c("y", "o"), ordered=T)


# Plot parameters
a1 = 0.25
a2 = 0.6
s = 0.4
ageColors <- c("deepskyblue", "firebrick")
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
					symbols = c("****", "***", "**", "*", "ns"))

p <- ggplot(data=dplyr::filter(ifn_data, celltype != "Neurons"), aes(x=age, y=ifng_response)) +
		    geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.35, h = 0)) +
		    geom_violin(aes(fill=age), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
		    scale_fill_manual(values=ageColors) +
		    facet_wrap(~celltype, scales = "free", nrow =2) +
		    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
			theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
			theme(strip.text.x = element_text(size = 7)) +
			theme(legend.position="none") +
			stat_compare_means(size = 1.7)
ggsave(paste0("plots/Fig3b_Stats_IFNgResponse_10celltypes_", Sys.Date(), ".pdf"), p, height=2.5, width=6)


#=====================================================================================
# Inspect receptor signature. 
# SupFig 4b

# Subset data to signature genes
colnames(d) <- tolower(colnames(d))
ifnr_data <- d[, colnames(d) %in% c("ifngr1", "ifngr2")]
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype")]
ifnr_data$ifngr <- rowSums(ifnr_data)
ifnr_data <- cbind(meta, ifnr_data)

# Reorder factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
ifnr_data$celltype <- factor(ifnr_data$celltype,  levels=CELLS, ordered=T)
ifnr_data$age <- factor(ifnr_data$age, levels=c("y", "o"), ordered=T)


# Plot parameters
a1 = 0.25
a2 = 0.6
s = 0.4
ageColors <- c("deepskyblue", "firebrick")

p <- ggplot(data=dplyr::filter(ifnr_data, celltype != "Neurons"), aes(x=age, y=ifngr)) +
		    geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.35, h = 0)) +
		    geom_violin(aes(fill=age), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
		    scale_fill_manual(values=ageColors) +
		    facet_wrap(~celltype, scales = "free", nrow =2) +
		    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
			theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
			theme(strip.text.x = element_text(size = 7)) +
			theme(legend.position="none") +
			stat_compare_means()
ggsave(paste0("plots/SupFig4b_IFNgReceptors_10celltypes_", Sys.Date(), ".pdf"), p, height=2.25, width=6)

#=====================================================================================
# Bst2
# Sup Figure 4i

d <- svz_alldata
# Reorder factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
d$celltype <- factor(ifnr_data$celltype,  levels=CELLS, ordered=T)
d$age <- factor(ifnr_data$age, levels=c("y", "o"), ordered=T)
p <- ggplot(data=dplyr::filter(d, celltype != "Neurons"), aes(x=age, y=Bst2)) +
		    geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.35, h = 0)) +
		    geom_violin(aes(fill=age), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
		    scale_fill_manual(values=ageColors) +
		    facet_wrap(~celltype, nrow =2) +
		    scale_y_continuous(limits = c(0,NA)) +
		    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
			theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
			theme(strip.text.x = element_text(size = 7)) +
			theme(legend.position="none")
			#stat_compare_means()
ggsave(paste0("plots/SupFig4i_Bst2_10celltypes_", Sys.Date(), ".pdf"), p, height=2.25, width=6)


