# SMART-Seq V5 Single-Cell T-cell Naive/Memory/Effect Characterization Analysis - tSNE
# Matthew Buckley

remove(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)

setwd("~/Dropbox/DBN2019/smartseq/visuals")
load("../seurat/data/tcl_combo.rda")

#==================================================================================================
# Custom tSNE
meta <- tcl@meta.data
d <- data.frame(id = colnames(tcl@data), Location = meta$Tissue,
                Age = meta$Age, Mouse = meta$Mouse, Replicate = meta$Replicate,
                tcl@dr$tsne@cell.embeddings)

pdf("plots/Fig2a_tsne_combinedOld_Salah.pdf")
p <- ggplot(data =  filter(d, Age == "Old"), aes(x = tSNE_1, y = tSNE_2, color = Location)) +
            geom_point(size = 5, alpha = .6) +
            ggtitle("tSNE Projection of T-Cell Transcriptomes") +
            scale_color_manual(values = c("#C70039", "#188ad0"))
p
dev.off()
#===========================
pdf("plots/SupFig3a_tsne_replicate_Salah.pdf")
p <- ggplot(data = filter(d, Age == "Old"), aes(x = tSNE_1, y = tSNE_2, color = Replicate)) +
            geom_point(size = 5, alpha = .6) +
            ggtitle("tSNE Projection of T-Cell Transcriptomes") +
            scale_color_manual(values = c("#C70039", "#188ad0")) +
            scale_shape_manual(values = c(16, 1))
p
dev.off()

#===========================
pdf("plots/SupFig3a_tsne_mouse_Salah.pdf")
p <- ggplot(data = filter(d, Age == "Old"), aes(x = tSNE_1, y = tSNE_2, color = Mouse)) +
            geom_point(size = 5, alpha = .6) +
            ggtitle("tSNE Projection of T-Cell Transcriptomes") +
            scale_color_brewer(palette = "Set1", direction = 1,
  aesthetics = "colour")
p
dev.off()




