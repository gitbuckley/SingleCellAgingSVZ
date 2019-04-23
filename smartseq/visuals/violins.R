# Make plotting dataframe from seurat object
# Matthew Buckley

#==================================================================================================
rm(list = ls())
library(tidyverse)
library(reshape2)
library(viridis)
library(Seurat)

setwd("~/Desktop/Dropbox/DBN2019/smartseq/visuals")
load("../seurat/data/tcl_combo.rda")

#==================================================================================================
# Load Seurat Object and make tidy dataframe for plotting.

data <- as.data.frame(t(as.matrix(tcl@data)))
data <- tibble::rownames_to_column(data, "cell")

# Long format
data.melt <- melt(data, id.vars = c("cell"), value.name = "expression")
colnames(data.melt)[2] <- "gene"

meta <- tcl@meta.data
meta <- tibble::rownames_to_column(meta, "cell")
full <- merge(meta, data.melt)

plot.df <- full

#==================================================================================================
# Step: Specify Genes of Interest


genes3 <- c("Cx3cr1")
# Notes: Havcr2 = Tim-3;   Spn = Cd43,  Prdm1 = Blimp-1
exhaustion <- c("Havcr2", "Pdcd1", "Lag3", "Cd44", "Cd43", "Spn", "Cd69", "Prdm1", "Ctla4")
ex.inhib <- c("Pdcd1", "Lag3", "Ctla4", "Havcr2", "Cd160" )
effector <- c("Gzmb", "Ifng", "Tnf", "Il2")
genes_new <- c("Cd3e", "Cd8a", "Cd4", "Cd62L", "Cd44", "Cd69",
             "Xcl1", "Itgal", "Itga4",  "Lamp1", "Gzmb", "Ifng")
genes_CD39 <-  c("Pdcd1", "Ifng", "Entpd1") # CD39
fig2genes <- c("Ifng", "Pdcd1")


#==================================================================================================
# Plot violins.
# Subset to genes of interest
dat <- dplyr::filter(plot.df, gene %in% fig2genes, Age == "Old")
dat$gene <- factor(dat$gene, levels = rev(fig2genes))

p <- ggplot(data = dat, aes(x = gene, y = expression, fill = Location)) +
            geom_violin(alpha = .75, scale = "width", position=position_dodge(.8)) +
            geom_point(alpha = 0.35, size = .7, stroke = 0, position = position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = .8)) +
            theme_classic() +
            labs(y = "Expression") +
            ggtitle("T-cell Phenotype") +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text.x=element_text(size=12)) +
            theme(axis.text.y=element_text(size=12)) +
            theme(axis.title.x = element_text(size = 12))
ggsave(paste0("plots/Fig2c_Pd_Ifn_", Sys.Date(), ".pdf"), p, height = 1.8, width = 3.1, bg = "transparent")

