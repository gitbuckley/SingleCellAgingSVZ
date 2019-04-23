# Manhattan-like plots highlighting DEGs by Cell type using MAST DE values.
# Cell type proportions by sample.
# Figure 1C
# Matthew Buckley 1/8/2019
 
rm(list = ls())
library(tidyverse)
library(tibble)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(NCmisc) #v1.1.5

# See sessionInfo() at end of script for version information.

# Modify path below. All subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/mast")

#==================================================================================================
df <- read_csv("data/mast_df.csv")

# Add z-score based on two sided null hypothesis.
df$z <- p.to.Z(df$p_val) * sign(df$avg_logFC)
df$z.adj <- p.to.Z(df$p_val_adj) * sign(df$avg_logFC)

# Randomize rows to reduce overplotting issues.
df <- df[sample(nrow(df)), ]

# Adjust factors
CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
           "OPC", "Oligodendrocytes", "Endothelial", "Mural_cells",
           "Microglia", "Macrophages", "T_cells")
df$celltype_factor <- factor(df$celltype,  levels=CELLS, ordered=T)

#==================================================================================================

# Adjust colors
new_colors <- c("#03c03c", "#0054b4", "#966fd6", "#A52A2A",
                "#ffdf00", "#aec6cf", "#ffb347", "#e5aa70",
                "#db7093", "#2F4F4F", "#e8000d")
names(new_colors) <- levels(df$celltype_factor)

# Make custom color column to facilitate grey coloring by threshold.
col <- new_colors[df$celltype_factor]
col[df$p_val_adj > 0.05] <- "#D3D3D3" # grey
df$col <- as.factor(col)

q <- ggplot(df, aes(x = celltype_factor, y = z, color = col)) +
  geom_jitter(width = 0.40, alpha = .55, size = 1) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("Changes in gene expression with age") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
q

#==================================================================================================
# Celltype age enrichment for each replicate
#==================================================================================================

load(file = "../seurat/data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata

# Group by cell type and sample, then count cells.
# Then divide sample cell type counts by total sample cell counts. 
keep = c("Celltype", "orig.ident")
counts <- d %>% select(one_of(keep)) %>% group_by(orig.ident, Celltype) %>% dplyr::mutate(count = n()) %>% distinct()
totals <- aggregate(counts$count, by=list(orig.ident=counts$orig.ident), FUN=sum)
colnames(totals) <- c("orig.ident", "sample.sum")
prop_df <- merge(counts, totals)
prop_df$prop <- prop_df$count / prop_df$sample.sum
write.csv(prop_df, file = "data/SupTable3Sheet3_CelltypeCounts.csv")

enrich <- c(); rep <- c(); cell <- c(); y_prop <- c(); o_prop <- c()
for (celltype in unique(prop_df$Celltype)) {
    for (r in c(1:3)) {
        old_prop <- prop_df[prop_df$orig.ident==paste0("o", r) & prop_df$Celltype==celltype, "prop"]
        young_prop <- prop_df[prop_df$orig.ident==paste0("y", r) & prop_df$Celltype==celltype, "prop"]
        y_prop <- c(y_prop, young_prop)
        o_prop <- c(o_prop, old_prop)
        e <- old_prop/young_prop
        enrich <- c(enrich, e)
        rep <- c(rep, r)
        cell <- c(cell, celltype)
    }
}

enrich_df <- data.frame("Rep" = rep, "Celltype" = cell, "Y_prop" = y_prop, "O_prop" = o_prop, "Log2Enrich" = log2(enrich))

# Supplementary Table 3, sheet 2:
write.csv(enrich_df, file = "data/SupTable3Sheet2_ProportionEnrichments.csv" )

#==================================================================================================

library(plyr)
# Calculate standard error
sumr <- plyr::ddply(enrich_df, c("Celltype"), summarise,
               N    = length(Log2Enrich),
               mean = mean(Log2Enrich),
               sd   = sd(Log2Enrich),
               se   = sd / sqrt(N)
)
enrich_df <- merge(enrich_df, sumr, by = "Celltype")
head(enrich_df)

# Adjust colors
enrich_df$celltype_factor <- factor(enrich_df$Celltype,  levels=CELLS, ordered=T)
names(new_colors) <- levels(enrich_df$celltype_factor)

p <- ggplot(enrich_df, aes(x = celltype_factor, y = Log2Enrich)) +
  geom_bar(aes(fill = celltype_factor), stat = "summary", fun.y = "mean",  alpha = .5) +
  geom_point(aes(shape = as.factor(Rep), color = celltype_factor), size = 2.5, alpha = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_errorbar(aes(x=celltype_factor, ymin=mean-se, ymax=mean+se), width=.15) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20, face = "plain", margin = margin(r = 20))) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Log2 Enrichment") +
  theme(legend.position="none") +
  scale_fill_manual(values = new_colors) +
  scale_color_manual(values = new_colors) +
  ggtitle("Changes in cell number with age") +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
p


#==================================================================================================
# Combine into one plot

pq <- plot_grid(p, NULL, q, ncol = 1, rel_heights = c(1, 0.15, 1.5))
ggsave("plots/Fig1C_ProportionsPlusDE.png", pq, height = 9, width = 8.4, dpi = 700)

#==================================================================================================
sessionInfo()

# R version 3.4.3 (2017-11-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
#  [1] plyr_1.8.4                 NCmisc_1.1.5              
#  [3] scales_1.0.0               viridis_0.5.1             
#  [5] viridisLite_0.3.0          RColorBrewer_1.1-2        
#  [7] bindrcpp_0.2.2             forcats_0.3.0             
#  [9] stringr_1.3.1              purrr_0.2.5               
# [11] readr_1.1.1                tidyr_0.8.1               
# [13] tibble_2.0.1               tidyverse_1.2.1           
# [15] MAST_1.4.1                 SummarizedExperiment_1.8.1
# [17] DelayedArray_0.4.1         matrixStats_0.54.0        
# [19] Biobase_2.38.0             GenomicRanges_1.30.3      
# [21] GenomeInfoDb_1.14.0        IRanges_2.12.0            
# [23] S4Vectors_0.16.0           BiocGenerics_0.24.0       
# [25] dplyr_0.7.6                Seurat_2.3.4              
# [27] Matrix_1.2-14              cowplot_0.9.3             
# [29] ggplot2_3.1.0   

# End of script
