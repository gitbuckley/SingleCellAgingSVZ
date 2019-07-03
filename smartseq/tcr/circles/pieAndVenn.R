# Pie Charts from T Cell TCR-B Chain Clonotypes

rm(list=ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(eulerr)

# Modify path if needed. All sebsequent paths are relative. 
setwd("~/Desktop/Dropbox/DBN2019/smartseq/tcr/circles")

# Read in tcr output from "tracer summarise" command from sbatchSummarize scripts
data_r1 <- read.table("../tracer_cluster/data/TracerOut_Mar15/filtered_TCRAB_summary/recombinants.txt", header=TRUE)
data_r2 <- read.table("../tracer_cluster/data/TracerOut_Oct19/filtered_TCRAB_summary/recombinants.txt", header=TRUE)
data_comb <- rbind(data_r1, data_r2) # 590 by 5

# Use sample name to create addition age and tissue variables
data2 <- separate(data_comb, cell_name, c("age_rep", "tissue", "other"), sep = "-", remove = F, extra = "drop")
df <- separate(data2, age_rep, c("age", "rep"), sep = -1, remove = F) # 590 by 10
df_prod <- filter(df, productive == "True") # 459  10

# Fix one label
df2 <- data.frame(lapply(df_prod, function(x) {gsub("BloodCD8", "Blood", x)}))

# Get counts of recombinant ID's by age and tissue AND MOUSE REPLICATE group:
dfc <- dplyr::count(df2, age, tissue, rep, locus, recombinant_id) # dim(dfc): 292 6
dfc <- dfc[order(dfc$n, decreasing=TRUE),]
# dfc$recomb_id_factor <- factor(dfc$recombinant_id, levels = rev(unique(dfc$recombinant_id)))
# dfc$n_factor <- factor(dfc$n)
dfc <- unite(dfc, sample, age, rep, sep = " ", remove = F)

#=======================================================================================
# Pie Charts

tcells <- dfc
tcells$tissue <- as.character(tcells$tissue)
tcells$locus <- as.character(tcells$locus)
tcells$tissue <- as.character(tcells$tissue)
old_tcells <- tcells %>% filter(age == "Old")

# All Old SVZ
d <-  old_tcells %>% filter(locus == "B", tissue == "SVZ")
large <- sum(d[d$n>2,]$n)
small <- sum(d[d$n==2,]$n)
single <- sum(d[d$n==1,]$n)
n <- large + small + single
plot_df <- melt(data.frame("large"=large/n, "small"=small/n, "single"=single/n))
colnames(plot_df) <- c("ClonotypeSize","Value")

p <- ggplot(plot_df, aes(x="", y=Value, fill=ClonotypeSize)) + geom_bar(width = 1, stat = "identity", color = "white")
p <- p + coord_polar("y", start=0)
p <- p + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
p <- p + theme_void()
p <- p + geom_text(aes(label=as.character(round(plot_df$Value,3))), position = position_stack(vjust = 0.5), size = 6)
ggsave("plots/Fig2e_Pie_Old_SVZ_clonotype.pdf", p, width=4, height=4)

# All Old Blood
d <-  old_tcells %>% filter(locus == "B", tissue == "Blood")
large <- sum(d[d$n>2,]$n)
small <- sum(d[d$n==2,]$n)
single <- sum(d[d$n==1,]$n)
n <- large + small + single
plot_df <- melt(data.frame("large"=large/n, "small"=small/n, "single"=single/n))
colnames(plot_df) <- c("ClonotypeSize","Value")

p <- ggplot(plot_df, aes(x="", y=Value, fill=ClonotypeSize)) + geom_bar(width = 1, stat = "identity", color = "white")
p <- p + coord_polar("y", start=0)
p <- p + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
p <- p + theme_void()
p <- p + geom_text(aes(label=as.character(round(plot_df$Value, 3))), position = position_stack(vjust = 0.5), size = 6)
ggsave("plots/Fig2e_Pie_Old_Blood_clonotype.pdf", p, width=4, height=4)

#=======================================================================================
# Venn Diagram Fig 2
d <-  old_tcells %>% filter(locus == "B", tissue == "SVZ")
svz_clones <- unique(d$recombinant_id)
length(svz_clones)# 76
d <-  old_tcells %>% filter(locus == "B", tissue == "Blood")
blood_clones <- unique(d$recombinant_id)
length(blood_clones) # 50
length(union(svz_clones, blood_clones)) # 123
length(intersect(svz_clones, blood_clones)) # 3

# Pie Fig 2
p <- euler(c(A = 73, B = 47, "A&B" = 3))
pdf("plots/Fig2f_Venn.pdf", width=4, height=4)
plot(p)
dev.off()

#=======================================================================================
# Venn Diagram Ext Fig 3

d <-  old_tcells %>% filter(locus == "B", rep == 1)
clones_1 <- unique(d$recombinant_id); length(clones_1) # 38
d <-  old_tcells %>% filter(locus == "B", rep == 2)
clones_2 <- unique(d$recombinant_id); length(clones_2) # 29
d <-  old_tcells %>% filter(locus == "B", rep == 3)
clones_3 <- unique(d$recombinant_id); length(clones_3) # 16
d <-  old_tcells %>% filter(locus == "B", rep == 4)
clones_4 <- unique(d$recombinant_id); length(clones_4) # 40

length(c(clones_1, clones_2, clones_3, clones_4)) # 123
length(union(c(clones_1, clones_2), c(clones_3, clones_4))) #123 --> no overlap

# no public clones
length(intersect(clones_1, clones_2)) # 0
length(intersect(clones_3, clones_4)) # 0
length(intersect(c(clones_1, clones_2), c(clones_3, clones_4))) # 0 --> no overlap confirmed

# Non-overlapping circles with proportional areas created in Illustrator. 

sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14

# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods  
# [8] base     

# other attached packages:
#  [1] eulerr_5.1.0          reshape2_1.4.3        bindrcpp_0.2.2       
#  [4] randomcoloR_1.1.0     RColorBrewer_1.1-2    ComplexHeatmap_1.17.1
#  [7] forcats_0.3.0         stringr_1.3.1         dplyr_0.7.6          
# [10] purrr_0.2.5           readr_1.1.1           tidyr_0.8.1          
# [13] tibble_2.0.1          ggplot2_3.1.0         tidyverse_1.2.1  


