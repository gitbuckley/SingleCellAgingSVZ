# PCA aNSC IFN Response Old 

rm(list=ls())
library(Seurat)
library(tidyverse)
library(dplyr)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/pca")

# Load Data
load(file = "../seurat/data/svz_celltypes_2019-01-31.rda")
svz__data <- as.matrix(svz@data)
load(file = "../data/mouse.gene.hallmark.kegg.reactome.rda")
genelist <- tolower(hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE)

#============================================================================
# Old cells, all  3 replicates
cells_loc <- svz@meta.data$Celltype %in% c("aNSCs_NPCs", "Astrocytes_qNSCs") & svz@meta.data$Age == "o"

# Subset to ifn genes and cells of interest
rownames(svz__data) <- tolower(rownames(svz__data))
ifn_data <- svz__data[rownames(svz__data) %in% genelist, cells_loc ]

# Ensure no rows have zero sum. PCA won't work with zero sum columns.
ifn_data <- ifn_data[rowSums(ifn_data)>0, ] 

PCA_int <- prcomp(t(ifn_data), scale = T, center = T, retx=T)
PCA_results <- PCA_int$x; dim(PCA_results)
summa <- summary(PCA_int); summa

# Obtain cell types specific labels
old1_counts <- sum(grepl("o1",colnames(ifn_data))) # 162
old2_counts <- sum(grepl("o2",colnames(ifn_data))) # 315
old3_counts <- sum(grepl("o3",colnames(ifn_data))) # 85

age_colors <- c(rep("#EF4136", old1_counts),
                rep("#EF4136", old2_counts),
                rep("#EF4136", old3_counts))

age_ifn_colors <- age_colors

# High ifn cells
ifn_data_sum <- apply(ifn_data,2,mean)
ifn_high <- colnames(ifn_data)[ifn_data_sum > (mean(ifn_data_sum) + 1.65*sd(ifn_data_sum))]
ifn_high_loc <- colnames(ifn_data) %in% ifn_high
age_ifn_colors[ifn_high_loc] <- "black"
age_ifn_fill <- age_ifn_colors
age_ifn_fill[ifn_high_loc] <- "darkred"

# Plotting dataframe
# Flip X coordinates (does not impact interpretation)
data_1 <- data.frame("PC1"=PCA_results[,1] * -1, "PC2"=PCA_results[,2], "Age_ifn_fill"=age_ifn_fill,
 "Age_color"=as.factor(age_colors), "Age_ifn_colors"=as.factor(age_ifn_colors), "Response" = ifn_data_sum)
p <- ggplot(data_1)
p<- p+geom_point(aes(x=PC1, y=PC2), fill=data_1$Age_ifn_fill, color=data_1$Age_ifn_colors, size=1.75, alpha=.4, shape=21, stroke=.6)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2  (", round(summa$importance[2,2], digits = 2)*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[2,1],digits = 2)*100,"% of Variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=8))
p<-p+theme(axis.title.x=element_text(size=10))
p<-p+theme(axis.text.y=element_text(size=8))
p<-p+theme(axis.title.y=element_text(size=10))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))
pdf(paste("plots/Fig4a_PCA_Rep123_sd1.65_aNSC_qNSC_", Sys.Date(),".pdf",sep=""),height=3.5,width=3)
print(p)
dev.off()

#============================================================================
# Old cells, Replicate 1

# Get cells of interest column locations
cells_loc <- svz@meta.data$Celltype %in% c("aNSCs_NPCs", "Astrocytes_qNSCs") & svz@meta.data$Age == "o" & svz@meta.data$Replicate == 1

# Subset to ifn genes and cells of interest
ifn_data <- svz__data[rownames(svz__data) %in% genelist, cells_loc ] # dim: 310 162

# Filter unexpressed genes
ifn_data <- ifn_data[rowSums(ifn_data)>0, ] # 184 162

PCA_int <- prcomp(t(ifn_data), scale = T, center = T, retx=T)
PCA_results <- PCA_int$x; dim(PCA_results)
summa <- summary(PCA_int); summa


# Obtain cell types specific labels
old1_counts <- sum(grepl("o1",colnames(ifn_data)))
#old2_counts <- sum(grepl("o2",colnames(ifn_data)))
#old3_counts <- sum(grepl("o3",colnames(ifn_data)))

age_colors <- c(rep("#EF4136", old1_counts))

age_ifn_colors <- age_colors

# High ifn cells
ifn_data_sum <- apply(ifn_data,2,mean)
ifn_high <- colnames(ifn_data)[ifn_data_sum > (mean(ifn_data_sum) + 1.65*sd(ifn_data_sum))]
ifn_high_loc <- colnames(ifn_data) %in% ifn_high
age_ifn_colors[ifn_high_loc] <- "#000000"

# Plotting dataframe
data_1 <- data.frame("PC1"=PCA_results[,1], "PC2"=PCA_results[,2], "Age_color"=as.factor(age_colors), "Age_ifn_colors"=as.factor(age_ifn_colors), "Response" = ifn_data_sum)
p <- ggplot(data_1)
p<- p+geom_point(aes(x=PC1, y=PC2), fill=data_1$Age_color, color=data_1$Age_ifn_colors, size=1.5, alpha=0.4, shape=21, stroke=1.2)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2  (", round(summa$importance[2,2], digits = 2)*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[2,1],digits = 2)*100,"% of Variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=14))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=14))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))

pdf(paste("plots/SupFig5a_Rep1_sd1.65_old_aNSC_qNSC_",Sys.Date(),".pdf",sep=""),height=5,width=5)
print(p)
dev.off()

#========================
# Old cells, Replicate 2

# Get cells of interest column locations
cells_loc <- svz@meta.data$Celltype %in% c("aNSCs_NPCs", "Astrocytes_qNSCs") & svz@meta.data$Age == "o" & svz@meta.data$Replicate == 2

# Subset to ifn genes and cells of interest
ifn_data <- svz__data[rownames(svz__data) %in% genelist, cells_loc ] # dim: 310 315

# Filter unexpressed genes
ifn_data <- ifn_data[rowSums(ifn_data)>0, ] # 221 315

PCA_int <- prcomp(t(ifn_data), scale = T, center = T, retx=T)
PCA_results <- PCA_int$x; dim(PCA_results)
summa <- summary(PCA_int); summa

# Obtain cell types specific labels
#old1_counts <- sum(grepl("o1",colnames(ifn_data)))
old2_counts <- sum(grepl("o2",colnames(ifn_data)))
#old3_counts <- sum(grepl("o3",colnames(ifn_data)))

age_colors <- c(rep("#EF4136", old2_counts))

age_ifn_colors <- age_colors

# High ifn cells
ifn_data_sum <- apply(ifn_data,2,mean)
ifn_high <- colnames(ifn_data)[ifn_data_sum > (mean(ifn_data_sum) + 1.65*sd(ifn_data_sum))]
ifn_high_loc <- colnames(ifn_data) %in% ifn_high
age_ifn_colors[ifn_high_loc] <- "#000000"

# Plotting dataframe
data_1 <- data.frame("PC1"=PCA_results[,1], "PC2"=PCA_results[,2], "Age_color"=as.factor(age_colors), "Age_ifn_colors"=as.factor(age_ifn_colors), "Response" = ifn_data_sum)
p <- ggplot(data_1)
p<- p+geom_point(aes(x=PC1, y=PC2), fill=data_1$Age_color, color=data_1$Age_ifn_colors, size=1.5, alpha=0.4, shape=21, stroke=1.2)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2  (", round(summa$importance[2,2], digits = 2)*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[2,1],digits = 2)*100,"% of Variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=14))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=14))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))

pdf(paste("plots/SupFig5a_Rep2_sd1.65_old_aNSC_qNSC_",Sys.Date(),".pdf",sep=""),height=5,width=5)
print(p)
dev.off()

#========================
# Old cells, Replicate 3

# Get cells of interest column locations
cells_loc <- svz@meta.data$Celltype %in% c("aNSCs_NPCs", "Astrocytes_qNSCs") & svz@meta.data$Age == "o" & svz@meta.data$Replicate == 3

# Subset to ifn genes and cells of interest
ifn_data <- svz__data[rownames(svz__data) %in% genelist, cells_loc ] # dim: 310  85

# Filter unexpressed genes
ifn_data <- ifn_data[rowSums(ifn_data)>0, ] # 182  85

PCA_int <- prcomp(t(ifn_data), scale = T, center = T, retx=T)
PCA_results <- PCA_int$x; dim(PCA_results)
summa <- summary(PCA_int); summa


# Obtain cell types specific labels
#old1_counts <- sum(grepl("o1",colnames(ifn_data)))
#old2_counts <- sum(grepl("o2",colnames(ifn_data)))
old3_counts <- sum(grepl("o3",colnames(ifn_data)))

age_colors <- c(rep("#EF4136", old3_counts))

age_ifn_colors <- age_colors

# High ifn cells
ifn_data_sum <- apply(ifn_data,2,mean)
ifn_high <- colnames(ifn_data)[ifn_data_sum > (mean(ifn_data_sum) + 1.65*sd(ifn_data_sum))]
ifn_high_loc <- colnames(ifn_data) %in% ifn_high
age_ifn_colors[ifn_high_loc] <- "#000000"

# Plotting dataframe
data_1 <- data.frame("PC1"=PCA_results[,1], "PC2"=PCA_results[,2], "Age_color"=as.factor(age_colors), "Age_ifn_colors"=as.factor(age_ifn_colors), "Response" = ifn_data_sum)
p <- ggplot(data_1)
p<- p+geom_point(aes(x=PC1, y=PC2), fill=data_1$Age_color, color=data_1$Age_ifn_colors, size=1.5, alpha=0.4, shape=21, stroke=1.2)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2  (", round(summa$importance[2,2], digits = 2)*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[2,1],digits = 2)*100,"% of Variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=14))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=14))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))

pdf(paste("plots/SupFig5a_Rep3_Refactored_sd1.65_old_aNSC_qNSC_",Sys.Date(),".pdf",sep=""),height=5,width=5)
print(p)
dev.off()

#============================================================================
#============================================================================
# FLUIDIGM DATA
# PCA aNSC IFN Response Old Only

data <- read.csv("../../C1/FluidigmC1Expression.csv", header=TRUE)
colnames(data)[1] <- "gene"

#============================================================================
# Get old cells only
olds <- grepl("O.*", colnames(data))
olds[1] <- TRUE # keep gene names
oldCells <- data[,olds]

# Subset to ifn genes and cells of interest
oldCells$gene <- tolower(oldCells$gene)
ifn_data <- oldCells[oldCells$gene %in% genelist, ] # dim: 96 138
ifn_data <- as.matrix(ifn_data[,2:ncol(ifn_data)]) # remove gene column

# Ensure no rows have zero sum.
ifn_data <- ifn_data[rowSums(ifn_data)>0, ] # 96 137

PCA_int <- prcomp(t(ifn_data), scale = T, center = T, retx=T)
PCA_results <- PCA_int$x; dim(PCA_results)
summa <- summary(PCA_int); summa

# Obtain cell types specific labels
old_counts <- sum(grepl("O",colnames(ifn_data)))
age_colors <- c(rep("#EF4136", old_counts))
age_ifn_colors <- age_colors

# High ifn cells
ifn_data_sum <- apply(ifn_data,2,mean)
ifn_high <- colnames(ifn_data)[ifn_data_sum > (mean(ifn_data_sum) + 1.65*sd(ifn_data_sum))]
ifn_high_loc <- colnames(ifn_data) %in% ifn_high
age_ifn_colors[ifn_high_loc] <- "black"
#age_ifn_fill <- age_ifn_colors
#age_ifn_fill[ifn_high_loc] <- "#000000"

# Plotting dataframe
data_1 <- data.frame("PC1"=PCA_results[,1], "PC2"=PCA_results[,2], "Age_color"=as.factor(age_colors), "Age_ifn_colors"=as.factor(age_ifn_colors), "Response" = ifn_data_sum)
p <- ggplot(data_1)
p<- p+geom_point(aes(x=PC1, y=PC2), fill=data_1$Age_color, color=data_1$Age_ifn_colors, size=1.5, alpha=0.4, shape=21, stroke=1.2)
p<-p+theme_classic()
p<- p+ labs(y = paste("PC2  (", round(summa$importance[2,2], digits = 2)*100,"% of Variance)", sep = ""), x =paste("PC1  (",round(summa$importance[2,1],digits = 2)*100,"% of Variance)", sep = ""))
p<-p+theme(axis.text.x=element_text(size=20))
p<-p+theme(axis.title.x=element_text(size=14))
p<-p+theme(axis.text.y=element_text(size=20))
p<-p+theme(axis.title.y=element_text(size=14))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))
pdf(paste("plots/SupFig5a_Fluidigm_sd1.65_aNSC_qNSC_", Sys.Date(),".pdf",sep=""),height=5,width=5)
print(p)
dev.off()

