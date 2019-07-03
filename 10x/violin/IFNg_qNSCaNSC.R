rm(list=ls())
library(Seurat)
library(tidyverse)
library(dplyr)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/violin")

# Load Data
load(file = "../seurat/data/svz_celltypes_2019-01-31.rda")
svz__data <- as.matrix(svz@data)

#===================================================================================================
# Ifng response genes

load(file = "../data/mouse.gene.hallmark.kegg.reactome.rda")
genelist <- tolower(hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE)

# Get cells of interest column locations
cells_loc <- svz@meta.data$Celltype %in% c("aNSCs_NPCs", "Astrocytes_qNSCs") # 2463 are TRUE

# Subset expression data to ifn genes and cells of interest
rownames(svz__data) <- tolower(rownames(svz__data))
ifn_data <- svz__data[rownames(svz__data) %in% genelist, cells_loc ] # dim: 310 x 2463

# Ensure no rows have zero sum.
ifn_data <- ifn_data[rowSums(ifn_data)>0, ] # 252 2463

# Subset metadata to cells of interest
meta_data <- svz@meta.data[cells_loc,]

# High ifn cells
num_sd <- 1.65
ifn_data_sum <- apply(ifn_data,2,mean)
ifn_high <- colnames(ifn_data)[ifn_data_sum > (mean(ifn_data_sum) + num_sd*sd(ifn_data_sum))]
ifn_high_loc <- colnames(ifn_data) %in% ifn_high
#age_ifn_colors[ifn_high_loc] <- "#000000"


# Convert row names to column
data <- add_rownames(meta_data, "Cell")
data$Response <- ifn_data_sum
data$High <- ifn_high_loc
data$Age <- relevel(as.factor(data$Age), "y")
data$Age_High <- as.character(data$Age)
data$Age_High[data$High] <- "High"
data$Age_High <- as.factor(data$Age_High)
head(data)
# Add all ifn genes 
ifn_data_t <- as.data.frame(t(ifn_data))
ifn_data_t <- add_rownames(ifn_data_t, "Cell")
# Combine
data <- merge(data, ifn_data_t) # dim: 2463  267

# Signature Violin
p <- ggplot(data)
p <- p + geom_violin(aes(x=Age, y=Response, fill=Age), scale="width", draw_quantiles = c(0.5), alpha = 0.6)
p <- p + geom_jitter(aes(x=Age, y=Response, color=Age_High), size = 1.3, position = position_jitter(width = .25,height=0),alpha=0.7)
p <- p + facet_grid(cols=vars(Replicate))
p<-p+ theme_classic() + theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
p<-p + theme(legend.position="none")
p<-p+theme(axis.text.y=element_text(size=15))
p<-p+theme(axis.text.x=element_text(size=20), axis.ticks.x = element_blank())
p<-p+theme(axis.title=element_text(size=20))
p<-p+theme(plot.title=element_text(size=20))
p<-p+theme(axis.title.y=element_text(vjust=1))
p<-p+theme(axis.title.x=element_text(vjust=-0.10))
p<-p+labs(y="IFN Response Gene Expression")
p<-p+theme(axis.title.x=element_blank())
p<-p+labs(title=NULL)
p<-p+scale_fill_manual(values=c("#40BBEC","#EF4136"))
p <- p + scale_color_manual(values=c("black", "#EF4136", "#40BBEC"))
ggsave(paste("plots/SupFig5c_Violin_aNSC_qNSC_sd_",num_sd,".pdf",sep=""),p, height=5,width=6)

#===================================================================================================
# PROPORTIONS
y1 <- data %>% filter(orig.ident=="y1")
y1_high <- sum(y1$High) / dim(y1)[1]
o1 <- data %>% filter(orig.ident=="o1")
o1_high <- sum(o1$High) / dim(o1)[1]

y2 <- data %>% filter(orig.ident=="y2")
y2_high <- sum(y2$High) / dim(y2)[1]
o2 <- data %>% filter(orig.ident=="o2")
o2_high <- sum(o2$High) / dim(o2)[1]

y3 <- data %>% filter(orig.ident=="y3")
y3_high <- sum(y3$High) / dim(y3)[1]
o3 <- data %>% filter(orig.ident=="o3")
o3_high <- sum(o3$High) / dim(o3)[1]

prop_df <- data.frame("Age" = c("y","o","y","o","y","o"),
					  "Bst2Hi" = c(y1_high, o1_high, y2_high, o2_high, y3_high, o3_high))
prop_df

#   Age     Bst2Hi
# 1   y 0.00000000
# 2   o 0.06172840
# 3   y 0.01547988
# 4   o 0.11111111
# 5   y 0.07225434
# 6   o 0.12941176
