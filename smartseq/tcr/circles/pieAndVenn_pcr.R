# Venn Charts from T Cell TCR-B  Clonotypes
# Nested PCR Data

rm(list=ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(eulerr)

# Modify path if needed. All sebsequent paths are relative. 
setwd("~/Desktop/Dropbox/DBN2019/smartseq/tcr/circles")

# Preprocess/clean/Exclude hippocampus 
tcells<-data.frame(read.csv("../data/TCR_sequencing_NestedPCR.csv", header=TRUE, skip=1))
tcells$X <- NULL
cytokines<-as.matrix(tcells[,23:31])
cytokine_sum<-apply(cytokines!="",1,sum)
tcells<-tcells[cytokine_sum>0,]
tcells<-tcells[tcells$Region!="Hipp",]
tcells <- tcells[tcells$B_seq != "",]

#Old TCells
old5_tcells<-tcells[as.vector(tcells$Mouse_ID)=="Old_5",]
old6_tcells<-tcells[as.vector(tcells$Mouse_ID)=="Old_6",]

old5_tcells_blood <- dplyr::filter(old5_tcells, Region=="Blood")
old6_tcells_blood <- dplyr::filter(old6_tcells, Region=="Blood")

old5_tcells_svz <- dplyr::filter(old5_tcells, Region=="SVZ")
old6_tcells_svz <- dplyr::filter(old6_tcells, Region=="SVZ")

# T Cell Clonality analysis 
# (Low yeilds from old mice 1 and 2)

# Old 3
length(old5_tcells$B_seq) # 161
length(unique(old5_tcells$B_seq)) # 69
length(old5_tcells_blood$B_seq) # 83
length(unique(old5_tcells_blood$B_seq)) # 36
length(old5_tcells_svz$B_seq) # 78
length(unique(old5_tcells_svz$B_seq)) # 40

# Old 4
length(old6_tcells$B_seq) # 148
length(unique(old6_tcells$B_seq)) # 86
length(old6_tcells_blood$B_seq) # 82
length(unique(old6_tcells_blood$B_seq)) # 60
length(old6_tcells_svz$B_seq) # 66
length(unique(old6_tcells_svz$B_seq)) # 30

# Ext Data 3f
tcells_blood_seq <- union(old5_tcells_blood$B_seq, old6_tcells_blood$B_seq)
length(tcells_blood_seq) # 96
tcells_svz_seq <- union(old5_tcells_svz$B_seq, old6_tcells_svz$B_seq)
length(tcells_svz_seq) # 70
length(intersect(tcells_blood_seq, tcells_svz_seq)) # 11 : overlap to subtract from 96 and 70.

p <- euler(c(A = 85, B = 69, "A&B" = 11))
pdf("plots/SupFig3f_Venn.pdf", width=4, height=4)
plot(p)
dev.off()

# Ext Data 3f
length(unique(old5_tcells$B_seq)) # 69
length(unique(old6_tcells$B_seq)) # 86
length(intersect(old6_tcells$B_seq, old5_tcells$B_seq)) # 0 / no overlap to subtract

# Non intersecting circles produced in Illustrator


















