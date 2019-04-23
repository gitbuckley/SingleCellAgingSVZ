# Colorful barplots scRNAseq
# Based off Ben Dulken's script.

rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(randomcoloR)
library(dplyr)

############################################################################################
# Combine and process TraCeR output files.

# Modify path as needed. All subsequent are relative. 
setwd("~/Desktop/Dropbox/DBN2019/smartseq/tcr/barplots")

# Read in tcr output from "tracer summarise" command from sbatchSummarize scripts
data <- read.table("../tracer_cluster/data/TracerOut_Mar15/filtered_TCRAB_summary/recombinants.txt", header=TRUE)
data_r2 <- read.table("../tracer_cluster/data/TracerOut_Oct19/filtered_TCRAB_summary/recombinants.txt", header=TRUE)

# Combine datasets
data_comb <- rbind(data, data_r2)

# Use sample name to create addition age and tissue variables
data2 <- separate(data_comb, cell_name, c("age_rep", "tissue", "other"), sep = "-", remove = F, extra = "drop")
df <- separate(data2, age_rep, c("age", "rep"), sep = -1, remove = F)
df_prod <- filter(df, productive == "True")

# Correct mislabeled tissue
df2 <- data.frame(lapply(df_prod, function(x) {
	gsub("BloodCD8", "Blood", x)}))

# Get counts of recombinant ID's by age and tissue AND MOUSE REPLICATE group:
dfc <- dplyr::count(df2, age, tissue, rep, locus, recombinant_id)
dfc <- dfc[order(dfc$n, decreasing=TRUE),]
# dfc$recomb_id_factor <- factor(dfc$recombinant_id, levels = rev(unique(dfc$recombinant_id)))
# dfc$n_factor <- factor(dfc$n)
dfc <- unite(dfc, sample, age, rep, sep = " ", remove = F)
write.table(dfc, file = "../data/TCRcounts_prod_all4.txt", sep = "\t", quote = F, row.names = F)

############################################################################################
#T Cell Clonality analysis using TCRB

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

tcells<-data.frame(read.table("../data/TCRcounts_prod_all4.txt",sep="\t",fill=T,header=T))
tcells<-tcells[tcells$locus=="B", ]

#Old TCells
old1_tcells<-tcells[as.vector(tcells$sample)=="Old 1",]
old2_tcells<-tcells[as.vector(tcells$sample)=="Old 2",]
old3_tcells<-tcells[as.vector(tcells$sample)=="Old 3",]
old4_tcells<-tcells[as.vector(tcells$sample)=="Old 4",]
young1_tcells<-tcells[as.vector(tcells$sample)=="Young 1",]
young2_tcells<-tcells[as.vector(tcells$sample)=="Young 2",]

#Clonality
#TCRB
for(j in 1:4){
  
  temp<-get(paste0("old",j,"_tcells"))
  tcells_seq<-temp[as.vector(temp$recombinant_id)!="",]
  tcells_seq<-tcells_seq[order(tcells_seq$recombinant_id),]
  abundance<-aggregate(tcells_seq$n,by=list(as.vector(tcells_seq$recombinant_id)),sum)
  abundance<-abundance[order(abundance[,2],decreasing=F),]
  
  regionvec<-c()
  for(i in 1:length(abundance[,1])){
    curr<-tcells_seq[grepl(abundance[i,1],as.vector(tcells_seq$recombinant_id)),]
    region<-as.vector(curr$tissue)
    regionvec<-c(regionvec,rep("orange",max(0,curr$n[grepl("SVZ",curr$tissue)])),
                 rep("red",max(0,curr$n[grepl("Blood",curr$tissue)])))
  }
  regionvec<-rev(regionvec)

  distinctCols <- distinctColorPalette(60)
  
  d<-data.frame(tcr=abundance[,1],prevalence=abundance[,2],id=rep(2,length(abundance[,1])))
  d2<-data.frame(id=rep(1,length(regionvec)),pos=rep(1,length(regionvec)),fac=regionvec)
  d$tcr<-factor(d$tcr,levels=abundance[,1],ordered=T)
  d$id<-factor(d$id,levels=c("1","2"),ordered=T)
  d2$id<-factor(d2$id,levels=c("1","2"),ordered=T)
  d2$fac<-factor(d2$fac,ordered=F)
  
  # Free axis
  p<-ggplot(d)
  p<-p+geom_bar(aes(x=d$id,y=d$prevalence,fill=d$tcr),stat="identity", width=1.4)
  p<-p+geom_bar(data=d2,aes(x=d2$id,y=d2$pos),  width=.5, fill=as.vector(d2$fac),stat="identity")
  p<-p+theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'),
    axis.ticks.length=unit(.25, "cm"))
  p <- p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  p <- p + theme(axis.title.x=element_blank(), axis.text.x = element_text(size = 24))
  p<-p+scale_fill_manual(values=distinctCols)
  p<-p + theme(legend.position="none") 
  p<-p+coord_flip()
  #p<-p+ylim(0,160)
  pdf(paste0("plots/Fig2D_OldMouse",j,"_TCRB_scRNAseq_", Sys.Date(), ".pdf"),width=10,height=3)
  print(p)
  dev.off()

  # Fixed axis
  p<-ggplot(d)
  p<-p+geom_bar(aes(x=d$id,y=d$prevalence,fill=d$tcr),stat="identity", width=1.4)
  p<-p+geom_bar(data=d2,aes(x=d2$id,y=d2$pos),  width=.5, fill=as.vector(d2$fac),stat="identity")
  p<-p+theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p+scale_fill_manual(values=distinctCols)
  p <- p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  p <- p + theme(axis.title.x=element_blank(), axis.text.x = element_text(size = 12))
  p<-p + theme(legend.position="none") 
  p<-p+coord_flip()
  p<-p+ylim(0,160)
  pdf(paste0("plots/SupFig3C_OldMouse",j,"_TCRB_scRNAseq_fixedaxis_", Sys.Date(), ".pdf"),width=8,height=3)
  print(p)
  dev.off()

}
