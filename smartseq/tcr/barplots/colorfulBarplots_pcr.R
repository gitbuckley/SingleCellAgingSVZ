# T Cell Clonality analysis - Mice 3 and 4 only for clonality analysis.
# (Low yeilds from old mice 1 and 2).
# Based off Ben Dulken's script.

rm(list=ls())
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(randomcoloR)

# Modify path as needed. All subsequent are relative.
setwd("~/Desktop/Dropbox/DBN2019/smartseq/tcr/barplots")
 
tcells<-data.frame(read.csv("../data/TCR_sequencing_NestedPCR.csv", header=TRUE, skip=1))
cytokines<-as.matrix(tcells[,24:32])
cytokine_sum<-apply(cytokines!="",1,sum)
tcells<-tcells[cytokine_sum>0,]
tcells<-tcells[tcells$Region!="Hipp",]

#Old TCells
old5_tcells<-tcells[as.vector(tcells$Mouse_ID)=="Old_5",]
old6_tcells<-tcells[as.vector(tcells$Mouse_ID)=="Old_6",]

young1_tcells<-tcells[as.vector(tcells$Mouse_ID)=="Young_1",]
young2_tcells<-tcells[as.vector(tcells$Mouse_ID)=="Young_2",]

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#===========================================================================================================================
#Clonality
#TCRB

for(j in c(5,6)){
  
  temp<-get(paste0("old",j,"_tcells"))
  tcells_seq<-temp[as.vector(temp$B_seq)!="",]
  tcells_seq<-tcells_seq[order(tcells_seq$B_seq),]
  aseq_rep<-as.vector(unique(tcells_seq$A_seq))
  a_col_vec<-col_vector[sample(1:length(col_vector),length(aseq_rep))]
  a_col_vec[aseq_rep==""]<-"#000000"
  abundance<-aggregate(c(rep(1,length(tcells_seq$B_seq))),by=list(as.vector(tcells_seq$B_seq)),sum)
  abundance<-abundance[order(abundance[,2],decreasing=F),]
  
  regionvec<-c()
  avec<-c()
  for(i in 1:length(abundance[,1])){
    curr<-tcells_seq[grepl(abundance[i,1],as.vector(tcells_seq$B_seq)),]
    region<-as.vector(curr$Region)
    avec<-c(avec,a_col_vec[match(as.vector(curr$A_seq),aseq_rep)])
    regionvec<-c(regionvec,rep("orange",sum(grepl("SVZ",region))),
                 rep("#a4c639",sum(grepl("Hipp",region))),
                 rep("red",sum(grepl("Blood",region))))
  }
  regionvec<-rev(regionvec)
  avec<-rev(avec)
  
  d<-data.frame(tcr=abundance[,1],prevalence=abundance[,2],id=rep(3,length(abundance[,1])))
  d2<-data.frame(id=rep(1,length(regionvec)),pos=rep(1,length(regionvec)),fac=regionvec)
  d3<-data.frame(id=rep(2,length(avec)),pos=rep(1,length(avec)),fac=avec)
  d$tcr<-factor(d$tcr,levels=abundance[,1],ordered=T)
  d$id<-factor(d$id,levels=c("1","2","3"),ordered=T)
  d2$id<-factor(d2$id,levels=c("1","2","3"),ordered=T)
  d2$fac<-factor(d2$fac,ordered=F)
  d3$id<-factor(d3$id,levels=c("1","2","3"),ordered=T)
  d3$fac<-factor(d3$fac,ordered=F)

  distinctCols <- distinctColorPalette(100)
  
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
  pdf(paste0("plots/Fig2E_PCR_OldMouse",j,"_TCRB_AddedCells_", Sys.Date(), ".pdf"),width=10,height=3)
  print(p)
  dev.off()

  p<-ggplot(d)
  p<-p+geom_bar(aes(x=d$id,y=d$prevalence,fill=d$tcr),stat="identity", width = 1.4)
  p<-p+geom_bar(data=d2,aes(x=d2$id,y=d2$pos),fill=as.vector(d2$fac),stat="identity", width = .5)
  p<-p+theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  p<-p+scale_fill_manual(values=distinctCols)
  p <- p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  p <- p + theme(axis.title.x=element_blank(), axis.text.x = element_text(size = 12))
  p<-p + theme(legend.position="none") 
  p<-p+coord_flip()
  p<-p+ylim(0,160)
  pdf(paste0("plots/SupFig3D_PCR_OldMouse",j,"_TCRB_AddedCells_fixedaxis_", Sys.Date(), ".pdf"),width=8,height=3)
  print(p)
  dev.off()
}
