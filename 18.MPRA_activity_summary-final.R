#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

##### install MPRAnalyze
## through github
# install_github("YosefLab/MPRAnalyze")
## or biocManager
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MPRAnalyze")

library(ggplot2)
library(splitstackshape)
library(grid)
library(gridExtra)
library(dplyr)
library(ggpmisc)
library(ggbeeswarm)
library(MPRAnalyze)

#################################### Show MER11 as example for detail seq analysis
Summary_Table2 = read.csv("input/Summary_Table2_2023_1_5.subfamilyInfo.csv")

##### designed library
Annotation = data.frame(read.delim("input/MPRA_combined_final_2022_8_25.tsv",sep="",header=T))

##### load the retrieved inserts by association analysis
inserts = read.delim("input/TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle.out.gz",header=F,sep="")
colnames(inserts)[1:2] = c("unique_MPRA","bc_asso")
Annotation = merge(Annotation,inserts[,c("unique_MPRA","bc_asso")],by="unique_MPRA",all.x=T)
Annotation$bc_asso_group = ifelse(!is.na(Annotation$bc_asso) & Annotation$bc_asso>=10,"bc_asso>=10","bc_asso<10")

color_Group = c("Consensus" = "#4d4d4d","Instance"="#4d4d4d","Positive"="#33a02c","Negative"="#1f78b4")
alpha_Group = c("Consensus" = 0.3,"Instance"= 0.3,"Positive"= 1,"Negative"= 1)
bc_min = 10

Date = "2023_10_17"
### read count tables
######## iPSC
### each replicate
iPSC_1 = read.csv("input/iPSC_1_counts.tsv",sep="\t")
iPSC_1 = iPSC_1[iPSC_1$n_obs_bc>=bc_min,]
colnames(iPSC_1) = paste(colnames(iPSC_1),".rep1",sep="")
iPSC_2 = read.csv("input/iPSC_2_counts.tsv",sep="\t")
iPSC_2 = iPSC_2[iPSC_2$n_obs_bc>=bc_min,]
colnames(iPSC_2) = paste(colnames(iPSC_2),".rep2",sep="")
iPSC_3 = read.csv("input/iPSC_3_counts.tsv",sep="\t")
iPSC_3 = iPSC_3[iPSC_3$n_obs_bc>=bc_min,]
colnames(iPSC_3) = paste(colnames(iPSC_3),".rep3",sep="")
iPSC =merge(iPSC_1,iPSC_2,by.x="name.rep1",by.y="name.rep2",all=T)
iPSC =merge(iPSC,iPSC_3,by.x="name.rep1",by.y="name.rep3",all=T)
colnames(iPSC)[1] = "unique_MPRA"
iPSC$replicates = ifelse(!is.na(iPSC$dna_count.rep1) & !is.na(iPSC$dna_count.rep2) & !is.na(iPSC$dna_count.rep3),3,NA)
iPSC$replicates = ifelse((!is.na(iPSC$dna_count.rep1) & !is.na(iPSC$dna_count.rep2) & is.na(iPSC$dna_count.rep3)) |
                           (!is.na(iPSC$dna_count.rep1) & is.na(iPSC$dna_count.rep2) & !is.na(iPSC$dna_count.rep3)) |
                           (is.na(iPSC$dna_count.rep1) & !is.na(iPSC$dna_count.rep2) & !is.na(iPSC$dna_count.rep3)),2,iPSC$replicates)
iPSC$replicates = ifelse((!is.na(iPSC$dna_count.rep1) & is.na(iPSC$dna_count.rep2) & is.na(iPSC$dna_count.rep3)) |
                           (is.na(iPSC$dna_count.rep1) & is.na(iPSC$dna_count.rep2) & !is.na(iPSC$dna_count.rep3)) |
                           (is.na(iPSC$dna_count.rep1) & !is.na(iPSC$dna_count.rep2) & is.na(iPSC$dna_count.rep3)),1,iPSC$replicates)
iPSC$ratio.ave = apply(iPSC[,c("ratio.rep1","ratio.rep2","ratio.rep3")],1,mean,na.rm=T)
iPSC$bc_group = ifelse(!is.na(iPSC$n_obs_bc.rep1) & iPSC$n_obs_bc.rep1>=10 & !is.na(iPSC$n_obs_bc.rep2 & iPSC$n_obs_bc.rep2>=10 & !is.na(iPSC$n_obs_bc.rep3) & iPSC$n_obs_bc.rep3>=10) |
                         (!is.na(iPSC$n_obs_bc.rep1) & iPSC$n_obs_bc.rep1>=10 & !is.na(iPSC$n_obs_bc.rep2) & iPSC$n_obs_bc.rep2>=10) |
                         (!is.na(iPSC$n_obs_bc.rep1) & iPSC$n_obs_bc.rep1>=10 & !is.na(iPSC$n_obs_bc.rep3) & iPSC$n_obs_bc.rep3>=10) |
                         (!is.na(iPSC$n_obs_bc.rep3) & iPSC$n_obs_bc.rep3>=10 & !is.na(iPSC$n_obs_bc.rep2) & iPSC$n_obs_bc.rep2>=10),
                       "obs_bc>=10","bc<10") 

######## NPC
### each replicate
NPC_1 = read.csv("input/NPC_1_counts.tsv",sep="\t")
NPC_1 = NPC_1[NPC_1$n_obs_bc>=bc_min,]
colnames(NPC_1) = paste(colnames(NPC_1),".rep1",sep="")
NPC_2 = read.csv("input/NPC_2_counts.tsv",sep="\t")
NPC_2 = NPC_2[NPC_2$n_obs_bc>=bc_min,]
colnames(NPC_2) = paste(colnames(NPC_2),".rep2",sep="")
NPC_3 = read.csv("input/NPC_3_counts.tsv",sep="\t")
NPC_3 = NPC_3[NPC_3$n_obs_bc>=bc_min,]
colnames(NPC_3) = paste(colnames(NPC_3),".rep3",sep="")
NPC =merge(NPC_1,NPC_2,by.x="name.rep1",by.y="name.rep2",all=T)
NPC =merge(NPC,NPC_3,by.x="name.rep1",by.y="name.rep3",all=T)
colnames(NPC)[1] = "unique_MPRA"

NPC$replicates = ifelse(!is.na(NPC$dna_count.rep1) & !is.na(NPC$dna_count.rep2) & !is.na(NPC$dna_count.rep3),3,NA)
NPC$replicates = ifelse((!is.na(NPC$dna_count.rep1) & !is.na(NPC$dna_count.rep2) & is.na(NPC$dna_count.rep3)) |
                          (!is.na(NPC$dna_count.rep1) & is.na(NPC$dna_count.rep2) & !is.na(NPC$dna_count.rep3)) |
                          (is.na(NPC$dna_count.rep1) & !is.na(NPC$dna_count.rep2) & !is.na(NPC$dna_count.rep3)),2,NPC$replicates)
NPC$replicates = ifelse((!is.na(NPC$dna_count.rep1) & is.na(NPC$dna_count.rep2) & is.na(NPC$dna_count.rep3)) |
                          (is.na(NPC$dna_count.rep1) & is.na(NPC$dna_count.rep2) & !is.na(NPC$dna_count.rep3)) |
                          (is.na(NPC$dna_count.rep1) & !is.na(NPC$dna_count.rep2) & is.na(NPC$dna_count.rep3)),1,NPC$replicates)
NPC$ratio.ave = apply(NPC[,c("ratio.rep1","ratio.rep2","ratio.rep3")],1,mean,na.rm=T)
NPC$bc_group = ifelse(!is.na(NPC$n_obs_bc.rep1) & NPC$n_obs_bc.rep1>=10 & !is.na(NPC$n_obs_bc.rep2 & NPC$n_obs_bc.rep2>=10 & !is.na(NPC$n_obs_bc.rep3) & NPC$n_obs_bc.rep3>=10) |
                        (!is.na(NPC$n_obs_bc.rep1) & NPC$n_obs_bc.rep1>=10 & !is.na(NPC$n_obs_bc.rep2) & NPC$n_obs_bc.rep2>=10) |
                        (!is.na(NPC$n_obs_bc.rep1) & NPC$n_obs_bc.rep1>=10 & !is.na(NPC$n_obs_bc.rep3) & NPC$n_obs_bc.rep3>=10) |
                        (!is.na(NPC$n_obs_bc.rep3) & NPC$n_obs_bc.rep3>=10 & !is.na(NPC$n_obs_bc.rep2) & NPC$n_obs_bc.rep2>=10),
                      "obs_bc>=10","bc<10") 
rm(iPSC_1,iPSC_2,iPSC_3)
rm(NPC_1,NPC_2,NPC_3)
Annotation_iPSC = merge(Annotation,iPSC,by="unique_MPRA",all.x=T)
Annotation_NPC = merge(Annotation,NPC,by="unique_MPRA",all.x=T)

Cell = "iPSC"
#################### Summary plot
for (Cell in c("iPSC","NPC")){
  if (Cell == "iPSC"){
    Annotation_plot = Annotation_iPSC
  } else {
    Annotation_plot = Annotation_NPC
  }
  Annotation_plot$bc_group = ifelse(is.na(Annotation_plot$bc_group),"obs_bc_no",Annotation_plot$bc_group)
  Annotation_plot$Group_bc = paste(Annotation_plot$bc_asso_group,Annotation_plot$bc_group)
  
  ###### plot 1: how many inserts were retrieved
  Annotation_plot$Family = ifelse(Annotation_plot$Group=="Consensus","",Annotation_plot$Family)
  Annotation_plot$Species = ifelse(Annotation_plot$Group=="Consensus","",Annotation_plot$Species)
  Annotation_plot$Frame = ifelse(Annotation_plot$Group=="Consensus","",Annotation_plot$Frame)
  Annotation_plot_sum = data.frame(Annotation_plot[Annotation_plot$Is_included=="Included_MPRA",] %>% group_by(Group,Family,Species,Frame,Group_bc) %>% dplyr::count())
  Annotation_plot_sum$ID = ifelse(Annotation_plot_sum$Group=="Consensus","Consensus",paste(Annotation_plot_sum$Species,Annotation_plot_sum$Family,Annotation_plot_sum$Frame,sep="-"))
  Annotation_plot_sum$ID = ifelse(Annotation_plot_sum$Group=="Negative" | Annotation_plot_sum$Group=="Positive",Annotation_plot_sum$Group,Annotation_plot_sum$ID)
  Annotation_plot_sum1 = Annotation_plot_sum[Annotation_plot_sum$Group == "Instance",]
  Annotation_plot_sum1 = Annotation_plot_sum1[order(paste(Annotation_plot_sum1$Family,Annotation_plot_sum1$Frame)),]
  Annotation_plot_sum2 = Annotation_plot_sum[Annotation_plot_sum$Group != "Instance",]
  Annotation_plot_sum = rbind(Annotation_plot_sum1,Annotation_plot_sum2)
  list_IDs = unique(Annotation_plot_sum$ID)
  
  Annotation_plot_sum_tmp = data.frame(Annotation_plot[Annotation_plot$Is_included=="Included_MPRA",] %>% group_by(Group,Family,Species,Frame) %>% dplyr::count())
  Annotation_plot_sum_tmp$ID = ifelse(Annotation_plot_sum_tmp$Group=="Consensus","Consensus",paste(Annotation_plot_sum_tmp$Species,Annotation_plot_sum_tmp$Family,Annotation_plot_sum_tmp$Frame,sep="-"))
  Annotation_plot_sum_tmp$ID = ifelse(Annotation_plot_sum_tmp$Group=="Negative" | Annotation_plot_sum_tmp$Group=="Positive",Annotation_plot_sum_tmp$Group,Annotation_plot_sum_tmp$ID)
  
  Annotation_plot_sum = merge(Annotation_plot_sum,Annotation_plot_sum_tmp[,c("ID","n")],by="ID",all.x=T)
  Annotation_plot_sum$freq = Annotation_plot_sum$n.x/Annotation_plot_sum$n.y
  Annotation_plot_sum$ID = factor(Annotation_plot_sum$ID,levels=list_IDs)
  
  Annotation_plot_label = Annotation_plot_sum[!duplicated(Annotation_plot_sum$ID),]
  Annotation_plot_label$ID = factor(Annotation_plot_label$ID,levels=list_IDs)
  
  # p = ggplot(data=Annotation_plot_sum, aes(x=ID, y= freq,group=Group_bc)) +
  #   geom_col(aes(fill=Group_bc))+
  #   scale_fill_brewer(palette = "GnBu")+
  #   ylim(0,1.2)+
  #   geom_text(data = Annotation_plot_label,
  #             aes(x = ID, group=1, y = 1 + 0.08,
  #                 label = n.y),
  #             color="black", position=position_dodge(.9), angle=75,hjust=.5)+
  #   ylab("Proportion of inserts")+
  #   xlab("Frame")+
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = rel(1)),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(),
  #     axis.line = element_line(colour = "black"),
  #     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.text=element_text(colour="black",size=rel(1),angle = 0),
  #     axis.text.x=element_text(colour="black",vjust=.5,hjust = 0.95,angle = 90),
  #     axis.title=element_text(colour="black",size=rel(1)),
  #     legend.key = element_rect(colour = "transparent", fill = "white"),
  #     legend.position="right",
  #     #legend.position = "none",
  #     legend.background = element_blank(),
  #     legend.text=element_text(size=rel(1)))
  # 
  # pdf(paste("Step5-1-",Cell,"_",Date,".pdf",sep=""),     
  #     width = 12,        # 5 x 300 pixels
  #     height = 5,
  #     pointsize = 10)
  # grid.draw(p)
  # dev.off()
  
  ########################## violin plot
  shape_species = c("Consensus" = 8,"hg19"=17,"macFas5"=15,"panTro4"=16,"Control"=16)
  color_species = c("Consensus" = "#6a3d9a","hg19"="#e31a1c","macFas5"="#1f78b4","panTro4"="#33a02c","Control"="grey44")
  alpha_species = c("Consensus" = 1,"hg19"=0.2,"macFas5"=0.2,"panTro4"=0.2,"Control" = 0.2)
  size_species = c("Consensus" = 1.5,"hg19"=1,"macFas5"=1,"panTro4"=1,"Control" = 1)
  
  if (Cell == "iPSC"){
    Annotation_plot = Annotation_iPSC
  } else {
    Annotation_plot = Annotation_NPC
  }
  Annotation_plot$bc_group = ifelse(is.na(Annotation_plot$bc_group),"obs_bc_no",Annotation_plot$bc_group)
  Annotation_plot$Group_bc = paste(Annotation_plot$bc_asso_group,Annotation_plot$bc_group)
  
  Annotation_plot2 = Annotation_plot[Annotation_plot$Is_included=="Included_MPRA",]
  Annotation_plot2$Species = ifelse(Annotation_plot2$Group == "Consensus","Consensus",Annotation_plot2$Species)
  Annotation_plot2$Species = ifelse(Annotation_plot2$Group == "Negative" | Annotation_plot2$Group == "Positive","Control",Annotation_plot2$Species)
  Annotation_plot2$ID2 = ifelse(Annotation_plot2$Group == "Instance" | Annotation_plot2$Group=="Consensus",paste(Annotation_plot2$Family,Annotation_plot2$Frame,sep="-"),NA)
  Annotation_plot2$ID2 = ifelse(Annotation_plot2$Group=="Negative" | Annotation_plot2$Group=="Positive",Annotation_plot2$Group,Annotation_plot2$ID2)
  
  # p = ggplot(Annotation_plot2, aes(x = ID2, y = log2(as.numeric(as.character(ratio.ave))),fill=Species)) +
  #   geom_boxplot(position = position_dodge(width = 0.9)) +
  #   #geom_quasirandom(aes(fill = Species,color= Species, alpha=Species,shape=Species),dodge.width = 0.9, varwidth = TRUE)+
  #   #geom_point(aes(fill = Species,color= Species, alpha=Species,shape=Species),position = position_jitterdodge(seed = 1, dodge.width = 0.9))+
  #   #geom_jitter(aes(fill = Species,color= Species, alpha=Species,shape=Species),width = 0.3,na.rm=TRUE)+
  #   scale_fill_manual(values=color_species)+
  #   scale_color_manual(values=color_species)+
  #   scale_shape_manual(values=shape_species)+
  #   scale_alpha_manual(values=alpha_species)+
  #   xlab("Family")+
  #   ylab("log2 ratio")+
  #   ggtitle(Cell)+
  #   #geom_hline(yintercept=c(0,1), linetype="dashed", size=1)+
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = rel(1)),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(), 
  #     axis.line = element_line(colour = "black"),
  #     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
  #     axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
  #     axis.title=element_text(colour="black",size=rel(1.4)),
  #     legend.key = element_rect(colour = "transparent", fill = "white"),
  #     legend.position="right",
  #     #legend.position = "none",
  #     legend.background = element_blank(),
  #     legend.text=element_text(size=rel(1)))
  # pdf(paste("Step5-2-",Cell,"_",Date,".pdf",sep=""),     
  #     width = 12,        # 5 x 300 pixels
  #     height = 5,
  #     pointsize = 10)
  # grid.draw(p)
  # dev.off()
  
  p = ggplot(Annotation_plot2, aes(x = ID2, y = log2(as.numeric(as.character(ratio.ave))))) +
    geom_violin(aes(),color = "black",fill=NA,scale="width") + 
    geom_jitter(aes(fill = Species,color= Species, alpha=Species,shape=Species),width = 0.3,na.rm=TRUE)+
    scale_fill_manual(values=color_species)+
    scale_color_manual(values=color_species)+
    scale_shape_manual(values=shape_species)+
    scale_alpha_manual(values=alpha_species)+
    xlab("Family")+
    ylab("log2 ratio")+
    ggtitle(Cell)+
    geom_hline(yintercept=c(0,1), linetype="dashed", size=1)+
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
      axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
      axis.title=element_text(colour="black",size=rel(1.4)),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      legend.position="right",
      #legend.position = "none",
      legend.background = element_blank(),
      legend.text=element_text(size=rel(1)))
  
  pdf(paste("Figure_S5F-",Cell,"_",Date,".pdf",sep=""),     
      width = 9,        # 5 x 300 pixels
      height = 5,
      pointsize = 10)
  grid.draw(p)
  dev.off()
  
  Summary_Table2$Species2 = ifelse(Summary_Table2$Group == "Consensus","Consensus",Summary_Table2$Species)
  Summary_Table2$Species2 = ifelse(Summary_Table2$Group %in% c("Positive","Negative"),"Control",Summary_Table2$Species2)
  Summary_Table2$Family2 = paste(Summary_Table2$Family,Summary_Table2$Frame,sep="_")
  if (Cell == "iPSC"){
    p = ggplot(Summary_Table2, aes(x = Family2, y = log2(iPSC.alpha))) +
      geom_violin(aes(),color = "black",fill=NA,scale="width") + 
      geom_jitter(aes(fill = Species2,color= Species2, alpha=Species2,shape=Species2,size=Species2),width = 0.3,na.rm=TRUE)+
      scale_fill_manual(values=color_species)+
      scale_color_manual(values=color_species)+
      scale_shape_manual(values=shape_species)+
      scale_alpha_manual(values=alpha_species)+
      scale_size_manual(values=size_species)+
      xlab("Family")+
      ylab("log2 alpha in iPSC")+
      ggtitle(Cell)+
      geom_hline(yintercept=c(0,1), linetype="dashed", size=1)+
      theme(
        plot.title = element_text(hjust = 0.5, size = rel(1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
        axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
        axis.title=element_text(colour="black",size=rel(1.4)),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.position="right",
        #legend.position = "none",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))
  } else {
    p = ggplot(Summary_Table2, aes(x = Family2, y = log2(NPC.alpha))) +
      geom_violin(aes(),color = "black",fill=NA,scale="width") + 
      geom_jitter(aes(fill = Species2,color= Species2, alpha=Species2,shape=Species2,size=Species2),width = 0.3,na.rm=TRUE)+
      scale_fill_manual(values=color_species)+
      scale_color_manual(values=color_species)+
      scale_shape_manual(values=shape_species)+
      scale_alpha_manual(values=alpha_species)+
      scale_size_manual(values=size_species)+
      xlab("Family")+
      ylab("log2 alpha in NPC")+
      ggtitle(Cell)+
      geom_hline(yintercept=c(0,1), linetype="dashed", size=1)+
      theme(
        plot.title = element_text(hjust = 0.5, size = rel(1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
        axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
        axis.title=element_text(colour="black",size=rel(1.4)),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))
  }
  pdf(paste("Figure_S5F-",Cell,"_",Date,".pdf",sep=""),     
      width = 9,        # 5 x 300 pixels
      height = 5,
      pointsize = 10)
  grid.draw(p)
  dev.off()
  
  ##### DNA
  my.formula = y~x
  # qlist = list()
  # i = 1
  # p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(dna_count.rep1)),as.numeric(as.character(dna_count.rep2)))) +
  #   geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
  #   scale_color_manual(values=color_Group)+
  #   scale_alpha_manual(values=alpha_Group)+
  #   xlab("rep1 (normalized DNA count)")+
  #   ylab("rep2 (normalized DNA count)")+
  #   geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  #   stat_poly_eq(formula = my.formula,
  #                eq.with.lhs = "italic(hat(y))~`=`~",
  #                p.digits = 3,
  #                aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
  #                parse = TRUE,size = 4.5,) +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = rel(1)),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(), 
  #     axis.line = element_line(colour = "black"),
  #     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
  #     axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
  #     axis.title=element_text(colour="black",size=rel(1.4)),
  #     legend.key = element_rect(colour = "transparent", fill = "white"),
  #     #legend.position="right",
  #     legend.position = "none",
  #     legend.background = element_blank(),
  #     legend.text=element_text(size=rel(1)))
  # qlist[[i]] <- ggplotGrob(p1)
  # i = i +1
  # p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(dna_count.rep1)),as.numeric(as.character(dna_count.rep3)))) +
  #   geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
  #   scale_color_manual(values=color_Group)+
  #   scale_alpha_manual(values=alpha_Group)+
  #   xlab("rep1 (normalized DNA count)")+
  #   ylab("rep3 (normalized DNA count)")+
  #   geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  #   stat_poly_eq(formula = my.formula,
  #                eq.with.lhs = "italic(hat(y))~`=`~",
  #                p.digits = 3,
  #                aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
  #                parse = TRUE,size = 4.5,) +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = rel(1)),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(), 
  #     axis.line = element_line(colour = "black"),
  #     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
  #     axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
  #     axis.title=element_text(colour="black",size=rel(1.4)),
  #     legend.key = element_rect(colour = "transparent", fill = "white"),
  #     #legend.position="right",
  #     legend.position = "none",
  #     legend.background = element_blank(),
  #     legend.text=element_text(size=rel(1)))
  # qlist[[i]] <- ggplotGrob(p1)
  # i = i +1
  # p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(dna_count.rep2)),as.numeric(as.character(dna_count.rep3)))) +
  #   geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
  #   scale_color_manual(values=color_Group)+
  #   scale_alpha_manual(values=alpha_Group)+
  #   xlab("rep2 (normalized DNA count)")+
  #   ylab("rep3 (normalized DNA count)")+
  #   geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  #   stat_poly_eq(formula = my.formula,
  #                eq.with.lhs = "italic(hat(y))~`=`~",
  #                p.digits = 3,
  #                aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
  #                parse = TRUE,size = 4.5,) +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = rel(1)),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(), 
  #     axis.line = element_line(colour = "black"),
  #     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
  #     axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
  #     axis.title=element_text(colour="black",size=rel(1.4)),
  #     legend.key = element_rect(colour = "transparent", fill = "white"),
  #     #legend.position="right",
  #     legend.position = "none",
  #     legend.background = element_blank(),
  #     legend.text=element_text(size=rel(1)))
  # qlist[[i]] <- ggplotGrob(p1)
  # i = i +1
  # pdf(paste("Step5-4-",Cell,"_",Date,".dna_count.pdf",sep=""),     
  #     width = 12,        # 5 x 300 pixels
  #     height = 4,
  #     pointsize = 10)
  # do.call("grid.arrange",c(qlist,ncol=3))
  # dev.off()
  # 
  # 
  # ##### rna count
  # qlist = list()
  # i = 1
  # p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(rna_count.rep1)),as.numeric(as.character(rna_count.rep2)))) +
  #   geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
  #   scale_color_manual(values=color_Group)+
  #   scale_alpha_manual(values=alpha_Group)+
  #   xlab("rep1 (normalized RNA count)")+
  #   ylab("rep2 (normalized RNA count)")+
  #   geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  #   stat_poly_eq(formula = my.formula,
  #                eq.with.lhs = "italic(hat(y))~`=`~",
  #                p.digits = 3,
  #                aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
  #                parse = TRUE,size = 4.5,) +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = rel(1)),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(), 
  #     axis.line = element_line(colour = "black"),
  #     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
  #     axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
  #     axis.title=element_text(colour="black",size=rel(1.4)),
  #     legend.key = element_rect(colour = "transparent", fill = "white"),
  #     #legend.position="right",
  #     legend.position = "none",
  #     legend.background = element_blank(),
  #     legend.text=element_text(size=rel(1)))
  # qlist[[i]] <- ggplotGrob(p1)
  # i = i +1
  # p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(rna_count.rep1)),as.numeric(as.character(rna_count.rep3)))) +
  #   geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
  #   scale_color_manual(values=color_Group)+
  #   scale_alpha_manual(values=alpha_Group)+
  #   xlab("rep1 (normalized RNA count)")+
  #   ylab("rep3 (normalized RNA count)")+
  #   geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  #   stat_poly_eq(formula = my.formula,
  #                eq.with.lhs = "italic(hat(y))~`=`~",
  #                p.digits = 3,
  #                aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
  #                parse = TRUE,size = 4.5,) +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = rel(1)),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(), 
  #     axis.line = element_line(colour = "black"),
  #     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
  #     axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
  #     axis.title=element_text(colour="black",size=rel(1.4)),
  #     legend.key = element_rect(colour = "transparent", fill = "white"),
  #     #legend.position="right",
  #     legend.position = "none",
  #     legend.background = element_blank(),
  #     legend.text=element_text(size=rel(1)))
  # qlist[[i]] <- ggplotGrob(p1)
  # i = i +1
  # p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(ratio.rep2)),as.numeric(as.character(ratio.rep3)))) +
  #   geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
  #   scale_color_manual(values=color_Group)+
  #   scale_alpha_manual(values=alpha_Group)+
  #   xlab("rep2 (normalized RNA count)")+
  #   ylab("rep3 (normalized RNA count)")+
  #   geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  #   stat_poly_eq(formula = my.formula,
  #                eq.with.lhs = "italic(hat(y))~`=`~",
  #                p.digits = 3,
  #                aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
  #                parse = TRUE,size = 4.5,) +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = rel(1)),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(), 
  #     axis.line = element_line(colour = "black"),
  #     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #     axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
  #     axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
  #     axis.title=element_text(colour="black",size=rel(1.4)),
  #     legend.key = element_rect(colour = "transparent", fill = "white"),
  #     #legend.position="right",
  #     legend.position = "none",
  #     legend.background = element_blank(),
  #     legend.text=element_text(size=rel(1)))
  # qlist[[i]] <- ggplotGrob(p1)
  # i = i +1
  # pdf(paste("Step5-5-",Cell,"_",Date,".RNA_count.pdf",sep=""),     
  #     width = 12,        # 5 x 300 pixels
  #     height = 4,
  #     pointsize = 10)
  # do.call("grid.arrange",c(qlist,ncol=3))
  # dev.off()
  # 
  ###### ratio
  qlist = list()
  i = 1
  p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(log2.rep1)),as.numeric(as.character(log2.rep2)))) +
    geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
    scale_color_manual(values=color_Group)+
    scale_alpha_manual(values=alpha_Group)+
    xlab("rep1 (normalized log2 ratio)")+
    ylab("rep2 (normalized log2 ratio)")+
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 p.digits = 3,
                 aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
                 parse = TRUE,size = 4.5,) +
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
      axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
      axis.title=element_text(colour="black",size=rel(1.4)),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      #legend.position="right",
      legend.position = "none",
      legend.background = element_blank(),
      legend.text=element_text(size=rel(1)))
  qlist[[i]] <- ggplotGrob(p1)
  i = i +1
  p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(log2.rep1)),as.numeric(as.character(log2.rep3)))) +
    geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
    scale_color_manual(values=color_Group)+
    scale_alpha_manual(values=alpha_Group)+
    xlab("rep1 (normalized log2 ratio)")+
    ylab("rep3 (normalized log2 ratio)")+
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 p.digits = 3,
                 aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
                 parse = TRUE,size = 4.5,) +
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
      axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
      axis.title=element_text(colour="black",size=rel(1.4)),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      #legend.position="right",
      legend.position = "none",
      legend.background = element_blank(),
      legend.text=element_text(size=rel(1)))
  qlist[[i]] <- ggplotGrob(p1)
  i = i +1
  p1 = ggplot(Annotation_plot2, aes(as.numeric(as.character(log2.rep2)),as.numeric(as.character(log2.rep3)))) +
    geom_point(aes(color=Group,alpha=Group),size=2,shape=19) + 
    scale_color_manual(values=color_Group)+
    scale_alpha_manual(values=alpha_Group)+
    xlab("rep2 (normalized log2 ratio)")+
    ylab("rep3 (normalized log2 ratio)")+
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 p.digits = 3,
                 aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
                 parse = TRUE,size = 4.5,) +
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
      axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
      axis.title=element_text(colour="black",size=rel(1.4)),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      #legend.position="right",
      legend.position = "none",
      legend.background = element_blank(),
      legend.text=element_text(size=rel(1)))
  qlist[[i]] <- ggplotGrob(p1)
  i = i +1
  pdf(paste("Figure_S5E-",Cell,"_",Date,".log2ratio.pdf",sep=""),     
      width = 12,        # 5 x 300 pixels
      height = 4,
      pointsize = 10)
  do.call("grid.arrange",c(qlist,ncol=3))
  dev.off()
}