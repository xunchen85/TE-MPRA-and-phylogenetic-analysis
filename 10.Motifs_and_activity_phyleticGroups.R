#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

#library(phyloseq)
library(treeio)
library(ggtree)
library(ggplot2)
library(aplot)
library(ape)
library(dplyr)
library(grid)
library(gridExtra)
library(ggstance)
library(splitstackshape)
library(RColorBrewer)
library(ggtreeExtra)
library(ggnewscale)
library(Biostrings)
library(randomcoloR)
library(ggbeeswarm)
library(ggpubr)

Date = "2023_9_13"

###################################################### motif
##### step 1: subfamily and phyletic group
MER11_hg191 = data.frame(species = "hg19",group="FG1",subfamily = c("11A_i","11A_k","11A_q","11A_m","11A_n","11A_p","11A_a","11A_o","11A_r","11A_l"))
MER11_hg192 = data.frame(species = "hg19",group="FG2",subfamily = c("11A_b","11A_c","11C_w","11A_j","11A_d","11B_i","11A_e","11A_f"))
MER11_hg193 = data.frame(species = "hg19",group="FG3",subfamily = c("11A_g","11B_h","11B_g","11C_v","11A_h","11B_k","11B_f","11B_n","11B_r","11B_o","11C_o","11B_q","11C_r"))
MER11_hg194 = data.frame(species = "hg19",group="FG4",subfamily = c("11B_e","11B_c","11B_b","11B_p","11B_m","11C_x","11B_d","11B_l","11C_q","11B_j","11C_a","11C_y","11B_a","11C_z","11C_t",
                                                                    "11C_a2","11C_s","11C_b2","11C_c","11C_b","11C_d","11C_e","11C_f","11C_g","11C_h","11C_p","11C_i","11C_j",
                                                                    "11C_n","11C_m","11C_k","11C_l"))

MER11_macque1 = data.frame(species = "macFas5",group="FG1",subfamily = c("MER11A_g17","MER11A_g18","MER11A_g15","MER11A_g2","MER11A_g22","MER11A_g33","MER11A_g19"))
MER11_macque2 = data.frame(species = "macFas5",group="FG2",subfamily = c("MER11A_g13","MER11A_g8","MER11A_g16","MER11A_g21"))
MER11_macque3 = data.frame(species = "macFas5",group="FG3",subfamily = c("MER11A_g5","MER11A_g24","MER11A_g25","MER11A_g26","MER11A_g23","MER11A_g30","MER11A_g11","MER11A_g6","MER11A_g9"))
MER11_macque4 = data.frame(species = "macFas5",group="FG4",subfamily = c("MER11A_g7","MER11A_g29","MER11A_g1","MER11A_g28","MER11A_g27","MER11A_g20","MER11A_g32","MER11A_g14",
                                                                         "MER11A_g3","MER11A_g31","MER11A_g10"))
MER11_both = rbind(MER11_hg191,MER11_hg192,MER11_hg193,MER11_hg194,
                   MER11_macque1,MER11_macque2,MER11_macque3,MER11_macque4)
MER11_both$Order = NA
MER11_both[MER11_both$species == "hg19",]$Order = 1:nrow(MER11_both[MER11_both$species == "hg19",])
MER11_both[MER11_both$species == "macFas5",]$Order = 1:nrow(MER11_both[MER11_both$species == "macFas5",])

rm(MER11_hg191,MER11_hg192,MER11_hg193,MER11_hg194,
   MER11_macque1,MER11_macque2,MER11_macque3,MER11_macque4)
head(MER11_both)

MER11_both = MER11_both[!duplicated(MER11_both$subfamily),c("group","subfamily","Order")]

##### step2: summary table 2
Summary_Table2 = read.csv("input/Summary_Table2_2023_1_5.subfamilyInfo.csv")

# achieved subfamily info of macFas5
macaque.tree = read.csv("input/macFas5_rmsk_TE_0bp_MER11A.family.rename_200bp_10_0.02_all_95_group_final.csv")

# load the macaque subfamily information
Summary_Table2.kept = merge(Summary_Table2,macaque.tree[,c("label","label.final.correctedName")],by.x="instanceID.renamed",by.y="label",all.x=T)

Summary_Table2.kept$label.final.correctedName = ifelse(Summary_Table2.kept$Family == "MER11A" &
                                                         Summary_Table2.kept$species == "macFas5",
                                                       Summary_Table2.kept$label.final.correctedName,
                                                       Summary_Table2.kept$subfamily.name.final)
Summary_Table2.kept$label.final.correctedName = ifelse(Summary_Table2.kept$Group != "Instance",NA,Summary_Table2.kept$label.final.correctedName)
Summary_Table2.kept[Summary_Table2.kept$Species == "panTro4",]$label.final.correctedName

# kept MER11, frames with high quality activity and instances with group info 
Summary_Table2.fimo = Summary_Table2.kept[!is.na(Summary_Table2.kept$is_kept) & 
                                            grepl("highQuality",Summary_Table2.kept$is_kept) &
                                            Summary_Table2.kept$Group == "Instance" &
                                            Summary_Table2.kept$Family %in% c("MER11A","MER11B","MER11C"),
                                          c("uniqueID_new","instanceID.renamed","label.final.correctedName","Species","Family","iPSC.alpha.Zscore","NPC.alpha.Zscore")]

Summary_Table2.fimo = merge(Summary_Table2.fimo,MER11_both,by.x="label.final.correctedName",by.y="subfamily",all.x=T)
Summary_Table2.fimo = Summary_Table2.fimo[!is.na(Summary_Table2.fimo$group),]

###### step3: load the fimo motif sets
# load fimo results
Fimo.all = read.delim("input/MER11_34_52_frame_fimo.tsv_2023_1_7.gz",header=F,sep="")
colnames(Fimo.all) = c("Family.frame","DB","motif_id","motif_name","uniqueID.new","start","stop","strand","score","p_value","q_value","matched_sequence")
Fimo.all$uniqueID.new.DB.motif_name = paste(Fimo.all$Family.frame,Fimo.all$uniqueID.new,Fimo.all$DB,Fimo.all$motif_name)

# kept unique motif per sequence
Fimo.all.sum = data.frame(Fimo.all %>% group_by(Family.frame,DB,uniqueID.new,motif_name) %>% dplyr::count())
Fimo.all.sum = Fimo.all.sum[Fimo.all.sum$DB == "JASPAR2022" & !grepl("C_",Fimo.all.sum$uniqueID.new),]
Fimo.all.sum$motif_name = toupper(Fimo.all.sum$motif_name)
Fimo.all.sum = merge(Fimo.all.sum,Summary_Table2.fimo,by.x="uniqueID.new",by.y="uniqueID_new",all.x=T)
Fimo.all.sum = Fimo.all.sum[Fimo.all.sum$Family.frame == "MER11_combined_F2",]
head(Fimo.all.sum)

####### step 4: summarize by motifs
### 4.0
# not used
#MER11_F2.motifs = c("RARA::RXRG","INSM1","DMBX1","CRX","NR5A1","FOXF2","RFX6","ZNF331","ZNF136",
#                    "ZIC3","TEAD4","TEAD1","TFCP2","GATA5","POU2F1::SOX2","SOX17","IKZF1","ZNF317",
#                    "SOX15","ZIM3","TFAP2A","TFAP2B","TFAP2C","RORA","CDX2","PPARG::RXRA","RORA","IRF6","ZBED2","HSF1","HSF4")
# final set
MER11_F2.motifs = c("RARA::RXRG","INSM1","DMBX1","CRX","NR5A1","FOXF2","RFX6","ZNF331","ZNF136",
                    "ZIC3","TEAD4","TEAD1","TFCP2","GATA5","POU2F1::SOX2","SOX17","IKZF1","ZNF317")
MER11_F2.motifs2 = c("RARA::RXRG","INSM1","DMBX1","CRX","NR5A1","FOXF2","RFX6","ZNF331","ZNF136",
                     "ZIC3","TEAD4","TEAD1","TFCP2","GATA5","POU2F1::SOX2","SOX17","IKZF1","ZNF317","SOX15","IRF6","ZBED2","HSF1","HSF4")

#Fimo.all.sum.plot = Fimo.all.sum[!is.na(Fimo.all.sum$group),]
Fimo.all.sum.plot = Fimo.all.sum

# kept candidate motifs
Fimo.all.sum.plot = Fimo.all.sum.plot[Fimo.all.sum.plot$motif_name %in% MER11_F2.motifs2,]

# count for each motif
Fimo.all.sum.plot.sum = data.frame(Fimo.all.sum.plot[!is.na(Fimo.all.sum.plot$group),] %>% group_by(label.final.correctedName,DB,motif_name,Order,Species) %>% dplyr::count())

# total count
Fimo.all.sum.plot.sum.total = data.frame(Summary_Table2.fimo %>% group_by(Family,label.final.correctedName,Order,Species) %>% dplyr::count())

# combined
Fimo.all.sum.plot.sum = merge(Fimo.all.sum.plot.sum,Fimo.all.sum.plot.sum.total[,c("label.final.correctedName","n")],by="label.final.correctedName",all.x=T)
Fimo.all.sum.plot.sum$perC = Fimo.all.sum.plot.sum$n.x/Fimo.all.sum.plot.sum$n.y
Fimo.all.sum.plot.sum$perC = ifelse(Fimo.all.sum.plot.sum$n.y>=10,Fimo.all.sum.plot.sum$perC,NA)
Fimo.all.sum.plot.sum$unique_ID = paste(Fimo.all.sum.plot.sum$Species,Fimo.all.sum.plot.sum$label.final.correctedName)


### 4.1 count for each motif per subfamily
Fimo.all.sum.plot.sum.subfamily = data.frame(Fimo.all.sum.plot[!is.na(Fimo.all.sum.plot$group),] %>% 
                                               group_by(Family,DB,motif_name,Species) %>% 
                                               dplyr::count())

Fimo.all.sum.plot.sum.subfamily$uniqueID.group = paste(Fimo.all.sum.plot.sum.subfamily$Species,
                                                       Fimo.all.sum.plot.sum.subfamily$Family,
                                                       Fimo.all.sum.plot.sum.subfamily$motif_name)
# total count
Fimo.all.sum.plot.sum.subfamily.total.tmp1 = data.frame(Summary_Table2.fimo %>% group_by(Family,Species) %>% dplyr::count())
Fimo.all.sum.plot.sum.subfamily.total = data.frame("tmp1"=NA,"tmp2"=NA)
for (Motif in MER11_F2.motifs2){
  Fimo.all.sum.plot.sum.subfamily.total.tmp2 = Fimo.all.sum.plot.sum.subfamily.total.tmp1
  Fimo.all.sum.plot.sum.subfamily.total.tmp2$motif = Motif
  if (nrow(Fimo.all.sum.plot.sum.subfamily.total) == 1){
    Fimo.all.sum.plot.sum.subfamily.total = Fimo.all.sum.plot.sum.subfamily.total.tmp2
  } else {
    Fimo.all.sum.plot.sum.subfamily.total = rbind(Fimo.all.sum.plot.sum.subfamily.total,Fimo.all.sum.plot.sum.subfamily.total.tmp2)
  }
}
rm (Fimo.all.sum.plot.sum.subfamily.total.tmp2,Fimo.all.sum.plot.sum.subfamily.total.tmp1)
Fimo.all.sum.plot.sum.subfamily.total$uniqueID.group = paste(Fimo.all.sum.plot.sum.subfamily.total$Species,
                                                             Fimo.all.sum.plot.sum.subfamily.total$Family,
                                                             Fimo.all.sum.plot.sum.subfamily.total$motif)
Fimo.all.sum.plot.sum.subfamily.total$Type = "subfamily"
Fimo.all.sum.plot.sum.subfamily.total = merge(Fimo.all.sum.plot.sum.subfamily.total,Fimo.all.sum.plot.sum.subfamily[,c("uniqueID.group","n")],by="uniqueID.group",all.x=T)
colnames(Fimo.all.sum.plot.sum.subfamily.total)[2] = "group"

### 4.2 count for each motif by each group
Fimo.all.sum.plot.sum.group = data.frame(Fimo.all.sum.plot[!is.na(Fimo.all.sum.plot$group),] %>% 
                                           group_by(group,DB,motif_name,Species) %>% 
                                           dplyr::count())

Fimo.all.sum.plot.sum.group$uniqueID.group = paste(Fimo.all.sum.plot.sum.group$Species,
                                                   Fimo.all.sum.plot.sum.group$group,
                                                   Fimo.all.sum.plot.sum.group$motif_name)
# total count per phyletic group
Fimo.all.sum.plot.sum.group.total.tmp1 = data.frame(Summary_Table2.fimo %>% group_by(group,Species) %>% dplyr::count())

Fimo.all.sum.plot.sum.group.total = data.frame("tmp1"=NA,"tmp2"=NA)
for (Motif in MER11_F2.motifs2){
  Fimo.all.sum.plot.sum.group.total.tmp2 = Fimo.all.sum.plot.sum.group.total.tmp1
  Fimo.all.sum.plot.sum.group.total.tmp2$motif = Motif
  if (nrow(Fimo.all.sum.plot.sum.group.total) == 1){
    Fimo.all.sum.plot.sum.group.total = Fimo.all.sum.plot.sum.group.total.tmp2
  } else {
    Fimo.all.sum.plot.sum.group.total = rbind(Fimo.all.sum.plot.sum.group.total,Fimo.all.sum.plot.sum.group.total.tmp2)
  }
}
rm (Fimo.all.sum.plot.sum.group.total.tmp2,Fimo.all.sum.plot.sum.group.total.tmp1)
Fimo.all.sum.plot.sum.group.total$uniqueID.group = paste(Fimo.all.sum.plot.sum.group.total$Species,
                                                         Fimo.all.sum.plot.sum.group.total$group,Fimo.all.sum.plot.sum.group.total$motif)
Fimo.all.sum.plot.sum.group.total$Type = "group"
Fimo.all.sum.plot.sum.group.total = merge(Fimo.all.sum.plot.sum.group.total,Fimo.all.sum.plot.sum.group[,c("uniqueID.group","n")],by="uniqueID.group",all.x=T)

Fimo.all.sum.plot.sum.group.total.both = rbind(Fimo.all.sum.plot.sum.group.total,Fimo.all.sum.plot.sum.subfamily.total)

# combined
Fimo.all.sum.plot.sum.group.total.both$perC = Fimo.all.sum.plot.sum.group.total.both$n.y/Fimo.all.sum.plot.sum.group.total.both$n.x
Fimo.all.sum.plot.sum.group.total.both$perC = ifelse(Fimo.all.sum.plot.sum.group.total.both$n.x>=10,Fimo.all.sum.plot.sum.group.total.both$perC,NA)
Fimo.all.sum.plot.sum.group.total.both$uniqueID2 = paste(Fimo.all.sum.plot.sum.group.total.both$Species,Fimo.all.sum.plot.sum.group.total.both$Type,Fimo.all.sum.plot.sum.group.total.both$motif)
############ plot 1: between phyletic groups and subfamilies

max_ones = data.frame(Fimo.all.sum.plot.sum.group.total.both %>%
                        group_by(uniqueID2,motif) %>%
                        dplyr::slice(which.max(perC)))
min_ones = data.frame(Fimo.all.sum.plot.sum.group.total.both %>%
                        group_by(uniqueID2,motif) %>%
                        dplyr::slice(which.min(perC)))

Fimo.all.sum.plot.sum.group.total.both.plot = merge(max_ones[,c("uniqueID2","group","Species","n.x","motif","Type","perC")],min_ones[,c("uniqueID2","perC")],by="uniqueID2",all=T)
Fimo.all.sum.plot.sum.group.total.both.plot$perC.diff = Fimo.all.sum.plot.sum.group.total.both.plot$perC.x-Fimo.all.sum.plot.sum.group.total.both.plot$perC.y
Fimo.all.sum.plot.sum.group.total.both.plot.hg19 = Fimo.all.sum.plot.sum.group.total.both.plot[Fimo.all.sum.plot.sum.group.total.both.plot$Species == "hg19" & 
                                                                                                 Fimo.all.sum.plot.sum.group.total.both.plot$motif %in% MER11_F2.motifs,]
Fimo.all.sum.plot.sum.group.total.both.plot.hg19$motif = factor(Fimo.all.sum.plot.sum.group.total.both.plot.hg19$motif,levels=MER11_F2.motifs)
glist = list()
Order = 1

p1 = ggplot(Fimo.all.sum.plot.sum.group.total.both.plot.hg19, aes(x=motif,y=perC.diff*100,group=Type)) +
  geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(perC.diff*100,1)),  color="black",
            position = position_dodge(0.9), size=4,vjust=-0.2, hjust=-0.2,angle=45)+
  scale_fill_manual(values =c("#016A70","#A2C579"))+
  ggtitle("MER11 subfamily group")+
  #geom_hline(yintercept=c(20,40), linetype="dashed", color = "black")+
  #ylim(0,100)+
  #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
  scale_x_discrete(drop = FALSE)+
  ylab("specificity")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.line.x = element_blank(),
        #axis.line.y = element_blank(),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.95),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1))) 
glist[[Order]] <- ggplotGrob(p1)
Order = Order + 1
p1 = ggplot(Fimo.all.sum.plot.sum.group.total.both.plot.hg19, aes(x=motif,y=perC.x*100,group=Type)) +
  geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(perC.x*100,1)),  color="black",
            position = position_dodge(0.9), size=4,vjust=-0.2, hjust=-0.2,angle=45)+
  scale_fill_manual(values =c("#016A70","#A2C579"))+
  ggtitle("MER11 subfamily group")+
  scale_x_discrete(drop = FALSE)+
  ylab("max")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.line.x = element_blank(),
        #axis.line.y = element_blank(),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.95),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1))) 
glist[[Order]] <- ggplotGrob(p1)
Order = Order + 1
p1 = ggplot(Fimo.all.sum.plot.sum.group.total.both.plot.hg19, aes(x=motif,y=perC.y*100,group=Type)) +
  geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(perC.y*100,1)),  color="black",
            position = position_dodge(0.9), size=4,vjust=-0.2, hjust=-0.2,angle=45)+
  scale_fill_manual(values =c("#016A70","#A2C579"))+
  ggtitle("MER11 subfamily group")+
  scale_x_discrete(drop = FALSE)+
  ylab("min")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.line.x = element_blank(),
        #axis.line.y = element_blank(),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.95),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1))) 
glist[[Order]] <- ggplotGrob(p1)
Order = Order + 1

pdf(paste("Figure_3G_S5H.pdf",sep=""),    # create PNG for the heat map        
    width = 12,        # 5 x 300 pixels
    height = 12,
    pointsize = 10)        # smaller font size
do.call("grid.arrange",c(glist,ncol=1))
dev.off()

############ plot 2: between between species
glist2.group = list()
Order1 = 1
motifName = "SOX15"
for(motifName in unique(MER11_F2.motifs2)){
  Fimo.all.sum.plot.sum.group.total2 = Fimo.all.sum.plot.sum.group.total.both[Fimo.all.sum.plot.sum.group.total.both$Type == "group",]
  # group
  Fimo.all.sum.plot.sum.group.total2 = Fimo.all.sum.plot.sum.group.total2[Fimo.all.sum.plot.sum.group.total2$motif == motifName,]
  if (nrow(Fimo.all.sum.plot.sum.group.total2) > 0) {
    m1 = ggplot(Fimo.all.sum.plot.sum.group.total2, aes(x=group,y=perC*100,group=Species)) +
      geom_line(aes(color=Species),size=1) + 
      geom_point(aes(color=Species,fill=Species),size=3)+
      ylab("% frames")+
      xlab("")+
      expand_limits(y = 0)+
      scale_x_discrete(drop = FALSE)+
      scale_color_manual(values=c("macFas5"="#2166ac","hg19"="#66c2a5","panTro4" = "#b2df8a"))+
      ggtitle(motifName) + 
      theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "#f0f0f0"),
            axis.line = element_line(colour = "black"),
            #axis.line.x = element_blank(),
            #axis.line.y = element_blank(),
            axis.text.y=element_text(colour="black",size=rel(1)),
            axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
            axis.title=element_text(colour="black",size=rel(1)),
            legend.position="right",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1))) 
    glist2.group[[Order1]] <- ggplotGrob(m1)
    Order1 = Order1 + 1
  }
}
pdf(paste("Figure_S10A.pdf",sep=""),    # create PNG for the heat map
    width = 14,        # 5 x 300 pixels
    height = 3.5*5,
    pointsize = 10)        # smaller font size
do.call("grid.arrange",c(glist2.group,ncol=5))
dev.off()

### plot 2: activity between species containing sox and pou::sox motifs
Fimo.all.sum.plot.FG4 = Fimo.all.sum.plot[Fimo.all.sum.plot$group == "FG4",]
Fimo.all.sum.plot.FG4 = Fimo.all.sum.plot.FG4[Fimo.all.sum.plot.FG4$motif_name %in% c("SOX15","SOX17","POU2F1::SOX2"),]
Fimo.all.sum.plot.FG4.no = Fimo.all.sum.plot.FG4[!(Fimo.all.sum.plot.FG4$motif_name %in% c("SOX15","SOX17","POU2F1::SOX2")),]
Fimo.all.sum.plot.FG4$x.axis = paste(Fimo.all.sum.plot.FG4$Species,Fimo.all.sum.plot.FG4$motif_name)
Fimo.all.sum.plot.FG4 = Fimo.all.sum.plot.FG4[!Fimo.all.sum.plot.FG4$x.axis %in% c("macFas5 SOX15","macFas5 SOX17"),]
Fimo.all.sum.plot.FG4$x.axis = factor(Fimo.all.sum.plot.FG4$x.axis,levels=c("macFas5 POU2F1::SOX2","hg19 POU2F1::SOX2","panTro4 POU2F1::SOX2",
                                                                            "hg19 SOX15","panTro4 SOX15","hg19 SOX17","panTro4 SOX17"))
Fimo.all.sum.plot.FG4.SOX15.id = Fimo.all.sum.plot.FG4[Fimo.all.sum.plot.FG4$motif_name == "SOX15",]$uniqueID.new
Fimo.all.sum.plot.FG4.SOX17.id = Fimo.all.sum.plot.FG4[Fimo.all.sum.plot.FG4$motif_name == "SOX17",]$uniqueID.new

Fimo.all.sum.plot.FG4 = Fimo.all.sum.plot.FG4[!(Fimo.all.sum.plot.FG4$motif_name == "POU2F1::SOX2" & 
                                                  Fimo.all.sum.plot.FG4$uniqueID.new %in% c(Fimo.all.sum.plot.FG4.SOX15.id,Fimo.all.sum.plot.FG4.SOX17.id)),]

my_comparisons = list(c("macFas5 POU2F1::SOX2","hg19 POU2F1::SOX2"),c("macFas5 POU2F1::SOX2","panTro4 POU2F1::SOX2"),
                      c("hg19 POU2F1::SOX2","hg19 SOX15"),c("hg19 POU2F1::SOX2","hg19 SOX17"),
                      c("panTro4 POU2F1::SOX2","panTro4 SOX15"),c("panTro4 POU2F1::SOX2","panTro4 SOX17"))
my_comparisons.hg19 = list(c("macFas5 POU2F1::SOX2","hg19 POU2F1::SOX2"),
                           c("hg19 POU2F1::SOX2","hg19 SOX15"),c("hg19 POU2F1::SOX2","hg19 SOX17"))
my_comparisons.panTro4 = list(c("macFas5 POU2F1::SOX2","panTro4 POU2F1::SOX2"),
                              c("panTro4 POU2F1::SOX2","panTro4 SOX15"),c("panTro4 POU2F1::SOX2","panTro4 SOX17"))

# 
p1_0 = ggplot(Fimo.all.sum.plot.FG4, aes(x = x.axis, y = iPSC.alpha.Zscore,color=Species)) +
  geom_quasirandom(aes(color = Species),method = "quasirandom",alpha=.7) + 
  stat_summary(fun = "mean", 
               fun.min = function(x)mean(x)-sd(x), 
               fun.max = function(x)mean(x) + sd(x), 
               geom = "pointrange", 
               position = position_dodge(width = .9),shape=95,color="black",size=.5)+
  scale_color_manual(values=c("macFas5"="#2166ac","hg19"="#66c2a5","panTro4" = "#b2df8a"))+
  stat_compare_means(method = "t.test",comparisons = my_comparisons,aes(label = paste0("", ..p.format..))) +
  xlab("") + 
  ylab("Insert length (bp)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        
        axis.text.x=element_text(colour="black",size=rel(1.2),angle = 90,hjust = 0.95,vjust = .2),
        axis.title=element_text(colour="black",size=rel(1.5)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        legend.position="none",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1.2)))

p1_1 = ggplot(Fimo.all.sum.plot.FG4[Fimo.all.sum.plot.FG4$Species != "panTro4",], aes(x = x.axis, y = iPSC.alpha.Zscore,color=Species)) +
  geom_quasirandom(aes(color = Species),method = "quasirandom",alpha=.7) + 
  stat_summary(fun = "mean", 
               fun.min = function(x)mean(x)-sd(x), 
               fun.max = function(x)mean(x) + sd(x), 
               geom = "pointrange", 
               position = position_dodge(width = .9),shape=95,color="black",size=.5)+
  scale_color_manual(values=c("macFas5"="#2166ac","hg19"="#66c2a5","panTro4" = "#b2df8a"))+
  stat_compare_means(method = "t.test",comparisons = my_comparisons.hg19,aes(label = paste0("", ..p.format..))) +
  xlab("") + 
  ylab("Insert length (bp)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        
        axis.text.x=element_text(colour="black",size=rel(1.2),angle = 90,hjust = 0.95,vjust = .2),
        axis.title=element_text(colour="black",size=rel(1.5)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        legend.position="none",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1.2)))

p1_2 = ggplot(Fimo.all.sum.plot.FG4[Fimo.all.sum.plot.FG4$Species != "hg19",], aes(x = x.axis, y = iPSC.alpha.Zscore,color=Species)) +
  geom_quasirandom(aes(color = Species),method = "quasirandom",alpha=.7) + 
  stat_summary(fun = "mean", 
               fun.min = function(x)mean(x)-sd(x), 
               fun.max = function(x)mean(x) + sd(x), 
               geom = "pointrange", 
               position = position_dodge(width = .9),shape=95,color="black",size=.5)+
  scale_color_manual(values=c("macFas5"="#2166ac","hg19"="#66c2a5","panTro4" = "#b2df8a"))+
  stat_compare_means(method = "t.test",comparisons = my_comparisons.panTro4,aes(label = paste0("", ..p.format..))) +
  xlab("") + 
  ylab("Insert length (bp)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        
        axis.text.x=element_text(colour="black",size=rel(1.2),angle = 90,hjust = 0.95,vjust = .2),
        axis.title=element_text(colour="black",size=rel(1.5)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        legend.position="none",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1.2)))

p1_0.2 = ggplot(Fimo.all.sum.plot.FG4[Fimo.all.sum.plot.FG4$iPSC.alpha.Zscore<=200,], aes(x = x.axis, y = iPSC.alpha.Zscore,color=Species)) +
  geom_quasirandom(aes(color = Species),method = "quasirandom",alpha=.7) + 
  stat_summary(fun = "mean", 
               fun.min = function(x)mean(x)-sd(x), 
               fun.max = function(x)mean(x) + sd(x), 
               geom = "pointrange", 
               position = position_dodge(width = .9),shape=95,color="black",size=.5)+
  scale_color_manual(values=c("macFas5"="#2166ac","hg19"="#66c2a5","panTro4" = "#b2df8a"))+
  stat_compare_means(method = "t.test",comparisons = my_comparisons,aes(label = paste0("", ..p.format..))) +
  xlab("") + 
  ylab("Insert length (bp)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        
        axis.text.x=element_text(colour="black",size=rel(1.2),angle = 90,hjust = 0.95,vjust = .2),
        axis.title=element_text(colour="black",size=rel(1.5)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        legend.position="none",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1.2)))

p1_1.2 = ggplot(Fimo.all.sum.plot.FG4[Fimo.all.sum.plot.FG4$Species != "panTro4"& Fimo.all.sum.plot.FG4$iPSC.alpha.Zscore<=200,], aes(x = x.axis, y = iPSC.alpha.Zscore,color=Species)) +
  geom_quasirandom(aes(color = Species),method = "quasirandom",alpha=.7) + 
  stat_summary(fun = "mean", 
               fun.min = function(x)mean(x)-sd(x), 
               fun.max = function(x)mean(x) + sd(x), 
               geom = "pointrange", 
               position = position_dodge(width = .9),shape=95,color="black",size=.5)+
  scale_color_manual(values=c("macFas5"="#2166ac","hg19"="#66c2a5","panTro4" = "#b2df8a"))+
  stat_compare_means(method = "t.test",comparisons = my_comparisons.hg19,aes(label = paste0("", ..p.format..))) +
  xlab("") + 
  ylab("Insert length (bp)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        
        axis.text.x=element_text(colour="black",size=rel(1.2),angle = 90,hjust = 0.95,vjust = .2),
        axis.title=element_text(colour="black",size=rel(1.5)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        legend.position="none",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1.2)))

p1_2.2 = ggplot(Fimo.all.sum.plot.FG4[Fimo.all.sum.plot.FG4$Species != "hg19" & Fimo.all.sum.plot.FG4$iPSC.alpha.Zscore<=200,], aes(x = x.axis, y = iPSC.alpha.Zscore,color=Species)) +
  geom_quasirandom(aes(color = Species),method = "quasirandom",alpha=.7) + 
  stat_summary(fun = "mean", 
               fun.min = function(x)mean(x)-sd(x), 
               fun.max = function(x)mean(x) + sd(x), 
               geom = "pointrange", 
               position = position_dodge(width = .9),shape=95,color="black",size=.5)+
  scale_color_manual(values=c("macFas5"="#2166ac","hg19"="#66c2a5","panTro4" = "#b2df8a"))+
  stat_compare_means(method = "t.test",comparisons = my_comparisons.panTro4,aes(label = paste0("", ..p.format..))) +
  xlab("") + 
  ylab("Insert length (bp)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        
        axis.text.x=element_text(colour="black",size=rel(1.2),angle = 90,hjust = 0.95,vjust = .2),
        axis.title=element_text(colour="black",size=rel(1.5)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        legend.position="none",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1.2)))

gA1 <- ggplotGrob(p1_0)
gA2 <- ggplotGrob(p1_0.2)
gA3 <- ggplotGrob(p1_1)
gA4 <- ggplotGrob(p1_1.2)
gA5 <- ggplotGrob(p1_2)
gA6 <- ggplotGrob(p1_2.2)

g = cbind(gA1,gA2,gA3,gA4,gA5,gA6,size = "last")
pdf(paste("Figure_5F_comparison-activity-motif_",Date,".pdf",sep=""),    # create PNG for the heat map
    width = 36,        # 5 x 300 pixels
    height = 7.5,
    pointsize = 10)        # smaller font size
grid.draw(g)
dev.off()


