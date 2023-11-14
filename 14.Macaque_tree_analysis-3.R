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

######################### step 1 achieve the list of tree files
# group info
MER11_hg191 = data.frame(species = "hg19",group="FG1",subfamily = c("11A_p","11A_m","11A_n","11A_q","11A_i","11A_k","11A_a","11A_o","11A_r","11A_l"))
MER11_hg192 = data.frame(species = "hg19",group="FG2",subfamily = c("11A_b","11A_c","11C_w","11A_j","11A_d","11A_e","11B_i","11A_f"))
MER11_hg193 = data.frame(species = "hg19",group="FG3",subfamily = c("11A_g","11B_h","11B_g","11A_h","11C_v","11B_k","11B_f","11B_n","11B_r","11B_o","11C_o","11B_q","11C_r"))
MER11_hg194 = data.frame(species = "hg19",group="FG4",subfamily = c("11B_e","11B_c","11B_b","11B_p","11B_m","11C_x","11B_d","11B_l","11C_q","11B_j","11C_a","11C_y","11B_a","11C_t",
                                                                    "11C_z","11C_s","11C_a2","11C_b2","11C_b","11C_c","11C_d","11C_e","11C_f","11C_g","11C_h","11C_j","11C_i","11C_p",
                                                                    "11C_n","11C_m","11C_k","11C_l"))

MER11_macque1 = data.frame(species = "macFas5",group="FG1",subfamily = c("MER11A_g17","MER11A_g18","MER11A_g15","MER11A_g2","MER11A_g22","MER11A_g19","MER11A_g33"))
MER11_macque2 = data.frame(species = "macFas5",group="FG2",subfamily = c("MER11A_g8","MER11A_g13","MER11A_g16","MER11A_g21"))
MER11_macque3 = data.frame(species = "macFas5",group="FG3",subfamily = c("MER11A_g5","MER11A_g24","MER11A_g25","MER11A_g26","MER11A_g23","MER11A_g30","MER11A_g11","MER11A_g6","MER11A_g9"))
MER11_macque4 = data.frame(species = "macFas5",group="FG4",subfamily = c("MER11A_g7","MER11A_g29","MER11A_g1","MER11A_g28","MER11A_g27","MER11A_g20","MER11A_g32","MER11A_g14",
                                                                         "MER11A_g31","MER11A_g10","MER11A_g3"))
MER11_both = rbind(MER11_hg191,MER11_hg192,MER11_hg193,MER11_hg194,MER11_macque1,MER11_macque2,MER11_macque3,MER11_macque4)

MER11_both$Order = NA
MER11_both[MER11_both$species == "hg19",]$Order = 1:nrow(MER11_both[MER11_both$species == "hg19",])
MER11_both[MER11_both$species == "macFas5",]$Order = 1:nrow(MER11_both[MER11_both$species == "macFas5",])
rm(MER11_hg191,MER11_hg192,MER11_hg193,MER11_hg194,MER11_macque1,MER11_macque2,MER11_macque3,MER11_macque4)

# subfamily info of human
# obtain hg19 subfamily info
Summary_Table1 = read.csv("input/Summary_Table1_2023_9_5.subfamilyInfo.csv")
Summary_Table2 = read.csv("input/Summary_Table2_2023_1_5.subfamilyInfo.csv")

# info subfamily
subfamily_info = Summary_Table1[Summary_Table1$species %in% c("hg19","macFas5") & Summary_Table1$TEfamily.x %in% c("MER11A","MER11B","MER11C"),]
subfamily_info$subfamily.name.final.species = NA

# achieved subfamily info of macFas5
subfamily_info.macFas5 = read.csv("input_trees/macFas5_rmsk_TE_0bp_MER11A.family.rename_200bp_10_0.02_all_95_group_final.csv")
subfamily_info = merge(subfamily_info,subfamily_info.macFas5[,c("label","consensus.name")],by.x="instanceID.renamed",by.y="label",all.x=T)
# combined
subfamily_info$subfamily.name.final.species = ifelse(subfamily_info$species == "hg19",subfamily_info$subfamily.name.final,subfamily_info$consensus.name.y)

#### intersect table between human and macaque (same step as another script)
intersect_Table = read.csv("input/Summary_Table1.intersect_reciprocal_hu_ch_ma_2023_2_20.csv")

intersect_Table.both = intersect_Table[intersect_Table$species.x %in% c("macFas5","hg19") & intersect_Table$liftover.x %in% c("hg19ToMacFas5","macFas5ToHg19"),]
intersect_Table.both$is_sameFamily = ifelse(intersect_Table.both$is_intersect_reciprocal == "no_ortholog","no_ortholog",NA)
intersect_Table.both$is_sameFamily = ifelse(intersect_Table.both$is_intersect_reciprocal != "no_ortholog" & intersect_Table.both$TEfamily.x == intersect_Table.both$TEfamily.y & !is.na(intersect_Table.both$TEfamily.y),"Consistent",intersect_Table.both$is_sameFamily)
intersect_Table.both$is_sameFamily = ifelse(intersect_Table.both$is_intersect_reciprocal != "no_ortholog" & (intersect_Table.both$TEfamily.x != intersect_Table.both$TEfamily.y | is.na(intersect_Table.both$TEfamily.y)),"Inconsistent",intersect_Table.both$is_sameFamily)
intersect_Table.both$is_sameFamily = ifelse(intersect_Table.both$is_sameFamily == "Inconsistent" & !is.na(intersect_Table.both$TEfamily.y),intersect_Table.both$TEfamily.y,intersect_Table.both$is_sameFamily)

# combine the intersected tables
subfamily_info.plot = merge(subfamily_info,intersect_Table.both[,c("instanceID.renamed","is_intersect_reciprocal","liftover.instance","is_sameFamily")],by="instanceID.renamed",all.x=T)

# exclude consensus sequences
subfamily_info.plot = subfamily_info.plot[grepl("^T_",subfamily_info.plot$instanceID.renamed),]

# set non ortholog group
subfamily_info.plot$subfamily.name.final = ifelse(subfamily_info.plot$is_sameFamily == "no_ortholog","no",subfamily_info.plot$subfamily.name.final)

# set non candidate families
subfamily_info.plot$subfamily.name.final = ifelse(is.na(subfamily_info.plot$subfamily.name.final) | grepl("_U$",subfamily_info.plot$subfamily.name.final),"other",subfamily_info.plot$subfamily.name.final)
subfamily_info.plot$is_existed = ifelse(subfamily_info.plot$subfamily.name.final == "no","no","yes")

# perc each subfamily
subfamily_info.plot.sum = data.frame(subfamily_info.plot %>% 
                                       group_by(species,subfamily.name.final.species) %>% 
                                       dplyr::summarise(n = n()) %>% mutate(total = sum(n)) %>% mutate(freq = n/ sum(n)))
# perc.liftover
subfamily_info.plot.tmp = data.frame(subfamily_info.plot[subfamily_info.plot$is_existed == "yes",] %>% 
                                       group_by(subfamily.name.final.species) %>% 
                                       dplyr::summarise(n = n()) %>% mutate(total = sum(n)) %>% mutate(freq = n/ sum(n)))

subfamily_info.plot.sum = merge(subfamily_info.plot.sum,subfamily_info.plot.tmp[,c("subfamily.name.final.species","n")],by="subfamily.name.final.species",all.x=T)

subfamily_info.plot.sum$intersect_perc = 0
subfamily_info.plot.sum$intersect_perc = ifelse(!is.na(subfamily_info.plot.sum$n.y),subfamily_info.plot.sum$n.y/subfamily_info.plot.sum$n.x,0)
subfamily_info.plot.sum = subfamily_info.plot.sum[!is.na(subfamily_info.plot.sum$subfamily.name.final.species) & !grepl("_U$",subfamily_info.plot.sum$subfamily.name.final.species),]

subfamily_info.plot.sum = merge(subfamily_info.plot.sum,MER11_both,by.x="subfamily.name.final.species",by.y="subfamily",all=T)
subfamily_info.plot.sum$uniqueID = paste(subfamily_info.plot.sum$group,subfamily_info.plot.sum$species.x)

my_comparisons = list(c("FG1 hg19","FG1 macFas5"),c("FG2 hg19","FG2 macFas5"),
                      c("FG3 hg19","FG3 macFas5"),c("FG4 hg19","FG4 macFas5"))

############################
# liftover rate
head(subfamily_info.plot.sum)
p1 = ggplot(subfamily_info.plot.sum[!is.na(subfamily_info.plot.sum$group),],aes(x=uniqueID,intersect_perc,group=uniqueID))+
  geom_quasirandom(aes(fill = species.x),method = "quasirandom",alpha=.7,size=2,shape=21) + 
  stat_summary(fun = "mean",
               fun.min = function(x)mean(x)-sd(x),
               fun.max = function(x)mean(x) + sd(x),
               geom = "pointrange",
               position = position_dodge(width = .9),shape=95)+
  ylab("liftover proportion")+
  xlab("Group")+
  scale_x_discrete(drop = FALSE)+
  scale_fill_manual(values=c("hg19"="#2166ac","macFas5"="#66c2a5"))+
  #scale_fill_manual(values=cluster.color)+
  #scale_color_manual(values=cluster.color)+
  stat_compare_means(method = "t.test",comparisons = my_comparisons,aes(label = paste0("", ..p.format..))) +
  #ggtitle(paste(Family))+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    #legend.position="right",
    legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))

# counts
p2 = ggplot(subfamily_info.plot.sum[!is.na(subfamily_info.plot.sum$group),],aes(x=uniqueID,n.x,intersect_perc,group=uniqueID))+
  #geom_jitter(aes(colour = species,fill=species,group=species),alpha=0.5)+
  geom_quasirandom(aes(fill = species.x),method = "quasirandom",alpha=.7,size=2,shape=21) + 
  #geom_violin(aes(group=uniqueID),fill=NA) + #add points colored by significance
  stat_summary(fun = "mean",
               fun.min = function(x)mean(x)-sd(x),
               fun.max = function(x)mean(x) + sd(x),
               geom = "pointrange",
               position = position_dodge(width = .9),shape=95)+
  ylab("liftover count")+
  xlab("Group")+
  scale_x_discrete(drop = FALSE)+
  scale_fill_manual(values=c("hg19"="#2166ac","macFas5"="#66c2a5"))+
  #scale_fill_manual(values=cluster.color)+
  #scale_color_manual(values=cluster.color)+
  stat_compare_means(method = "t.test",comparisons = my_comparisons,aes(label = paste0("", ..p.format..))) +
  #ggtitle(paste(Family))+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    #legend.position="right",
    legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))

###
subfamily_info.plot.sum.tmp = data.frame(subfamily_info.plot.sum %>% 
                                           group_by(species.x,group,uniqueID) %>% 
                                           dplyr::summarise(n = sum(n.x)))

p3 = ggplot(subfamily_info.plot.sum.tmp[!is.na(subfamily_info.plot.sum.tmp$group),], aes(uniqueID, n)) + 
  geom_col(aes(fill=species.x)) + 
  scale_fill_manual(values=c("hg19"="#2166ac","macFas5"="#66c2a5"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1))) 

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
gC <- ggplotGrob(p3)
g2 = cbind(gA,gB,gC,size = "first")
pdf(paste("Figure_4D_4E",".pdf",sep=""),    # create PNG for the heat map
    width = 12,        # 5 x 300 pixels
    height = 4,
    pointsize = 10)        # smaller font size
grid.draw(g2)
dev.off()

#################################
#################################
######################## activity
Summary_Table2.MER11 = merge(Summary_Table2,subfamily_info[,c("instanceID.renamed","subfamily.name.final.species")],by.x="instanceID.renamed",all.x=T)
Summary_Table2.MER11 = Summary_Table2.MER11[!is.na(Summary_Table2.MER11$subfamily.name.final.species),]
Summary_Table2.MER11 = merge(Summary_Table2.MER11,MER11_both,by.x="subfamily.name.final.species",by.y="subfamily",all.x=T)
Summary_Table2.MER11$x_axis = paste(Summary_Table2.MER11$group,Summary_Table2.MER11$species.x)
Summary_Table2.MER11 = merge(Summary_Table2.MER11,subfamily_info.plot[,c("instanceID.renamed","is_existed")],by="instanceID.renamed",all.x=T)

##################### barplot
#############################################
Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.MER11$iPSC.alpha.Zscore.isActive == "None","No activity",NA)
Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.MER11$iPSC.alpha.Zscore.isActive != "None" & 
                                                                 Summary_Table2.MER11$iPSC.alpha.Zscore < 2,"<2", Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.MER11$iPSC.alpha.Zscore.isActive != "None" & 
                                                                 Summary_Table2.MER11$iPSC.alpha.Zscore >= 2 & Summary_Table2.MER11$iPSC.alpha.Zscore < 4 ,"2-4" ,Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.MER11$iPSC.alpha.Zscore.isActive != "None" & 
                                                                 Summary_Table2.MER11$iPSC.alpha.Zscore >= 4 & Summary_Table2.MER11$iPSC.alpha.Zscore < 6 ,"4-6", Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.MER11$iPSC.alpha.Zscore.isActive != "None" & 
                                                                 Summary_Table2.MER11$iPSC.alpha.Zscore >= 6 & Summary_Table2.MER11$iPSC.alpha.Zscore < 8 ,"6-8", Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.MER11$iPSC.alpha.Zscore.isActive != "None" & 
                                                                 Summary_Table2.MER11$iPSC.alpha.Zscore >= 8 , ">=8", Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group = factor(Summary_Table2.MER11$iPSC.alpha.Zscore.isActive.group,levels=c(">=8","6-8","4-6","2-4","<2","No activity"))

### summary
Summary_Table2.MER11.highQuality = Summary_Table2.MER11[grepl("highQuality",Summary_Table2.MER11$is_kept),]
# excluded NA values that refer to low quality in each cell
Summary_Table2.MER11.highQuality = Summary_Table2.MER11.highQuality[!is.na(Summary_Table2.MER11.highQuality$iPSC.alpha.Zscore.isActive.group),]

Alpha.Zscore.isActive.group.color = c(">=8" = "#67001f","6-8"="#b2182b","4-6"="#d6604d","2-4"="#f4a582","<2"="#fddbc7","No activity"="#f0f0f0")
color_alpha.group = c("alpha>=4"="#bd0026","alpha>=2"="#fd8d3c","low"="#f0f0f0")

##################### plot 3: activity distribution (summary)
# iPSC
Summary_Table2.MER11.highQuality_sum = data.frame(Summary_Table2.MER11.highQuality[!is.na(Summary_Table2.MER11.highQuality$group),] %>% group_by(x_axis,Frame,iPSC.alpha.Zscore.isActive.group) %>% dplyr::summarise(n = n()))
Summary_Table2.MER11.highQuality_sum$uniqueName = paste(Summary_Table2.MER11.highQuality_sum$x_axis,Summary_Table2.MER11.highQuality_sum$Frame)
Summary_Table2.MER11.highQuality_sum_tmp = data.frame(Summary_Table2.MER11.highQuality[!is.na(Summary_Table2.MER11.highQuality$group),] %>% group_by(x_axis,Frame) %>% dplyr::summarise(total = n()))
Summary_Table2.MER11.highQuality_sum_tmp$uniqueName = paste(Summary_Table2.MER11.highQuality_sum_tmp$x_axis,Summary_Table2.MER11.highQuality_sum_tmp$Frame)
Summary_Table2.MER11.highQuality_sum = merge(Summary_Table2.MER11.highQuality_sum,Summary_Table2.MER11.highQuality_sum_tmp,by="uniqueName",all=T)
Summary_Table2.MER11.highQuality_sum$perC = Summary_Table2.MER11.highQuality_sum$n/Summary_Table2.MER11.highQuality_sum$total
Summary_Table2.MER11.highQuality_sum$iPSC.alpha.Zscore.isActive.group = factor(Summary_Table2.MER11.highQuality_sum$iPSC.alpha.Zscore.isActive.group,levels=rev(c(">=8","6-8","4-6","2-4","<2","No activity")))
Summary_Table2.MER11.highQuality_sum$x_axis.x = factor(Summary_Table2.MER11.highQuality_sum$x_axis.x)
Summary_Table2.MER11.highQuality_sum$perC = ifelse(Summary_Table2.MER11.highQuality_sum$total<10,NA,Summary_Table2.MER11.highQuality_sum$perC)
Summary_Table2.MER11.highQuality_sum_noActivity = Summary_Table2.MER11.highQuality_sum[Summary_Table2.MER11.highQuality_sum$iPSC.alpha.Zscore.isActive.group == "No activity",]
colnames(Summary_Table2.MER11.highQuality_sum_noActivity)[which(colnames(Summary_Table2.MER11.highQuality_sum_noActivity) == "n")] = "n.noActivity"
Summary_Table2.MER11.highQuality_sum = merge(Summary_Table2.MER11.highQuality_sum,Summary_Table2.MER11.highQuality_sum_noActivity[,c("uniqueName","n.noActivity")],by="uniqueName",all.x=T)
Summary_Table2.MER11.highQuality_sum$n.Activity = Summary_Table2.MER11.highQuality_sum$total - Summary_Table2.MER11.highQuality_sum$n.noActivity
p1<-ggplot(Summary_Table2.MER11.highQuality_sum[Summary_Table2.MER11.highQuality_sum$Frame.x == "F1",],aes(x=x_axis.x,y=perC, fill=iPSC.alpha.Zscore.isActive.group))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = total,y=1))+
  geom_text(aes(label = n.Activity,y=0.5))+
  geom_hline(yintercept = c(0.5), linetype="dotted",color = "white", linewidth=1)+
  scale_fill_manual(values=Alpha.Zscore.isActive.group.color) + 
  scale_x_discrete(drop = FALSE)+
  ggtitle("iPS cells")+
  ylab("Proportion")+
  xlab("")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
p2<-ggplot(Summary_Table2.MER11.highQuality_sum[Summary_Table2.MER11.highQuality_sum$Frame.x == "F2",],aes(x=x_axis.x,y=perC, fill=iPSC.alpha.Zscore.isActive.group))+
  geom_bar(position="stack", stat="identity")+
  geom_hline(yintercept = c(0.5), linetype="dotted",color = "white", linewidth=1)+
  geom_text(aes(label = total,y=1))+
  geom_text(aes(label = n.Activity,y=0.5))+
  scale_fill_manual(values=Alpha.Zscore.isActive.group.color) + 
  scale_x_discrete(drop = FALSE)+
  ggtitle("iPS cells")+
  ylab("Proportion")+
  xlab("")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
g2 = cbind(gA,gB,size = "first")
pdf(paste("Figure_4F",".pdf",sep=""),    # create PNG for the heat map
    width = 12,        # 5 x 300 pixels
    height = 4,
    pointsize = 10)        # smaller font size
grid.draw(g2)
dev.off()