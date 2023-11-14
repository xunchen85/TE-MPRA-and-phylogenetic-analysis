#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

library(splitstackshape)
library(dplyr)
library(ggplot2)
library(gplots)
library(grid)
library(gridExtra)
library(ggridges)
library(ggpmisc)
library(ggbeeswarm)
library(ggpubr)

options(scipen=999)
Date = "2023_4_21"
########################### Table 1: instance info
Summary_Table1 = read.csv(file = "input/Summary_Table1_2022_8_9.csv")

###### how many instances were analyzed
Summary_Table1.sum = data.frame(Summary_Table1[Summary_Table1$species != "Ancient",] %>% group_by(TEfamily,species,frames_tested,F1_tested,F2_tested) 
                                %>% dplyr::summarise(n = n()))
Summary_Table1.sum.len = data.frame(Summary_Table1[Summary_Table1$species != "Ancient" & 
                                                     ((grepl("MER34",Summary_Table1$TEfamily) & Summary_Table1$Instance_len>=250)|
                                                        (!grepl("MER34",Summary_Table1$TEfamily) & Summary_Table1$Instance_len>=500)),] %>% group_by(TEfamily,species,frames_tested,F1_tested,F2_tested) 
                                    %>% dplyr::summarise(n = n()))

# no filtering
Summary_Table1.sum2 = Summary_Table1.sum
Summary_Table1.sum2$Group = paste(Summary_Table1.sum2$F1_tested,Summary_Table1.sum2$F2_tested)
Summary_Table1.sum2$TEfamily2 = paste(Summary_Table1.sum2$TEfamily,Summary_Table1.sum2$species)
Summary_Table1.sum2$Group = ifelse(Summary_Table1.sum2$Group == "NA NA","Not tested",Summary_Table1.sum2$Group)
Summary_Table1.sum2$Group = ifelse(Summary_Table1.sum2$Group == "Tested NA","F1",Summary_Table1.sum2$Group)
Summary_Table1.sum2$Group = ifelse(Summary_Table1.sum2$Group == "NA Tested","F2",Summary_Table1.sum2$Group)
Summary_Table1.sum2$Group = ifelse(Summary_Table1.sum2$Group == "Tested Tested","F1 & F2",Summary_Table1.sum2$Group)
Summary_Table1.sum2$Group = factor(Summary_Table1.sum2$Group,levels = rev(c("F1 & F2", "F1","F2","Not tested")))
#write.csv(Summary_Table1.sum2,file = "Summary_Table1_2022_8_9.supplementary.csv")

######################### plots
#Summary_Table1.sum2 = read.csv("Summary_Table1_2022_8_9.supplementary.csv")
colCombination = c("F1 & F2"="#225ea8","F1"="#7fcdbb","F2"="#edf8b1","Not tested"="#f0f0f0")

# filter by length
Summary_Table1.sum2 = Summary_Table1.sum.len
Summary_Table1.sum2$Group = paste(Summary_Table1.sum2$F1_tested,Summary_Table1.sum2$F2_tested)
Summary_Table1.sum2$TEfamily2 = paste(Summary_Table1.sum2$TEfamily,Summary_Table1.sum2$species)
Summary_Table1.sum2$Group = ifelse(Summary_Table1.sum2$Group == "NA NA","Not tested",Summary_Table1.sum2$Group)
Summary_Table1.sum2$Group = ifelse(Summary_Table1.sum2$Group == "Tested NA","F1",Summary_Table1.sum2$Group)
Summary_Table1.sum2$Group = ifelse(Summary_Table1.sum2$Group == "NA Tested","F2",Summary_Table1.sum2$Group)
Summary_Table1.sum2$Group = ifelse(Summary_Table1.sum2$Group == "Tested Tested","F1 & F2",Summary_Table1.sum2$Group)
Summary_Table1.sum2$Group = factor(Summary_Table1.sum2$Group,levels = rev(c("F1 & F2", "F1","F2","Not tested")))
colCombination = c("F1 & F2"="#225ea8","F1"="#7fcdbb","F2"="#edf8b1","Not tested"="#f0f0f0")

p1 = ggplot(Summary_Table1.sum2[Summary_Table1.sum2$TEfamily !="MER11D",],aes(x=TEfamily2,n,fill=Group))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=colCombination)+
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
    #legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
pdf("Figure_S5C.pdf",    # create PNG for the heat map        
    height = 4,        # 5 x 300 pixels
    width = 8,
    pointsize = 10 )        # smaller font size
grid.draw(p1)
dev.off()

########################## Table 2:
## load table 2 with scaled alpha value
Summary_Table2.final = read.csv("input/Summary_Table2_2022_8_9.csv")
head(Summary_Table2.final)

color_retrieved2=c("highQuality" = "#0868ac","lowQuality"="#a8ddb5","Not retrieved"="#bdbdbd")

################################### plot 3, retrieved proportion
# barplot
Summary_Table2.final.plot = Summary_Table2.final[Summary_Table2.final$Is_included == "Included_MPRA",]
Summary_Table2.final.plot$uniqueName = ifelse(Summary_Table2.final.plot$Group == "Positive","Positive",NA)
Summary_Table2.final.plot$uniqueName = ifelse(Summary_Table2.final.plot$Group == "Negative","Negative",Summary_Table2.final.plot$uniqueName)
Summary_Table2.final.plot$uniqueName = ifelse(Summary_Table2.final.plot$Group == "Consensus","Consensus",Summary_Table2.final.plot$uniqueName)
Summary_Table2.final.plot$uniqueName = ifelse(Summary_Table2.final.plot$Group == "Instance",paste(Summary_Table2.final.plot$Family,Summary_Table2.final.plot$Frame,sep="_"),Summary_Table2.final.plot$uniqueName)
Summary_Table2.final.plot.sum = data.frame(Summary_Table2.final.plot %>% group_by(uniqueName,is_kept) 
                                           %>% dplyr::summarise(n = n()))
Summary_Table2.final.plot.sum$is_kept = as.character(Summary_Table2.final.plot.sum$is_kept)
Summary_Table2.final.plot.sum$is_kept = ifelse(is.na(Summary_Table2.final.plot.sum$is_kept),"Not retrieved",Summary_Table2.final.plot.sum$is_kept)
Summary_Table2.final.plot.sum$is_kept = factor(Summary_Table2.final.plot.sum$is_kept,levels=rev(c("highQuality.iPSC/NPC","highQuality.iPSC","highQuality.NPC","lowQuality","Not retrieved")))

Summary_Table2.final.plot.sum$uniqueName = factor(Summary_Table2.final.plot.sum$uniqueName,levels=c("MER11A_F1","MER11B_F1","MER11C_F1","MER11A_F2","MER11B_F2","MER11C_F2",
                                                                                                    "MER34_F1","MER34A_F1","MER34A1_F1","MER34B_F1","MER34C__F1","MER34C_F1","MER34C2_F1","MER34D_F1",
                                                                                                    "MER52A_F1","MER52C_F1","MER52A_F2","MER52C_F2","MER52D_F2","Consensus","Negative","Positive"))

Summary_Table2.final.plot.sum2 = data.frame(Summary_Table2.final.plot %>% group_by(uniqueName) 
                                            %>% dplyr::summarise(n = n()))

head(Summary_Table2.final.plot)
color_retrieved=c("highQuality.iPSC/NPC" = "#0868ac","highQuality.iPSC" = "#2b8cbe","highQuality.NPC"="#4eb3d3","lowQuality"="#a8ddb5","Not retrieved"="#bdbdbd")
p1 = ggplot(data = Summary_Table2.final.plot.sum, aes(x=uniqueName, y=n)) +
  geom_bar(position="fill", stat="identity",aes(fill=is_kept))+
  scale_fill_manual(values=color_retrieved)+
  geom_text(data = Summary_Table2.final.plot.sum2,
            aes(x = uniqueName, group=1, y = 1 - 0.08,
                label = n),position=position_dodge(.9), angle=90,hjust=.5)+
  ylab("Proportion of inserts")+
  xlab("Frame")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.x.top = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y.right = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5,hjust = 0.95,angle = 90),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))

pdf("Figure_S5D.pdf",    # create PNG for the heat map        
    height = 4,        # 5 x 300 pixels
    width = 8,
    pointsize = 10 )        # smaller font size
grid.draw(p1)
dev.off()

################################### plot 4, activity score
## iPSC
Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive == "None","No activity",NA)
Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive != "None" & 
                                                                      Summary_Table2.final.plot$iPSC.alpha.Zscore < 2,"<2", Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive != "None" & 
                                                                      Summary_Table2.final.plot$iPSC.alpha.Zscore >= 2 & Summary_Table2.final.plot$iPSC.alpha.Zscore < 4 ,"2-4" ,Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive != "None" & 
                                                                      Summary_Table2.final.plot$iPSC.alpha.Zscore >= 4 & Summary_Table2.final.plot$iPSC.alpha.Zscore < 6 ,"4-6", Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive != "None" & 
                                                                      Summary_Table2.final.plot$iPSC.alpha.Zscore >= 6 & Summary_Table2.final.plot$iPSC.alpha.Zscore < 8 ,"6-8", Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive != "None" & 
                                                                      Summary_Table2.final.plot$iPSC.alpha.Zscore >= 8 , ">=8", Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group = factor(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group,levels=c(">=8","6-8","4-6","2-4","<2","No activity"))
head(Summary_Table2.final.plot[is.na(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group),],20)
unique(Summary_Table2.final.plot$iPSC.alpha.Zscore.isActive.group)
## NPC
Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$NPC.alpha.Zscore.isActive == "None","No activity",NA)
Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$NPC.alpha.Zscore.isActive != "None" & 
                                                                     Summary_Table2.final.plot$NPC.alpha.Zscore < 2,"<2", Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$NPC.alpha.Zscore.isActive != "None" & 
                                                                     Summary_Table2.final.plot$NPC.alpha.Zscore >= 2 & Summary_Table2.final.plot$NPC.alpha.Zscore < 4 ,"2-4" ,Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$NPC.alpha.Zscore.isActive != "None" & 
                                                                     Summary_Table2.final.plot$NPC.alpha.Zscore >= 4 & Summary_Table2.final.plot$NPC.alpha.Zscore < 6 ,"4-6", Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$NPC.alpha.Zscore.isActive != "None" & 
                                                                     Summary_Table2.final.plot$NPC.alpha.Zscore >= 6 & Summary_Table2.final.plot$NPC.alpha.Zscore < 8 ,"6-8", Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group = ifelse(Summary_Table2.final.plot$NPC.alpha.Zscore.isActive != "None" & 
                                                                     Summary_Table2.final.plot$NPC.alpha.Zscore >= 8 , ">=8", Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group = factor(Summary_Table2.final.plot$NPC.alpha.Zscore.isActive.group,levels=c(">=8","6-8","4-6","2-4","<2","No activity"))
#
Alpha.Zscore.isActive.group.color = c(">=8" = "#67001f","6-8"="#b2182b","4-6"="#d6604d","2-4"="#f4a582","<2"="#fddbc7","No activity"="#f0f0f0")
color_alpha.group = c("alpha>=4"="#bd0026","alpha>=2"="#fd8d3c","low"="#f0f0f0")

## only kept the high quality ones in the summary plot
# excluded both low quality in iPSC and NPC
Summary_Table2.final.plot.highQuality = Summary_Table2.final.plot[grepl("highQuality",Summary_Table2.final.plot$is_kept),]
# excluded NA values that refer to low quality in each cell
Summary_Table2.final.plot.highQuality = Summary_Table2.final.plot.highQuality[!is.na(Summary_Table2.final.plot.highQuality$iPSC.alpha.Zscore.isActive),]

Summary_Table2.final.plot.highQuality$iPSC.alpha.Zscore.isActive.group = factor(Summary_Table2.final.plot.highQuality$iPSC.alpha.Zscore.isActive.group)
Summary_Table2.final.plot.highQuality$NPC.alpha.Zscore.isActive.group = factor(Summary_Table2.final.plot.highQuality$NPC.alpha.Zscore.isActive.group)

# iPSC
Summary_Table2.final.plot.highQuality_sum = data.frame(Summary_Table2.final.plot.highQuality %>% group_by(uniqueName,iPSC.alpha.Zscore.isActive.group) %>% dplyr::summarise(n = n()))
Summary_Table2.final.plot.highQuality_sum_tmp = data.frame(Summary_Table2.final.plot.highQuality %>% group_by(uniqueName) %>% dplyr::summarise(total = n()))
Summary_Table2.final.plot.highQuality_sum= merge(Summary_Table2.final.plot.highQuality_sum,Summary_Table2.final.plot.highQuality_sum_tmp,by="uniqueName",all.x=T)
Summary_Table2.final.plot.highQuality_sum$perC = Summary_Table2.final.plot.highQuality_sum$n/Summary_Table2.final.plot.highQuality_sum$total
Summary_Table2.final.plot.highQuality_sum = Summary_Table2.final.plot.highQuality_sum[!(Summary_Table2.final.plot.highQuality_sum$uniqueName %in% c("Consensus","Negative")),]
Summary_Table2.final.plot.highQuality_sum$uniqueName = factor(Summary_Table2.final.plot.highQuality_sum$uniqueName,levels=c("MER11A_F1","MER11B_F1","MER11C_F1","MER11A_F2","MER11B_F2","MER11C_F2",
                                                                                                                            "MER34_F1","MER34A_F1","MER34A1_F1","MER34B_F1","MER34C__F1","MER34C_F1","MER34C2_F1","MER34D_F1",
                                                                                                                            "MER52A_F1","MER52C_F1","MER52A_F2","MER52C_F2","MER52D_F2","Positive"))
Summary_Table2.final.plot.highQuality_sum$iPSC.alpha.Zscore.isActive.group = factor(Summary_Table2.final.plot.highQuality_sum$iPSC.alpha.Zscore.isActive.group,levels=rev(c(">=8","6-8","4-6","2-4","<2","No activity")))
p1<-ggplot(Summary_Table2.final.plot.highQuality_sum,aes(x=uniqueName,y=perC, fill=iPSC.alpha.Zscore.isActive.group))+
  geom_bar(position="stack", stat="identity")+
  geom_hline(yintercept = c(0.5), linetype="dotted",color = "white", linewidth=1)+
  scale_fill_manual(values=Alpha.Zscore.isActive.group.color) + 
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

Summary_Table2.final.plot.highQuality_sum = data.frame(Summary_Table2.final.plot.highQuality %>% group_by(uniqueName,NPC.alpha.Zscore.isActive.group) %>% dplyr::summarise(n = n()))
Summary_Table2.final.plot.highQuality_sum_tmp = data.frame(Summary_Table2.final.plot.highQuality %>% group_by(uniqueName) %>% dplyr::summarise(total = n()))
Summary_Table2.final.plot.highQuality_sum= merge(Summary_Table2.final.plot.highQuality_sum,Summary_Table2.final.plot.highQuality_sum_tmp,by="uniqueName",all.x=T)
Summary_Table2.final.plot.highQuality_sum$perC = Summary_Table2.final.plot.highQuality_sum$n/Summary_Table2.final.plot.highQuality_sum$total

# exclude negative and consensus group
Summary_Table2.final.plot.highQuality_sum = Summary_Table2.final.plot.highQuality_sum[!(Summary_Table2.final.plot.highQuality_sum$uniqueName %in% c("Consensus","Negative")),]
Summary_Table2.final.plot.highQuality_sum$uniqueName = factor(Summary_Table2.final.plot.highQuality_sum$uniqueName,levels=c("MER11A_F1","MER11B_F1","MER11C_F1","MER11A_F2","MER11B_F2","MER11C_F2",
                                                                                                                            "MER34_F1","MER34A_F1","MER34A1_F1","MER34B_F1","MER34C__F1","MER34C_F1","MER34C2_F1","MER34D_F1",
                                                                                                                            "MER52A_F1","MER52C_F1","MER52A_F2","MER52C_F2","MER52D_F2","Positive"))
Summary_Table2.final.plot.highQuality_sum$NPC.alpha.Zscore.isActive.group = factor(Summary_Table2.final.plot.highQuality_sum$NPC.alpha.Zscore.isActive.group,levels=rev(c(">=8","6-8","4-6","2-4","<2","No activity")))
p2<-ggplot(Summary_Table2.final.plot.highQuality_sum,aes(x=uniqueName,y=perC, fill=NPC.alpha.Zscore.isActive.group))+
  geom_bar(position="stack", stat="identity")+
  geom_hline(yintercept = c(0.5), linetype="dotted",color = "white", linewidth=1)+
  scale_fill_manual(values=Alpha.Zscore.isActive.group.color) + 
  ggtitle("NP cells")+
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
g = cbind(gA,gB, size = "first")

pdf("Figure_3D.pdf",    # create PNG for the heat map        
    width = 16,
    height = 5,
    pointsize = 10)        # smaller font size
grid.draw(g)
dev.off()
