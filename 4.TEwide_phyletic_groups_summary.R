#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

library(ggplot2)
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

Date = "2023_9_11"
##### step 1: subfamily and functional group
phyletic_group_table = read.csv("input/TEwide_group.list_2023_9_4.functional_group_edited.csv")
phyletic_group_table$uniqueID = paste(phyletic_group_table$tree,phyletic_group_table$family.cluster,phyletic_group_table$subfamily2_1)

##### step 2: summary by each subfamily and phyletic group
phyletic_group_table.subfamily.total = data.frame(phyletic_group_table %>% 
                                                    group_by(tree,family.cluster,subfamily2_1,family.count_200bp) %>% 
                                                    dplyr::summarise(family.count_200bp.clustered = sum(subfamily.count_200bp,na.rm=T)))
phyletic_group_table.subfamily.total$uniqueID = paste(phyletic_group_table.subfamily.total$tree,phyletic_group_table.subfamily.total$family.cluster,phyletic_group_table.subfamily.total$subfamily2_1)

phyletic_group_table.subfamily = data.frame(phyletic_group_table %>% 
                                              group_by(tree,family.cluster,subfamily2_1,family.count_200bp,functional.group) %>% 
                                              dplyr::summarise(subfamily.count_200bp.per_FG = sum(subfamily.count_200bp,na.rm=T)))
phyletic_group_table.subfamily$uniqueID = paste(phyletic_group_table.subfamily$tree,phyletic_group_table.subfamily$family.cluster,phyletic_group_table.subfamily$subfamily2_1)
phyletic_group_table.subfamily = merge(phyletic_group_table.subfamily,phyletic_group_table.subfamily.total[,c("uniqueID","family.count_200bp.clustered")],by="uniqueID",all.x=T)
phyletic_group_table.subfamily$perC.per_FG = phyletic_group_table.subfamily$subfamily.count_200bp.per_FG/phyletic_group_table.subfamily$family.count_200bp.clustered
phyletic_group_table.subfamily = phyletic_group_table.subfamily[order(-phyletic_group_table.subfamily$perC.per_FG),]
phyletic_group_table.subfamily$uniqueID.topFG = paste(phyletic_group_table.subfamily$uniqueID,phyletic_group_table.subfamily$functional.group)
phyletic_group_table.subfamily.top = phyletic_group_table.subfamily[!duplicated(phyletic_group_table.subfamily$uniqueID),]
phyletic_group_table.subfamily.top$uniqueID.type2 = paste(phyletic_group_table.subfamily.top$tree,phyletic_group_table.subfamily.top$functional.group)

phyletic_group_table.subfamily$is_top = ifelse(phyletic_group_table.subfamily$uniqueID.topFG %in% phyletic_group_table.subfamily.top$uniqueID.topFG,"yes1","type1")

#####
group_colors = c("PG1"="#0692C9","PG2"="#8638F4","PG3"="#06894A","PG4"="#608E06")

p1 = ggplot(phyletic_group_table.subfamily[phyletic_group_table.subfamily$family.cluster == "MER11A",], aes(x=subfamily2_1,
                                                                                                            y=subfamily.count_200bp.per_FG,group=functional.group)) +
  geom_col(aes(fill=functional.group),alpha=0.7) + 
  geom_text(aes(label = subfamily.count_200bp.per_FG), size = 5, hjust = 0.5, vjust = 0.5, position = "stack")+
  scale_fill_manual(values =group_colors)+
  scale_x_discrete(drop = FALSE)+
  ylab("# instances")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.line.x = element_blank(),
        #axis.line.y = element_blank(),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.95),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1))) 

pdf(paste("Figure_1G",".pdf",sep=""),    # create PNG for the heat map
    width = 4,        # 5 x 300 pixels
    height = 4,
    pointsize = 10)        # smaller font size
grid.draw(p1)
dev.off()



#####
# type 1: proportion of sequences annotated as the top phyetic and non-top phyletic group
phyletic_group_table.type1 = data.frame(phyletic_group_table.subfamily %>% 
                                          group_by(uniqueID,tree,family.cluster,subfamily2_1,family.count_200bp,family.count_200bp.clustered,is_top) %>% 
                                          dplyr::summarise(subfamily.count_200bp.type1_2 = sum(subfamily.count_200bp.per_FG,na.rm=T)))
head(phyletic_group_table.type1)

# type 2: proportion of sequences from other subfamilies annotated as the the this top phyletic group
phyletic_group_table.group = data.frame(phyletic_group_table %>% 
                                          group_by(tree,family.cluster,functional.group,subfamily2_1) %>% 
                                          dplyr::summarise(subfamily.count_200bp.FG = sum(subfamily.count_200bp,na.rm=T)))
phyletic_group_table.group$uniqueID.topFG = paste(phyletic_group_table.group$tree,phyletic_group_table.group$family.cluster,phyletic_group_table.group$subfamily2_1,phyletic_group_table.group$functional.group)
phyletic_group_table.group$uniqueID = paste(phyletic_group_table.group$tree,phyletic_group_table.group$family.cluster,phyletic_group_table.group$subfamily2_1)
phyletic_group_table.group$is_top = ifelse(phyletic_group_table.group$uniqueID.topFG %in% phyletic_group_table.subfamily[phyletic_group_table.subfamily$is_top == "yes1",]$uniqueID.topFG,"yes2","type2")
phyletic_group_table.type2 = data.frame(phyletic_group_table.group %>% 
                                          group_by(tree,family.cluster,functional.group,is_top) %>% 
                                          dplyr::summarise(subfamily.count_200bp.type1_2 = sum(subfamily.count_200bp.FG,na.rm=T)))
# retrieve the corresponding subfamily
phyletic_group_table.type2$uniqueID.type2 = paste(phyletic_group_table.type2$tree,phyletic_group_table.type2$functional.group)
phyletic_group_table.type2 = phyletic_group_table.type2[phyletic_group_table.type2$uniqueID.type2 %in% phyletic_group_table.subfamily.top$uniqueID.type2,]
phyletic_group_table.type2 = merge(phyletic_group_table.type2,phyletic_group_table.subfamily.top[,c("uniqueID.type2","subfamily2_1")],by="uniqueID.type2",all.x=T)
phyletic_group_table.type2$uniqueID = paste(phyletic_group_table.type2$tree,phyletic_group_table.type2$family.cluster,phyletic_group_table.type2$subfamily2_1)
##### combined the table
phyletic_group_table.error.total1 = phyletic_group_table[!duplicated(phyletic_group_table$uniqueID),c("uniqueID","tree","family.cluster","subfamily2_1")]
phyletic_group_table.error.total1$error.type = "type1"
phyletic_group_table.error.total2 = phyletic_group_table[!duplicated(phyletic_group_table$uniqueID),c("uniqueID","tree","family.cluster","subfamily2_1")]
phyletic_group_table.error.total2$error.type = "type2"
phyletic_group_table.error.total = rbind(phyletic_group_table.error.total1,phyletic_group_table.error.total2)
phyletic_group_table.error.total$uniqueID.error.type = paste(phyletic_group_table.error.total$uniqueID,phyletic_group_table.error.total$error.type)
phyletic_group_table.error = rbind(phyletic_group_table.type1[,c("uniqueID","is_top","subfamily.count_200bp.type1_2")],
                                   phyletic_group_table.type2[,c("uniqueID","is_top","subfamily.count_200bp.type1_2")])
phyletic_group_table.error$uniqueID.error.type = paste(phyletic_group_table.error$uniqueID,phyletic_group_table.error$is_top)
phyletic_group_table.error = phyletic_group_table.error[phyletic_group_table.error$is_top %in% c("type1","type2"),]
phyletic_group_table.error.total = merge(phyletic_group_table.error.total,phyletic_group_table.error[,c("uniqueID.error.type","subfamily.count_200bp.type1_2")],by.x="uniqueID.error.type",all=T)
phyletic_group_table.error.total$subfamily.count_200bp.type1_2.value = ifelse(is.na(phyletic_group_table.error.total$subfamily.count_200bp.type1_2),0,phyletic_group_table.error.total$subfamily.count_200bp.type1_2)

##### ordered by the number of phyletic groups per cluster
phyletic_group_table.sum = data.frame(phyletic_group_table[!duplicated(phyletic_group_table[,c("tree","functional.group")]),] %>% 
                                        group_by(tree) %>% 
                                        dplyr::summarise(FG.count = n()))

phyletic_group_table.error.total = merge(phyletic_group_table.error.total,phyletic_group_table.sum,by="tree",all.x=T)
head(phyletic_group_table.error.total)
phyletic_group_table.error.total = phyletic_group_table.error.total[with(phyletic_group_table.error.total, order(-FG.count,family.cluster, -subfamily.count_200bp.type1_2.value)),]
phyletic_group_table.error.total$subfamily2_1 = factor(phyletic_group_table.error.total$subfamily2_1,levels=unique(phyletic_group_table.error.total$subfamily2_1))

##### plot
p2 = ggplot(phyletic_group_table.error.total, aes(x=subfamily2_1,y=subfamily.count_200bp.type1_2.value,group=error.type)) +
  geom_bar(aes(fill=error.type),stat="identity", position=position_dodge()) + 
  geom_text(aes(label=subfamily.count_200bp.type1_2), vjust=-.2, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values =c("#016A70","#A2C579"))+
  #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
  scale_x_discrete(drop = FALSE)+
  ylab("# instances")+
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

pdf(paste("Figure_S11A",".pdf",sep=""),    # create PNG for the heat map
    width = 14,        # 5 x 300 pixels
    height = 5,
    pointsize = 10)        # smaller font size
grid.draw(p2)
dev.off()

##### summary 
# total analyzed
sum(phyletic_group_table.subfamily.top$family.count_200bp)
# total clustered
sum(phyletic_group_table.subfamily.top$family.count_200bp.clustered)

# type 1
sum(phyletic_group_table.error.total[phyletic_group_table.error.total$error.type == "type1",]$subfamily.count_200bp.type1_2.value)

# 26 subfamilies list
subfamlies = phyletic_group_table.error.total[phyletic_group_table.error.total$error.type == "type1" & phyletic_group_table.error.total$subfamily.count_200bp.type1_2.value>0,]

# type 2 (minus the redundant number)
sum(phyletic_group_table.error.total[phyletic_group_table.error.total$error.type == "type2",]$subfamily.count_200bp.type1_2.value)-294-30

# proportion across the 26 subfamilies
sum(phyletic_group_table.error.total[phyletic_group_table.error.total$error.type == "type1",]$subfamily.count_200bp.type1_2.value)/sum(phyletic_group_table.subfamily.top[phyletic_group_table.subfamily.top$subfamily2_1 %in% subfamlies$subfamily2_1,]$family.count_200bp)

