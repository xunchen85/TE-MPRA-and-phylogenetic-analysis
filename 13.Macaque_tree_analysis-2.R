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

Date = "2023_9_9"

######################### step 1 achieve the list of tree files
# load the fasta
fasta.plot.all = readDNAStringSet("input_trees_consensus/macFas5_MER11ABC_v1.fa")
subfamily.macaque.new = names(fasta.plot.all)
sequence = paste(fasta.plot.all)
fasta.df.all <- data.frame(subfamily.macaque.new, sequence)

# subfamily info of human
info_subfamilies = read.csv("input/MER11_subtree_info_2022_12_27.csv")

# intersect table between human and macaque (same step as another script)
intersect_Table = read.csv("input/Summary_Table1.intersect_reciprocal_hu_ch_ma_2023_2_20.csv")
intersect_Table.macFas5 = intersect_Table[intersect_Table$species.x == "macFas5",]
intersect_Table.macFas5$is_sameFamily = ifelse(intersect_Table.macFas5$is_intersect_reciprocal == "no_ortholog","no_ortholog",NA)
intersect_Table.macFas5$is_sameFamily = ifelse(intersect_Table.macFas5$is_intersect_reciprocal != "no_ortholog" & intersect_Table.macFas5$TEfamily.x == intersect_Table.macFas5$TEfamily.y & !is.na(intersect_Table.macFas5$TEfamily.y),"Consistent",intersect_Table.macFas5$is_sameFamily)
intersect_Table.macFas5$is_sameFamily = ifelse(intersect_Table.macFas5$is_intersect_reciprocal != "no_ortholog" & (intersect_Table.macFas5$TEfamily.x != intersect_Table.macFas5$TEfamily.y | is.na(intersect_Table.macFas5$TEfamily.y)),"Inconsistent",intersect_Table.macFas5$is_sameFamily)
intersect_Table.macFas5$is_sameFamily = ifelse(intersect_Table.macFas5$is_sameFamily == "Inconsistent" & !is.na(intersect_Table.macFas5$TEfamily.y),intersect_Table.macFas5$TEfamily.y,intersect_Table.macFas5$is_sameFamily)

# obtain macFas5 subfamily info
subfamily_info.macFas5 = read.csv("input_trees/macFas5_rmsk_TE_MER11A_0bp_10_0.02_all_95_group_final.csv")
subfamily_info.macFas5.unique = subfamily_info.macFas5[!duplicated(subfamily_info.macFas5$label.final.correctedName),]

# obtain hg19 subfamily info
Summary_Table1 = read.csv("input/Summary_Table1_2023_9_5.subfamilyInfo.csv")
Summary_Table2 = read.csv("input/Summary_Table2_2023_1_5.subfamilyInfo.csv")

# obtain hg19 consensus tree order
# updated on 20293/9/9
#subfamily_info.hg19.order = read.csv("../Final_subfamilies_2023_2_17/Step14_1_hg19_MER11ABC_combined.consensus2.mafft.prank.best.fas.nonrev_dna_rerooted.csv")
subfamily_info.hg19.order.tmp = read.csv("input/TEwide_group.list_2023_9_4.functional_group_edited.csv")
subfamily_info.hg19.order.tmp = subfamily_info.hg19.order.tmp[subfamily_info.hg19.order.tmp$subfamily2_1 %in% c("MER11A","MER11B","MER11C","MER11D"),]
subfamily_info.hg19.order = merge(info_subfamilies[info_subfamilies$Family %in% c("11A","11B","11C"),],subfamily_info.hg19.order.tmp[,c("subfamily.species","Order","branch.length.ToKeptNode")],by.x="subfamily.name.final",by.y="subfamily.species",all.x=T)
subfamily_info.hg19.order$subfamily.name.final2 = subfamily_info.hg19.order$subfamily.name.final
subfamily_info.hg19.order = subfamily_info.hg19.order[order(subfamily_info.hg19.order$Order),]

# combine the intersected tables
intersect_Table.macFas5 = merge(intersect_Table.macFas5,Summary_Table1[,c("instanceID.renamed","subfamily.name.final")],by.y="instanceID.renamed",by.x="liftover.instance",all.x=T)

## load the divergent rates table computed by Zicong
#Div_paired = read.delim("div/RAxML_distances.MER11_3_2023_4_25_originalName",sep="",header=F)
# updated on 2023/9/9
Div_paired = read.delim("input_trees_consensus/RAxML_distances.hg19_macFas5_MER11ABC_v2",header=F,sep="")
Div_paired$is_kept = ifelse((grepl("MER11A_g",Div_paired$V1) & !grepl("MER11A_g",Div_paired$V2)) | 
                              (grepl("MER11A_g",Div_paired$V2) & !grepl("MER11A_g",Div_paired$V1)),"kept_h_m",NA)
Div_paired$is_kept = ifelse((grepl("MER11A_g",Div_paired$V1) & grepl("MER11A_g",Div_paired$V2)),"kept_m",Div_paired$is_kept)
Div_paired_macaque = Div_paired[!is.na(Div_paired$is_kept) & Div_paired$is_kept == "kept_m",]

Div_paired = Div_paired[!is.na(Div_paired$is_kept) & Div_paired$is_kept == "kept_h_m",]
Div_paired$subfamily.macaque = ifelse(grepl("MER11A_g",Div_paired$V1),Div_paired$V1,Div_paired$V2)
Div_paired$subfamily.human = ifelse(!grepl("MER11A_g",Div_paired$V1),Div_paired$V1,Div_paired$V2)
Div_paired$subfamily.human = gsub("^h\\|","",Div_paired$subfamily.human)

# correct the subfamily name
fasta.plot.original = readDNAStringSet("input_trees_consensus/macFas5_MER11ABC_v1.fa")
seq_name.original = names(fasta.plot.original)
sequence.original = paste(fasta.plot.original)
fasta.df.original <- data.frame(seq_name.original, sequence.original)
fasta.df.original = merge(fasta.df.original,fasta.df.all,by.y="sequence",by.x="sequence.original",all.x=T)
Div_paired = merge(Div_paired,fasta.df.original[,c("seq_name.original","subfamily.macaque.new")],by.x="subfamily.macaque",by.y="seq_name.original",all.x=T)

Div_paired = merge(Div_paired,info_subfamilies[,c("consensus.name","subfamily.name.final")],by.x="subfamily.human",by.y="consensus.name",all.x=T)
Div_paired$subfamily.name.final = ifelse(is.na(Div_paired$subfamily.name.final),Div_paired$subfamily.human,Div_paired$subfamily.name.final)

# add 11D as well
# Div_paired$subfamily.name.final = factor(Div_paired$subfamily.name.final,levels=c(subfamily_info.hg19.order$subfamily.name.final2,c("11D_a","11D_b","11D_c","MER11D")))


#############################################
##################### plot 1
# load the non rerooted tree including all sequences 
Tree = "macFas5_MER11ABC_v2.mafft.prank.best.fas.rev_dna"
tr <- read.tree(paste("input_trees_consensus/",Tree,".contree",sep=""))
tr$node.label = paste("Node",1:tr$Nnode,"/",tr$node.label,sep="")

### plot
# tree plot
g1 <- ggtree(tr,branch.length = "branch.length",root.position = 0) + 
  geom_tiplab(align=TRUE) + hexpand(.01) + 
  theme_tree2()

g2 = ggtree(tr, layout = 'equal_angle',size=0.03)+geom_tiplab()+geom_tippoint()

pdf(paste("Figure_S9B-",Tree,"_unrooted",Date,".pdf",sep=""),    # create PNG for the heat map
    width = 12,        # 5 x 300 pixels
    height = 12,
    pointsize = 10)        # smaller font size
grid.draw(g2)
dev.off()

#############################################
##################### plot 2
####
dat = read.csv("input_trees_consensus/macFas5_MER11ABC_v1.mafft.prank.best.fas.nonrev_dna.roottest.csv",skip = 12,header=T)
dat$rank = 1:nrow(dat)
tr.all <- read.tree(paste("input_trees_consensus/macFas5_MER11ABC_v1.mafft.prank.best.fas.nonrev_dna.roottest.trees",sep=""))
#for (k in dat$rank) {
  k=16
  tr = tr.all[[k]]
  tr$node.label = paste("Node",1:tr$Nnode,"/",tr$node.label,sep="")
  
  # actual order
  tr2 <- ladderize(tr, right = FALSE)
  is_tip <- tr2$edge[,2] <= length(tr2$tip.label)
  ordered_tips <- tr2$edge[is_tip, 2]
  ordered_subfamilies.table.tmp = data.frame(tr2$tip.label[ordered_tips])
  colnames(ordered_subfamilies.table.tmp) = "subfamily"
  
  #
  g <- ggtree(tr,branch.length = "branch.length",root.position = 0) + 
    #geom_tiplab() + 
    geom_tiplab(align=TRUE) + 
    geom_text(aes(x=branch, label=round(branch.length,4), color="red",vjust=-.5))+
    ggtitle(paste("#:",nrow(dat),"(",k,");ID:",dat[dat$rank == k,]$ID,";",dat[dat$rank == k,]$p.AU,sep=""))+
    #geom_text(aes(x=branch,label=round(branch.length,5)), hjust=-.1, size=2, color="black") +
    hexpand(.01) + 
    geom_tippoint(size=2)+
    theme_tree2()+
    theme(title=element_text(size=20,face="bold"))
  
  # divergent rate
  Div_paired_macaque2 = Div_paired_macaque[!duplicated(Div_paired_macaque$V2),]
  Div_paired_macaque2$V1 = Div_paired_macaque2$V2
  Div_paired_macaque2$V3 = 0
  Div_paired_macaque2 = data.frame("V1"=tr$tip.label,"V2"=tr$tip.label,"V3"=0,"is_kept"="kept_m")
  Div_paired_macaque3 = Div_paired_macaque[,c("V2","V1","V3","is_kept")]
  colnames(Div_paired_macaque3) = colnames(Div_paired_macaque)
  Div_paired_macaque.plot = rbind(Div_paired_macaque,Div_paired_macaque2,Div_paired_macaque3)
  Div_paired_macaque.plot = Div_paired_macaque.plot[Div_paired_macaque.plot$V1 %in% tr$tip.label & 
                                                      Div_paired_macaque.plot$V2 %in% tr$tip.label,]
  Div_paired_macaque.plot$V1 = factor(Div_paired_macaque.plot$V1,levels=(ordered_subfamilies.table.tmp$subfamily))
  Div_paired_macaque.plot$V2 = factor(Div_paired_macaque.plot$V2,levels=(ordered_subfamilies.table.tmp$subfamily))
  
  Div_paired_macaque.plot[Div_paired_macaque.plot$V1 == Div_paired_macaque.plot$V2,]
  ##
  g2 <- ggplot(Div_paired_macaque.plot, aes(y=V2, x=V1)) + 
    geom_tile(aes(fill=V3)) + 
    scale_x_discrete(drop = FALSE)+
    scale_y_discrete(drop = FALSE)+
    scale_fill_distiller("My Scale", palette = "Spectral", trans = "reverse")+
    xlab(NULL) + ylab(NULL) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "#f0f0f0"),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text.y=element_text(colour="black",size=rel(1)),
          axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
          axis.title=element_text(colour="black",size=rel(1)),
          legend.position="right",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1)))
  p <- g2 %>% insert_left(g)
  pdf(paste("Figure_S9C-","macFas5","-",k,"_similarity-1.pdf",sep=""),    # create PNG for the heat map
      width = 20,        # 5 x 300 pixels
      height = 10,
      pointsize = 10)        # smaller font size
  grid.draw(p)
  dev.off()
#}

# load the tree
#for (k in 16) {
  k = 16
  tr = tr.all[[k]]
  tr$node.label = paste("Node",1:tr$Nnode,"/",tr$node.label,sep="")
  
  ## 2.3 organize into the table
  label_tabel_originalOrder = data.frame(tr$edge)
  label_tabel_originalOrder$edgeOrder = 1:nrow(label_tabel_originalOrder)
  label_tabel_originalOrder$edgeID = paste(label_tabel_originalOrder$X1,label_tabel_originalOrder$X2)
  label_table = as_tibble(tr)
  label_table$Order = 1:nrow(label_table)
  label_table$cluster = NA
  label_table$collapse = NA
  label_table$edgeID = paste(label_table$parent,label_table$node)
  label_table = merge(label_table,label_tabel_originalOrder[,c("edgeID","edgeOrder")],by="edgeID",all.x=T)
  label_table = label_table[order(label_table$Order),]
  rm(label_tabel_originalOrder)
  label_table = data.frame(label_table)
  
  ###### tree order and distance
  # distance
  tr.dis = data.frame(dist.nodes(tr))
  colnames(tr.dis) = label_table$label
  row.names(tr.dis) = label_table$label
  tr.dis$subfamily = rownames(tr.dis)
  tr.dis = tr.dis[rownames(tr.dis) %in% label_table$label,c("subfamily","Node1/")]
  tr.dis = tr.dis[!grepl("^Node",tr.dis$subfamily),]
  
  # actual order
  tr2 <- ladderize(tr, right = FALSE)
  is_tip <- tr2$edge[,2] <= length(tr2$tip.label)
  ordered_tips <- tr2$edge[is_tip, 2]
  ordered_subfamilies = tr2$tip.label[ordered_tips]
  
  # retrieve the subfamily names and liftover table
  subfamily_info.macFas5.plot = merge(subfamily_info.macFas5,intersect_Table.macFas5[,c("instanceID.renamed","is_intersect_reciprocal","liftover.instance","is_sameFamily","subfamily.name.final")],by.x="label",by.y="instanceID.renamed",all.x=T)
  # exclude consensus sequences
  subfamily_info.macFas5.plot = subfamily_info.macFas5.plot[grepl("^T_",subfamily_info.macFas5.plot$label),]
  # set non ortholog group
  subfamily_info.macFas5.plot$subfamily.name.final = ifelse(subfamily_info.macFas5.plot$is_sameFamily == "no_ortholog","no",subfamily_info.macFas5.plot$subfamily.name.final)
  # set non candidate families
  subfamily_info.macFas5.plot$subfamily.name.final = ifelse(is.na(subfamily_info.macFas5.plot$subfamily.name.final) | grepl("_U$",subfamily_info.macFas5.plot$subfamily.name.final),"other",subfamily_info.macFas5.plot$subfamily.name.final)
  subfamily_info.macFas5.plot$is_existed = ifelse(subfamily_info.macFas5.plot$subfamily.name.final == "no","no","yes")
  
  # perc each subfamily
  subfamily_info.macFas5.plot.sum = data.frame(subfamily_info.macFas5.plot[subfamily_info.macFas5.plot$subfamily.name.final!="no",] %>% 
                                                 group_by(label.final.correctedName,subfamily.name.final) %>% 
                                                 dplyr::summarise(n = n()) %>% mutate(total = sum(n)) %>% mutate(freq = n/ sum(n)))
  # perc.liftover
  subfamily_info.macFas5.plot.sum.tmp = data.frame(subfamily_info.macFas5.plot %>% 
                                                     group_by(label.final.correctedName,is_existed) %>% 
                                                     dplyr::summarise(n = n()) %>% mutate(total = sum(n)) %>% mutate(freq = n/ sum(n)))
  subfamily_info.macFas5.plot.sum.tmp2 = subfamily_info.macFas5.plot.sum.tmp[subfamily_info.macFas5.plot.sum.tmp$label.final.correctedName %in% ordered_subfamilies,]
  
  ### plot
  # tree plot
  g <- ggtree(tr,branch.length = "branch.length",root.position = 0) + 
    geom_tiplab(align=TRUE) + hexpand(.01) + 
    theme_tree2()
  # pdf(paste("Step19_2_",Tree,"_rerooted.pdf",sep=""),    # create PNG for the heat map
  #     width = 3,        # 5 x 300 pixels
  #     height = 6,
  #     pointsize = 10)        # smaller font size
  # grid.draw(g)
  # dev.off()
  
  # barplot
  subfamily_info.macFas5.plot.sum.tmp$label.final.correctedName = factor(subfamily_info.macFas5.plot.sum.tmp$label.final.correctedName,levels=c(ordered_subfamilies,c("MER11A_g4","MER11A_g12","MER11A_U")))
  # p1 <- ggplot(subfamily_info.macFas5.plot.sum.tmp, aes(label.final.correctedName, freq)) + 
  #   geom_col(aes(fill=is_existed)) + 
  #   scale_fill_manual(values=c("no"="#e0e0e0","yes"="#66c2a5"))+
  #   coord_flip() +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(),
  #         axis.line = element_line(colour = "black"),
  #         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         axis.text.y=element_text(colour="black",size=rel(1)),
  #         axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
  #         axis.title=element_text(colour="black",size=rel(1)),
  #         legend.position="none",
  #         legend.background = element_blank(),
  #         legend.text=element_text(size=rel(1))) 
  
  subfamily_info.macFas5.plot.sum.tmp2 = subfamily_info.macFas5.plot.sum.tmp[subfamily_info.macFas5.plot.sum.tmp$label.final.correctedName %in% ordered_subfamilies,]
  subfamily_info.macFas5.plot.sum.tmp2$label.final.correctedName = factor(subfamily_info.macFas5.plot.sum.tmp2$label.final.correctedName,levels=c(ordered_subfamilies))
  p1_2 <- ggplot(subfamily_info.macFas5.plot.sum.tmp2, aes(label.final.correctedName, freq)) + 
    geom_col(aes(fill=is_existed)) + 
    scale_fill_manual(values=c("no"="#e0e0e0","yes"="#66c2a5"))+
    coord_flip() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text.y=element_text(colour="black",size=rel(1)),
          axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
          axis.title=element_text(colour="black",size=rel(1)),
          legend.position="none",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1))) 
  
  ## heatmap of divergent rate
  # order the divergent rate table
  
  Div_paired.plot = Div_paired
  Div_paired.plot = Div_paired.plot[!grepl("_U",Div_paired.plot$subfamily.macaque.new),]
  Div_paired.plot = Div_paired.plot[!is.na(Div_paired.plot$subfamily.name.final),]
  Div_paired.plot2 = Div_paired.plot[Div_paired.plot$subfamily.name.final %in% subfamily_info.hg19.order$subfamily.name.final2 & 
                                       Div_paired.plot$subfamily.macaque.new %in% ordered_subfamilies,]
  Div_paired.plot$subfamily.macaque.new = factor(Div_paired.plot$subfamily.macaque.new,levels=c(ordered_subfamilies,c("MER11A_g4","MER11A_g12")))
  Div_paired.plot$subfamily.name.final = factor(Div_paired.plot$subfamily.name.final,levels=c(unique(subfamily_info.hg19.order$subfamily.name.final2,c("11D_a","11D_b","11D_c","MER11D"))))
  Div_paired.plot2$subfamily.macaque.new = factor(Div_paired.plot2$subfamily.macaque.new,levels=c(ordered_subfamilies))
  Div_paired.plot2$subfamily.name.final = factor(Div_paired.plot2$subfamily.name.final,levels=c(subfamily_info.hg19.order[!grepl("_U$",subfamily_info.hg19.order$subfamily.name.final2),]$subfamily.name.final2))
  
  subfamily_info.macFas5.plot.sum.plot = subfamily_info.macFas5.plot.sum
  subfamily_info.macFas5.plot.sum.plot = subfamily_info.macFas5.plot.sum.plot[!grepl("_U",subfamily_info.macFas5.plot.sum.plot$label.final.correctedName),]
  subfamily_info.macFas5.plot.sum.plot = subfamily_info.macFas5.plot.sum.plot[!is.na(subfamily_info.macFas5.plot.sum.plot$subfamily.name.final) & subfamily_info.macFas5.plot.sum.plot$subfamily.name.final != "other",]
  subfamily_info.macFas5.plot.sum.plot2 = subfamily_info.macFas5.plot.sum.plot[subfamily_info.macFas5.plot.sum.plot$subfamily.name.final %in% subfamily_info.hg19.order[!grepl("_U$",subfamily_info.hg19.order$subfamily.name.final2),]$subfamily.name.final2 & 
                                                                                 subfamily_info.macFas5.plot.sum.plot$label.final.correctedName %in% ordered_subfamilies,]
  subfamily_info.macFas5.plot.sum.plot$label.final.correctedName = factor(subfamily_info.macFas5.plot.sum.plot$label.final.correctedName,levels=c(ordered_subfamilies,c("MER11A_g4","MER11A_g12")))
  subfamily_info.macFas5.plot.sum.plot$subfamily.name.final = factor(subfamily_info.macFas5.plot.sum.plot$subfamily.name.final,levels=c(unique(subfamily_info.hg19.order$subfamily.name.final2,c("11D_a","11D_b","11D_c","MER11D"))))
  subfamily_info.macFas5.plot.sum.plot2$label.final.correctedName = factor(subfamily_info.macFas5.plot.sum.plot2$label.final.correctedName,levels=c(ordered_subfamilies))
  subfamily_info.macFas5.plot.sum.plot2$subfamily.name.final = factor(subfamily_info.macFas5.plot.sum.plot2$subfamily.name.final,levels=c(subfamily_info.hg19.order[!grepl("_U$",subfamily_info.hg19.order$subfamily.name.final2),]$subfamily.name.final2))
  
  # 
  subfamily_info.macFas5.unique.plot = subfamily_info.macFas5.unique[!grepl("_U",subfamily_info.macFas5.unique$label.final.correctedName),]
  subfamily_info.macFas5.unique.plot$label.final.correctedName = factor(subfamily_info.macFas5.unique.plot$label.final.correctedName,levels=c(ordered_subfamilies,c("MER11A_g4","MER11A_g12")))
  subfamily_info.macFas5.unique.plot2 = subfamily_info.macFas5.unique.plot[subfamily_info.macFas5.unique.plot$label.final.correctedName %in% ordered_subfamilies,]
  subfamily_info.macFas5.unique.plot2$label.final.correctedName = factor(subfamily_info.macFas5.unique.plot2$label.final.correctedName,levels=ordered_subfamilies)
  
  head(Div_paired.plot)
  # p2 <- ggplot(Div_paired.plot, aes(y=subfamily.macaque.new, x=subfamily.name.final)) + 
  #   geom_tile(aes(fill=V3)) + 
  #   geom_point(data=subfamily_info.macFas5.plot.sum.plot,aes(y=label.final.correctedName,x=subfamily.name.final,size=freq),
  #              shape=21,fill=NA)+
  #   scale_x_discrete(drop = FALSE)+
  #   scale_y_discrete(drop = FALSE)+
  #   scale_size(limits=c(0,1))+
  #   #scale_fill_viridis_c()+
  #   scale_fill_gradient2(low="#2166ac",mid = "white", high = "#b2182b",midpoint = 0.1)+
  #   xlab(NULL) + ylab(NULL) + 
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_rect(fill = "#f0f0f0"),
  #         axis.line = element_line(colour = "black"),
  #         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         axis.text.y=element_text(colour="black",size=rel(1)),
  #         axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
  #         axis.title=element_text(colour="black",size=rel(1)),
  #         legend.position="right",
  #         legend.background = element_blank(),
  #         legend.text=element_text(size=rel(1)))
  p2_2<- ggplot(Div_paired.plot2, aes(y=subfamily.macaque.new, x=subfamily.name.final)) + 
    geom_tile(aes(fill=V3)) + 
    geom_point(data=subfamily_info.macFas5.plot.sum.plot2,aes(y=label.final.correctedName,x=subfamily.name.final,size=freq),
               shape=21,color="#b2182b",fill=NA)+
    scale_x_discrete(drop = FALSE)+
    scale_y_discrete(drop = FALSE)+
    scale_size(limits=c(0,1))+
    scale_fill_viridis_c()+
    #scale_fill_gradient2(low="#2166ac",mid = "white", high = "#b2182b")+
    #scale_fill_gradient2(low="navy", mid="white", high="red") +
    xlab(NULL) + ylab(NULL) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "#f0f0f0"),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text.y=element_text(colour="black",size=rel(1)),
          axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
          axis.title=element_text(colour="black",size=rel(1)),
          legend.position="right",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1)))
  
  ## p3, actual counts
  # p3 <- ggplot(subfamily_info.macFas5.unique.plot, aes(label.final.correctedName,n)) + 
  #   geom_col(aes(),fill = "#737373") + 
  #   geom_text(aes(y=n,x=label.final.correctedName,label=n))+
  #   coord_flip() +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(),
  #         axis.line = element_line(colour = "black"),
  #         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  #         axis.text.y=element_text(colour="black",size=rel(1)),
  #         axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
  #         axis.title=element_text(colour="black",size=rel(1)),
  #         legend.position="none",
  #         legend.background = element_blank(),
  #         legend.text=element_text(size=rel(1))) 
  p3_2 <- ggplot(subfamily_info.macFas5.unique.plot2, aes(label.final.correctedName,n)) + 
    geom_col(aes(),fill = "#737373") + 
    geom_text(aes(y=n,x=label.final.correctedName,label=n))+
    coord_flip() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text.y=element_text(colour="black",size=rel(1)),
          axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
          axis.title=element_text(colour="black",size=rel(1)),
          legend.position="none",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1))) 
  
  # gA <- ggplotGrob(p1)
  # gB <- ggplotGrob(p2)
  # gC <- ggplotGrob(p3)
  # g2 = cbind(gA,gB,gC,size = "first")
  # pdf(paste("Step19_2_",Tree,"_rerooted_heatmap.pdf",sep=""),    # create PNG for the heat map
  #     width = 36,        # 5 x 300 pixels
  #     height = 6,
  #     pointsize = 10)        # smaller font size
  # grid.draw(g2)
  # dev.off()
  
  gA <- ggplotGrob(p1_2)
  gB <- ggplotGrob(p2_2)
  gC <- ggplotGrob(p3_2)
  g2 = cbind(gA,gB,gC,size = "first")
  pdf(paste("Figure_4C-","macFas5","-",k,"_similarity-2.pdf",sep=""),    # create PNG for the heat map
      width = 36,        # 5 x 300 pixels
      height = 6,
      pointsize = 10)        # smaller font size
  grid.draw(g2)
  dev.off()
#}

# merge
tr.dis2 = merge(tr.dis,data.frame(label_table),by.x="subfamily",by.y="label",all.x=T)
tr.dis2 = tr.dis2[match(ordered_subfamilies,tr.dis2$subfamily),]
tr.dis2$Order.rerooted.tree = 1:nrow(tr.dis2)
#write.csv(tr.dis2,file = paste("Figure_4C-",Tree,"_rerooted.csv",sep=""))

## 
# load FASTA
fasta.plot = readDNAStringSet(paste("input_trees_consensus/",gsub(".nonrev_dna|.rev_dna","",Tree),sep=""))
seq_name = names(fasta.plot)
sequence = paste(fasta.plot)
fasta.df <- data.frame(seq_name, sequence)
fasta.df = fasta.df[fasta.df$seq_name %in% tr.dis2$subfamily,]
fasta.df = fasta.df[match(tr.dis2$subfamily,fasta.df$seq_name),]

if (grepl("^hg19",Tree)){
  fasta.df$seq_name = tr.dis2$subfamily.name.final2
}
fasta.df = fasta.df[!is.na(fasta.df$seq_name),]
fasta.df$sequence = gsub("B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z","N",fasta.df$sequence)
fasta.df$seq_name = paste(">",fasta.df$seq_name,sep="")
# write to a re-ordered fasta file
write.table(fasta.df, sep="\n",row.names = FALSE, col.names = FALSE, quote = FALSE,file=paste(Tree,".reordered",sep=""))

# p_iris_msa = msaplot(p=ggtree(tr,aes(),size=0.2)+geom_tiplab(align=TRUE) + hexpand(.01) + 
#                        theme_tree2(), 
#                      fasta=paste("input/",Tree,".reordered",sep=""),height=1,offset=.1,width=3)
# pdf(paste("Step19_2_",Tree,"_rerooted_msa.pdf",sep=""),    # create PNG for the heat map
#     width = 24,        # 5 x 300 pixels
#     height = 6,
#     pointsize = 10)        # smaller font size
# grid.draw(p_iris_msa)
# dev.off() 

#############################################
##################### plot 3: activity distribution
# ordered_subfamilies = subfamily.macaque.new
# Summary_Table2 = Summary_Table2[,!colnames(Summary_Table2) %in% c("X.x","X.y")]
# subfamily_info.macFas5.activity = merge(subfamily_info.macFas5,Summary_Table2[,!grepl("Seq",colnames(Summary_Table2))],by.x="label",by.y="instanceID.renamed",all.y=T)
# 
# # plot 3.1 activity plot
# subfamily_info.macFas5.activity.plot = subfamily_info.macFas5.activity[subfamily_info.macFas5.activity$species == "macFas5",]
# subfamily_info.macFas5.activity.plot_MER11 = subfamily_info.macFas5.activity.plot[subfamily_info.macFas5.activity.plot$Family == "MER11A",]
# subfamily_info.macFas5.activity.plot_MER11$label.final.correctedName
# subfamily_info.macFas5.activity.plot_MER11 = subfamily_info.macFas5.activity.plot_MER11[subfamily_info.macFas5.activity.plot_MER11$label.final.correctedName %in% c(ordered_subfamilies,c("MER11A_g4","MER11A_g12")),]
# subfamily_info.macFas5.activity.plot_MER11$label.final.correctedName = factor(subfamily_info.macFas5.activity.plot_MER11$label.final.correctedName,levels=c(ordered_subfamilies,c("MER11A_g4","MER11A_g12")))
# 
# ##### 
# subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group = ifelse(subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive == "None","No activity",NA)
# subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group = ifelse(subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive != "None" & 
#                                                                                        subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore < 2,"<2", subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group)
# subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group = ifelse(subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive != "None" & 
#                                                                                        subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore >= 2 & subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore < 4 ,"2-4" ,subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group)
# subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group = ifelse(subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive != "None" & 
#                                                                                        subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore >= 4 & subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore < 6 ,"4-6", subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group)
# subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group = ifelse(subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive != "None" & 
#                                                                                        subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore >= 6 & subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore < 8 ,"6-8", subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group)
# subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group = ifelse(subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive != "None" & 
#                                                                                        subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore >= 8 , ">=8", subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group)
# subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group = factor(subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group,levels=c(">=8","6-8","4-6","2-4","<2","No activity"))
# head(subfamily_info.macFas5.activity.plot_MER11[is.na(subfamily_info.macFas5.activity.plot_MER11$iPSC.alpha.Zscore.isActive.group),],20)
# 
# Alpha.Zscore.isActive.group.color = c(">=8" = "#67001f","6-8"="#b2182b","4-6"="#d6604d","2-4"="#f4a582","<2"="#fddbc7","No activity"="#f0f0f0")
# color_alpha.group = c("alpha>=4"="#bd0026","alpha>=2"="#fd8d3c","low"="#f0f0f0")
# 
# 
# ### summary
# subfamily_info.macFas5.activity.plot_MER11.highQuality = subfamily_info.macFas5.activity.plot_MER11[grepl("highQuality",subfamily_info.macFas5.activity.plot_MER11$is_kept),]
# # excluded NA values that refer to low quality in each cell
# subfamily_info.macFas5.activity.plot_MER11.highQuality = subfamily_info.macFas5.activity.plot_MER11.highQuality[!is.na(subfamily_info.macFas5.activity.plot_MER11.highQuality$iPSC.alpha.Zscore.isActive.group),]
# 
# # iPSC
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum = data.frame(subfamily_info.macFas5.activity.plot_MER11.highQuality %>% group_by(label.final.correctedName,Frame,iPSC.alpha.Zscore.isActive.group) %>% dplyr::summarise(n = n()))
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$uniqueName = paste(subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$label.final.correctedName,subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$Frame)
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum_tmp = data.frame(subfamily_info.macFas5.activity.plot_MER11.highQuality %>% group_by(label.final.correctedName,Frame) %>% dplyr::summarise(total = n()))
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum_tmp$uniqueName = paste(subfamily_info.macFas5.activity.plot_MER11.highQuality_sum_tmp$label.final.correctedName,subfamily_info.macFas5.activity.plot_MER11.highQuality_sum_tmp$Frame)
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum = merge(subfamily_info.macFas5.activity.plot_MER11.highQuality_sum,subfamily_info.macFas5.activity.plot_MER11.highQuality_sum_tmp,by="uniqueName",all=T)
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$perC = subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$n/subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$total
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$iPSC.alpha.Zscore.isActive.group = factor(subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$iPSC.alpha.Zscore.isActive.group,levels=rev(c(">=8","6-8","4-6","2-4","<2","No activity")))
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$label.final.correctedName.x = factor(subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$label.final.correctedName.x,levels=c(ordered_subfamilies,c("MER11A_g4","MER11A_g12")))
# subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$perC = ifelse(subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$total<10,NA,subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$perC)
# p1<-ggplot(subfamily_info.macFas5.activity.plot_MER11.highQuality_sum[subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$Frame.x == "F1",],aes(x=label.final.correctedName.x,y=perC, fill=iPSC.alpha.Zscore.isActive.group))+
#   geom_bar(position="stack", stat="identity")+
#   geom_hline(yintercept = c(0.5), linetype="dotted",color = "white", linewidth=1)+
#   scale_fill_manual(values=Alpha.Zscore.isActive.group.color) + 
#   scale_x_discrete(drop = FALSE)+
#   ggtitle("iPS cells")+
#   ylab("Proportion")+
#   xlab("")+
#   theme(
#     plot.title = element_text(hjust = 0.5, size = rel(1)),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#     axis.text=element_text(colour="black",size=rel(1),angle = 0),
#     axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
#     axis.title=element_text(colour="black",size=rel(1)),
#     legend.key = element_rect(colour = "transparent", fill = "white"),
#     legend.position="right",
#     legend.background = element_blank(),
#     legend.text=element_text(size=rel(1)))
# p2<-ggplot(subfamily_info.macFas5.activity.plot_MER11.highQuality_sum[subfamily_info.macFas5.activity.plot_MER11.highQuality_sum$Frame.x == "F2",],aes(x=label.final.correctedName.x,y=perC, fill=iPSC.alpha.Zscore.isActive.group))+
#   geom_bar(position="stack", stat="identity")+
#   geom_hline(yintercept = c(0.5), linetype="dotted",color = "white", linewidth=1)+
#   scale_fill_manual(values=Alpha.Zscore.isActive.group.color) + 
#   scale_x_discrete(drop = FALSE)+
#   ggtitle("iPS cells")+
#   ylab("Proportion")+
#   xlab("")+
#   theme(
#     plot.title = element_text(hjust = 0.5, size = rel(1)),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#     axis.text=element_text(colour="black",size=rel(1),angle = 0),
#     axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
#     axis.title=element_text(colour="black",size=rel(1)),
#     legend.key = element_rect(colour = "transparent", fill = "white"),
#     legend.position="right",
#     legend.background = element_blank(),
#     legend.text=element_text(size=rel(1)))
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# g2 = rbind(gA,gB,size = "first")
# pdf(paste("Step19_2_",Tree,"_barplot.pdf",sep=""),    # create PNG for the heat map
#     width = 8,        # 5 x 300 pixels
#     height = 8,
#     pointsize = 10)        # smaller font size
# grid.draw(g2)
# dev.off()

##### box plot
# p1 = ggplot(subfamily_info.macFas5.activity.plot_MER11,aes(x=label.final.correctedName,iPSC.alpha.Zscore,group=label.final.correctedName))+
#   geom_jitter(aes(colour = iPSC.alpha.Zscore.isActive,fill=iPSC.alpha.Zscore.isActive,group=iPSC.alpha.Zscore.isActive),alpha=0.5)+
#   geom_violin(aes(),fill=NA) + #add points colored by significance
#   
#   stat_summary(fun = "mean",
#                fun.min = function(x)mean(x)-sd(x),
#                fun.max = function(x)mean(x) + sd(x),
#                geom = "pointrange",
#                position = position_dodge(width = .9),shape=95)+
#   ylab("Z-scaled activity in iPS cells")+
#   xlab("")+
#   scale_x_discrete(drop = FALSE)+
#   #scale_fill_manual(values=cluster.color)+
#   #scale_color_manual(values=cluster.color)+
#   #stat_compare_means(method = "t.test",paired=TRUE,aes(label = paste0("p = ", ..p.format..)))+
#   #ggtitle(paste(Family))+
#   theme(
#     plot.title = element_text(hjust = 0.5, size = rel(1)),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#     axis.text=element_text(colour="black",size=rel(1),angle = 0),
#     axis.text.x=element_text(colour="black",vjust=.5,hjust=0.95,angle = 90),
#     axis.title=element_text(colour="black",size=rel(1)),
#     legend.key = element_rect(colour = "transparent", fill = "white"),
#     #legend.position="right",
#     legend.position = "none",
#     legend.background = element_blank(),
#     legend.text=element_text(size=rel(1)))
# 
# ##### motifs
# Summary_Table2.fimoPlot = Summary_Table2[!is.na(Summary_Table2$iPSC.alpha.range),]   ## exclude positive and negative
# # 
# Fimo.all = read.delim("input/MER11_34_52_frame_fimo.tsv_2023_1_7.gz",header=F,sep="")
# colnames(Fimo.all) = c("Family.frame","DB","motif_id","motif_name","uniqueID.new","start","stop","strand","score","p_value","q_value","matched_sequence")
# Fimo.all$uniqueID.new.DB.motif_name = paste(Fimo.all$Family.frame,Fimo.all$uniqueID.new,Fimo.all$DB,Fimo.all$motif_name)
# Fimo.all.sum = data.frame(Fimo.all %>% group_by(Family.frame,DB,uniqueID.new,motif_name) %>% dplyr::count())
# 
# Fimo.all.sum.plot = Fimo.all.sum[Fimo.all.sum$uniqueID.new %in% subfamily_info.macFas5.activity.plot_MER11$uniqueID_new,]
# Fimo.all.sum.plot = merge(Fimo.all.sum.plot,subfamily_info.macFas5.activity.plot_MER11[,c("uniqueID_new","label.final.correctedName")],by.x="uniqueID.new",by.y="uniqueID_new",all.x=T)
# #
# 
# Fimo.all.sum.plot.sum = data.frame(Fimo.all.sum.plot %>% group_by(label.final.correctedName,DB,motif_name) %>% dplyr::count())
# Fimo.all.sum.plot.sum.total = data.frame(subfamily_info.macFas5.activity.plot_MER11 %>% group_by(label.final.correctedName) %>% dplyr::count())
# Fimo.all.sum.plot.sum = merge(Fimo.all.sum.plot.sum,Fimo.all.sum.plot.sum.total,by="label.final.correctedName",all.x=T)
# head(Fimo.all.sum.plot.sum)
# Fimo.all.sum.plot.sum$perC = Fimo.all.sum.plot.sum$n.x/Fimo.all.sum.plot.sum$n.y
# Fimo.all.sum.plot.sum = Fimo.all.sum.plot.sum[Fimo.all.sum.plot.sum$motif_name %in% c("CRX","DMBX1","ETV2","FOXF2","GATA5","GSC","GSC2","INSM1","MECOM","PITX2","PITX1","PITX3","OTX2",
#                                                                                       "RARA::RXRG","RFX6","TEAD1","TEAD2","ZIC1::ZIC2","ZIC2","ZIC3","ZNF136","ZNF317"),]
# 
# Fimo.all.sum.plot.sum$label.final.correctedName = factor(Fimo.all.sum.plot.sum$label.final.correctedName,levels=ordered_subfamilies)
# Fimo.all.sum.plot.sum$motif_name = factor(Fimo.all.sum.plot.sum$motif_name,levels = c("CRX","DMBX1","ETV2","FOXF2","GATA5","GSC","GSC2","INSM1","MECOM","PITX2","PITX1","PITX3","OTX2",
#                                                                                       "RARA::RXRG","RFX6","TEAD1","TEAD2","ZIC1::ZIC2","ZIC2","ZIC3","ZNF136","ZNF317"))
# 
# ### plot
# p1 = ggplot(Fimo.all.sum.plot.sum, aes(x=label.final.correctedName,y=motif_name)) +
#   geom_tile(aes(fill=perC*100)) + 
#   #scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 50,na.value = 'black')+
#   scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
#   scale_x_discrete(drop = FALSE)+
#   xlab("divergent rate")+
#   ylab("% instances containing the motif")+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = "#f0f0f0"),
#         axis.line = element_line(colour = "black"),
#         #axis.line.x = element_blank(),
#         #axis.line.y = element_blank(),
#         axis.text.y=element_text(colour="black",size=rel(1)),
#         axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
#         axis.title=element_text(colour="black",size=rel(1)),
#         legend.position="right",
#         legend.background = element_blank(),
#         legend.text=element_text(size=rel(1))) 
# 
# pdf(paste("Step19_2_",Tree,"_motif.pdf",sep=""),    # create PNG for the heat map
#     width = 10,        # 5 x 300 pixels
#     height = 5.5,
#     pointsize = 10)        # smaller font size
# grid.draw(p1)
# dev.off()
