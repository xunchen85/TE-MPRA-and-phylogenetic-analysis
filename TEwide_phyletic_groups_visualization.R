#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

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

## divergent rate
Div_paired = read.delim("input/RAxML_distances.2023_6_28",header=F,sep="")
colnames(Div_paired) = c("Group","from","to","score")
# add the same sequences 
Div_paired.same1 = Div_paired[!duplicated(Div_paired[,c("Group","from")]),]
Div_paired.same1$to = Div_paired.same1$from
Div_paired.same2 = Div_paired[!duplicated(Div_paired[,c("Group","to")]),]
Div_paired.same2$from = Div_paired.same2$to
Div_paired.same = rbind(Div_paired.same1,Div_paired.same2)
Div_paired.same = Div_paired.same[!duplicated(Div_paired.same[,c("Group","from")]),]
Div_paired.same$score = 0
Div_paired = rbind(Div_paired,Div_paired.same)
rm(Div_paired.same1,Div_paired.same2,Div_paired.same)

### divergent rate for MER11
# MER11 order of human and macaque
MER11_hg191 = data.frame(species = "hg19",group="FG1",subfamily = c("11A_p","MER11A","11A_m","11A_n","11A_q","11A_i","11A_k","11A_a","11A_o","11A_r","11A_l"))
MER11_hg192 = data.frame(species = "hg19",group="FG2",subfamily = c("11A_b","11A_c","11C_w","11A_j","11A_d","11A_e","11B_i","11A_f"))
MER11_hg193 = data.frame(species = "hg19",group="FG3",subfamily = c("11A_g","11B_h","11B_g","11A_h","11C_v","11B_k","11B_f","11B_n","11B_r","11B_o","11C_o","11B_q","11C_r"))
MER11_hg194 = data.frame(species = "hg19",group="FG4",subfamily = c("11B_e","11B_c","11B_b","11B_p","11B_m","11C_x","11B_d","11B_l","11C_q","11B_j","11C_a","11C_y","11B_a","11C_t",
                                                                    "11C_z","11C_s","11C_a2","MER11B","11C_b2","11C_b","11C_c","11C_d","11C_e","11C_f","11C_g","11C_h","11C_j","11C_i","11C_p",
                                                                    "MER11C","11C_n","11C_m","11C_k","11C_l"))

MER11_macque0 = data.frame(species = "macFas5",group="FG2",subfamily = c("MER11A_g8","MER11A_g13"))
MER11_macque1 = data.frame(species = "macFas5",group="FG1",subfamily = c("MER11A_g17","MER11A_g18","MER11A_g15","MER11A_g2","MER11A_g22","MER11A_g19","MER11A_g33"))
MER11_macque2 = data.frame(species = "macFas5",group="FG2",subfamily = c("MER11A_g16","MER11A_g21"))
MER11_macque3 = data.frame(species = "macFas5",group="FG3",subfamily = c("MER11A_g5","MER11A_g24","MER11A_g25","MER11A_g26","MER11A_g23","MER11A_g30","MER11A_g11","MER11A_g6","MER11A_g9"))
MER11_macque4 = data.frame(species = "macFas5",group="FG4",subfamily = c("MER11A_g7","MER11A_g29","MER11A_g1","MER11A_g28","MER11A_g27","MER11A_g20","MER11A_g32","MER11A_g14",
                                                                         "MER11A_g31","MER11A_g10","MER11A_g3"))
MER11_both = rbind(MER11_hg191,MER11_hg192,MER11_hg193,MER11_hg194,
                   MER11_macque0,MER11_macque1,MER11_macque2,MER11_macque3,MER11_macque4)
MER11_both$Order = NA
MER11_both[MER11_both$species == "hg19",]$Order = 1:nrow(MER11_both[MER11_both$species == "hg19",])
MER11_both[MER11_both$species == "macFas5",]$Order = 1:nrow(MER11_both[MER11_both$species == "macFas5",])
rm(MER11_hg191,MER11_hg192,MER11_hg193,MER11_hg194,MER11_macque1,MER11_macque2,MER11_macque3,MER11_macque4)

# rename table
Div_paired_MER11.renamed.table = read.csv("input/MER11_subtree_info_2022_12_27.csv")
Div_paired_MER11.renamed.table$consensus.name = paste("h|",Div_paired_MER11.renamed.table$consensus.name,sep="")
Div_paired_MER11.renamed.table[nrow(Div_paired_MER11.renamed.table)+1,] = NA
Div_paired_MER11.renamed.table[nrow(Div_paired_MER11.renamed.table),c("subfamily.name.final","consensus.name")] = c("MER11A","MER11A")
Div_paired_MER11.renamed.table[nrow(Div_paired_MER11.renamed.table)+1,] = NA
Div_paired_MER11.renamed.table[nrow(Div_paired_MER11.renamed.table),c("subfamily.name.final","consensus.name")] = c("MER11B","MER11B")
Div_paired_MER11.renamed.table[nrow(Div_paired_MER11.renamed.table)+1,] = NA
Div_paired_MER11.renamed.table[nrow(Div_paired_MER11.renamed.table),c("subfamily.name.final","consensus.name")] = c("MER11C","MER11C")
Div_paired_MER11.renamed.table[nrow(Div_paired_MER11.renamed.table)+1,] = NA
Div_paired_MER11.renamed.table[nrow(Div_paired_MER11.renamed.table),c("subfamily.name.final","consensus.name")] = c("MER11D","MER11D")

# human vs human
Div_paired_MER11.hh = Div_paired[Div_paired$Group == "hg19_MER11ABC",]
Div_paired_MER11.hh$X.axis = Div_paired_MER11.hh$from
Div_paired_MER11.hh$Y.axis = Div_paired_MER11.hh$to
Div_paired_MER11.hh2 = Div_paired[Div_paired$Group == "hg19_MER11ABC",]
Div_paired_MER11.hh2$X.axis = Div_paired_MER11.hh$to
Div_paired_MER11.hh2$Y.axis = Div_paired_MER11.hh$from
Div_paired_MER11.hh = rbind(Div_paired_MER11.hh,Div_paired_MER11.hh2)
Div_paired_MER11.hh = merge(Div_paired_MER11.hh,Div_paired_MER11.renamed.table[,c("subfamily.name.final","consensus.name")],by.x="X.axis",by.y="consensus.name",all.x=T)
Div_paired_MER11.hh = merge(Div_paired_MER11.hh,Div_paired_MER11.renamed.table[,c("subfamily.name.final","consensus.name")],by.x="Y.axis",by.y="consensus.name",all.x=T)
Div_paired_MER11.hh$subfamily.name.final.x = factor(Div_paired_MER11.hh$subfamily.name.final.x,levels = factor(MER11_both[MER11_both$species == "hg19",]$subfamily))
Div_paired_MER11.hh$subfamily.name.final.y = factor(Div_paired_MER11.hh$subfamily.name.final.y,levels = factor(MER11_both[MER11_both$species == "hg19",]$subfamily))

# human vs macaque
Div_paired_MER11 = Div_paired[Div_paired$Group == "hg19_macFas5_MER11ABC",]
Div_paired_MER11.mh = Div_paired_MER11[(grepl("MER11A_g",Div_paired_MER11$from) & !grepl("MER11A_g",Div_paired_MER11$to)) |
                                         (!grepl("MER11A_g",Div_paired_MER11$from) & grepl("MER11A_g",Div_paired_MER11$to)),]
Div_paired_MER11.mh$Y.axis = ifelse(grepl("MER11A_g",Div_paired_MER11.mh$from),Div_paired_MER11.mh$from,Div_paired_MER11.mh$to)
Div_paired_MER11.mh$X.axis = ifelse(!grepl("MER11A_g",Div_paired_MER11.mh$from),Div_paired_MER11.mh$from,Div_paired_MER11.mh$to)
Div_paired_MER11.mh = merge(Div_paired_MER11.mh,Div_paired_MER11.renamed.table[,c("subfamily.name.final","consensus.name")],by.x="X.axis",by.y="consensus.name")
Div_paired_MER11.mh$subfamily.name.final = factor(Div_paired_MER11.mh$subfamily.name.final,levels = factor(MER11_both[MER11_both$species == "hg19",]$subfamily))
Div_paired_MER11.mh$Y.axis = factor(Div_paired_MER11.mh$Y.axis,levels = factor(MER11_both[MER11_both$species == "macFas5",]$subfamily))
Div_paired_MER11.mh = Div_paired_MER11.mh[!is.na(Div_paired_MER11.mh$Y.axis),]

### LTR consensus sequence name from Carter et al.
Consensus_LTR7 = read.delim("input/LTR7_consensus_sequences_cedrid_2023_6_19.fa",header=F,sep="")
Consensus_LTR7 = Consensus_LTR7[grepl(">",Consensus_LTR7$V1),]
Consensus_LTR7$V1 = gsub(">","",Consensus_LTR7$V1)
colnames(Consensus_LTR7) = c("sequence.ID","sequence.name")

### liftover of subfamily level
#liftover_subfamily = read.delim("../Final_liftover_intersect_2022_11_23/TEwide_liftover.subfamily.summary_2023_6_26",header=F,sep="")
liftover_subfamily = read.delim("input/TEwide_liftover.phyleticGroup.summary_2023_8_29",header=F,sep="")

colnames(liftover_subfamily) = c("species.liftover","Family","Family_Count_analyzed","Subfamily","Node_ID","Subfamily_Count","liftover_count","liftover_count_len200","Subfamily_color")
liftover_subfamily$perC = liftover_subfamily$liftover_count_len200/liftover_subfamily$Subfamily_Count
liftover_subfamily$species = ifelse(grepl("^hg19",liftover_subfamily$species.liftover),"hg19","macFas5")

# rename 11C
liftover_subfamily$Subfamily = ifelse(grepl("^11|^34|^52",liftover_subfamily$Subfamily),paste("h|",liftover_subfamily$Subfamily,sep=""),liftover_subfamily$Subfamily)

### instance information based on tree analysis
#summary_table <- read.delim("../Final_TEwide_tree_2023_4_23/TEwide_group.list_2023_6_26",header=F,sep="")
#summary_table_25c = summary_table[summary_table$V25 !="",]
#summary_table_24c = summary_table[summary_table$V25 =="",]
#summary_table_24c$V25 = "-"
#summary_table_24c = summary_table_24c[,c(1:19,25,20,21,22,23,24)]
#colnames(summary_table_24c) = colnames(summary_table_25c)
#summary_table.missed <- read.delim("../Final_TEwide_tree_2023_4_23/TEwide_group.list_2023_6_26.missed_subfamilies",header=F,sep="")
#summary_table.missed$V25 = "-"
#summary_table.missed = summary_table.missed[,c(1:19,25,20,21,22,23,24)]
#colnames(summary_table.missed) = colnames(summary_table_25c)
#summary_table.MER11_34_52 <- read.delim("../Final_TEwide_tree_2023_4_23/Step4_hg19_rmsk_TE_0bp_MER11_34_52_group_final_2022_11_9.list",header=F,sep="")
#summary_table.MER11_34_52$V8 = summary_table.MER11_34_52$V4
#summary_table = rbind(summary_table_25c,summary_table_24c,summary_table.missed,summary_table.MER11_34_52)
#summary_table = summary_table[!grepl("^MER34",summary_table$V1),]
#summary_table = summary_table[!grepl("^MER52",summary_table$V1),]
#summary_table = summary_table[!grepl("macFas5_rmsk",summary_table$V25),]
#write.table(summary_table, sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE,file="TEwide_group.list_final_2023_8_29.hg19")
summary_table <- read.delim("input/TEwide_group.list_final_2023_8_29.hg19",header=F,sep="")

### selected trees
rooted.trees.summary <- read.csv("input/root_summary_2023_8_29.csv")
rooted.trees.summary$Combination = ifelse(rooted.trees.summary$Combination == "hg19_G5_v1","hgg19_MER11ABC_v1",rooted.trees.summary$Combination)

### kept species
Kept_species = c("hg19ToHg38","hg19ToPanTro6","hg19ToGorGor3","hg19ToPonAbe2","hg19ToNomLeu3",
                 "hg19ToMacFas5","hg19ToPapAnu2","hg19ToCalJac3","hg19ToMicMur1","hg19ToMm10")

#############################################
##################### plot 1
# re-rooted tree 
tree_list = list.files("input_trees_consensus/")
tree_list = tree_list[grepl("roottest.trees$",tree_list)]
tree_list = tree_list[grepl("hg19_",tree_list)]
tree_list = gsub(".mafft.prank.best.fas.nonrev_dna.roottest.trees","",tree_list)

tree_IDs_orders = c("hg19_G1_v1","hg19_G10_v1","hg19_G11_v1","hg19_G12_v1","hg19_G13_v1","hg19_G14_v1","hg19_G15_v1",
                    "hg19_G16_v1","hg19_G17_v1","hg19_G18_v1","hg19_G19_v1","hg19_G2_v1","hg19_G3_v1","hg19_G4_v1",
                    "hg19_G6_v1","hg19_G8_v1","hg19_G9_v1","hg19_MER11ABC_v1")

tree_IDs_kept = c(4,2,3,18,1,2,1,4,2,7,1,1,3,2,1,2,5,23)

names(tree_IDs_kept) = tree_list

Tree = "hg19_MER11ABC_v1"
ordered_subfamilies.table = data.frame("subfamily"=NA,"species"=NA,"subfamily.species"=NA,"Order"=NA,"family.group"=NA,"treeOrder" = NA)
label_table_final = data.frame("tmp1"=NA,"tmp2"=NA)
for (Tree in tree_list) {
  # 1.1 load the tree
  dat = rooted.trees.summary[rooted.trees.summary$Combination == Tree,]

  # 1.2 load trees
  tr.all <- read.tree(paste("input_trees_consensus/",Tree,".mafft.prank.best.fas.nonrev_dna.roottest.trees",sep=""))

  # 1.3
  k=23
  for (k in dat$rank) {
    if(k != tree_IDs_kept[Tree]){
      next
    }
    tr = tr.all[[k]]
    tr$node.label = paste("Node",1:tr$Nnode,"/",tr$node.label,sep="")
    # load the summary info
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
    label_table$tree = Tree
    label_table$treeOrder = k
    
    ###### tree order and distance to root
    # distance
    tr.dis = data.frame(dist.nodes(tr))
    colnames(tr.dis) = label_table$label
    row.names(tr.dis) = label_table$label
    tr.dis$subfamily = rownames(tr.dis)
    tr.dis = tr.dis[rownames(tr.dis) %in% label_table$label,c("subfamily","Node1/")]
    tr.dis = tr.dis[!grepl("^Node",tr.dis$subfamily),]
    colnames(tr.dis) = c("label","branch.length.ToKeptNode")
    label_table = merge(label_table,tr.dis,by="label",all.x=T)
    
    # actual order
    tr2 <- ladderize(tr, right = FALSE)
    is_tip <- tr2$edge[,2] <= length(tr2$tip.label)
    ordered_tips <- tr2$edge[is_tip, 2]
    ordered_subfamilies.table.tmp = data.frame(tr2$tip.label[ordered_tips])
    colnames(ordered_subfamilies.table.tmp) = "subfamily"
    
    if(grepl("hg19",Tree)){
      Species = "hg19"
      ordered_subfamilies.table.tmp$species = Species
      ordered_subfamilies.table.tmp$subfamily.species = ifelse(ordered_subfamilies.table.tmp$species == "hg19" & grepl("_g",ordered_subfamilies.table.tmp$subfamily),paste("h",ordered_subfamilies.table.tmp$subfamily,sep="_"),
                                                               ordered_subfamilies.table.tmp$subfamily)
    } else {
      Species = "macFas5"
      ordered_subfamilies.table.tmp$species = Species
      ordered_subfamilies.table.tmp$subfamily.species = ifelse(ordered_subfamilies.table.tmp$species == "macFas5" & grepl("_g",ordered_subfamilies.table.tmp$subfamily),paste("m",ordered_subfamilies.table.tmp$subfamily,sep="_"),
                                                               ordered_subfamilies.table.tmp$subfamily)
    }
    ordered_subfamilies.table.tmp$Order = 1:nrow(ordered_subfamilies.table.tmp)
    ordered_subfamilies.table.tmp$family.group = gsub("Step21_|.mafft.prank.best.fas.nonrev_dna","",Tree)
    ordered_subfamilies.table.tmp$treeOrder = k
    label_table = merge(label_table,ordered_subfamilies.table.tmp[,c("subfamily","Order")],by.x="label",by.y="subfamily",all.x=T)
    colnames(label_table)[which(colnames(label_table) == "Order.x")] = "Order"
    colnames(label_table)[which(colnames(label_table) == "Order.y")] = "Order_subfamily"
    
    ### plot 1.1 tree plot
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
    ### plot 1.2 actual count plot
    if(unique(ordered_subfamilies.table.tmp$species) == "hg19") {
      Kept_species.plot = Kept_species
      liftover_subfamily.plot = liftover_subfamily[liftover_subfamily$Subfamily %in% ordered_subfamilies.table.tmp$subfamily & liftover_subfamily$species == "hg19",]
      #liftover_subfamily.plot = liftover_subfamily.plot[liftover_subfamily.plot$Subfamily_Count>=10,] 
      liftover_subfamily.plot = liftover_subfamily.plot[liftover_subfamily.plot$species.liftover %in% Kept_species.plot,]
      liftover_subfamily.plot$species.liftover = factor(liftover_subfamily.plot$species.liftover,levels=Kept_species.plot)
      liftover_subfamily.plot$Subfamily = factor(liftover_subfamily.plot$Subfamily,levels=rev(tr$tip.label))
      liftover_subfamily.plot1 = liftover_subfamily.plot[!duplicated(liftover_subfamily.plot[,c("Subfamily","Subfamily_Count")]),]
      liftover_subfamily.plot2 = liftover_subfamily.plot[liftover_subfamily.plot$species.liftover == "hg19ToMacFas5",]
    } else {
      Kept_species.plot = "macFas5ToHg19"
      liftover_subfamily.plot = liftover_subfamily[liftover_subfamily$Subfamily %in% ordered_subfamilies.table.tmp$subfamily & liftover_subfamily$species == "macFas5",]
      #liftover_subfamily.plot = liftover_subfamily.plot[liftover_subfamily.plot$Subfamily_Count>=10,] 
      liftover_subfamily.plot = liftover_subfamily.plot[liftover_subfamily.plot$species.liftover %in% Kept_species.plot,]
      liftover_subfamily.plot$species.liftover = factor(liftover_subfamily.plot$species.liftover,levels=Kept_species.plot)
      liftover_subfamily.plot$Subfamily = factor(liftover_subfamily.plot$Subfamily,levels=rev(tr$tip.label))
      liftover_subfamily.plot1 = liftover_subfamily.plot[!duplicated(liftover_subfamily.plot[,c("Subfamily","Subfamily_Count")]),]
      liftover_subfamily.plot2 = liftover_subfamily.plot[liftover_subfamily.plot$species.liftover == "macFas5ToHg19",]
    }
    
    ###
    label_table.final.tmp = merge(label_table,liftover_subfamily.plot2,by.x="label",by.y="Subfamily",all.x=T)
    if(k == tree_IDs_kept[Tree]){
      if (nrow(label_table_final) == 1) {
        label_table_final = label_table.final.tmp
      } else {
        label_table_final = rbind(label_table_final,label_table.final.tmp)
      }
      ordered_subfamilies.table = rbind(ordered_subfamilies.table,ordered_subfamilies.table.tmp)
    }
    
    # actual counts 
    p1.count <- ggplot(liftover_subfamily.plot1, aes(Subfamily, Subfamily_Count)) + 
      scale_x_discrete(drop = FALSE)+
      geom_col(aes(fill=Family)) + 
      geom_text(aes(label = Subfamily_Count))+
      scale_fill_brewer(palette="Accent")+
      coord_flip() + theme_tree2() + theme(legend.position='none')
    # liftover across species
    p1.liftover.species <- ggplot(liftover_subfamily.plot, aes(x=species.liftover, y=Subfamily)) + 
      geom_tile(aes(fill=perC)) + 
      geom_text(aes(label=round(perC,1)),color="white") + 
      scale_x_discrete(drop = FALSE)+
      scale_y_discrete(drop = FALSE)+
      scale_fill_viridis_c() + 
      theme_minimal() + xlab(NULL) + ylab(NULL)
    # liftover between hg19 and macFas5
    p1.liftover.human_macque <-ggplot(liftover_subfamily.plot2, aes(Subfamily, perC)) + 
      scale_x_discrete(drop = FALSE)+
      geom_col(aes(),fill = "#66c2a5") + 
      ylim(c(0,1))+
      coord_flip() +
      theme(#panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="none",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1))) 
    
    ## pair-wise similarity
    Div_paired.plot = Div_paired[Div_paired$from %in% ordered_subfamilies.table.tmp$subfamily & 
                                   Div_paired$to %in% ordered_subfamilies.table.tmp$subfamily,]
    if(Species == "hg19"){
      Div_paired.plot$from = gsub("^h_","",Div_paired.plot$from)
      Div_paired.plot$to = gsub("^h_","",Div_paired.plot$to)
    } else {
      Div_paired.plot$from = gsub("^m_","",Div_paired.plot$from)
      Div_paired.plot$to = gsub("^m_","",Div_paired.plot$to)
    }
    Div_paired.plot.order = merge(Div_paired.plot,ordered_subfamilies.table.tmp[,c("subfamily","Order")],all.x=T,by.x="from",by.y="subfamily")
    Div_paired.plot.order = merge(Div_paired.plot.order,ordered_subfamilies.table.tmp[,c("subfamily","Order")],all.x=T,by.x="to",by.y="subfamily")
    Div_paired.plot.order$X.axis = ifelse(Div_paired.plot.order$Order.x>=Div_paired.plot.order$Order.y,Div_paired.plot.order$from,Div_paired.plot.order$to)
    Div_paired.plot.order$Y.axis = ifelse(Div_paired.plot.order$Order.x>=Div_paired.plot.order$Order.y,Div_paired.plot.order$to,Div_paired.plot.order$from)
    Div_paired.plot.order2 = Div_paired.plot.order
    Div_paired.plot.order2$Y.axis = ifelse(Div_paired.plot.order2$Order.x>=Div_paired.plot.order2$Order.y,Div_paired.plot.order2$from,Div_paired.plot.order2$to)
    Div_paired.plot.order2$X.axis = ifelse(Div_paired.plot.order2$Order.x>=Div_paired.plot.order2$Order.y,Div_paired.plot.order2$to,Div_paired.plot.order2$from)
    Div_paired.plot.order = rbind(Div_paired.plot.order,Div_paired.plot.order2)
    #Div_paired.plot.order$X.axis = factor(Div_paired.plot.order$X.axis,levels=rev(ordered_subfamilies.table.tmp$subfamily))
    if (nrow(Div_paired.plot.order) == 0) {
      Div_paired.plot.order[1:length(ordered_subfamilies.table.tmp$subfamily),] = NA
      Div_paired.plot.order$X.axis = ordered_subfamilies.table.tmp$subfamily
      Div_paired.plot.order$Y.axis = ordered_subfamilies.table.tmp$subfamily
      Div_paired.plot.order$score = 1
    }
    Div_paired.plot.order$X.axis = factor(Div_paired.plot.order$X.axis,levels=(ordered_subfamilies.table.tmp$subfamily))
    Div_paired.plot.order$Y.axis = factor(Div_paired.plot.order$Y.axis,levels=ordered_subfamilies.table.tmp$subfamily)
    
    ##
    p1.similarity <- ggplot(Div_paired.plot.order, aes(x=X.axis, y=Y.axis)) + 
      geom_tile(aes(fill=score)) + 
      #geom_point(data=subfamily_info.macFas5.plot.sum.plot,aes(y=label.final.correctedName,x=subfamily.name.final,size=freq),
      #           shape=21,fill=NA)+
      #scale_size(limits=c(0,1))+
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
            axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.95),
            axis.title=element_text(colour="black",size=rel(1)),
            legend.position="right",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1)))
    ###
    p4 <- p1.similarity %>% insert_left(g) %>% insert_right(p1.liftover.species, width=.5) %>% insert_right(p1.liftover.human_macque, width=.5) %>% insert_right(p1.count, width=.5) 
    
    Width = 12
    if (nrow(ordered_subfamilies.table.tmp)<5){
      Width = 12
    } else if (nrow(ordered_subfamilies.table.tmp)>=5 & nrow(ordered_subfamilies.table.tmp)<=10){
      Width = 24
    } else {
      Width = 32
    }
    if(Tree == "hg19_MER11ABC_v1") {
      pdf("Figure_1E.pdf",    # create PNG for the heat map
          width = 40.5,        # 5 x 300 pixels
          height = 12,
          pointsize = 10)        # smaller font size
      grid.draw(p4)
      dev.off()
    } else {
      pdf(paste("Figure_1E_",Tree,".pdf",sep=""),    # create PNG for the heat map
          width = Width,        # 5 x 300 pixels
          height = 12,
          pointsize = 10)        # smaller font size
      grid.draw(p4)
      dev.off()
    }
  }
}

ordered_subfamilies.table =ordered_subfamilies.table[-1,]
# write.csv(label_table_final,file = "TEwide_group.list_2023_9_4.branchLength.csv")
# write.csv(ordered_subfamilies.table,file = "TEwide_group.list_2023_9_4.functional_group.csv")
# determined phyletic groups will be manually added to this table "TEwide_group.list_2023_9_4.functional_group_edited.csv"

######################### plot the liftover ratio between family and phyletic groups

ordered_subfamilies.table.final = read.csv("input/TEwide_group.list_2023_9_4.functional_group_edited.csv")

group_colors = c("PG1"="#0692C9","PG2"="#8638F4","PG3"="#06894A","PG4"="#608E06")

glist = list()
Order = 1
Tree = "hg19_G5_v1"
for (Tree in unique(ordered_subfamilies.table.final$tree)) {
  ordered_subfamilies.table.tmp = ordered_subfamilies.table.final[ordered_subfamilies.table.final$tree == Tree,]
  ordered_subfamilies.table.tmp.sum.g = data.frame(ordered_subfamilies.table.tmp %>% 
                                                     group_by(functional.group) %>% 
                                                     dplyr::summarise(intersect = sum(liftover_count.macque),total = sum(subfamily.count_200bp)))
  ordered_subfamilies.table.tmp.sum.g$perC = ordered_subfamilies.table.tmp.sum.g$intersect/ordered_subfamilies.table.tmp.sum.g$total
  ordered_subfamilies.table.tmp.sum.g$Type = "group"
  
  ordered_subfamilies.table.tmp.sum.f = data.frame(ordered_subfamilies.table.tmp %>% 
                                                     group_by(subfamily2_1) %>% 
                                                     dplyr::summarise(intersect = sum(liftover_count.macque),total = sum(subfamily.count_200bp)))
  ordered_subfamilies.table.tmp.sum.f$perC = ordered_subfamilies.table.tmp.sum.f$intersect/ordered_subfamilies.table.tmp.sum.f$total
  ordered_subfamilies.table.tmp.sum.f$Type = "subfamily"
  colnames(ordered_subfamilies.table.tmp.sum.f)[1] = "Group"
  colnames(ordered_subfamilies.table.tmp.sum.g)[1] = "Group"
  ordered_subfamilies.table.tmp.sum2 = rbind(ordered_subfamilies.table.tmp.sum.g,ordered_subfamilies.table.tmp.sum.f)
  head(ordered_subfamilies.table.tmp)
  # liftover between hg19 and macFas5
  p1_2 = ggplot(ordered_subfamilies.table.tmp, aes(x = subfamily2_1, y = perC.macque)) +
    geom_quasirandom(aes(fill = functional.group,size=subfamily.count_200bp),
                     method = "quasirandom",alpha=.7,shape=21) + 
    stat_summary(fun = "mean", 
                 fun.min = function(x)mean(x)-sd(x), 
                 fun.max = function(x)mean(x) + sd(x), 
                 geom = "pointrange", 
                 position = position_dodge(width = .9),shape=95,color="black",size=.5)+
    scale_fill_manual(values=group_colors)+
    xlab("") + 
    ggtitle(unique(ordered_subfamilies.table.tmp$subfamily2_1))+
    ylab("% instances shared with macaques") + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text=element_text(colour="black",size=rel(1),angle = 0),
          
          axis.text.x=element_text(colour="black",size=rel(1),angle = 90,hjust = 0.95,vjust = .2),
          axis.title=element_text(colour="black",size=rel(1)),
          #        legend.position=c(0.8,0.8),
          #        legend.position="bottom",
          legend.position="none",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1.2)))
  glist[[Order]] <-ggplotGrob(p1_2) 
  Order = Order + 1
  p1_2 = ggplot(ordered_subfamilies.table.tmp.sum2, aes(x=Group,y=perC*100,group=Type)) +
    geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
    geom_text(aes(label=round(perC*100,1)), vjust=-.2, color="black",
              position = position_dodge(0.9), size=3.5)+
    scale_fill_manual(values =c("#016A70","#A2C579"))+
    #geom_hline(yintercept=c(20,40), linetype="dashed", color = "black")+
    #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
    ggtitle(unique(ordered_subfamilies.table.tmp$subfamily2_1))+
    scale_x_discrete(drop = FALSE)+
    ylab("max % peaks-associated instances")+
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
          legend.position="none",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1))) 
  glist[[Order]] <-ggplotGrob(p1_2) 
  Order = Order + 1
}

pdf(paste("Figure_1F",".pdf",sep=""),    # create PNG for the heat map
    width = 12,        # 5 x 300 pixels
    height = 28,
    pointsize = 10)        # smaller font size
do.call("grid.arrange",c(glist,ncol=4))
dev.off()


###################### plot 2, similarity of LTR7 original families and our families
Tree = "hg19_G6_v2.mafft.prank.best.fas.nonrev_dna"
head(ordered_subfamilies.table.final)

ordered_subfamilies.table.LTR7 = ordered_subfamilies.table.final[ordered_subfamilies.table.final$tree == "Step21_hg19_G6_v1",]

## pair-wise similarity
Div_paired.plot = Div_paired[Div_paired$Group == "hg19_LTR7_dfam",]
Div_paired.plot = Div_paired.plot[(Div_paired.plot$from %in% c(ordered_subfamilies.table.LTR7$subfamily) & Div_paired.plot$to %in% Consensus_LTR7$sequence.ID) | 
                                    (Div_paired.plot$to %in% c(ordered_subfamilies.table.LTR7$subfamily) & Div_paired.plot$from %in% Consensus_LTR7$sequence.ID),]

Div_paired.plot$from = gsub("^h_","",Div_paired.plot$from)
Div_paired.plot$to = gsub("^h_","",Div_paired.plot$to)

Div_paired.plot$Y.axis = ifelse(Div_paired.plot$from %in% ordered_subfamilies.table.LTR7$subfamily,Div_paired.plot$from,Div_paired.plot$to)
Div_paired.plot$X.axis = ifelse(Div_paired.plot$from %in% ordered_subfamilies.table.LTR7$subfamily,Div_paired.plot$to,Div_paired.plot$from)
Div_paired.plot = merge(Div_paired.plot,Consensus_LTR7,by.x="X.axis",by.y="sequence.ID",all.x=T)

Div_paired.plot$Y.axis = factor(Div_paired.plot$Y.axis,levels = ordered_subfamilies.table.LTR7$subfamily)
Div_paired.plot = Div_paired.plot[!(Div_paired.plot$sequence.name %in% c("LTR7_Hera","LTR7_mMol","LTR7_mMyo","LTR7_pDis","LTR7_pKuh","LTR7_rAeg","LTR7_rFer")),]
Div_paired.plot$sequence.name = factor(Div_paired.plot$sequence.name,levels=c("LTR7C","LTR7B","LTR7BC","LTR7d2","LTR7d1","LTR7Y",
                                                                              "LTR7u2","LTR7u1","LTR7up1","LTR7up2",
                                                                              "LTR7a1","LTR7a2","LTR7B0","LTR7B2","LTR7B3","LTR7YY","LTR7","LTR7D3","LTR7up3"))
p1.similarity <- ggplot(Div_paired.plot, aes(x=sequence.name, y=Y.axis)) + 
  geom_tile(aes(fill=score)) + 
  scale_x_discrete(drop = FALSE)+
  scale_y_discrete(drop = FALSE)+
  scale_fill_distiller("score", palette = "Spectral", trans = "reverse")+
  xlab(NULL) + ylab(NULL) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#f0f0f0"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.95),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))


pdf(paste("Figure_S11C","_rerooted-LTR7-1.pdf",sep=""),    # create PNG for the heat map
    width = 6,        # 5 x 300 pixels
    height = 10,
    pointsize = 10)        # smaller font size
grid.draw(p1.similarity)
dev.off()

