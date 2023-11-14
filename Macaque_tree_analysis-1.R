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
library(ape)
library(dplyr)
library(grid)
library(gridExtra)
library(ggstance)
library(splitstackshape)
library(RColorBrewer)
library(ggtreeExtra)
library(ggnewscale)
library(castor)
#library(TreeTools)
library(Biostrings)
library(randomcoloR)

##########################
# 1.1 tree files
Families = list.files("input_trees/")
Families = Families[grepl("contree$",Families)]
Families = Families[grepl("^macFas5_rmsk",Families)]
Families = gsub(".mafft.prank.opt.gt99.fa.contree","",Families)

shape_species_cr = c("Ancient" = 18,"Consensus" = 18,"hg19"=NA,"macFas5"=NA,"panTro4"=NA)
color_species_cr = c("Ancient" = "black","Consensus" = "black","hg19"=NA,"macFas5"=NA,"panTro4"=NA)

######################### read annotation
Summary_Table1 = read.csv("input/Summary_Table1_2023_9_5.subfamilyInfo.csv")
Summary_Table2 = read.csv("input/Summary_Table2_2023_1_5.subfamilyInfo.csv")

# intersect group
intersect_Table = read.csv("input/Summary_Table1.intersect_reciprocal_hu_ch_ma_2023_2_20.csv")
intersect_Table.macFas5 = intersect_Table[intersect_Table$species.x == "macFas5",]
intersect_Table.macFas5$is_sameFamily = ifelse(intersect_Table.macFas5$is_intersect_reciprocal == "no_ortholog","no_ortholog",NA)
intersect_Table.macFas5$is_sameFamily = ifelse(intersect_Table.macFas5$is_intersect_reciprocal != "no_ortholog" & intersect_Table.macFas5$TEfamily.x == intersect_Table.macFas5$TEfamily.y & !is.na(intersect_Table.macFas5$TEfamily.y),"Consistent",intersect_Table.macFas5$is_sameFamily)
intersect_Table.macFas5$is_sameFamily = ifelse(intersect_Table.macFas5$is_intersect_reciprocal != "no_ortholog" & (intersect_Table.macFas5$TEfamily.x != intersect_Table.macFas5$TEfamily.y | is.na(intersect_Table.macFas5$TEfamily.y)),"Inconsistent",intersect_Table.macFas5$is_sameFamily)

######################### step 1 
### 1.1 variables
Family = "macFas5_rmsk_TE_MER11A_0bp"
N.leaf = 10
M.dist = 0.02
bootstrap.value.min = 95
min.dist.tip = "all"

### 1.2 determine subtrees of each family
rowID = 1
#for(Family in Families) {
for(Family in c("macFas5_rmsk_TE_MER11A_0bp")) {
    
  ########################## step 1: build the tree
  ######## 1.2.1.1 use the consensus tree
  tr <- read.tree(paste("input_trees/",Family,".mafft.prank.opt.gt99.fa.contree",sep=""))
  fasta.plot = readDNAStringSet(paste("input_trees/",Family,".mafft.prank.opt.gt99.fa",sep=""))
  tr$node.label = paste("Node",1:tr$Nnode,"/",tr$node.label,sep="")
  
  # 1.2.2 convert the FASTA files
  seq_name = names(fasta.plot)
  sequence = paste(fasta.plot)
  fasta.df <- data.frame(seq_name, sequence)
  
  # drop tips
  ######## 1.2.3 achieve the table      
  if(file.exists(paste("input_trees/",Family,"_",N.leaf,"_",M.dist,"_",min.dist.tip,"_",bootstrap.value.min,"_group_final.csv",sep=""))){
    label_table.final.tip = read.csv(file=paste("input_trees/",Family,"_",N.leaf,"_",M.dist,"_",min.dist.tip,"_",bootstrap.value.min,"_group_final.csv",sep="")) 
    ## group
    groupInfo <- split(label_table.final.tip$label, label_table.final.tip$label.final.correctedName)
    ###### alter the name if it is not correct
    label_table.final.tip.color = label_table.final.tip[!duplicated(label_table.final.tip$cluster.color),]
    cluster.color = label_table.final.tip.color$cluster.color
    names(cluster.color) = label_table.final.tip.color$label.final.correctedName
  } else {
    label_table.final.tip.sum.combined = data.frame("tmp1"=NA,"tmp2"=NA)
    ## organize
    label_table_originalOrder = data.frame(tr$edge)
    label_table_originalOrder$edgeOrder = 1:nrow(label_table_originalOrder)
    label_table_originalOrder$edgeID = paste(label_table_originalOrder$X1,label_table_originalOrder$X2)
    # label
    label_table = as_tibble(tr)
    label_table$Order = 1:nrow(label_table)
    #label_table$cluster = NA
    #label_table$collapse = NA
    label_table$edgeID = paste(label_table$parent,label_table$node)
    label_table = merge(label_table,label_table_originalOrder[,c("edgeID","edgeOrder")],by="edgeID",all.x=T)
    label_table = label_table[order(label_table$Order),]
    label_table$node.info = label_table$label
    label_table = cSplit(label_table,"node.info",sep="/",type.convert = as.character)
    label_table$node.order = as.numeric(gsub("^Node","",label_table$node.info_1))
    label_table$node.bootstrap = as.numeric(as.character(label_table$node.info_2))
    rm(label_table_originalOrder)
    label_table.full = data.frame(label_table)
    if (nrow(label_table)<=20){
      next
    }
    ######## 1.2.3 select clusters with the thresholds
    # exclude tips in
    # label_table = label_table.full[label_table.full$branch.length<= min.dist.tip | grepl("^Node",label_table.full$label),]
    # tr.clean = drop.tip(tr,label_table.full[label_table.full$branch.length>min.dist.tip & !grepl("^Node",label_table.full$label),]$label)
    # compute the distances
    tr.dis = data.frame(dist.nodes(tr))
    colnames(tr.dis) = label_table.full$label
    row.names(tr.dis) = label_table.full$label
    
    tr.dis = tr.dis[rownames(tr.dis) %in% label_table$label,colnames(tr.dis) %in% label_table$label]
    # Inspect each node
    label_table.node = data.frame(label_table[grepl("^Node",label_table$label),])
    label_table.node = label_table.node[order(-label_table.node$node.order),]
    label_table.node$Nleaf.all = NA
    label_table.node$Nleaf = NA
    label_table.node$Sleaf = NA
    label_table.node$Mdist = NA
    label_table.node$Mdist.parentNode = NA
    label_table.node$isKept = NA
    
    ## obtain the qualified nodes 
    label_table.node[label_table.node$node.info_1 == "Node1",]$Mdist.parentNode = 1
    label_table.node[label_table.node$node.info_1 == "Node1",]$node.bootstrap = 100
    label_table.node[label_table.node$node.info_1 == "Node1",]$Mdist = 1
    label_table.node = label_table.node[!is.na(label_table.node$label),]
    
    Grouped.tips = c()
    UnGrouped.tips = c()
    
    # remove nodes with fewer than minimum tips
    for (Node.rowID in 1:nrow(label_table.node)) {
      tree_sub = get_subtree_at_node(tr, label_table.node[Node.rowID,]$label)
      list_labels = tree_sub$subtree$tip.label
      if (label_table.node[Node.rowID,]$node.info_1 != "Node1"){
        label_table.node[Node.rowID,]$Mdist = min(tr.dis[label_table.node[Node.rowID,]$label,!(colnames(tr.dis) %in% c(tree_sub$subtree$tip.label,label_table.node$label))])
      }
      label_table.node[Node.rowID,]$Nleaf.all = length(list_labels)
      label_table.node[Node.rowID,]$Mdist.parentNode = min(tr.dis[label_table.node[Node.rowID,]$label,(colnames(tr.dis) %in% label_table.node[((Node.rowID+1):nrow(label_table.node)),]$label)])
    }
    label_table.node = label_table.node[label_table.node$Nleaf.all>=N.leaf,]
    label_table.node = label_table.node[label_table.node$node.bootstrap>=bootstrap.value.min,]
    
    # determine the subtrees
    for (Node.rowID in 1:nrow(label_table.node)) {
      tree_sub = get_subtree_at_node(tr, label_table.node[Node.rowID,]$label)
      list_labels = tree_sub$subtree$tip.label
      if (label_table.node[Node.rowID,]$Mdist < M.dist |
          length(list_labels) < N.leaf |
          label_table.node[Node.rowID,]$node.bootstrap < bootstrap.value.min){  ## if the minimum distance to other instances is too small or have fewer tips
        next
      } else if (!any(list_labels %in% c(Grouped.tips,UnGrouped.tips))) {
        if (length(list_labels) < N.leaf) {
          next
        } else {
          label_table.node[Node.rowID,]$Nleaf = length(list_labels)
          label_table.node[Node.rowID,]$Sleaf = as.character(paste(list_labels,collapse = ' '))
          label_table.node[Node.rowID,]$isKept = "kept"
          Grouped.tips = c(Grouped.tips,list_labels)
        }
      } else if (any(list_labels %in% c(Grouped.tips,UnGrouped.tips))){  ## if it is larger subtree
        if (length(list_labels[!(list_labels %in% c(Grouped.tips,UnGrouped.tips))]) >= N.leaf){
          label_table.node[Node.rowID,]$Nleaf = length(list_labels[!(list_labels %in% c(Grouped.tips,UnGrouped.tips))])
          label_table.node[Node.rowID,]$Sleaf = as.character(paste(list_labels[!(list_labels %in% c(Grouped.tips,UnGrouped.tips))],collapse = ' '))
          label_table.node[Node.rowID,]$isKept = "kept"
          Grouped.tips = c(Grouped.tips,list_labels[!(list_labels %in% c(Grouped.tips,UnGrouped.tips))])
        } else {
          UnGrouped.tips = c(UnGrouped.tips,list_labels[!(list_labels %in% c(Grouped.tips,UnGrouped.tips))])
        }
      } 
    }
    # UnGrouped.tips = UnGrouped.tips[!UnGrouped.tips %in% Grouped.tips]
    
    # kept cluster
    label_table.node.kept = label_table.node[label_table.node$isKept == "kept" & !is.na(label_table.node$isKept),]
    label_table.node.kept  = label_table.node.kept[order(-label_table.node.kept$edgeOrder),]
    
    ########## 1.2.4 Assign clusters
    label_table.final = label_table
    label_table.final$cluster.final = "Ungrouped"
    Row = 1
    for (Row in 1:nrow(label_table.node.kept)){
      label_table.final$cluster.final = ifelse(label_table.final$label %in% as.list(strsplit(label_table.node.kept[Row,]$Sleaf, '\\s+'))[[1]] &
                                                 !(label_table.final$label %in% UnGrouped.tips),label_table.node.kept[Row,]$label,label_table.final$cluster.final)
    }
    
    # color and min.dist.tip
    label_table.final.tip = data.frame(label_table.final[!grepl("^Node",label_table.final$label),])
    
    # label_table.final.tip = merge(label_table.final.tip,Annotation_plot[,c("uniqueID_new","TEfamily","species")],by.x="label",by.y="uniqueID_new",all.x=T)
    
    label_table.final.tip$cluster.final = factor(label_table.final.tip$cluster.final)
    # summary
    label_table.final.tip.sum = data.frame(label_table.final.tip %>% 
                                             group_by(cluster.final) %>% 
                                             summarise(n = n()) %>%
                                             mutate(sum=sum(n),freq = n / sum(n)))
    label_table.final.tip.sum.nU = label_table.final.tip.sum[order(-label_table.final.tip.sum$freq),]
    label_table.final.tip.sum.U = label_table.final.tip.sum[grepl("Ungrouped",label_table.final.tip.sum$cluster.final),]
    label_table.final.tip.sum = rbind(label_table.final.tip.sum.nU[!grepl("Ungrouped",label_table.final.tip.sum.nU$cluster.final),],label_table.final.tip.sum.U)
    #label_table.final.tip.sum$label.final = paste(label_table.final.tip.sum$cluster.final,":","(",label_table.final.tip.sum$sum,",",round(label_table.final.tip.sum$freq,2),")",sep="")
    label_table.final.tip.sum$label.final = paste(label_table.final.tip.sum$cluster.final,":","(",label_table.final.tip.sum$sum,",",label_table.final.tip.sum$n,")",sep="")
    label_table.final.tip.sum$label.final.corrected = paste(gsub("hg19_rmsk_|_0bp|macFas5_rmsk_TE_0bp_|.family3.rename_200bp","",Family),"_g",1:nrow(label_table.final.tip.sum),":","(",label_table.final.tip.sum$sum,",",label_table.final.tip.sum$n,")",sep="")
    label_table.final.tip.sum$label.final.corrected = ifelse(grepl("Ungrouped",label_table.final.tip.sum$label.final),
                                                             paste(gsub("hg19_rmsk_|_0bp|macFas5_rmsk_TE_0bp_|.family3.rename_200bp","",Family),"_U",":","(",label_table.final.tip.sum$sum,",",label_table.final.tip.sum$n,")",sep=""),
                                                             label_table.final.tip.sum$label.final.corrected)
    label_table.final.tip.sum$label.final.correctedName = paste(gsub("hg19_rmsk_|_0bp|macFas5_rmsk_TE_0bp_|.family3.rename_200bp","",Family),"_g",1:nrow(label_table.final.tip.sum),sep="")
    label_table.final.tip.sum$label.final.correctedName = ifelse(grepl("Ungrouped",label_table.final.tip.sum$label.final),
                                                                 paste(gsub("hg19_rmsk_|_0bp|macFas5_rmsk_TE_0bp_|.family3.rename_200bp","",Family),"_U",sep=""),
                                                                 label_table.final.tip.sum$label.final.correctedName)
    label_table.final.tip.sum$combination = paste(Family,N.leaf,M.dist,bootstrap.value.min)
    if (nrow(label_table.final.tip.sum.combined) == 1) {
      label_table.final.tip.sum.combined = label_table.final.tip.sum
    } else {
      label_table.final.tip.sum.combined = rbind(label_table.final.tip.sum.combined,label_table.final.tip.sum)
    }
    label_table.final.tip.sum.uniq = label_table.final.tip.sum[!duplicated(label_table.final.tip.sum$cluster.final),]
    
    label_table.final.tip = merge(label_table.final.tip,label_table.final.tip.sum.uniq[,c("cluster.final","label.final","label.final.corrected","label.final.correctedName","n")],by="cluster.final",all.x=T)
    
    # group info
    groupInfo <- split(label_table.final.tip$label, label_table.final.tip$label.final.correctedName)
    ###### plot
    cluster.color = distinctColorPalette(length(unique(label_table.final.tip$label.final.correctedName)))
    names(cluster.color) = unique(label_table.final.tip$label.final.correctedName)
    cluster.color[grepl("_U",names(cluster.color))] = "black"      ### assign to grey color
    
    cluster.color.frame = data.frame(cluster.color)
    cluster.color.frame$label.final.correctedName = rownames(cluster.color.frame)
    label_table.final.tip = merge(label_table.final.tip,cluster.color.frame,by="label.final.correctedName",all.x=T)
    
    # assign
    label_table.final.tip$branch.length.ToKeptNode = NA
    rowTip = 1
    for(rowTip in 1:nrow(label_table.final.tip)){
      if (label_table.final.tip[rowTip,]$cluster.final %in% colnames(tr.dis)){
        label_table.final.tip[rowTip,]$branch.length.ToKeptNode = tr.dis[label_table.final.tip[rowTip,]$label,label_table.final.tip[rowTip,]$cluster.final]
      } else {
        label_table.final.tip[rowTip,]$branch.length.ToKeptNode = label_table.final.tip[rowTip,]$branch.length
      }
    }
    
    #for (min.dist.tip in c(0.5,1,100)){
    # exclude outliers in the figures
    if (min.dist.tip == "all"){
      min.dist.tip2 = 100
    }
    label_table.final.tip$is_removed = ifelse(label_table.final.tip$branch.length > min.dist.tip2,"removed","kept")
    label_table.tipsRemoved = data.frame(label_table.final.tip[label_table.final.tip$is_removed == "removed",])
    label_table.final.tip = label_table.final.tip[,!colnames(label_table.final.tip) %in% c("node.info_1","node.info_2"),]
    label_table.final.tip$node.info = label_table.final.tip$cluster.final
    label_table.final.tip = data.frame(cSplit(label_table.final.tip,"node.info",sep="/",type.convert = as.character()))
    # achieve the proportion
    label_table.final.tip.tmp = label_table.final.tip
    label_table.final.tip.tmp$label.final = gsub("\\(","\\,",label_table.final.tip.tmp$label.final)
    label_table.final.tip.tmp$label.final = gsub("\\:","\\,",label_table.final.tip.tmp$label.final)
    label_table.final.tip.tmp = data.frame(cSplit(label_table.final.tip.tmp,"label.final",sep=",",type.convert = as.character))
    label_table.final.tip$TEfamily.representative = label_table.final.tip.tmp$label.final.correctedName
    label_table.final.tip$count = label_table.final.tip.tmp$label.final_3
    
    label_table.final.tip$perC = gsub("\\)","",label_table.final.tip.tmp$label.final_4)
    label_table.final.tip$consensus.name = label_table.final.tip$TEfamily.representative
    # write to the table
    write.csv(label_table.final.tip,file=paste(Family,"_",N.leaf,"_",M.dist,"_",min.dist.tip,"_",bootstrap.value.min,"_group_final.csv",sep="")) 
  }
  
  ############# add liftover info
  # add the intersect group
  Summary_Table1.plot = Summary_Table1[Summary_Table1$species %in% c("macFas5","Ancient"),]
  Summary_Table1.plot = merge(Summary_Table1.plot,intersect_Table.macFas5[,c("instanceID.renamed","is_intersect_reciprocal","is_sameFamily")],all.x=T,by="instanceID.renamed")
  Summary_Table1.plot$is_sameFamily = ifelse(Summary_Table1.plot$species == "Ancient","consensus",Summary_Table1.plot$is_sameFamily)
  shape.intersect = c("no_ortholog"=NA,"Consistent"=16,"Inconsistent"=16,"consensus"=18)
  color.intersect = c("no_ortholog"=NA,"Consistent"="black","Inconsistent"="grey","consensus"="black")
  
  # remove the ones with long branch length
  #tr.clean = drop.tip(tr,label_table.tipsRemoved$label)
  # not remove it
  tr.clean = tr
  
  # 1. standard tree with barplot 
  tr.clean = groupOTU(tr.clean, groupInfo)
  # p_iris.sd <- ggtree(tr.clean,size=0.1,aes(color=group)) +
  #   #geom_label(aes(label=node),size=1)+
  #   scale_color_manual(values= cluster.color) +
  #   guides(color = guide_legend(override.aes = list(size = 8)))
  # p_iris.sd = p_iris.sd %<+% Summary_Table1.plot[,c("instanceID.renamed","is_sameFamily")] + 
  #   new_scale_color() + 
  #   geom_tippoint(aes(shape=is_sameFamily,color=is_sameFamily),size=1)+
  #   scale_color_manual(values=color.intersect)+
  #   scale_shape_manual(values = shape.intersect)
  
  # 3. unrooted tree
  p_iris.un <- ggtree(tr.clean, layout = 'equal_angle',size=0.03,aes(color=group))+
    scale_color_manual(values= cluster.color)+
    guides(color = guide_legend(override.aes = list(size = 8)))
  # add annotation of liftover
  p_iris.un = p_iris.un %<+% Summary_Table1.plot[,c("instanceID.renamed","is_sameFamily")] + 
    new_scale_color() + 
    geom_tippoint(aes(shape=is_sameFamily,color=is_sameFamily),size=1)+
    scale_color_manual(values=color.intersect)+
    scale_shape_manual(values = shape.intersect)
  
  # 4. circular tree
  # p_iris.cr <- ggtree(tr.clean, layout = 'circular',size=0.03,aes(color=group))+
  #   scale_color_manual(values= cluster.color)+
  #   guides(color = guide_legend(override.aes = list(size = 8)))
  
  
  pdf(paste("Figure_4B-",Family,"_",N.leaf,"_",M.dist,"_",min.dist.tip,"_",bootstrap.value.min,"_un_final.pdf",sep=""),    # create PNG for the heat map
      width = 24,        # 5 x 300 pixels
      height = 24,
      pointsize = 10)        # smaller font size
  grid.draw(p_iris.un)
  dev.off() 
  
  ########################## step 2: generate the consensus sequences
  ### remove the original consensus if they have
  label_table.final.tip.noConsensus = label_table.final.tip[grepl("T_",label_table.final.tip$label),]
  
  ### consensus sequence
  consensus_table = data.frame("sequence_name"=NA,"sequence"=NA)
  Order_consensus = 1
  
  for(labal.final.each in unique(label_table.final.tip.noConsensus$consensus.name)){
    sequence_name = label_table.final.tip.noConsensus[label_table.final.tip.noConsensus$consensus.name == labal.final.each,]$label
    if(length(sequence_name) == 0){
      next
    }
    sequence = fasta.df[fasta.df$seq_name %in% sequence_name,]$sequence
    consensus_table[Order_consensus,]$sequence_name = paste(">",labal.final.each,sep="")
    consensus_table[Order_consensus,]$sequence = consensusString(sequence, ambiguityMap="N", threshold=0.51)
    Order_consensus = Order_consensus + 1
  }
  # assign them to "N"
  ## kept ambiguous DNAs
  #write.table(consensus_table, sep="\n",row.names = FALSE, col.names = FALSE, quote = FALSE,
  #            file=paste("Step19_",Family,"_",N.leaf,"_",M.dist,"_",min.dist.tip,"_",bootstrap.value.min,"_consensus.fa",sep=""))
  ## remove Gaps
  consensus_table$sequence = gsub("-","",consensus_table$sequence)
  #write.table(consensus_table, sep="\n",row.names = FALSE, col.names = FALSE, quote = FALSE,
  #            file=paste("Step19_",Family,"_",N.leaf,"_",M.dist,"_",min.dist.tip,"_",bootstrap.value.min,"_consensus_noGap.fa",sep=""))
  consensus_table = consensus_table[!grepl("_U$",consensus_table$sequence_name),]
  write.table(consensus_table, sep="\n",row.names = FALSE, col.names = FALSE, quote = FALSE,
              file=paste("Figure4B-",Family,"_",N.leaf,"_",M.dist,"_",min.dist.tip,"_",bootstrap.value.min,"_consensus_noGap_noUngrouped.fa",sep=""))
  
  
}


