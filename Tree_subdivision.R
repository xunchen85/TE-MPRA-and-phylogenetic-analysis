
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
library(Biostrings)
library(randomcoloR)

#########################
setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/Project_Neurogenesis/Final_edited_version_2022_11_25/Final_scripts/")

# 1.1 tree files
Families = list.files("input_trees/")
Families = Families[grepl("contree$",Families)]
Families = gsub(".mafft.prank.opt.gt99.fa.contree","",Families)

# 1.2 summary table
Summary_Table1 = read.csv("input/Summary_Table1_2022_8_9.csv")
Annotation_plot = Summary_Table1
Annotation_plot$uniqueID_new = Annotation_plot$instanceID.renamed
shape_species = c("Ancient" = 18,"Consensus" = 18,"hg19"=NA,"macFas5"=NA,"panTro4"=NA)

######################### step 1 
### 1.1 variables
Family = "hg19_rmsk_TE_MER11A_0bp"
N.leaf = 10
M.dist = 0.02
bootstrap.value.min = 95
min.dist.tip = 100

### 1.2 determine subtrees of each family
rowID = 1
for(Family in Families[!grepl("hg19",Families)]) {
  label_table.final.tip.sum.combined = data.frame("tmp1"=NA,"tmp2"=NA)
  
  ######## 1.2.1.2 or use the tree file
  tr <- read.tree(paste("input_trees//",Family,".mafft.prank.opt.gt99.fa.contree",sep=""))
  fasta.plot = readDNAStringSet(paste("input_trees/",Family,".mafft.prank.opt.gt99.fa",sep=""))
  tr$node.label = paste("Node",1:tr$Nnode,"/",tr$node.label,sep="")

  # 1.2.2 convert the FASTA files
  seq_name = names(fasta.plot)
  sequence = paste(fasta.plot)
  fasta.df <- data.frame(seq_name, sequence)
  
  # drop tips
  ######## 1.2.3 achieve the table                 
  ## organize
  label_table_originalOrder = data.frame(tr$edge)
  label_table_originalOrder$edgeOrder = 1:nrow(label_table_originalOrder)
  label_table_originalOrder$edgeID = paste(label_table_originalOrder$X1,label_table_originalOrder$X2)
  # label
  label_table = as_tibble(tr)
  label_table$Order = 1:nrow(label_table)
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
  
  label_table.final.tip = merge(label_table.final.tip,Annotation_plot[,c("uniqueID_new","TEfamily","species")],by.x="label",by.y="uniqueID_new",all.x=T)
  
  label_table.final.tip$cluster.final = factor(label_table.final.tip$cluster.final)
  # summary
  label_table.final.tip.sum = data.frame(label_table.final.tip %>% 
                                           group_by(cluster.final,TEfamily) %>% 
                                           summarise(n = n()) %>%
                                           mutate(sum=sum(n),freq = n / sum(n)))
  label_table.final.tip.sum = label_table.final.tip.sum[order(-label_table.final.tip.sum$freq),]
  label_table.final.tip.sum$label.final = paste(label_table.final.tip.sum$cluster.final,":",label_table.final.tip.sum$TEfamily,"(",label_table.final.tip.sum$sum,",",round(label_table.final.tip.sum$freq,2),")",sep="")
  label_table.final.tip.sum$combination = paste(Family,N.leaf,M.dist,bootstrap.value.min)
  if (nrow(label_table.final.tip.sum.combined) == 1) {
    label_table.final.tip.sum.combined = label_table.final.tip.sum
  } else {
    label_table.final.tip.sum.combined = rbind(label_table.final.tip.sum.combined,label_table.final.tip.sum)
  }
  label_table.final.tip.sum.uniq = label_table.final.tip.sum[!duplicated(label_table.final.tip.sum$cluster.final),]
  
  label_table.final.tip = merge(label_table.final.tip,label_table.final.tip.sum.uniq[,c("cluster.final","label.final")],by="cluster.final",all.x=T)
  
  # group info
  groupInfo <- split(label_table.final.tip$label, label_table.final.tip$label.final)
  ###### plot
  #cluster.color = colorRampPalette(brewer.pal(8, "Paired"))(n = length(unique(label_table.final.tip$label.final)))
  cluster.color = distinctColorPalette(length(unique(label_table.final.tip$label.final)))
  names(cluster.color) = unique(label_table.final.tip$label.final)
  #cluster.color[grepl("Ungrouped:",names(cluster.color))] = "#bdbdbd"      ### assign to grey color
  cluster.color[grepl("Ungrouped:",names(cluster.color))] = "black"      ### assign to grey color
  
  cluster.color.frame = data.frame(cluster.color)
  cluster.color.frame$label.final = rownames(cluster.color.frame)
  label_table.final.tip = merge(label_table.final.tip,cluster.color.frame,by="label.final",all.x=T)
  
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
  label_table.final.tip$is_removed = ifelse(label_table.final.tip$branch.length > min.dist.tip,"removed","kept")
  label_table.tipsRemoved = data.frame(label_table.final.tip[label_table.final.tip$is_removed == "removed",])
  label_table.final.tip = label_table.final.tip[,!colnames(label_table.final.tip) %in% c("node.info_1","node.info_2"),]
  label_table.final.tip$node.info = label_table.final.tip$cluster.final
  label_table.final.tip = data.frame(cSplit(label_table.final.tip,"node.info",sep="/",type.convert = as.character()))
  
  
  # remove the ones with long branch length
  tr.clean = drop.tip(tr,label_table.tipsRemoved$label)
  # not remove it
  tr.clean = tr
  # 1. standard tree with barplot 
  tr.clean = groupOTU(tr.clean, groupInfo)
  
  # 3. unrooted tree
  p_iris.un <- ggtree(tr.clean, layout = 'equal_angle',size=0.03,aes(color=group))+
    scale_color_manual(values= cluster.color)+
    guides(color = guide_legend(override.aes = list(size = 8)))
  p_iris.un = p_iris.un %<+% Annotation_plot[,c("uniqueID_new","species","TEfamily")]+
    new_scale_color() + 
    geom_tippoint(aes(shape=species,color=TEfamily),size=0.5)+
    scale_color_brewer(palette = "Set1")+
    scale_shape_manual(values = shape_species)
  
  if (Family == "hg19_rmsk_TE_MER11B_0bp") {
    pdf(paste("Figure_1C-",Family,"_un_final.pdf",sep=""),    # create PNG for the heat map
        width = 24,        # 5 x 300 pixels
        height = 24,
        pointsize = 10)        # smaller font size
    grid.draw(p_iris.un)
    dev.off() 
  } else if (Family %in% c("hg19_rmsk_TE_MER11A_0bp","hg19_rmsk_TE_MER11C_0bp","hg19_rmsk_TE_MER11D_0bp")){
    pdf(paste("Figure_S3C-",Family,"_un_final.pdf",sep=""),    # create PNG for the heat map
        width = 24,        # 5 x 300 pixels
        height = 24,
        pointsize = 10)        # smaller font size
    grid.draw(p_iris.un)
    dev.off() 
  } else {
    pdf(paste(Family,"_un_final.pdf",sep=""),    # create PNG for the heat map
        width = 24,        # 5 x 300 pixels
        height = 24,
        pointsize = 10)        # smaller font size
    grid.draw(p_iris.un)
    dev.off() 
  }
  
  # achieve the proportion
  label_table.final.tip.tmp = label_table.final.tip
  label_table.final.tip.tmp$label.final = gsub("\\(","\\,",label_table.final.tip.tmp$label.final)
  label_table.final.tip.tmp$label.final = gsub("\\:","\\,",label_table.final.tip.tmp$label.final)
  label_table.final.tip.tmp = data.frame(cSplit(label_table.final.tip.tmp,"label.final",sep=",",type.convert = as.character))
  label_table.final.tip$TEfamily.representative = label_table.final.tip.tmp$label.final_2
  label_table.final.tip$count = label_table.final.tip.tmp$label.final_3
  
  label_table.final.tip$perC = gsub("\\)","",label_table.final.tip.tmp$label.final_3)
  label_table.final.tip$consensus.name = paste(gsub("^MER","",label_table.final.tip$TEfamily.representative),gsub("Node","N",gsub("Ungrouped","U",label_table.final.tip$node.info_1)),label_table.final.tip$count,sep="|")
  label_table.final.tip$is_consensus = ifelse(grepl("MER|LTR|ERV",label_table.final.tip$label),"consensus","instance")
  
  write.csv(label_table.final.tip,file=paste(Family,"_group_final.csv",sep="")) 
  rm(label_table.final.tip.tmp)
  
  ### consensus sequence
  consensus_table = data.frame("sequence_name"=NA,"sequence"=NA)
  Order_consensus = 1
  
  for(labal.final.each in unique(label_table.final.tip$consensus.name)){
    sequence_name = label_table.final.tip[label_table.final.tip$consensus.name == labal.final.each & label_table.final.tip$is_consensus == "instance",]$label
    if(length(sequence_name) == 0){
      next
    }
    sequence = fasta.df[fasta.df$seq_name %in% sequence_name,]$sequence
    consensus_table[Order_consensus,]$sequence_name = paste(">",labal.final.each,sep="")
    consensus_table[Order_consensus,]$sequence = consensusString(sequence, ambiguityMap="N", threshold=0.51)
    Order_consensus = Order_consensus + 1
  }
  consensus_table$sequence = gsub("-","",consensus_table$sequence)
  consensus_table = consensus_table[!grepl("\\|U\\|",consensus_table$sequence_name),]
  write.table(consensus_table, sep="\n",row.names = FALSE, col.names = FALSE, quote = FALSE,
              file=paste(Family,"_consensus_noGap_noUnGrouped.fa",sep=""))
}

