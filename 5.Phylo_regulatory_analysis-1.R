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
library(ape)
library(dplyr)
library(grid)
library(gridExtra)
library(ggstance)
library(splitstackshape)
library(RColorBrewer)
library(ggtreeExtra)
library(ggnewscale)
library("Biostrings")
library(randomcoloR)

#########################
Family="MER11A_F1"
color_species = c("Ancient" = "#6a3d9a","Consensus" = "#6a3d9a","hg19"="#e31a1c","macFas5"="#1f78b4","panTro4"="#33a02c")
shape_species = c("Ancient" = 18,"Consensus" = 18,"hg19"=16,"macFas5"=16,"panTro4"=16)

shape_species_cr = c("Ancient" = 18,"Consensus" = 18,"hg19"=NA,"macFas5"=NA,"panTro4"=NA)
color_species_cr = c("Ancient" = "black","Consensus" = "black","hg19"=NA,"macFas5"=NA,"panTro4"=NA)

######################### read annotation
# updated on 2023/9/5
Summary_Table1 = read.csv("input/Summary_Table1_2023_9_5.subfamilyInfo.csv")
Summary_Table1 = Summary_Table1[,!colnames(Summary_Table1) %in% c("X","X.x","X.y")]

Summary_Table2 = read.csv("input/Summary_Table2_2023_1_5.subfamilyInfo.csv")
Summary_Table2 = Summary_Table2[,!colnames(Summary_Table2) %in% c("X","X.x","X.y")]

colnames(Summary_Table2)[which(colnames(Summary_Table2) == "Instance_coordinate.x")] = "Instance_coordinate"

MER11_distance = Summary_Table1[Summary_Table1$TEfamily.x %in% c("MER11A","MER11B","MER11C"),c("subfamily.name.final","branch.length.ToKeptNode")]
colnames(MER11_distance)[2] = "distanceToNode1"

Summary_Table2$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2$iPSC.alpha.Zscore.isActive == "None","No activity",NA)
Summary_Table2$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2$iPSC.alpha.Zscore.isActive != "None" & 
                                                           Summary_Table2$iPSC.alpha.Zscore < 2,"<2", Summary_Table2$iPSC.alpha.Zscore.isActive.group)
Summary_Table2$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2$iPSC.alpha.Zscore.isActive != "None" & 
                                                           Summary_Table2$iPSC.alpha.Zscore >= 2 & Summary_Table2$iPSC.alpha.Zscore < 4 ,"2-4" ,Summary_Table2$iPSC.alpha.Zscore.isActive.group)
Summary_Table2$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2$iPSC.alpha.Zscore.isActive != "None" & 
                                                           Summary_Table2$iPSC.alpha.Zscore >= 4 & Summary_Table2$iPSC.alpha.Zscore < 6 ,"4-6", Summary_Table2$iPSC.alpha.Zscore.isActive.group)
Summary_Table2$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2$iPSC.alpha.Zscore.isActive != "None" & 
                                                           Summary_Table2$iPSC.alpha.Zscore >= 6 & Summary_Table2$iPSC.alpha.Zscore < 8 ,"6-8", Summary_Table2$iPSC.alpha.Zscore.isActive.group)
Summary_Table2$iPSC.alpha.Zscore.isActive.group = ifelse(Summary_Table2$iPSC.alpha.Zscore.isActive != "None" & 
                                                           Summary_Table2$iPSC.alpha.Zscore >= 8 , ">=8", Summary_Table2$iPSC.alpha.Zscore.isActive.group)

Summary_Table2$iPSC.alpha.Zscore.isActive.group = factor(Summary_Table2$iPSC.alpha.Zscore.isActive.group,levels=c(">=8","6-8","4-6","2-4","<2","No activity"))
iPSC.alpha.Zscore.isActive.group.color = c(">=8" = "#67001f","6-8"="#b2182b","4-6"="#d6604d","2-4"="#f4a582","<2"="#fddbc7","No activity"="#f0f0f0","Missing"="#e0ecf4")

######################### step 1 achieve the list of tree files
# 1.1 tree files
Families = list.files("input/")
Families = Families[grepl("treefile$",Families)]
Families = gsub(".mafft.prank.opt.gt99.fa.treefile","",Families)
Families = Families[!grepl("LTR7only|mafft.opt.gt99",Families)]         #### remove LTR7only tree file
Family = "hg19_MER11_combined.instance"

# 1.2 add ENCODE TF and H293K KRAB-ZNF profiles         ######
Roadmap = read.delim("input/hg19_instance_final_2022_12_28.combined.bed",header=F,sep="")
colnames(Roadmap) = c("Mark","Frame_chr","Start","End","uniqueID_new","uniqueID","strand")
Roadmap$Mark = gsub("_hg19","",Roadmap$Mark)
Roadmap$Mark = gsub("-rep",".",Roadmap$Mark)

# retrieve DNase (updated on 2022/12/29 using the new BEDtools parameters -e -f 0.5 -F 0.5)
Roadmap.original.instance = read.delim("input/hg19_instance_final_2022_12_29.roadmap.combined.bed",header=F,sep="")
Roadmap.original = Roadmap.original.instance[Roadmap.original.instance$V1 %in% c("DNase","H3K9me3"),]
colnames(Roadmap.original) = c("Mark","Frame_chr","Start","End","uniqueID_new","uniqueID","strand")

# combine
Roadmap = rbind(Roadmap,Roadmap.original)

# 
Kept_mark  = read.csv("input/ENCODE_kept_samples_2022_10_30.csv")
Kept_mark$Mark.cell = paste(Kept_mark$Experiment.target,Kept_mark$cell)
Roadmap = merge(Roadmap,Kept_mark[c("sampleID","Experiment.target","Mark.cell")],by.x="Mark",by.y="sampleID",all.x=T)
Roadmap = Roadmap[!is.na(Roadmap$Experiment.target) | !grepl("ENCFF",Roadmap$Mark),]
Roadmap$Experiment.target = ifelse(is.na(Roadmap$Experiment.target),Roadmap$Mark,Roadmap$Experiment.target)
Roadmap$Experiment.target = ifelse(is.na(Roadmap$Experiment.target),Roadmap$Mark,Roadmap$Experiment.target)

Roadmap$Mark.cell = ifelse(is.na(Roadmap$Mark.cell),Roadmap$Experiment.target,Roadmap$Mark.cell)
Roadmap$Mark.cell = gsub("H3K9me3 H1$","H3K9me3_H1",Roadmap$Mark.cell)
Roadmap$Mark.cell = gsub(" H1$","",Roadmap$Mark.cell)
Roadmap$Pos = Roadmap$Mark.cell

######################### step 2 plot the circular plot at the instance level
### plot
Families = c("hg19_rmsk_TE_MER11A_0bp","hg19_rmsk_TE_MER11B_0bp","hg19_rmsk_TE_MER11C_0bp")
Family = "hg19_rmsk_TE_MER11A_0bp"
bootstrap = 95
TreeType = "contree"
N.leaf = 10
M.dist = 0.02

for(Family in Families) {
  ######## 2.1.2 use the consensus tree
  tr <- read.tree(paste("input_trees/",Family,".mafft.prank.opt.gt99.fa.",TreeType,sep=""))
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
  if (nrow(label_table) <=20){
    next
  }

  ## add the grouping info
  Annotation_plot = Summary_Table1
  Annotation_plot$uniqueID_new = Annotation_plot$instanceID.renamed

  ######## format the annotation file
  ## step 3.1: labels, groups (instance level)
  # 1bp_coordinate
  Annotation_plot$Instance_coordinate_species_1bp = paste(Annotation_plot$species,":",Annotation_plot$chr,":",Annotation_plot$start+1,"-",Annotation_plot$end,sep="")
  Annotation_plot$Instance_coordinate_species_1bp = ifelse(Annotation_plot$species == "Ancient",Annotation_plot$Instance_coordinate_species,Annotation_plot$Instance_coordinate_species_1bp)
  
  # filtering
  Annotation_plot = Annotation_plot[Annotation_plot$instanceID.renamed %in% tr$tip.label,]
  #Annotation_plot = Annotation_plot[grepl(gsub("hg19_rmsk_TE_0bp_MER|.family3.rename_200bp","",Family),Annotation_plot$consensus.name),]
  
  #Annotation_plot = merge(Annotation_plot,Groupings_tree[,!(colnames(Groupings_tree) %in% c("TEfamily","species"))],by.x="uniqueID_new",by.y="label",all.x=T)
  Annotation_plot$clusterID = Annotation_plot$cluster.final
  Annotation_plot$collapse_col = Annotation_plot$cluster.color

  # color
  col_cluster = Annotation_plot[!duplicated(Annotation_plot$clusterID),]$collapse_col
  names(col_cluster) = Annotation_plot[!duplicated(Annotation_plot$clusterID),]$clusterID
  
  #col_family = colorRampPalette(brewer.pal(8, "Set1"))(n = length(unique(Annotation_plot$TEfamily)))
  col_family = distinctColorPalette(length(unique(Annotation_plot$TEfamily.x)))
  names(col_family) = unique(Annotation_plot$TEfamily)
  
  col_cluster_family = c(col_cluster,col_family)
  col_cluster_family['Ungrouped'] = "white"
  # species
  Annotation_plot$species = ifelse(Annotation_plot$species %in% c("hg19","panTro4","macFas5"),Annotation_plot$species,"Consensus")
  Annotation_plot$species = factor(Annotation_plot$species)
  
  Annotation_plot = Annotation_plot[match(tr$tip.label,Annotation_plot$uniqueID_new),]
  # exclude NA values
  Annotation_plot = Annotation_plot[!is.na(Annotation_plot$uniqueID_new),]
  rownames(Annotation_plot) = Annotation_plot$uniqueID_new
  
  ### Info of families
  # clusterID
  head(Annotation_plot)
  Annotation_plot_formatted_1 = Annotation_plot[,c("uniqueID_new","clusterID")]
  colnames(Annotation_plot_formatted_1) = c("uniqueID_new","Pos")
  Annotation_plot_formatted_1$Group = "cluster"
  
  # family
  Annotation_plot_formatted_2 = Annotation_plot[,c("uniqueID_new","TEfamily.x")]
  colnames(Annotation_plot_formatted_2) = c("uniqueID_new","Pos")
  Annotation_plot_formatted_2$Group = "family"
  
  # combined family and cluster
  Annotation_plot_formatted_3 = Annotation_plot_formatted_1
  Annotation_plot_formatted_3$Group = "cluster2"
  Annotation_plot_formatted_species = rbind(Annotation_plot_formatted_1,Annotation_plot_formatted_3)
  
  # div
  if (grepl("MER11",Family)){
    Annotation_plot2 = merge(Annotation_plot,MER11_distance[,c("subfamily.name.final","distanceToNode1")],by="subfamily.name.final",all.x=T)
    Annotation_plot_formatted_3 = Annotation_plot2[,c("uniqueID_new","distanceToNode1")]
    Annotation_plot_formatted_3 = Annotation_plot_formatted_3[!grepl("MER11",Annotation_plot_formatted_3$uniqueID_new),]
  } else {
    Annotation_plot_formatted_3 = Annotation_plot[,c("uniqueID_new","div_rate")]
  }
  colnames(Annotation_plot_formatted_3) = c("uniqueID_new","divergent_rate")
  Annotation_plot_formatted_3$Group = "div"
  Annotation_plot_formatted_div = Annotation_plot_formatted_3
  
  # make a duplicated layer
  Annotation_plot_formatted_div2 = Annotation_plot_formatted_div
  Annotation_plot_formatted_div2$Group = "div2"
  Annotation_plot_formatted_div = rbind(Annotation_plot_formatted_div,Annotation_plot_formatted_div2)
  
  ## grouping and color
  groupInfo <- split(Annotation_plot$uniqueID_new, Annotation_plot$clusterID)
  
  ############ plot
  # grouping
  tr = groupOTU(tr, groupInfo)
  
  # circular
  p_iris <- ggtree(tr, layout = 'fan',size=0.1,aes(color=group))+
    #p_iris <- ggtree(tr, layout = 'fan',size=0.1,aes(color=group),open.angle = 10)+
    scale_color_manual(values=col_cluster)+
    guides(color = guide_legend(override.aes = list(size = 8)))
  colnames(Annotation_plot)
  p_iris = p_iris %<+% Annotation_plot[,c("uniqueID_new","species","TEfamily.x")] + 
    new_scale_color() + 
    geom_tippoint(aes(shape=species,color=species),size=0.8)+
    scale_color_manual(values=color_species_cr) +
    scale_shape_manual(values = shape_species_cr)
  
  # family and species
  p_iris2 = p_iris +  
    new_scale_fill() + 
    geom_fruit(
      data=Annotation_plot_formatted_species,
      geom=geom_tile,
      mapping=aes(y=uniqueID_new,x=Group,fill=Pos),
      pwidth=0.02, # width of the external layer, default is 0.2 times of x range of tree.
      offset=0.02,
      axis.params=list(
        text.size=3,
        axis="none", # add axis text of the layer.
        text.angle=90, # the text angle of x-axis.
        hjust=0.5,  # adjust the horizontal position of text of axis.
        vjust=0
      )
    ) + 
    scale_fill_manual(values=col_cluster_family)
  #scale_fill_manual(values=col_family)
  
  # divergent rate
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  p_iris3 = p_iris +      # remove the family info
    new_scale_fill() + 
    geom_fruit(
      data=Annotation_plot_formatted_div,
      geom=geom_tile,
      mapping=aes(y=uniqueID_new,x=Group,fill=divergent_rate),
      pwidth=0.02, # width of the external layer, default is 0.2 times of x range of tree.
      offset=0.02,
      axis.params=list(
        text.size=3,
        axis="none", # add axis text of the layer.
        text.angle=90, # the text angle of x-axis.
        hjust=0.5,  # adjust the horizontal position of text of axis.
        vjust=0
      )
    ) + scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = median(unique(Annotation_plot_formatted_div$divergent_rate),na.rm=T),na.value = 'black')
  
  
  # epigenetic mark
  # select epigenetic marks
  if (grepl("MER11",Family)){
    Roadmap.plot = Roadmap[Roadmap$Mark.cell %in% c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","H3K9me3","H3K9me3_H1","H3K27me3",
                                                    "BCL11A","RXRA","POU5F1","EP300","TCF12","SP1","TEAD4","YY1","CEBPB","USF1","NANOG","ZNF143",
                                                    "KAP1","ZNF525","ZNF808","ZNF440","ZNF433","ZNF33A.2","ZNF468","ZNF611"),]
    Roadmap.plot$Pos = factor(Roadmap.plot$Pos,levels=c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","H3K9me3","H3K9me3_H1","H3K27me3",
                                                        "BCL11A","RXRA","POU5F1","NANOG","EP300","YY1","TCF12","SP1","TEAD4","CEBPB","USF1","ZNF143",
                                                        "ZNF525","ZNF808","ZNF440","ZNF433","ZNF33A.2","ZNF468","ZNF611","KAP1"))
    colMark = c("DNase" = "#b2df8a","H3K27ac" = "#33a02c","H3K4me1" = "#33a02c","H3K4me2" = "#33a02c","H3K4me3" = "#33a02c","H3K9me3" = "#bf812d","H3K9me3_H1" = "#bf812d", "H3K27me3" = "#bf812d",
                "BCL11A"="0868ac","RXRA"="0868ac","POU5F1"="0868ac","EP300" = "#0868ac","TCF12" = "#0868ac","SP1" = "#0868ac","TEAD4" = "#0868ac","YY1" = "#0868ac","CEBPB" = "#0868ac","USF1" = "#0868ac","NANOG" = "#0868ac","ZNF143" = "#0868ac",
                "ZNF525"="#4eb3d3","ZNF808"="#4eb3d3","ZNF440"="#4eb3d3","ZNF433"="#4eb3d3","ZNF33A.2"="#4eb3d3","ZNF468"="#4eb3d3","ZNF611"="#4eb3d3","KAP1"="#4eb3d3")
  } else if (grepl("MER34",Family)){
    Roadmap.plot = Roadmap[Roadmap$Mark.cell %in% c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","H3K9me3","H3K9me3_H1","H3K27me3","H2BK5ac",
                                                    "POU5F1","CHD7","H2AFZ",
                                                    "ZNF211","ZNF211.2"),]
    Roadmap.plot$Pos = factor(Roadmap.plot$Pos,levels=c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","H2BK5ac","H3K9me3","H3K9me3_H1","H3K27me3",
                                                        "POU5F1","CHD7","H2AFZ",
                                                        "ZNF211","ZNF211.2"))
    colMark = c("DNase" = "#b2df8a","H3K27ac" = "#33a02c","H3K4me1" = "#33a02c","H3K4me2" = "#33a02c","H3K4me3" = "#33a02c","H3K9me3" = "#bf812d", "H3K9me3_H1" = "#bf812d","H3K27me3" = "#bf812d","H2BK5ac" = "#33a02c",
                "POU5F1"="0868ac","CHD7" = "#0868ac","H2AFZ" = "#0868ac",
                "ZNF211"="#4eb3d3","ZNF211.2"="#4eb3d3")
  } else if (grepl("MER52",Family)){
    Roadmap.plot = Roadmap[Roadmap$Mark.cell %in% c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","H4K8ac","H2BK5ac","H4K91ac","H3K9me3","H3K9me3_H1","H3K27me3",
                                                    "REST","H2AFZ","RAD21","TBP","MAX","ZNF143",
                                                    "ZNF736","ZNF273","ZNF680","ZNF680.2","ZNF519","ZNF793","ZNF28","ZNF479"),]
    Roadmap.plot$Pos = factor(Roadmap.plot$Pos,levels=c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","H4K8ac","H2BK5ac","H4K91ac","H3K9me3","H3K9me3_H1","H3K27me3",
                                                        "REST","H2AFZ","RAD21","TBP","MAX","ZNF143",
                                                        "ZNF736","ZNF273","ZNF680","ZNF680.2","ZNF519","ZNF793","ZNF28","ZNF479"))
    colMark = c("DNase" = "#b2df8a","H3K27ac" = "#33a02c","H3K4me1" = "#33a02c","H3K4me2" = "#33a02c","H3K4me3" = "#33a02c","H4K8ac" = "#33a02c","H2BK5ac" = "#33a02c","H4K91ac" = "#33a02c","H3K9me3" = "#bf812d", "H3K9me3_H1" = "#bf812d","H3K27me3" = "#bf812d",
                "REST"="0868ac","H2AFZ"="0868ac","RAD21" = "#0868ac","TBP" = "#0868ac","MAX" = "#0868ac","ZNF143" = "#0868ac",
                "ZNF736"="#4eb3d3","ZNF273"="#4eb3d3","ZNF680"="#4eb3d3","ZNF680.2"="#4eb3d3","ZNF519"="#4eb3d3","ZNF793"="#4eb3d3","ZNF28"="#4eb3d3","ZNF479"="#4eb3d3")
  } 
  Roadmap.plot.active = Roadmap.plot[!grepl("ZNF|H3K9me3|H3K27me3|KAP1",Roadmap.plot$Pos) | Roadmap.plot$Pos == "ZNF143",]  ## added ZNF143 because it may be activator 2022/12/28
  Roadmap.plot.active$Pos = factor(Roadmap.plot.active$Pos,levels(Roadmap.plot$Pos)[levels(Roadmap.plot$Pos) %in% Roadmap.plot.active$Pos])

  p_iris2.active = p_iris +  
    new_scale_fill() + 
    geom_fruit(
      data=Annotation_plot_formatted_species,
      geom=geom_tile,
      mapping=aes(y=uniqueID_new,x=Group,fill=Pos),
      pwidth=0.04, # width of the external layer, default is 0.2 times of x range of tree.
      offset=0.03,
      axis.params=list(
        text.size=3,
        axis="none", # add axis text of the layer.
        text.angle=90, # the text angle of x-axis.
        hjust=0.5,  # adjust the horizontal position of text of axis.
        vjust=0
      )
    ) + 
    scale_fill_manual(values=col_cluster_family)
  p_iris3.active = p_iris + 
    new_scale_fill() + 
    geom_fruit(
      data=Annotation_plot_formatted_div,
      geom=geom_tile,
      mapping=aes(y=uniqueID_new,x=Group,fill=divergent_rate),
      pwidth=0.04, # width of the external layer, default is 0.2 times of x range of tree.
      offset=0.03,
      axis.params=list(
        text.size=3,
        axis="none", # add axis text of the layer.
        text.angle=90, # the text angle of x-axis.
        hjust=0.5,  # adjust the horizontal position of text of axis.
        vjust=0
      )
    ) + scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = median(unique(Annotation_plot_formatted_div$divergent_rate),na.rm=T),na.value = 'black')
  
    p_iris4.active = p_iris3.active + 
      new_scale_fill() + 
      geom_fruit(
        data=Roadmap.plot.active,
        geom=geom_tile,
        mapping=aes(y=uniqueID_new,x=Pos,fill=Pos),
        pwidth=0.9, # width of the external layer, default is 0.2 times of x range of tree.
        offset=0.04,
        axis.params=list(
          text.size=3,
          #axis="x", # add axis text of the layer.
          axis="none", # add axis text of the layer.
          text.angle=90, # the text angle of x-axis.
          hjust=0.5,  # adjust the horizontal position of text of axis.
          vjust=0
        )
      ) + scale_fill_manual(values=colMark) 
    
   
    pdf(paste("Figure_2A-",Family,"_circular.active.pdf",sep=""),    # create PNG for the heat map
        width = 5,        # 5 x 300 pixels
        height = 5,
        pointsize = 10)        # smaller font size
    grid.draw(p_iris4.active)
    dev.off()
  
  ## 
  #Annotation_plot = Annotation_plot[order(Annotation_plot$Order),]
  write.table(Annotation_plot,file=paste(Family,"-instance_annotation",".csv",sep=""),sep="\t")
  
  # load FASTA
  fasta.plot = readDNAStringSet(paste("input_trees/",Family,".mafft.prank.opt.gt99.fa",sep=""))
  seq_name = names(fasta.plot)
  sequence = paste(fasta.plot)
  fasta.df <- data.frame(seq_name, sequence)
  fasta.df = fasta.df[fasta.df$seq_name %in% Annotation_plot$uniqueID_new,]
  fasta.df = fasta.df[match(Annotation_plot$uniqueID_new,fasta.df$seq_name),]
  fasta.df$sequence = gsub("B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z","N",fasta.df$sequence)
  fasta.df$seq_name = paste(">",fasta.df$seq_name,sep="")
  # write to a re-ordered fasta file
  write.table(fasta.df, sep="\n",row.names = FALSE, col.names = FALSE, quote = FALSE,file=paste(Family,"-mafft.prank.opt.gt99.reordered.fa",sep=""))
  ## create the sequence MSA
  p_iris_msa = msaplot(p=ggtree(tr,aes(color=group),size=0.2) %<+% Annotation_plot[,c("uniqueID_new","species","TEfamily.x")] + 
                         scale_color_manual(values= col_cluster)+
                         new_scale_color() + 
                         geom_tippoint(aes(shape=species,color=TEfamily.x),size=0.8,show.legend=FALSE)+
                         scale_shape_manual(values = shape_species),
                       fasta=paste(Family,"-mafft.prank.opt.gt99.reordered.fa",sep=""),height=1,offset=0,width=3)
  pdf(paste("Figure2_2A-",Family,"_msa_",".pdf",sep=""),    # create PNG for the heat map
      width = 24,        # 5 x 300 pixels
      height = 12,
      pointsize = 10)        # smaller font size
  grid.draw(p_iris_msa)
  dev.off()
}
