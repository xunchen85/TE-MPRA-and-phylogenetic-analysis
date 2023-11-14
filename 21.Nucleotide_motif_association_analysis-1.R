#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

library(treeio)
library(ggtree)
library(ggmsa)
library(ggplot2)
library(gtable)    
library(ggrepel)
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
library(ComplexHeatmap)
library(pheatmap)
library(corrplot)
library(ggbeeswarm)
library(FSA)
library(ggpubr)
library(ComplexUpset)
library(ggplot2movies)
library(ComplexUpset)

Date = "2023_9_8"

#########################
Family="MER11_combined_F1"
color_species = c("Ancient" = "#6a3d9a","Consensus" = "#6a3d9a","hg19"="#e31a1c","macFas5"="#1f78b4","panTro4"="#33a02c")
shape_species = c("Ancient" = 18,"Consensus" = 18,"hg19"=NA,"macFas5"=NA,"panTro4"=NA)

########################################################
######################### step 1 read annotation and subfamilyinfo previously generated
Summary_Table2 = read.csv("input/Summary_Table2_2023_1_5.subfamilyInfo.csv")
Summary_Table2 =Summary_Table2[!duplicated(Summary_Table2$uniqueID_new),]

# 
Summary_Table2$instanceID.renamed = ifelse(is.na(Summary_Table2$instanceID.renamed),Summary_Table2$uniqueID_new,Summary_Table2$instanceID.renamed)
# 
Summary_Table2$subfamily.name.final = ifelse(Summary_Table2$Group != "Instance",NA,Summary_Table2$subfamily.name.final)

# exclude positive and negative control
Summary_Table2 = Summary_Table2[!(Summary_Table2$Family %in% c("Positive","Negative")),]

# exclude redundant frames
Summary_Table2 = Summary_Table2[Summary_Table2$Is_included == "Included_MPRA",]

# keep high quality ones
Summary_Table2 = Summary_Table2[Summary_Table2$is_kept %in% c("highQuality.iPSC/NPC","highQuality.iPSC","highQuality.NPC"),]

# MPRA activity 
#Summary_Table2 = data.frame(Summary_Table2[,c("Instance_coordinate","uniqueID_new","Family","uniqueID_MPRA","Group","Frame","Species","iPSC.activity.mean","NPC.activity.mean","iPSC.alpha","NPC.alpha","logFC","is_kept")])
Summary_Table2$log2.iPSC.alpha = log2(Summary_Table2$iPSC.alpha)

Summary_Table2$iPSC.alpha.range = ifelse(Summary_Table2$iPSC.alpha >=4 & !is.na(Summary_Table2$log2.iPSC.alpha) & Summary_Table2$iPSC.alpha.Zscore.isActive == "Active","alpha>=4","None")
Summary_Table2$iPSC.alpha.range = ifelse(Summary_Table2$iPSC.alpha >=2 & Summary_Table2$iPSC.alpha < 4 &
                                           !is.na(Summary_Table2$iPSC.alpha) & Summary_Table2$iPSC.alpha.Zscore.isActive == "Active","alpha>=2",Summary_Table2$iPSC.alpha.range)
Summary_Table2$iPSC.alpha.range = ifelse(Summary_Table2$iPSC.alpha <2 &
                                           !is.na(Summary_Table2$iPSC.alpha) & Summary_Table2$iPSC.alpha.Zscore.isActive == "Active","alpha<2",Summary_Table2$iPSC.alpha.range)

Summary_Table2$iPSC.alpha.range = ifelse(is.na(Summary_Table2$iPSC.alpha),"Missing",Summary_Table2$iPSC.alpha.range)

color_alpha.group = c("alpha>=4"="#bd0026","alpha>=2"="#fd8d3c","alpha<2"="#fdbd3c","None"="#f0f0f0","Missing"="#ffffff")

Summary_Table2 = Summary_Table2[order(-Summary_Table2$Order.final),]
Summary_Table2$Family.frame = ifelse(grepl("MER11",Summary_Table2$Family),paste("MER11_combined_",Summary_Table2$Frame,sep=""),NA)
Summary_Table2$Family.frame = ifelse(grepl("MER34",Summary_Table2$Family),paste("MER34_combined_",Summary_Table2$Frame,sep=""),Summary_Table2$Family.frame)
Summary_Table2$Family.frame = ifelse(grepl("MER52",Summary_Table2$Family),paste("MER52_combined_",Summary_Table2$Frame,sep=""),Summary_Table2$Family.frame)

# iPSC alpha.Zscore group
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
iPSC.alpha.Zscore.isActive.group.color = c(">=8" = "#67001f","6-8"="#b2182b","4-6"="#d6604d","2-4"="#f4a582","<2"="#fddbc7","No activity"="#f0f0f0")

# load the TF expression in Taka's data
Motif_exp = read.csv("input/rsem_fpkmTable.csv")
Motif_exp$logFC =log2(as.numeric(as.character(Motif_exp$X72h))+0.01)-log2(as.numeric(as.character(Motif_exp$X0h))+0.01)
Motif_exp$Gene.Symbol = gsub("POU5F1_variant2","POU5F1",Motif_exp$Gene.Symbol)
Motif_exp$Gene.Symbol = gsub("FIGLA_variant1","FIGLA",Motif_exp$Gene.Symbol)
Motif_exp$Gene.Symbol = gsub("KLF14_variant1","KLF14",Motif_exp$Gene.Symbol)
Motif_exp$Gene.Symbol = gsub("ZBTB12_variant1","ZBTB12",Motif_exp$Gene.Symbol)
Motif_exp$Gene.Symbol = gsub("_variant1","",Motif_exp$Gene.Symbol)
Motif_exp[grepl("ZFP",Motif_exp$Gene.Symbol),]
Motif_exp$X0h = as.numeric(as.character(Motif_exp$X0h))
Motif_exp$X72h = as.numeric(as.character(Motif_exp$X72h))

#################### motif cluster and order
# # motif and tree order
load("input/JASPAR.tree.RData")
tree_labels = tree$labels[tree$order]
JASPAR_tree = data.frame("order" = 1:length(tree_labels),"label" = tree_labels)
JASPAR_tree$motifIDs = JASPAR_tree$label
JASPAR_tree$motifIDs = gsub("_MA"," MA",JASPAR_tree$motifIDs)
JASPAR_tree$motifIDs = gsub("_MA"," ",JASPAR_tree$motifIDs)

JASPAR_tree = cSplit(JASPAR_tree,"motifIDs",sep=" ",type.convert = as.character)
JASPAR_tree$motifIDs_uniq = JASPAR_tree$motifIDs_2
JASPAR_tree = cSplit(JASPAR_tree,"motifIDs_uniq",sep="_",type.convert = as.character)
JASPAR_tree$motifIDs_2 = gsub("_",".",JASPAR_tree$motifIDs_2)
JASPAR_tree = JASPAR_tree[order(JASPAR_tree$order),]
# motif cluster
JASPAR_tree.color <- read.csv(file="input/JASPAR.parsed_tree_cluster_1.json.csv")
JASPAR_tree.color <- as.data.frame(JASPAR_tree.color)

# combine
JASPAR_tree.final = data.frame(merge(JASPAR_tree,JASPAR_tree.color,by.x="motifIDs_2",by.y="motifID",all=T))
JASPAR_tree.final$motif.name = toupper(JASPAR_tree.final$motif.name)
JASPAR_tree.final$gene.name = JASPAR_tree.final$motif.name
JASPAR_tree.final$gene.name = gsub("-",":",JASPAR_tree.final$gene.name)      ## update on 2023/2/21
JASPAR_tree.final = data.frame(cSplit(JASPAR_tree.final,"gene.name",sep=":"))
JASPAR_tree.final = merge(JASPAR_tree.final,Motif_exp[,c("Gene.Symbol","Biotype","X0h","X72h")],by.x="gene.name_1",by.y="Gene.Symbol",all.x=T)

# some did not find matched genes
JASPAR_tree.final = JASPAR_tree.final[order(JASPAR_tree.final$order),]     ### order by the tree order
JASPAR_tree.final = JASPAR_tree.final[!duplicated(JASPAR_tree.final$motif.name),]          ### remove redundant motifs
JASPAR_tree.final$X0h.modified = ifelse(is.na(JASPAR_tree.final$X0h),0,JASPAR_tree.final$X0h)

# updated on 2023/9/8
list_MER11.distance = read.csv(file="input/TEwide_group.list_2023_9_4.functional_group_edited.csv")
list_MER11.distance = list_MER11.distance[list_MER11.distance$subfamily2_1 %in% c("MER11A","MER11B","MER11C"),]
list_MER11.distance = list_MER11.distance[!is.na(list_MER11.distance$subfamily.species),c("functional.group","subfamily.species","branch.length.ToKeptNode","subfamily2_1","Order","functional.group")]
colnames(list_MER11.distance) = c("group","subfamily","distance","family","order.tree","FG")
list_MER11.distance$group = gsub("PG","G",list_MER11.distance$group)
list_MER11.distance$FG = gsub("PG","G",list_MER11.distance$FG)
head(list_MER11.distance)

# kept the order based on the combined tree
list_MER11.distance = list_MER11.distance[order(list_MER11.distance$order.tree),]
list_MER11.distance = list_MER11.distance[!is.na(list_MER11.distance$subfamily),]


# subfamily of macaque
# list_MER11.macFas5.distance = read.csv(file="input/macFas5_MER11A_subfamily.mafft.prank.best.fas.nonrev_dna_rerooted.csv")
list_MER11.macFas5.distance = read.csv(file="input_trees_consensus/macFas5_MER11ABC_v1.mafft.prank.best.fas.nonrev_dna.keptroot.csv")
list_MER11.macFas5.distance$FG = NA
# list_MER11.macFas5.distance[1:2,]$FG = "G2"
# list_MER11.macFas5.distance[3:9,]$FG = "G1"
# list_MER11.macFas5.distance[10:11,]$FG = "G2"
# list_MER11.macFas5.distance[12:20,]$FG = "G3"
list_MER11.macFas5.distance[1:7,]$FG = "G1"
list_MER11.macFas5.distance[8:11,]$FG = "G2"
list_MER11.macFas5.distance[12:20,]$FG = "G3"
list_MER11.macFas5.distance[21:nrow(list_MER11.macFas5.distance),]$FG = "G4"

## new order, revised on 2023 Mar. 27
list_MER11 = list_MER11.distance

######################### add the macaque subfamiliy information to the summary table
macaque.tree = read.csv("input_trees/macFas5_rmsk_TE_0bp_MER11A.family.rename_200bp_10_0.02_all_95_group_final.csv")

# load the macaque subfamily information
Summary_Table2.kept = merge(Summary_Table2,macaque.tree[,c("label","label.final.correctedName")],by.x="instanceID.renamed",by.y="label",all.x=T)

Summary_Table2.kept$label.final.correctedName = ifelse(Summary_Table2.kept$Family == "MER11A" &
                                                         Summary_Table2.kept$species == "macFas5",
                                                       Summary_Table2.kept$label.final.correctedName,
                                                       Summary_Table2.kept$subfamily.name.final)

Summary_Table2.kept$label.final.correctedName = ifelse(Summary_Table2.kept$Group != "Instance",NA,Summary_Table2.kept$label.final.correctedName)

################################################## write fasta files for FG4 between human and macaque
Species = "macFas5"
for(Species in c("hg19","macFas5","panTro4")){
  for(Frame in c("F1","F2")){
    if (Species %in% c("hg19","panTro4")){
      Summary_Table2.fasta.tmp = Summary_Table2.kept[Summary_Table2.kept$Species == Species & 
                                                       Summary_Table2.kept$Group == "Instance" &
                                                       Summary_Table2.kept$Frame == Frame &
                                                       Summary_Table2.kept$label.final.correctedName %in% list_MER11.distance[list_MER11.distance$FG == "G4",]$subfamily,]
    } else {
      Summary_Table2.fasta.tmp = Summary_Table2.kept[Summary_Table2.kept$Species == Species & 
                                                       Summary_Table2.kept$Group == "Instance" &
                                                       Summary_Table2.kept$Frame == Frame &
                                                       Summary_Table2.kept$label.final.correctedName %in% list_MER11.macFas5.distance[list_MER11.macFas5.distance$FG == "G4",]$subfamily,]
    }
    Summary_Table2.fasta.tmp$Sequence.name = paste(">",Summary_Table2.fasta.tmp$instanceID.renamed," ",Summary_Table2.fasta.tmp$iPSC.alpha.Zscore.isActive,sep="")
    Summary_Table2.fasta.tmp = Summary_Table2.fasta.tmp[,c("Sequence.name","Sequence")]
    write.table(Summary_Table2.fasta.tmp, sep="\n",row.names = FALSE, col.names = FALSE, quote = FALSE,file=paste("MER11_G4_",Species,"_",Frame,"_",Date,".fa",sep=""))
  }
}

##################################################
######################### step 2: prepare the input files for the motif and nucleotide association analysis 
Summary_Table2.fimoPlot = Summary_Table2.kept[!is.na(Summary_Table2.kept$iPSC.alpha.range),]   ## exclude positive and negative

# load the fimo results
Fimo.all = read.delim("input/MER11_34_52_frame_fimo.tsv_2023_1_7.gz",header=F,sep="")
# JASPAR
Fimo.all = Fimo.all[Fimo.all$V2 != "KZNF",]
colnames(Fimo.all) = c("Family.frame","DB","motif_id","motif_name","uniqueID.new","start","stop","strand","score","p_value","q_value","matched_sequence")
Fimo.all = Fimo.all[,!(colnames(Fimo.all) == "motif_id")]

Fimo.all$uniqueID.new.DB.motif_name = paste(Fimo.all$Family.frame,Fimo.all$uniqueID.new,Fimo.all$DB,Fimo.all$motif_name)
# JASPAR for now
Fimo.all.sum = data.frame(Fimo.all[Fimo.all$DB == "JASPAR2022",] %>% group_by(Family.frame,DB,uniqueID.new,motif_name) %>% dplyr::count())
#
Family = "MER11_combined_F2"
Families = unique(Fimo.all$Family.frame)

for(Family in Families[grepl("_combined_F",Families)]) {
  ### 2.1 
  # load the annotation
  Annotation_plot.group = Summary_Table2.fimoPlot[Summary_Table2.fimoPlot$Family.frame == Family & 
                                                    !is.na(Summary_Table2.fimoPlot$Family.frame) & 
                                                    !(Summary_Table2.fimoPlot$Family %in% "MER11D"),]    ## there is no MER11D frame was designed
  
  
  Annotation_plot.group$subfamily.name.final = ifelse(is.na(Annotation_plot.group$subfamily.name.final),paste(gsub("MER","",Annotation_plot.group$Family),"_N",sep=""),Annotation_plot.group$subfamily.name.final)
  Annotation_plot.group$label.final.correctedName = ifelse(is.na(Annotation_plot.group$label.final.correctedName),paste(gsub("MER","",Annotation_plot.group$Family),"_N",sep=""),Annotation_plot.group$label.final.correctedName)
  # exclude consensus sequences
  Annotation_plot.group = Annotation_plot.group[Annotation_plot.group$Group == "Instance",]
  #write.csv(Annotation_plot.group,file="../test.csv")
  
  #Annotation_plot.group = Summary_Table2.fimoPlot
  # load the fasta file
  fasta.plot = readDNAStringSet(paste("input_trees/",Family,".mafft.opt.gt99.reordered.fa",sep=""))
  seq_name = names(fasta.plot)
  sequence = paste(fasta.plot)
  fasta.df <- data.frame(seq_name, sequence)
  # add the grouping info
  fasta.df = merge(fasta.df,Annotation_plot.group,by.x="seq_name",by.y="uniqueID_new",all.y=T)  ## only keep the group information
  # split by each nucleotide
  fasta.df.rate = data.frame(Reduce(rbind, strsplit(fasta.df$sequence, "")))
  rownames(fasta.df.rate) = fasta.df$seq_name
  fasta.df.rate.t = data.frame(t(fasta.df.rate))
  
  ### 2.2 summary
  fasta.df.rate.t.sum.final = data.frame("Group1" = "First","Group" = "First","nt.all"=0,"nt.A"=0,"nt.T"=0,"nt.C"=0,"nt.G"=0,"nt.M"=0,"nt.N"=0,
                                         "nt.all.r"=0,"nt.A.r"=0,"nt.T.r"=0,"nt.C.r"=0,"nt.G.r"=0,"nt.M.r"=0,"nt.N.r"=0,"nt.pos"=0)
  if (grepl("MER11",Family)){
    Groupings = data.frame("Group2"=c("All","G1","G2","G3","G4"))
    Groupings$Group1 = ifelse(Groupings$Group2 == "All","All",NA)
    Groupings$Group1 = ifelse(grepl("^G",Groupings$Group2),"group",Groupings$Group1)
  } else {
    Groupings = data.frame("Group2"=c("All"))
    Groupings$Group1 = ifelse(Groupings$Group2 == "All","All",NA)
  }
  
  # unique families and subfamilies based on the hg19
  Annotation_plot.group.unique = Annotation_plot.group[Annotation_plot.group$species == "hg19",]
  Annotation_plot.group.unique = Annotation_plot.group.unique[!duplicated(paste(Annotation_plot.group.unique$TEfamily.y,Annotation_plot.group.unique$subfamily.name.final)),c("TEfamily.y","subfamily.name.final")]
  Annotation_plot.group.unique = Annotation_plot.group.unique[!is.na(Annotation_plot.group.unique$TEfamily.y),]
  
  Groupings = merge(Groupings,Annotation_plot.group.unique,by.x="Group2",by.y="subfamily.name.final",all.x=T)
  Groupings$TEfamily.y = ifelse(is.na(Groupings$TEfamily.y) & Groupings$Group1 == "group","MER11",Groupings$TEfamily.y)
  colnames(Groupings)[3] = "TEfamily"
  
  # not analyze it based on subfamilies
  Groupings = Groupings[Groupings$Group1 %in% c("All","group"),]
  
  ############################## 2.2.1 prepare nucleotide ped files; each nucleotide summary
  RowID = 5
  for (RowID in 1:nrow(Groupings)){
    kept.sample1 = ""
    if (Groupings[RowID,]$Group1 == "All"){
      kept.sample1 = colnames(fasta.df.rate.t)
    } else if (Groupings[RowID,]$Group1 == "group"){
      kept.sample1 = Annotation_plot.group[Annotation_plot.group$label.final.correctedName %in% c(list_MER11[list_MER11$FG == Groupings[RowID,]$Group2,]$subfamily,
                                                                                                  list_MER11.macFas5.distance[list_MER11.macFas5.distance$FG == Groupings[RowID,]$Group2,]$subfamily),]$uniqueID_new
    }
    
    fasta.df.rate.t.tmp = data.frame(fasta.df.rate.t[,kept.sample1])
    if (ncol(fasta.df.rate.t.tmp) < 2) {
      next
    } else {
      fasta.df.rate.t.tmp.2 = fasta.df.rate.t.tmp[,!grepl("^C_",colnames(fasta.df.rate.t.tmp))]
    }
    #
    if (ncol(fasta.df.rate.t.tmp) == 0) {
      fasta.df.rate.t.sum = data.frame("Group1" = Groupings[RowID,]$Group1,"Group" = Groupings[RowID,]$Group2,"nt.all"=0,"nt.A"=0,"nt.T"=0,"nt.C"=0,"nt.G"=0,"nt.M"=0,"nt.N"=0,
                                       "nt.all.r"=0,"nt.A.r"=0,"nt.T.r"=0,"nt.C.r"=0,"nt.G.r"=0,"nt.M.r"=0,"nt.N.r"=0,"nt.pos"=0)
    } else if (ncol(fasta.df.rate.t.tmp) > 0 & ncol(fasta.df.rate.t.tmp.2) == 0){
      fasta.df.rate.t.sum = fasta.df.rate.t.tmp[,c(1,1)]
      fasta.df.rate.t.sum$Group1 = Groupings[RowID,]$Group1
      fasta.df.rate.t.sum$Group = Groupings[RowID,]$Group2
      fasta.df.rate.t.tmp = data.frame()
      fasta.df.rate.t.sum$nt.all = ncol(fasta.df.rate.t.tmp)
      fasta.df.rate.t.sum$nt.A = rowSums(fasta.df.rate.t.tmp == "A",na.rm=T)
      fasta.df.rate.t.sum$nt.T = rowSums(fasta.df.rate.t.tmp == "T",na.rm=T)
      fasta.df.rate.t.sum$nt.C = rowSums(fasta.df.rate.t.tmp == "C",na.rm=T)
      fasta.df.rate.t.sum$nt.G = rowSums(fasta.df.rate.t.tmp == "G",na.rm=T)
      fasta.df.rate.t.sum$nt.M = rowSums(fasta.df.rate.t.tmp == "-",na.rm=T)
      fasta.df.rate.t.sum$nt.N = rowSums(fasta.df.rate.t.tmp == "N",na.rm=T)
      fasta.df.rate.t.sum$nt.all.r = ncol(fasta.df.rate.t.tmp.2)
      fasta.df.rate.t.sum$nt.A.r = 0
      fasta.df.rate.t.sum$nt.T.r = 0
      fasta.df.rate.t.sum$nt.C.r = 0
      fasta.df.rate.t.sum$nt.G.r = 0
      fasta.df.rate.t.sum$nt.M.r = 0
      fasta.df.rate.t.sum$nt.N.r = 0
      fasta.df.rate.t.sum$nt.pos = rownames(fasta.df.rate.t.tmp)
      fasta.df.rate.t.sum = fasta.df.rate.t.sum[,c(-1,-2)]
    } else if (ncol(fasta.df.rate.t.tmp) > 0 & ncol(fasta.df.rate.t.tmp.2) > 0){
      fasta.df.rate.t.sum = fasta.df.rate.t.tmp[,c(1,1)]
      fasta.df.rate.t.sum$Group1 = Groupings[RowID,]$Group1
      fasta.df.rate.t.sum$Group = Groupings[RowID,]$Group2
      fasta.df.rate.t.sum$nt.all = ncol(fasta.df.rate.t.tmp)
      fasta.df.rate.t.sum$nt.A = rowSums(fasta.df.rate.t.tmp == "A",na.rm=T)
      fasta.df.rate.t.sum$nt.T = rowSums(fasta.df.rate.t.tmp == "T",na.rm=T)
      fasta.df.rate.t.sum$nt.C = rowSums(fasta.df.rate.t.tmp == "C",na.rm=T)
      fasta.df.rate.t.sum$nt.G = rowSums(fasta.df.rate.t.tmp == "G",na.rm=T)
      fasta.df.rate.t.sum$nt.M = rowSums(fasta.df.rate.t.tmp == "-",na.rm=T)
      fasta.df.rate.t.sum$nt.N = rowSums(fasta.df.rate.t.tmp == "N",na.rm=T)
      fasta.df.rate.t.sum$nt.all.r = ncol(fasta.df.rate.t.tmp.2)
      fasta.df.rate.t.sum$nt.A.r = rowSums(fasta.df.rate.t.tmp.2 == "A",na.rm=T)
      fasta.df.rate.t.sum$nt.T.r = rowSums(fasta.df.rate.t.tmp.2 == "T",na.rm=T)
      fasta.df.rate.t.sum$nt.C.r = rowSums(fasta.df.rate.t.tmp.2 == "C",na.rm=T)
      fasta.df.rate.t.sum$nt.G.r = rowSums(fasta.df.rate.t.tmp.2 == "G",na.rm=T)
      fasta.df.rate.t.sum$nt.M.r = rowSums(fasta.df.rate.t.tmp.2 == "-",na.rm=T)
      fasta.df.rate.t.sum$nt.N.r = rowSums(fasta.df.rate.t.tmp.2 == "N",na.rm=T)
      fasta.df.rate.t.sum$nt.pos = rownames(fasta.df.rate.t.tmp)
      fasta.df.rate.t.sum = fasta.df.rate.t.sum[,c(-1,-2)]
    }
    fasta.df.rate.t.sum.final = rbind(fasta.df.rate.t.sum.final,fasta.df.rate.t.sum)
  }
  fasta.df.rate.t.sum.final = fasta.df.rate.t.sum.final[-1,]
  fasta.df.rate.t.sum.final = fasta.df.rate.t.sum.final[!is.na(fasta.df.rate.t.sum.final$Group),]
  
  # 2.2.2 each nucleotide majority
  fasta.df.rate.t.sum.final$majority = ifelse(fasta.df.rate.t.sum.final$nt.A >= fasta.df.rate.t.sum.final$nt.C &
                                                fasta.df.rate.t.sum.final$nt.A >= fasta.df.rate.t.sum.final$nt.T &
                                                fasta.df.rate.t.sum.final$nt.A >= fasta.df.rate.t.sum.final$nt.G,"A",NA)
  fasta.df.rate.t.sum.final$majority = ifelse(fasta.df.rate.t.sum.final$nt.T >= fasta.df.rate.t.sum.final$nt.C &
                                                fasta.df.rate.t.sum.final$nt.T >= fasta.df.rate.t.sum.final$nt.A &
                                                fasta.df.rate.t.sum.final$nt.T >= fasta.df.rate.t.sum.final$nt.G,"T",fasta.df.rate.t.sum.final$majority)
  fasta.df.rate.t.sum.final$majority = ifelse(fasta.df.rate.t.sum.final$nt.C >= fasta.df.rate.t.sum.final$nt.A &
                                                fasta.df.rate.t.sum.final$nt.C >= fasta.df.rate.t.sum.final$nt.T &
                                                fasta.df.rate.t.sum.final$nt.C >= fasta.df.rate.t.sum.final$nt.G,"C",fasta.df.rate.t.sum.final$majority)
  fasta.df.rate.t.sum.final$majority = ifelse(fasta.df.rate.t.sum.final$nt.G >= fasta.df.rate.t.sum.final$nt.C &
                                                fasta.df.rate.t.sum.final$nt.G >= fasta.df.rate.t.sum.final$nt.T &
                                                fasta.df.rate.t.sum.final$nt.G >= fasta.df.rate.t.sum.final$nt.A,"G",fasta.df.rate.t.sum.final$majority)
  fasta.df.rate.t.sum.final$majority.rate = ifelse(fasta.df.rate.t.sum.final$majority == "A",fasta.df.rate.t.sum.final$nt.A,NA)
  fasta.df.rate.t.sum.final$majority.rate = ifelse(fasta.df.rate.t.sum.final$majority == "T",fasta.df.rate.t.sum.final$nt.T,fasta.df.rate.t.sum.final$majority.rate)
  fasta.df.rate.t.sum.final$majority.rate = ifelse(fasta.df.rate.t.sum.final$majority == "C",fasta.df.rate.t.sum.final$nt.C,fasta.df.rate.t.sum.final$majority.rate)
  fasta.df.rate.t.sum.final$majority.rate = ifelse(fasta.df.rate.t.sum.final$majority == "G",fasta.df.rate.t.sum.final$nt.G,fasta.df.rate.t.sum.final$majority.rate)
  
  ### 2.2.3 convert it to .ped format
  fasta.df.rate.ped.tmp = data.frame(t(fasta.df.rate.t))
  fasta.df.rate.ped.final = data.frame(fasta.df.rate.ped.tmp[,c(1,1)])
  colnames(fasta.df.rate.ped.final) = c("tmp1","tmp1.2")
  fasta.df.rate.ped.final$sampleID = rownames(fasta.df.rate.ped.final)
  fasta.df.rate.ped.final$Order = 1:nrow(fasta.df.rate.ped.final)
  
  # achieve iPSC.alpha after the MAD normalization; iPSC.alpha.Zscore was used for the association analysis
  fasta.df.rate.ped.final = merge(fasta.df.rate.ped.final,Annotation_plot.group[,c("uniqueID_new","Family","iPSC.alpha.Zscore","logFC")],by.x="sampleID",by.y="uniqueID_new",all.x=T)
  colnames(fasta.df.rate.ped.final)[which(colnames(fasta.df.rate.ped.final)=="iPSC.alpha.Zscore")] = "iPSC.alpha"
  fasta.df.rate.ped.final = fasta.df.rate.ped.final[,c("Family","sampleID","Order","Order","Order","iPSC.alpha","logFC")]
  fasta.df.rate.ped.final = fasta.df.rate.ped.final[order(fasta.df.rate.ped.final$Order),]
  fasta.df.rate.ped.final$iPSC.alpha = ifelse(is.na(fasta.df.rate.ped.final$iPSC.alpha),999,fasta.df.rate.ped.final$iPSC.alpha)
  fasta.df.rate.ped.final$logFC = ifelse(is.na(fasta.df.rate.ped.final$logFC),999,fasta.df.rate.ped.final$logFC)
  fasta.df.rate.ped.final[,c(3,4)] = 0
  fasta.df.rate.ped.final[,5] = "other"
  fasta.df.rate.ped.final.kept = fasta.df.rate.ped.final
  
  # info table 1: compared to all families
  fasta.df.rate.ped.final = fasta.df.rate.ped.final.kept
  fasta.df.rate.ped.final.inDel = fasta.df.rate.ped.final
  
  fasta.df.rate.t.sum.ped = fasta.df.rate.t.sum.final[fasta.df.rate.t.sum.final$Group1 == "All",]
  rowID = 1
  fasta.df.rate.ped.final_BB = fasta.df.rate.ped.final
  for (rowID in 1:nrow(fasta.df.rate.t.sum.ped)){
    # first allele
    fasta.df.rate.ped.final[,(ncol(fasta.df.rate.ped.final)+1)] = ifelse(fasta.df.rate.ped.tmp[,rowID] == fasta.df.rate.t.sum.ped[rowID,]$majority &
                                                                           fasta.df.rate.ped.tmp[,rowID] != "-","A","0")
    fasta.df.rate.ped.final[,(ncol(fasta.df.rate.ped.final))] = ifelse(fasta.df.rate.ped.tmp[,rowID] != fasta.df.rate.t.sum.ped[rowID,]$majority &
                                                                         fasta.df.rate.ped.tmp[,rowID] != "-","B",fasta.df.rate.ped.final[,(ncol(fasta.df.rate.ped.final))])
    fasta.df.rate.ped.final.inDel[,(ncol(fasta.df.rate.ped.final.inDel)+1)] = ifelse(fasta.df.rate.ped.tmp[,rowID] != "-","A","B")
    colnames(fasta.df.rate.ped.final)[ncol(fasta.df.rate.ped.final)] = fasta.df.rate.t.sum.ped[rowID,]$nt.pos
    colnames(fasta.df.rate.ped.final.inDel)[ncol(fasta.df.rate.ped.final.inDel)] = fasta.df.rate.t.sum.ped[rowID,]$nt.pos
    # second allele
    fasta.df.rate.ped.final[,(ncol(fasta.df.rate.ped.final)+1)] = ifelse(fasta.df.rate.ped.tmp[,rowID] == fasta.df.rate.t.sum.ped[rowID,]$majority &
                                                                           fasta.df.rate.ped.tmp[,rowID] != "-","A","0")
    fasta.df.rate.ped.final[,(ncol(fasta.df.rate.ped.final))] = ifelse(fasta.df.rate.ped.tmp[,rowID] != fasta.df.rate.t.sum.ped[rowID,]$majority &
                                                                         fasta.df.rate.ped.tmp[,rowID] != "-","B",fasta.df.rate.ped.final[,(ncol(fasta.df.rate.ped.final))])
    fasta.df.rate.ped.final.inDel[,(ncol(fasta.df.rate.ped.final.inDel)+1)] = ifelse(fasta.df.rate.ped.tmp[,rowID] != "-","A","B")
    colnames(fasta.df.rate.ped.final)[ncol(fasta.df.rate.ped.final)] = paste(fasta.df.rate.t.sum.ped[rowID,]$nt.pos,"_2",sep="")
    colnames(fasta.df.rate.ped.final.inDel)[ncol(fasta.df.rate.ped.final.inDel)] = paste(fasta.df.rate.t.sum.ped[rowID,]$nt.pos,"_2",sep="")
  }
  
  # info table 2: compared to each family; I grouped all frames from a family (originally annotated) together;
  fasta.df.rate.t.sum.ped.families = fasta.df.rate.t.sum.final[fasta.df.rate.t.sum.final$Group1 == "family",]
  rowID = 1
  fasta.df.rate.ped.final.each = data.frame("tmp1"=NA,"tmp2"=NA)
  fasta.df.rate.ped.final.each.inDel = fasta.df.rate.ped.final.each
  for (eachfamily in unique(fasta.df.rate.t.sum.ped.families$Group)){
    print(eachfamily)
    fasta.df.rate.ped.final.each.sub = fasta.df.rate.ped.final.kept[fasta.df.rate.ped.final.kept$sampleID %in% Annotation_plot.group[Annotation_plot.group$Family == eachfamily & !is.na(Annotation_plot.group$Family),]$uniqueID_new,]
    fasta.df.rate.ped.final.each.sub.inDel = fasta.df.rate.ped.final.kept[fasta.df.rate.ped.final.kept$sampleID %in% Annotation_plot.group[Annotation_plot.group$Family == eachfamily,]$uniqueID_new,]
    print(nrow(fasta.df.rate.ped.final.each.sub))
    fasta.df.rate.t.sum.ped.each = fasta.df.rate.t.sum.ped.families[fasta.df.rate.t.sum.ped.families$Group == eachfamily,]
    fasta.df.rate.ped.tmp.each = fasta.df.rate.ped.tmp[rownames(fasta.df.rate.ped.tmp) %in% Annotation_plot.group[Annotation_plot.group$Family == eachfamily,]$uniqueID_new,]
    for (rowID in 1:nrow(fasta.df.rate.t.sum.ped)){
      # first allele
      fasta.df.rate.ped.final.each.sub[,(ncol(fasta.df.rate.ped.final.each.sub)+1)] = ifelse(fasta.df.rate.ped.tmp.each[,rowID] == fasta.df.rate.t.sum.ped.each[fasta.df.rate.t.sum.ped.each$nt.pos == fasta.df.rate.t.sum.ped[rowID,]$nt.pos,]$majority &
                                                                                               fasta.df.rate.ped.tmp.each[,rowID] != "-","A","0")
      fasta.df.rate.ped.final.each.sub[,(ncol(fasta.df.rate.ped.final.each.sub))] = ifelse(fasta.df.rate.ped.tmp.each[,rowID] != fasta.df.rate.t.sum.ped.each[fasta.df.rate.t.sum.ped.each$nt.pos == fasta.df.rate.t.sum.ped[rowID,]$nt.pos,]$majority &
                                                                                             fasta.df.rate.ped.tmp.each[,rowID] != "-","B",fasta.df.rate.ped.final.each.sub[,(ncol(fasta.df.rate.ped.final.each.sub))])
      fasta.df.rate.ped.final.each.sub.inDel[,(ncol(fasta.df.rate.ped.final.each.sub.inDel)+1)] = ifelse(fasta.df.rate.ped.tmp.each[,rowID] != "-","A","B")
      colnames(fasta.df.rate.ped.final.each.sub)[ncol(fasta.df.rate.ped.final.each.sub)] = fasta.df.rate.t.sum.ped[rowID,]$nt.pos
      colnames(fasta.df.rate.ped.final.each.sub.inDel)[ncol(fasta.df.rate.ped.final.each.sub.inDel)] = fasta.df.rate.t.sum.ped[rowID,]$nt.pos
      # second allele
      fasta.df.rate.ped.final.each.sub[,(ncol(fasta.df.rate.ped.final.each.sub)+1)] = ifelse(fasta.df.rate.ped.tmp.each[,rowID] == fasta.df.rate.t.sum.ped.each[fasta.df.rate.t.sum.ped.each$nt.pos == fasta.df.rate.t.sum.ped[rowID,]$nt.pos,]$majority &
                                                                                               fasta.df.rate.ped.tmp.each[,rowID] != "-","A","0")
      fasta.df.rate.ped.final.each.sub[,(ncol(fasta.df.rate.ped.final.each.sub))] = ifelse(fasta.df.rate.ped.tmp.each[,rowID] != fasta.df.rate.t.sum.ped.each[fasta.df.rate.t.sum.ped.each$nt.pos == fasta.df.rate.t.sum.ped[rowID,]$nt.pos,]$majority &
                                                                                             fasta.df.rate.ped.tmp.each[,rowID] != "-","B",fasta.df.rate.ped.final.each.sub[,(ncol(fasta.df.rate.ped.final.each.sub))])
      fasta.df.rate.ped.final.each.sub.inDel[,(ncol(fasta.df.rate.ped.final.each.sub.inDel)+1)] = ifelse(fasta.df.rate.ped.tmp.each[,rowID] != "-","A","B")
      colnames(fasta.df.rate.ped.final.each.sub)[ncol(fasta.df.rate.ped.final.each.sub)] = paste(fasta.df.rate.t.sum.ped[rowID,]$nt.pos,"_2",sep="")
      colnames(fasta.df.rate.ped.final.each.sub.inDel)[ncol(fasta.df.rate.ped.final.each.sub.inDel)] = paste(fasta.df.rate.t.sum.ped[rowID,]$nt.pos,"_2",sep="")
    }
    if (nrow(fasta.df.rate.ped.final.each) == 1){
      fasta.df.rate.ped.final.each = fasta.df.rate.ped.final.each.sub
      fasta.df.rate.ped.final.each.inDel = fasta.df.rate.ped.final.each.sub.inDel
    } else {
      fasta.df.rate.ped.final.each = rbind(fasta.df.rate.ped.final.each,fasta.df.rate.ped.final.each.sub)
      fasta.df.rate.ped.final.each.inDel = rbind(fasta.df.rate.ped.final.each.inDel,fasta.df.rate.ped.final.each.sub.inDel)
    }
  }
  
  ############################## 2.2.4 prepare motif ped files
  Fimo.all.sum.plot = Fimo.all.sum[Fimo.all.sum$Family.frame == Family,]
  list.frame = unique(Annotation_plot.group$uniqueID_new)
  list.motif = unique(Fimo.all.sum.plot$motif_name)
  motif.df.matrix.final = matrix(nrow = length(list.frame), ncol = length(list.motif)*2) 
  motif.df.matrix.final[is.na(motif.df.matrix.final)] <- "A"
  motif.df.matrix.final = data.frame(motif.df.matrix.final) 
  rownames(motif.df.matrix.final) = list.frame
  colnames(motif.df.matrix.final) = paste(rep(list.motif,each=2),"_",c(1,2),sep="")
  
  #
  motif.df.matrix.final_1 = motif.df.matrix.final[,c(1,2)]
  motif.df.matrix.final_1$sampleID = rownames(motif.df.matrix.final_1)
  motif.df.matrix.final_1$Order = 1:nrow(motif.df.matrix.final_1)
  motif.df.matrix.final_1$Order2 = 0
  motif.df.matrix.final_1$Order3 = "other"
  motif.df.matrix.final_1 = merge(motif.df.matrix.final_1,Annotation_plot.group[,c("uniqueID_new","Family","iPSC.alpha.Zscore","logFC")],by.x="sampleID",by.y="uniqueID_new",all.x=T)
  motif.df.matrix.final_1$iPSC.alpha = motif.df.matrix.final_1$iPSC.alpha.Zscore
  motif.df.matrix.final_1 = motif.df.matrix.final_1[,c("Family","sampleID","Order","Order2","Order3","iPSC.alpha","logFC")]
  motif.df.matrix.final_1 = motif.df.matrix.final_1[order(motif.df.matrix.final_1$Order),]
  motif.df.matrix.final_1$iPSC.alpha = ifelse(is.na(motif.df.matrix.final_1$iPSC.alpha),999,motif.df.matrix.final_1$iPSC.alpha)
  motif.df.matrix.final_1$logFC = ifelse(is.na(motif.df.matrix.final_1$logFC),999,motif.df.matrix.final_1$logFC)
  motif.df.matrix.final_1[,c(3,4)] = 0
  motif.df.matrix.final_1[,5] = "other"
  
  # assign by motif
  motif.df.matrix.final_BB = motif.df.matrix.final
  for (rowID in 1:nrow(Fimo.all.sum.plot)){
    motif.df.matrix.final[which(rownames(motif.df.matrix.final) == Fimo.all.sum.plot[rowID,]$uniqueID.new),
                          which(colnames(motif.df.matrix.final) == paste(Fimo.all.sum.plot[rowID,]$motif_name,"_2",sep=""))] = "B"
    motif.df.matrix.final_BB[which(rownames(motif.df.matrix.final_BB) == Fimo.all.sum.plot[rowID,]$uniqueID.new),
                             which(colnames(motif.df.matrix.final_BB) == paste(Fimo.all.sum.plot[rowID,]$motif_name,"_1",sep=""))] = "B"
    motif.df.matrix.final_BB[which(rownames(motif.df.matrix.final_BB) == Fimo.all.sum.plot[rowID,]$uniqueID.new),
                             which(colnames(motif.df.matrix.final_BB) == paste(Fimo.all.sum.plot[rowID,]$motif_name,"_2",sep=""))] = "B"
    if (Fimo.all.sum.plot[rowID,]$n >1){
      motif.df.matrix.final[which(rownames(motif.df.matrix.final) == Fimo.all.sum.plot[rowID,]$uniqueID.new),
                            which(colnames(motif.df.matrix.final) == paste(Fimo.all.sum.plot[rowID,]$motif_name,"_1",sep=""))] = "B"
    }
  }
  motif.df.matrix.final = cbind(motif.df.matrix.final_1,motif.df.matrix.final)
  motif.df.matrix.final_BB = cbind(motif.df.matrix.final_1,motif.df.matrix.final_BB)
  
  # convert to .map format for variant
  fasta.df.rate.t.sum.map = fasta.df.rate.t.sum.ped[,c(1,2,3,4)]
  fasta.df.rate.t.sum.map[,1] = 1
  #fasta.df.rate.t.sum.map[,2] = rownames(fasta.df.rate.t.sum.map)
  fasta.df.rate.t.sum.map[,2] = fasta.df.rate.t.sum.ped$nt.pos
  fasta.df.rate.t.sum.map[,3] = 0
  fasta.df.rate.t.sum.map[,4] = gsub("X","",fasta.df.rate.t.sum.map[,2])
  
  # convert to .map format for motif
  motif.df.matrix.final.map = data.frame("Group"=0,Motif=list.motif,"Other"=0,"Order"=1:length(list.motif))
  motif.df.matrix.final.map[,1] = 1
  motif.df.matrix.final.map[,3] = 0
  
  ### 2.2.5 write to the map and ped files
  RowID = 79
  Species = "hg19"
  RowID = 1
  for (Species in c("All","hg19","macFas5","panTro4","hg19_panTro4")){
    for (RowID in 1:nrow(Groupings)){
      ## candidate instances 
      kept.sample2 = ""
      if (Groupings[RowID,]$Group1 == "All" & Species != "All"){
        next
      } 
      if (Groupings[RowID,]$Group1 == "All"){    ## group across all families
        kept.sample2 = fasta.df.rate.ped.final$sampleID
      } else if (Groupings[RowID,]$Group1 == "group"){
        kept.sample2 = Annotation_plot.group[Annotation_plot.group$label.final.correctedName %in% c(list_MER11[list_MER11$FG == Groupings[RowID,]$Group2,]$subfamily,
                                                                                                    list_MER11.macFas5.distance[list_MER11.macFas5.distance$FG == Groupings[RowID,]$Group2,]$subfamily),]$uniqueID_new
      }
      
      # filter by species added on 2023_4_10
      if (Species == "All"){
        kept.sample = unique(Annotation_plot.group[Annotation_plot.group$uniqueID_new %in% kept.sample2,]$uniqueID_new)
      } else if (Species == "hg19_panTro4") {
        kept.sample = unique(Annotation_plot.group[Annotation_plot.group$uniqueID_new %in% kept.sample2 & Annotation_plot.group$Species %in% c("hg19","panTro4"),]$uniqueID_new)
      } else {
        kept.sample = unique(Annotation_plot.group[Annotation_plot.group$uniqueID_new %in% kept.sample2 & Annotation_plot.group$Species == Species,]$uniqueID_new)
      }
      ## filtering
      fasta.df.rate.ped.final.sub = fasta.df.rate.ped.final[fasta.df.rate.ped.final$sampleID %in% kept.sample,]
      fasta.df.rate.ped.final.inDel.sub = fasta.df.rate.ped.final.inDel[fasta.df.rate.ped.final.inDel$sampleID %in% kept.sample,]
      fasta.df.rate.ped.final.each.sub = fasta.df.rate.ped.final.each[fasta.df.rate.ped.final.each$sampleID %in% kept.sample,]
      fasta.df.rate.ped.final.each.inDel.sub = fasta.df.rate.ped.final.each.inDel[fasta.df.rate.ped.final.each.inDel$sampleID %in% kept.sample,]
      motif.df.matrix.final.sub = motif.df.matrix.final[motif.df.matrix.final$sampleID %in% kept.sample,]
      motif.df.matrix.final_BB.sub = motif.df.matrix.final_BB[motif.df.matrix.final_BB$sampleID %in% kept.sample,]
      
      # map files
      write.table(motif.df.matrix.final.map,file=paste("Step9_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".motif.map",sep=""),
                  sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      write.table(fasta.df.rate.t.sum.map,file=paste("Step9_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variant.map",sep=""),
                  sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      write.csv(Annotation_plot.group[Annotation_plot.group$uniqueID_new %in% kept.sample,c("uniqueID_new","Species","Family","Group","subfamily.name.final","label.final.correctedName")],file=paste("Step9_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".map.csv",sep=""),
                sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      # relative to all families
      write.table(fasta.df.rate.ped.final.sub[,-7],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variantAll.ped",sep=""),
                  sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      write.table(fasta.df.rate.ped.final.inDel.sub[,-7],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variantAll.indel.ped",sep=""),
                  sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      # write.table(fasta.df.rate.ped.final.sub[,-6],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variantAll.logFC.ped",sep=""),
      #             sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      # write.table(fasta.df.rate.ped.final.inDel.sub[,-6],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variantAll.indel.logFC.ped",sep=""),
      #             sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      # # relative to each family
      # write.table(fasta.df.rate.ped.final.each.sub[,-7],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variantEach.ped",sep=""),
      #             sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      # write.table(fasta.df.rate.ped.final.each.inDel.sub[,-7],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variantEach.indel.ped",sep=""),
      #             sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      # write.table(fasta.df.rate.ped.final.each.sub[,-6],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variantEach.logFC.ped",sep=""),
      #             sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      # write.table(fasta.df.rate.ped.final.each.inDel.sub[,-6],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".variantEach.indel.logFC.ped",sep=""),
      #             sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      # motif
      write.table(motif.df.matrix.final.sub[,-7],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".motif.ped",sep=""),
                  sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
      # write.table(motif.df.matrix.final.sub[,-6],file=paste("Association_",Species,"_",Family,"_",Groupings[RowID,]$Group2,"_",Date,".motif.logFC.ped",sep=""),
      #             sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  ### 2.2.6 output files
  write.csv(Annotation_plot.group,file=paste("Association_",Family,"_",Date,"_all.info.csv",sep=""))  
  write.csv(fasta.df.rate.t,file=paste("Association_",Family,"_",Date,"_all.variant.csv",sep=""))  
  write.csv(fasta.df.rate.t.sum.final,file=paste("Association_",Family,"_",Date,"_all.variant.sum.csv",sep="")) 
  write.csv(Fimo.all.sum.plot,file=paste("Association_",Family,"_",Date,"_all.motif.sum.csv",sep="")) 
}

##################################################
######################## performed the motif association analysis using plink
