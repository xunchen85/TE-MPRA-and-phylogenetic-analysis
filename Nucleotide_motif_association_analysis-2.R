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
#install.packages("ggplot2movies")
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
Families = c("MER11_combined_F1","MER11_combined_F2","MER34_combined_F1","MER52_combined_F1","MER52_combined_F2")

##################################################
########################we then performed the motif association analysis using plink with the default parameters

##################################################
######################## step 3: motif association plot; we then plot the p value
Kept_motifs = data.frame("Family"=NA,"Group"=NA,"Motif"=NA,"p-value"=NA)

linear.regression.output.kept = read.delim("input/MER11_34_52.glm.linear_2023_7_10",header=F,sep="")

colnames(linear.regression.output.kept) = c("Species","Family","Combination","Type","CHROM","POS","ID","REF","ALT","A1","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")  
linear.regression.output.kept$Combination.type = paste(linear.regression.output.kept$Combination,linear.regression.output.kept$Type)
linear.regression.output.kept$Combination.type = paste(linear.regression.output.kept$Combination,linear.regression.output.kept$Type)
linear.regression.output.kept$ID = toupper(linear.regression.output.kept$ID)
linear.regression.output.kept$BETA = ifelse(linear.regression.output.kept$REF == "B",-(linear.regression.output.kept$BETA),linear.regression.output.kept$BETA)
linear.regression.output.kept$BETA.group = ifelse(linear.regression.output.kept$BETA >0,"Positive","Negative")
linear.regression.output.kept$log10.P = -log10(linear.regression.output.kept$P)

# motif
linear.regression.output = linear.regression.output.kept[grepl("motif",linear.regression.output.kept$Type),]
linear.regression.output = merge(linear.regression.output,JASPAR_tree.final,by.x="ID",by.y="motif.name",all.x=T)
topN = 10
color_isActive = c("Active"="red","None"= "black")
color_isActive2 = c("Active"="red","None"= "white")

color_motif = c("Single"="#969696","Multiple"= "#000000")
shape_genotype = c("Positive" = 24,"Negative" = 25)
Family = "MER11_combined_F2"

kept_motifs_MER11_F2 = data.frame("Frame"=NA,"Family"=NA,"Type"=NA,"Motif"=NA,"Species"=NA,"thredhold"=NA,"is_top"=NA)
for(Family in Families[grepl("_combined_F",Families)]){
  #
  Summary_Table2.fimoAssic = Summary_Table2.kept[!is.na(Summary_Table2.kept$iPSC.alpha.range),]   ## exclude positive and negative
  Annotation_plot.group = Summary_Table2.fimoAssic[Summary_Table2.fimoAssic$Family.frame == Family & 
                                                     !is.na(Summary_Table2.fimoAssic$Family.frame) & 
                                                     !(Summary_Table2.fimoAssic$Family %in% "MER11D"),]    ## there is no MER11D frame was designed
  Annotation_plot.group$subfamily.name.final = ifelse(is.na(Annotation_plot.group$subfamily.name.final),paste(gsub("MER","",Annotation_plot.group$Family),"_N",sep=""),Annotation_plot.group$subfamily.name.final)
  Annotation_plot.group = Annotation_plot.group[order(-Annotation_plot.group$iPSC.alpha.Zscore),]     # order by the iPSC.alpha.Zscore
  if (grepl("^MER11",Family)){
    #Groupings = data.frame("Group2"=c("All","G1","G2","G3","G4"))
    Groupings = data.frame("Group2"=c("All","G4"))
    Groupings$Group1 = ifelse(Groupings$Group2 == "All","All",NA)
    Groupings$Group1 = ifelse(grepl("^G",Groupings$Group2),"group",Groupings$Group1)
  } else {
    Groupings = data.frame("Group2"=c("All"))
    Groupings$Group1 = ifelse(Groupings$Group2 == "All","All",NA)
  }
  #
  rowID = 1
  Species = "All"
  for (Species in c("All","hg19","hg19_panTro4","Comparison")){  ## species
    for(rowID in 1:nrow(Groupings)){
      Group = Groupings[rowID,]$Group2
      Type = "motif"
#      for (Type in c("motif","variantAll.indel","variantAll")){
      #for (Type in unique(linear.regression.output.kept$Type)){
        if (grepl("_BB",Type)){
          next
        }
        # 1. Annotation plot group
        if (Groupings[rowID,]$Group1 == "All"){
          Annotation_plot.group.sub = Annotation_plot.group
          linear.regression.output.sub = linear.regression.output[linear.regression.output$Combination.type == paste(Group,Type) & 
                                                                    linear.regression.output$Family == Family,]
        } else if (Groupings[rowID,]$Group1 == "group"){
          Annotation_plot.group.sub = Annotation_plot.group[Annotation_plot.group$label.final.correctedName %in% c(list_MER11[list_MER11$FG == Group,]$subfamily,
                                                                                                                   list_MER11.macFas5.distance[list_MER11.macFas5.distance$FG == Group,]$subfamily),]
          Annotation_plot.group.sub = Annotation_plot.group.sub[!is.na(Annotation_plot.group.sub$label.final.correctedName),]
          linear.regression.output.sub = linear.regression.output[linear.regression.output$Combination.type == paste(paste("MER11_F",Group,sep=""),Type) & 
                                                                    linear.regression.output$Family == Family,]
        } 
        
        # reverse order by the value
        Annotation_plot.group.sub = Annotation_plot.group.sub[order(-Annotation_plot.group.sub$iPSC.alpha.Zscore),] # order by the iPSC.alpha
        Annotation_plot.group.sub$uniqueID_new = factor(Annotation_plot.group.sub$uniqueID_new,levels=unique(Annotation_plot.group.sub$uniqueID_new))
        
        # 2. linear regression output
        # if (!grepl("logFC",Type)){
        #   linear.regression.output.sub.both = linear.regression.output[linear.regression.output$Combination.type %in% c(paste(paste("MER11_F",Group,sep=""),Type),paste(paste("MER11_F",Group,sep="")," ",Type,".logFC",sep="")) & 
        #                                                                  linear.regression.output$Family == Family,]
        # }  

        ### select based on the species 2023/7/7
        if (Species == "All") {
          Annotation_plot.group.sub = Annotation_plot.group.sub[Annotation_plot.group.sub$species %in% c("hg19","panTro4","macFas5"),]
        } else if (Species == "Comparison" & !(Groupings[rowID,]$Group1 %in% c("Combined","Combined2"))) {
          Annotation_plot.group.sub = Annotation_plot.group.sub[Annotation_plot.group.sub$species %in% c("hg19","panTro4","macFas5"),]
        } else {
          Annotation_plot.group.sub = Annotation_plot.group.sub[Annotation_plot.group.sub$species %in% Species,]
        }
        if (Species == "Comparison"){
          linear.regression.output.sub = linear.regression.output.sub[linear.regression.output.sub$Species %in% c("hg19","panTro4","macFas5"),]
          linear.regression.output.sub$Combination = paste(linear.regression.output.sub$Species,linear.regression.output.sub$Combination)
          linear.regression.output.sub$Combination.type = paste(linear.regression.output.sub$Species,linear.regression.output.sub$Combination.type)
          # if (!grepl("logFC",Type)){
          #   linear.regression.output.sub.both = linear.regression.output.sub.both[linear.regression.output.sub.both$Species %in% c("hg19","panTro4","macFas5"),]
          #   linear.regression.output.sub.both$Combination = paste(linear.regression.output.sub.both$Species,linear.regression.output.sub.both$Combination)
          #   linear.regression.output.sub.both$Combination.type = paste(linear.regression.output.sub.both$Species,linear.regression.output.sub.both$Combination.type)
          # }
        } else {
          linear.regression.output.sub = linear.regression.output.sub[linear.regression.output.sub$Species == Species,]
          # if (!grepl("logFC",Type)){
          #   linear.regression.output.sub.both = linear.regression.output.sub.both[linear.regression.output.sub.both$Species == Species,]
          # }
        }
        
        ### skip 
        if (nrow(linear.regression.output.sub) == 0 | (Species == "Comparison" & Groupings[rowID,]$Group1 %in% c("Combined","Combined2"))){
          next
        }
        
        ### prepare data for the correlation between ipsc activity and log2FC
        # if (!grepl("logFC",Type)){
        #   linear.regression.output.sub.both$uniqueID = paste(linear.regression.output.sub.both$Combination,linear.regression.output.sub.both$Family,linear.regression.output.sub.both$ID)
        #   linear.regression.output.sub.both1=linear.regression.output.sub.both[!grepl("logFC",linear.regression.output.sub.both$Type),]
        #   linear.regression.output.sub.both2=linear.regression.output.sub.both[grepl("logFC",linear.regression.output.sub.both$Type),]
        #   linear.regression.output.sub.both.final = merge(linear.regression.output.sub.both1,linear.regression.output.sub.both2,by="uniqueID",all=T)
        #   linear.regression.output.sub.both.final$Label.final = ifelse(linear.regression.output.sub.both.final$log10.P.x>=10 | linear.regression.output.sub.both.final$log10.P.y >=10,linear.regression.output.sub.both.final$ID.x,NA)
        # }
        
        linear.regression.output.sub = linear.regression.output.sub[order(linear.regression.output.sub$P),]
        linear.regression.output.sub$log10.P = -log10(linear.regression.output.sub$P)
        linear.regression.output.sub$is_sig_P10 = ifelse(linear.regression.output.sub$log10.P>=10,"*",NA)
        linear.regression.output.sub = linear.regression.output.sub[order(-linear.regression.output.sub$log10.P),]
        linear.regression.output.sub$ID = as.character(linear.regression.output.sub$ID)
        #linear.regression.output.sub$ID = factor(linear.regression.output.sub$ID,levels = unique(linear.regression.output.sub$ID))
        
        # if the motif allele is the majority
        linear.regression.output.sub$BETA = ifelse(linear.regression.output.sub$REF == "B",-(linear.regression.output.sub$BETA),linear.regression.output.sub$BETA)
        linear.regression.output.sub$BETA.group = ifelse(linear.regression.output.sub$BETA >0,"Positive","Negative")
        
        # get JASPAR tree and order
        JASPAR_tree.final.sub = JASPAR_tree.final[JASPAR_tree.final$motif.name %in% linear.regression.output.sub$ID,]
        linear.regression.output.sub$Label = NA
        linear.regression.output.sub$Label = ifelse(linear.regression.output.sub$log10.P>=5,as.character(linear.regression.output.sub$ID),linear.regression.output.sub$Label)
        linear.regression.output.sub$ID = factor(linear.regression.output.sub$ID,levels = unique(JASPAR_tree.final.sub$motif.name))
        
        # cluster and branch color
        linear.regression.output.sub.color = linear.regression.output.sub[!duplicated(linear.regression.output.sub$cluster_color),]
        linear.regression.output.sub$branch_color2 = linear.regression.output.sub$branch_color
        
        ## color lowly expressed TFs
        linear.regression.output.sub$branch_color2 = ifelse(linear.regression.output.sub$X0h < 1,"low.exp",linear.regression.output.sub$branch_color2)
        linear.regression.output.sub$branch_color = factor(linear.regression.output.sub$branch_color)
        linear.regression.output.sub$branch_color2 = factor(linear.regression.output.sub$branch_color2)
        linear.regression.output.sub$cluster_color = factor(linear.regression.output.sub$cluster_color)
        linear.regression.output.sub$BETA.group = factor(linear.regression.output.sub$BETA.group)
        
        # set the maximum expression and BETA color
        linear.regression.output.sub$X0h = ifelse(linear.regression.output.sub$X0h>=100,100,linear.regression.output.sub$X0h)    ### still set a maximum of 100
        linear.regression.output.sub$X0h.modified = ifelse(linear.regression.output.sub$X0h.modified>=100,100,linear.regression.output.sub$X0h)
        linear.regression.output.sub$X0h.modified = ifelse(linear.regression.output.sub$X0h.modified>1,linear.regression.output.sub$X0h.modified,NA)
        linear.regression.output.sub$BETA.modified = ifelse(linear.regression.output.sub$branch_color2 == "low.exp",NA,linear.regression.output.sub$BETA)
        
        # assign color
        branch.color = levels(linear.regression.output.sub$branch_color)
        names(branch.color) = levels(linear.regression.output.sub$branch_color)
        branch.color2 = levels(linear.regression.output.sub$branch_color2)
        names(branch.color2) = levels(linear.regression.output.sub$branch_color2)
        branch.color2['low.exp'] = "white"
        cluster.color = levels(linear.regression.output.sub$cluster_color)
        names(cluster.color) = levels(linear.regression.output.sub$cluster_color)
        
        # remove redundancy by branch and beta group
        linear.regression.output.sub.uniq = linear.regression.output.sub[order(-linear.regression.output.sub$log10.P),]
        linear.regression.output.sub.uniq = linear.regression.output.sub.uniq[!duplicated(linear.regression.output.sub.uniq[,c("branch_color","BETA.group","Combination.type")]),]
        linear.regression.output.sub.uniq$ID = factor(linear.regression.output.sub.uniq$ID,levels = levels(linear.regression.output.sub.uniq$ID)[levels(linear.regression.output.sub.uniq$ID) %in% linear.regression.output.sub.uniq$ID])
        
        # keep associated motifs log10 p value >=5
        linear.regression.output.sub.top = linear.regression.output.sub[linear.regression.output.sub$log10.P>=5,] 
        linear.regression.output.sub.top = linear.regression.output.sub.top[order(-linear.regression.output.sub.top$log10.P),]
        
        # exclude motifs from the same branch and the same association direction (updated 2023/7/11)
        # linear.regression.output.sub.top.uniq = linear.regression.output.sub.top[!duplicated(linear.regression.output.sub.top$branch_color),]
        linear.regression.output.sub.top.uniq = linear.regression.output.sub.top[!duplicated(linear.regression.output.sub.top[,c("branch_color","BETA.group")]),]
        
        # 3. Motifs 
        linear.regression.output.sub.top.order1 = linear.regression.output.sub.top[linear.regression.output.sub.top$BETA>0,]
        linear.regression.output.sub.top.order1 = linear.regression.output.sub.top.order1[order(-linear.regression.output.sub.top.order1$log10.P),]
        linear.regression.output.sub.top.order2 = linear.regression.output.sub.top[linear.regression.output.sub.top$BETA<=0,]
        linear.regression.output.sub.top.order2 = linear.regression.output.sub.top.order2[order(linear.regression.output.sub.top.order2$log10.P),]
        linear.regression.output.sub.top.order = rbind(linear.regression.output.sub.top.order1,linear.regression.output.sub.top.order2)
        rm(linear.regression.output.sub.top.order1,linear.regression.output.sub.top.order2)
        
        # Fimo dataframe, focus on the JASPAR2022 database
        Fimo.all.sum.sub = Fimo.all.sum[Fimo.all.sum$uniqueID.new %in% as.character(Annotation_plot.group.sub$uniqueID_new),]
        Fimo.all.sum.sub = Fimo.all.sum.sub[Fimo.all.sum.sub$DB == "JASPAR2022",]
        Fimo.all.sum.sub$motif_name = toupper(Fimo.all.sum.sub$motif_name)
        Fimo.all.sum.sub$uniqueID2 = paste(Fimo.all.sum.sub$DB,Fimo.all.sum.sub$uniqueID.new,Fimo.all.sum.sub$motif_name)
        Fimo.all.sum.sub = Fimo.all.sum.sub[order(-Fimo.all.sum.sub$n),]
        Fimo.all.sum.sub = Fimo.all.sum.sub[!duplicated(Fimo.all.sum.sub$uniqueID2),]
        Fimo.all.sum.sub = Fimo.all.sum.sub[Fimo.all.sum.sub$motif_name %in% as.character(linear.regression.output.sub.top$ID),]
        if (nrow(Fimo.all.sum.sub) == 0) {
          Fimo.all.sum.sub = Fimo.all.sum[Fimo.all.sum$motif_name %in% unique(as.character(linear.regression.output.sub[1:topN,]$ID)),]
          Group2 = paste(Group,"-log10(padj) <5")
          linear.regression.output.sub.top = linear.regression.output.sub[1:topN,]
        } else {
          Group2 = paste(Group,"-log10(padj) >=5")
          Fimo.all.sum.sub$motif_name = factor(Fimo.all.sum.sub$motif_name,levels=rev(unique(linear.regression.output.sub.top.order$ID)))
        }
        Fimo.all.sum.sub$uniqueID.new = factor(Fimo.all.sum.sub$uniqueID.new,levels=unique(as.character(Annotation_plot.group.sub$uniqueID_new)))
        Fimo.all.sum.sub$n.group = ifelse(Fimo.all.sum.sub$n >1,"Multiple","Single")
        
        ### 4. motif association plot
        #for (Branch in c("top","all")){
        for (Branch in c("all")){
          if (Branch == "top"){
            linear.regression.output.plot = linear.regression.output.sub.uniq
          } else {
            linear.regression.output.plot = linear.regression.output.sub
          }
          p3<-ggplot(linear.regression.output.plot, aes(ID,log10.P)) +
            geom_point(aes(shape = BETA.group,size=BETA,fill = (X0h)))+
            scale_shape_manual(values=shape_genotype)+
            scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 1,na.value = "grey50")+
            geom_hline(yintercept = c(5,10), linetype="dotted",color = "black", linewidth=.5)+
            geom_text_repel(aes(label = Label,color = Combination),max.overlaps = 100) +
            scale_x_discrete(drop = FALSE)+
            ggtitle(paste(Family,Group))+
            xlab("")+
            ylab("-log10(P value)")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  #panel.background = element_rect(fill = "#f0f0f0"),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.text.y=element_text(colour="black",size=rel(0.7)),
                  axis.text.x=element_blank(),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1)))
          # assign low expressed color to grey
          p3_beta<-ggplot(linear.regression.output.plot, aes(ID,log10.P)) +
            geom_point(aes(shape = BETA.group,fill=BETA,size = (X0h)))+
            scale_shape_manual(values=shape_genotype)+
            scale_size(limits = c(0,100))+
            scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 0,na.value = "grey50")+
            geom_hline(yintercept = c(5,10), linetype="dotted",color = "black", linewidth=.5)+
            geom_text_repel(aes(label = Label,color = Combination),max.overlaps = 100) +
            #scale_color_manual(values=color_genotype)+
            scale_x_discrete(drop = FALSE)+
            ggtitle(paste(Family,Group))+
            xlab("")+
            ylab("-log10(P value)")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  #panel.background = element_rect(fill = "#f0f0f0"),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.text.y=element_text(colour="black",size=rel(0.7)),
                  axis.text.x=element_blank(),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1))) 
          p3_branch<-ggplot(linear.regression.output.plot, aes(x=ID,y=1)) +
            geom_tile(aes(fill = branch_color),col = NA,size=0.1) + 
            scale_fill_manual(values =branch.color)+
            scale_x_discrete(drop = FALSE)+
            xlab("x")+
            ylab("brahch")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  #axis.line = element_line(colour = "black"),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(), 
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1)))
          
          p3_cluster<-ggplot(linear.regression.output.plot, aes(x=ID,y=1)) +
            geom_tile(aes(fill = cluster_color),col = NA,size=0.1) + 
            scale_fill_manual(values =cluster.color)+
            scale_x_discrete(drop = FALSE)+
            xlab("x")+
            ylab("cluster")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  #axis.line = element_line(colour = "black"),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(), 
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1)))
          gC <- ggplotGrob(p3)
          gC_beta <- ggplotGrob(p3_beta)
          gC_branch <- ggplotGrob(p3_branch)
          gC_cluster <- ggplotGrob(p3_cluster)
          
          # between species
          if (Species == "Comparison" & Groupings[rowID,]$Group2 == "G4"){
            linear.regression.output.plot.h = linear.regression.output.plot[linear.regression.output.plot$Species %in% c("hg19"),]
            linear.regression.output.plot.h$BETA.group = as.character(linear.regression.output.plot.h$BETA.group)
            linear.regression.output.plot.c = linear.regression.output.plot[linear.regression.output.plot$Species %in% c("panTro4"),]
            linear.regression.output.plot.c$BETA.group = as.character(linear.regression.output.plot.c$BETA.group)
            linear.regression.output.plot.m = linear.regression.output.plot[linear.regression.output.plot$Species %in% c("macFas5"),]
            linear.regression.output.plot.m$BETA.group = as.character(linear.regression.output.plot.m$BETA.group)
            linear.regression.output.plot.comparison = merge(linear.regression.output.plot.h[,c("ID","log10.P","BETA","BETA.group")],linear.regression.output.plot.c[,c("ID","log10.P","BETA","BETA.group")],by="ID",all=T)
            linear.regression.output.plot.comparison = merge(linear.regression.output.plot.comparison,linear.regression.output.plot.m[,c("ID","log10.P","BETA","BETA.group")],by="ID",all=T)
            linear.regression.output.plot.comparison$log10.P = ifelse(is.na(linear.regression.output.plot.comparison$log10.P),0,linear.regression.output.plot.comparison$log10.P)
            linear.regression.output.plot.comparison$log10.P.x = ifelse(is.na(linear.regression.output.plot.comparison$log10.P.x),0,linear.regression.output.plot.comparison$log10.P.x)
            linear.regression.output.plot.comparison$log10.P.y = ifelse(is.na(linear.regression.output.plot.comparison$log10.P.y),0,linear.regression.output.plot.comparison$log10.P.y)
            
            linear.regression.output.plot.comparison$label.final = ifelse(linear.regression.output.plot.comparison$log10.P>=5 | linear.regression.output.plot.comparison$log10.P.x >=5 | 
                                                                            linear.regression.output.plot.comparison$log10.P.y >=5,as.character(linear.regression.output.plot.comparison$ID),NA)
            linear.regression.output.plot.comparison$BETA.final.m = ifelse(linear.regression.output.plot.comparison$log10.P.x >= linear.regression.output.plot.comparison$log10.P,
                                                                           linear.regression.output.plot.comparison$BETA.x,
                                                                           linear.regression.output.plot.comparison$BETA)
            linear.regression.output.plot.comparison$BETA.final.mg = ifelse(linear.regression.output.plot.comparison$log10.P.x >= linear.regression.output.plot.comparison$log10.P,
                                                                            linear.regression.output.plot.comparison$BETA.group.x,
                                                                            linear.regression.output.plot.comparison$BETA.group)
            linear.regression.output.plot.comparison$BETA.final.c = ifelse(linear.regression.output.plot.comparison$log10.P.x >= linear.regression.output.plot.comparison$log10.P.y,
                                                                           linear.regression.output.plot.comparison$BETA.x,
                                                                           linear.regression.output.plot.comparison$BETA.y)
            linear.regression.output.plot.comparison$BETA.final.cg = ifelse(linear.regression.output.plot.comparison$log10.P.x >= linear.regression.output.plot.comparison$log10.P.y,
                                                                            linear.regression.output.plot.comparison$BETA.group.x,
                                                                            linear.regression.output.plot.comparison$BETA.group.y)
            linear.regression.output.plot.comparison = merge(linear.regression.output.plot.comparison,linear.regression.output.plot[!duplicated(linear.regression.output.plot$ID),c("ID","X0h")],by="ID",all.x=T)
            pS1 <- ggplot(linear.regression.output.plot.comparison, aes(log10.P.x,log10.P.y)) +
              geom_point(aes(shape = BETA.final.cg,fill=(X0h)))+
              scale_shape_manual(values=shape_genotype)+
              scale_size(limits = c(0,100))+
              scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",na.value = "grey50")+
              geom_hline(yintercept = c(5,10), linetype="dotted",color = "white", linewidth=.5)+
              geom_vline(xintercept = c(5,10), linetype="dotted",color = "white", linewidth=.5)+
              geom_text_repel(aes(label = label.final),max.overlaps = 100) +
              ggtitle("Association between motifs and activity change")+
              xlab("human (-log10 (pvalue))")+
              ylab("chimpanzee (-log10 (pvalue))")+
              ggtitle("human vs chimpanzee")+
              theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "#f0f0f0"),
                    #panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                    axis.text.y=element_text(colour="black",size=rel(1)),
                    axis.text.x=element_text(colour="black",size=rel(1)),
                    axis.title=element_text(colour="black",size=rel(1)),
                    legend.position="right",
                    legend.background = element_blank(),
                    legend.text=element_text(size=rel(1))) 
            pS2 <- ggplot(linear.regression.output.plot.comparison, aes(log10.P.x,log10.P)) +
              geom_point(aes(shape = BETA.final.mg,fill=(X0h)))+
              scale_shape_manual(values=shape_genotype)+
              scale_size(limits = c(0,100))+
              scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",na.value = "grey50")+
              geom_hline(yintercept = c(5,10), linetype="dotted",color = "white", linewidth=.5)+
              geom_vline(xintercept = c(5,10), linetype="dotted",color = "white", linewidth=.5)+
              geom_text_repel(aes(label = label.final),max.overlaps = 100) +
              ggtitle("Association between motifs and activity change")+
              xlab("human (-log10 (pvalue))")+
              ylab("macaque (-log10 (pvalue))")+
              ggtitle("human vs macaque")+
              theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "#f0f0f0"),
                    #panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                    axis.text.y=element_text(colour="black",size=rel(1)),
                    axis.text.x=element_text(colour="black",size=rel(1)),
                    axis.title=element_text(colour="black",size=rel(1)),
                    legend.position="right",
                    legend.background = element_blank(),
                    legend.text=element_text(size=rel(1))) 
            
            gS1 <- ggplotGrob(pS1)
            gS2 <- ggplotGrob(pS2)
            gSC = rbind(gS1,gS2,size = "first")
            
            pdf(paste("Figure_5C_",Family,"_",Species,"_",Group,"_",Type,"_",Date,"-species_comparison.pdf",sep=""),    # create PNG for the heat map
                width = 6,        # 5 x 300 pixels
                height = 10,
                pointsize = 10)        # smaller font size
            grid.draw(gSC)
            dev.off()
          }
          
          g2 = rbind(gC,gC_beta,gC_branch,gC_cluster, size = "first")
          if (Group == "All" & Species == "All") {
            Width = 24
            pdf(paste("Figure_5A_S6_",Family,"_",Species,"_",Group,"_",Type,"_",Branch,"_",Date,"-motif.pdf",sep=""),    # create PNG for the heat map
              width = Width,        # 5 x 300 pixels
              height = 16,
              pointsize = 10)        # smaller font size
            grid.draw(g2)
            dev.off()
          }
        }
        
        motifOrder.logP_5 = unique(as.character(linear.regression.output.sub.top$ID))
        motifOrder.logP_5.uniq = motifOrder.logP_5[motifOrder.logP_5 %in% linear.regression.output.sub.top.uniq$ID]
        motifOrder.logP_10 = unique(as.character(linear.regression.output.sub.top[linear.regression.output.sub.top$log10.P>=10,]$ID))
        motifOrder.logP_10.uniq = motifOrder.logP_10[motifOrder.logP_10 %in% as.character(linear.regression.output.sub.top.uniq$ID)]
        
        ###########################
        ### 5. correlation analysis
        if (Groupings[rowID,]$Group1 %in% c("group","All") & !grepl("logFC",Type) & length(unique(as.character(linear.regression.output.sub.top$ID)))>1){
          Candidate_motifs = unique(as.character(linear.regression.output.sub.top$ID))
          Candidate_instances = unique(as.character(Annotation_plot.group.sub$uniqueID_new))
          Fimo.all.sum.sub.convert.df = data.frame(matrix(ncol=length(Candidate_motifs),nrow=length(Candidate_instances)))
          colnames(Fimo.all.sum.sub.convert.df) = Candidate_motifs
          rownames(Fimo.all.sum.sub.convert.df) = Candidate_instances
          for (rowID.tmp in 1:nrow(Fimo.all.sum.sub.convert.df)){
            for (colID.tmp in 1:ncol(Fimo.all.sum.sub.convert.df)){
              Fimo.all.sum.sub.tmp = Fimo.all.sum.sub[Fimo.all.sum.sub$uniqueID.new == rownames(Fimo.all.sum.sub.convert.df)[rowID.tmp] & Fimo.all.sum.sub$motif_name == colnames(Fimo.all.sum.sub.convert.df)[colID.tmp],]
              if (nrow(Fimo.all.sum.sub.tmp)>0)
                Fimo.all.sum.sub.convert.df[rowID.tmp,colID.tmp] = 1
              #Fimo.all.sum.sub.convert.df[rowID.tmp,colID.tmp] = Fimo.all.sum.sub.tmp$n
            }
          }
          rm(rowID.tmp,colID.tmp)
          rm(Fimo.all.sum.sub.tmp)
          Fimo.all.sum.sub.convert.df[is.na(Fimo.all.sum.sub.convert.df)] <-0
          Fimo.all.sum.sub2.convert.df = Fimo.all.sum.sub.convert.df[,colnames(Fimo.all.sum.sub.convert.df) %in% motifOrder.logP_10]
          summary.tmp = data.frame(apply(Fimo.all.sum.sub.convert.df, 2, function(x) length(unique(x))))
          colnames(summary.tmp) = "uniq.count"
          summary.tmp$motif = rownames(summary.tmp)
          if (nrow(summary.tmp[summary.tmp$uniq.count>1,]) == 0) {
            next
          }
          
          ####################### upset plot
          Fimo.all.sum.sub.convert.df.upset = Fimo.all.sum.sub.convert.df
          Fimo.all.sum.sub.convert.df.upset$uniqueID_new = rownames(Fimo.all.sum.sub.convert.df.upset)
          Fimo.all.sum.sub.convert.df.upset = merge(Fimo.all.sum.sub.convert.df.upset,Annotation_plot.group.sub[,c("uniqueID_new","iPSC.alpha.Zscore.isActive.group","iPSC.alpha.Zscore")],by.x="uniqueID_new",all.x=T)
          Fimo.all.sum.sub.convert.df.upset$iPSC.alpha.Zscore.isActive.group = factor(Fimo.all.sum.sub.convert.df.upset$iPSC.alpha.Zscore.isActive.group,levels=rev(names(iPSC.alpha.Zscore.isActive.group.color)))
          Fimo.all.sum.sub.convert.df.upset$Group = Fimo.all.sum.sub.convert.df.upset$iPSC.alpha.Zscore.isActive.group
          Fimo.all.sum.sub.convert.df.upset$Group2 = ifelse(Fimo.all.sum.sub.convert.df.upset$Group == "No activity","None","Active")
          Fimo.all.sum.sub.convert.df.upset$Group2 = factor(Fimo.all.sum.sub.convert.df.upset$Group2,levels=c("None","Active"))
          color_isActive2  = c("Active"= "red","None"="white")

          # if (length(motifOrder.logP_5.uniq) > 1 & Species == "All" & Group == "G4" & Family == "MER11_combined_F2"){
          #   pdf(paste("Step9_3_",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-upset2.pdf",sep=""),                      
          #       width = 6,        # 5 x 300 pixels
          #       height = 12,
          #       pointsize = 10)        # smaller font size
          #   print(upset(
          #     Fimo.all.sum.sub.convert.df.upset,
          #     motifOrder.logP_5.uniq,min_size=10,width_ratio=0.1, sort_intersections_by=c('degree', 'cardinality'),
          #     base_annotations=list(
          #       'Intersection size'=intersection_size(
          #         counts=TRUE,
          #         text=list(check_overlap=TRUE),
          #         mapping=aes(fill=Group)
          #       ) + 
          #         scale_fill_manual(values=iPSC.alpha.Zscore.isActive.group.color) + 
          #         ggtitle("log10 P >=5 active") +
          #         theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
          #               panel.background = element_blank()) +
          #         ylab('Group')),
          #     set_sizes=(
          #       upset_set_size(
          #         geom=geom_bar(
          #           aes(fill=Group, x=group),
          #           width=0.8),
          #         position='right') +
          #         scale_fill_manual(values=iPSC.alpha.Zscore.isActive.group.color)
          #     ),
          #     # moves legends over the set sizes
          #     guides='over',
          #     annotations = list(
          #       'Length'=list(
          #         aes=aes(x=intersection, y=iPSC.alpha.Zscore),
          #         geom=geom_boxplot(na.rm=TRUE)
          #       ),
          #       'Rating'=list(
          #         aes=aes(x=intersection, y=iPSC.alpha.Zscore),
          #         geom=list(
          #           geom_jitter(aes(), na.rm=TRUE),
          #           geom_violin(alpha=0.5, na.rm=TRUE))
          #       ),
          #       'MPRA Rating'=list(
          #         aes=aes(x=intersection, fill=Group2),
          #         geom=list(
          #           geom_bar(stat='count', position='fill', na.rm=TRUE),
          #           geom_text(
          #             aes(
          #               label=!!aes_percentage(relative_to='intersection'),
          #               group=Group2
          #             ),
          #             color="black",
          #             stat='count',
          #             position=position_fill(vjust = .5)
          #           ),
          #           scale_y_continuous(labels=scales::percent_format())
          #         )
          #       )
          #     )
          #   )
          #   )
          #   dev.off()
          # }
          
          # if (length(motifOrder.logP_10) > 1 & Species == "All" & Group == "G4" & Family == "MER11_combined_F2"){
          #   Width = 6
          #   Height = 12
          #   pdf(paste("Step9_3_",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-upset3.pdf",sep=""),                      
          #       width = Width,        # 5 x 300 pixels
          #       height = Height,
          #       pointsize = 10)        # smaller font size
          #   print(upset(
          #     Fimo.all.sum.sub.convert.df.upset,
          #     motifOrder.logP_10,min_size=10,width_ratio=0.1,sort_intersections_by=c('degree', 'cardinality'), 
          #     base_annotations=list(
          #       'Intersection size'=intersection_size(
          #         counts=TRUE,
          #         text=list(check_overlap=TRUE),
          #         mapping=aes(fill=Group)
          #       ) + 
          #         scale_fill_manual(values=iPSC.alpha.Zscore.isActive.group.color) + 
          #         ggtitle("log10 P >=10 all") +
          #         theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
          #               panel.background = element_blank()) +
          #         ylab('Group')),
          #     set_sizes=(
          #       upset_set_size(
          #         geom=geom_bar(
          #           aes(fill=Group, x=group),
          #           width=0.8),
          #         position='right') +
          #         scale_fill_manual(values=iPSC.alpha.Zscore.isActive.group.color)
          #     ),
          #     # moves legends over the set sizes
          #     guides='over',
          #     annotations = list(
          #       'Length'=list(
          #         aes=aes(x=intersection, y=iPSC.alpha.Zscore),
          #         geom=geom_boxplot(na.rm=TRUE)
          #       ),
          #       'Rating'=list(
          #         aes=aes(x=intersection, y=iPSC.alpha.Zscore),
          #         geom=list(
          #           geom_jitter(aes(), na.rm=TRUE),
          #           geom_violin(alpha=0.5, na.rm=TRUE))
          #       ),
          #       'MPAA Rating'=list(
          #         aes=aes(x=intersection, fill=Group2),
          #         geom=list(
          #           geom_bar(stat='count', position='fill', na.rm=TRUE),
          #           geom_text(
          #             aes(
          #               label=!!aes_percentage(relative_to='intersection'),
          #               group=Group2
          #             ),
          #             color="black",
          #             stat='count',
          #             position=position_fill(vjust = .5)
          #           ),
          #           scale_y_continuous(labels=scales::percent_format())
          #         )
          #       )
          #     )
          #   )
          #   )
          #   dev.off()
          #   if (Family == "MER11_combined_F2" & Type == "motif" & Species == "All"){
          #     if (Group == "MER11A") {
          #       Height = 24
          #     }
          #     Width = 6
          #     pdf(paste("Step9_3_",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-upset3_2.pdf",sep=""),                      
          #         width = Width,        # 5 x 300 pixels
          #         height = Height,
          #         pointsize = 10)        # smaller font size
          #     print(upset(
          #       Fimo.all.sum.sub.convert.df.upset,
          #       motifOrder.logP_10,min_size=10,width_ratio=0.1,sort_intersections_by=c('degree', 'cardinality'), 
          #       base_annotations=list(
          #         'Intersection size'=intersection_size(
          #           counts=TRUE,
          #           text=list(check_overlap=TRUE),
          #           mapping=aes(fill=Group)
          #         ) + 
          #           scale_fill_manual(values=iPSC.alpha.Zscore.isActive.group.color) + 
          #           ggtitle("log10 P >=10 all") +
          #           theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
          #                 panel.background = element_blank()) +
          #           ylab('Group')),
          #       set_sizes=(
          #         upset_set_size(
          #           geom=geom_bar(
          #             aes(fill=Group, x=group),
          #             width=0.8),
          #           position='right') +
          #           scale_fill_manual(values=iPSC.alpha.Zscore.isActive.group.color)
          #       ),
          #       # moves legends over the set sizes
          #       guides='over',
          #       annotations = list(
          #         'Length'=list(
          #           aes=aes(x=intersection, y=iPSC.alpha.Zscore),
          #           geom=geom_boxplot(na.rm=TRUE)
          #         ),
          #         'Rating'=list(
          #           aes=aes(x=intersection, y=iPSC.alpha.Zscore),
          #           geom=list(
          #             geom_jitter(aes(), na.rm=TRUE),
          #             geom_violin(alpha=0.5, na.rm=TRUE))
          #         ),
          #         'MPRA Rating'=list(
          #           aes=aes(x=intersection, fill=Group2),
          #           geom=list(
          #             geom_bar(stat='count', position='fill', na.rm=TRUE),
          #             geom_text(
          #               aes(
          #                 label=!!aes_percentage(relative_to='intersection'),
          #                 group=Group2
          #               ),
          #               color="black",
          #               stat='count',
          #               position=position_fill(vjust = .5)
          #             ),
          #             scale_y_continuous(labels=scales::percent_format())
          #           )
          #         )
          #       )
          #     )
          #     )
          #     dev.off()
          #   }
          
            if (length(motifOrder.logP_10.uniq) > 1 & Species == "All" & Group == "G4" & Family == "MER11_combined_F2"){
              Width = 6
              Height = 12
              if (Family == "MER11_combined_F2" & Type == "motif"){
                Height = 16
              }
              pdf(paste("Figure_5B-",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-upset4.pdf",sep=""),
                  width = Width,        # 5 x 300 pixels
                  height = 12,
                  pointsize = 10)        # smaller font size
              print(upset(
                Fimo.all.sum.sub.convert.df.upset,
                motifOrder.logP_10.uniq,min_size=10,width_ratio=0.1,sort_intersections_by=c('degree', 'cardinality'),
                base_annotations=list(
                  'Intersection size'=intersection_size(
                    counts=TRUE,
                    text=list(check_overlap=TRUE),
                    mapping=aes(fill=Group)
                  ) +
                    scale_fill_manual(values=iPSC.alpha.Zscore.isActive.group.color) +
                    ggtitle("log10 P >=10 active") +
                    theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
                          panel.background = element_blank()) +
                    ylab('Group')),
                set_sizes=(
                  upset_set_size(
                    geom=geom_bar(
                      aes(fill=Group, x=group),
                      width=0.8),
                    position='right') +
                    scale_fill_manual(values=iPSC.alpha.Zscore.isActive.group.color)
                ),
                # moves legends over the set sizes
                guides='over',
                annotations = list(
                  'Length'=list(
                    aes=aes(x=intersection, y=iPSC.alpha.Zscore),
                    geom=geom_boxplot(na.rm=TRUE)
                  ),
                  'Rating'=list(
                    aes=aes(x=intersection, y=iPSC.alpha.Zscore),
                    geom=list(
                      geom_jitter(aes(), na.rm=TRUE),
                      geom_violin(alpha=0.5, na.rm=TRUE))
                  ),
                  'MPRA Rating'=list(
                    aes=aes(x=intersection, fill=Group2),
                    geom=list(
                      geom_bar(stat='count', position='fill', na.rm=TRUE),
                      geom_text(
                        aes(
                          label=!!aes_percentage(relative_to='intersection'),
                          group=Group2
                        ),
                        color="black",
                        stat='count',
                        position=position_fill(vjust = .5)
                      ),
                      scale_y_continuous(labels=scales::percent_format())
                    )
                  )
                )
              )
              )
              dev.off()
            }
          # } 
          
          # log10 P >= 5
          res.logP_5 <- cor(as.matrix(Fimo.all.sum.sub.convert.df),use = "complete.obs")
          # res.logP_5.uniq <- cor(as.matrix(Fimo.all.sum.sub.convert.df[,colnames(Fimo.all.sum.sub.convert.df) %in% as.character(linear.regression.output.sub.top.uniq$ID)]),use = "complete.obs")
          # kept_motifs_MER11_F2.tmp = data.frame("Frame"=rep(Family,nrow(res.logP_5)),
          #                                       "Family"=rep(Group,nrow(res.logP_5)),
          #                                       "Type"=rep(Type,nrow(res.logP_5)),
          #                                       "Motif"=c(colnames(res.logP_5)),
          #                                       "Species"=rep(Species,nrow(res.logP_5)),
          #                                       "thredhold" = "logP_5",
          #                                       "is_top"=NA)
          # if (nrow(res.logP_5.uniq)>=1){
          #   kept_motifs_MER11_F2.tmp$is_top = ifelse(kept_motifs_MER11_F2.tmp$Motif %in% colnames(res.logP_5.uniq),"yes",kept_motifs_MER11_F2.tmp$is_top)
          # }
          # kept_motifs_MER11_F2 = rbind(kept_motifs_MER11_F2,kept_motifs_MER11_F2.tmp)
          # 
          # 
          pdf(paste("Figure_S7-",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-corrplot.logP_5.pdf",sep=""),
              width = 6,
              height = 12,
              pointsize = 10)        # smaller font size
          p.corr1 = corrplot(res.logP_5,
                             addCoef.col = 'black',
                             #tl.srt = 45,
                             #p.mat = testRes$p,
                             p.mat = NULL,
                             #number.cex = 1.5,
                             method = 'square', ## or 'circle'
                             diag = TRUE,
                             order = 'hclust',
                             sig.level = 0.005,
                             #insig='blank',
                             type = 'upper',
                             bg="white",
                             hclust.method = "complete",
                             cl.pos = 'b', col = COL2('RdBu'))
          dev.off()
          # motifOrder.logP_5 = rownames(p.corr1$corr)
          # if (length(as.character(linear.regression.output.sub.top.uniq$ID))>1){
          #   pdf(paste("Step9_3_",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-corrplot.logP_5.uniq.pdf",sep=""),                      
          #       width = 12,
          #       height = 12,
          #       pointsize = 10)        # smaller font size
          #   p.corr1.uniq = corrplot(res.logP_5.uniq, 
          #                           addCoef.col = 'black', 
          #                           #tl.srt = 45,
          #                           #p.mat = testRes$p,
          #                           p.mat = NULL,
          #                           #number.cex = 1.5, 
          #                           method = 'square', ## or 'circle'
          #                           diag = TRUE,
          #                           order = 'hclust',
          #                           sig.level = 0.005,
          #                           #insig='blank',
          #                           type = 'upper',
          #                           bg="white",
          #                           hclust.method = "complete",
          #                           cl.pos = 'b', col = COL2('RdBu'))      
          #   dev.off()
          # }
          
          # # log10 P >= 10
          # if (length(motifOrder.logP_10) > 1) {
          #   summary.tmp = data.frame(apply(Fimo.all.sum.sub2.convert.df, 2, function(x) length(unique(x))))
          #   colnames(summary.tmp) = "uniq.count"
          #   summary.tmp$motif = rownames(summary.tmp)
          #   if (nrow(summary.tmp[summary.tmp$uniq.count>1,]) == 0){
          #     next
          #   }
          #   res.logP_10 <- cor(as.matrix(Fimo.all.sum.sub2.convert.df),use = "complete.obs")
          #   kept_motifs_MER11_F2.tmp = data.frame("Frame"=rep(Family,nrow(res.logP_10)),
          #                                         "Family"=rep(Group,nrow(res.logP_10)),
          #                                         "Type"=rep(Type,nrow(res.logP_10)),
          #                                         "Motif"=c(colnames(res.logP_10)),
          #                                         "Species"=rep(Species,nrow(res.logP_10)),
          #                                         "thredhold" = "logP_10",
          #                                         "is_top"=NA)
          #   pdf(paste("Step9_3_",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-corrplot.logP_10.pdf",sep=""),                      
          #       width = 12,
          #       height = 12,
          #       pointsize = 10)        # smaller font size
          #   p.corr2 = corrplot(res.logP_10,
          #                      addCoef.col = 'black', 
          #                      #tl.srt = 45,
          #                      #p.mat = testRes$p,
          #                      p.mat = NULL,
          #                      #number.cex = 1.5, 
          #                      method = 'square', ## or 'circle'
          #                      diag = TRUE,
          #                      order = 'hclust',
          #                      sig.level = 0.005,
          #                      #insig='blank',
          #                      type = 'upper',
          #                      bg="white",
          #                      hclust.method = "complete",
          #                      cl.pos = 'b', col = COL2('RdBu'))      
          #   dev.off()
          #   motifOrder.logP_10 = rownames(p.corr2$corr)
          #   
          #   # unique 
          #   Fimo.all.sum.sub2.convert.df.uniq = Fimo.all.sum.sub2.convert.df[,colnames(Fimo.all.sum.sub2.convert.df) %in% as.character(linear.regression.output.sub.top.uniq$ID)]
          #   if (length(colnames(Fimo.all.sum.sub2.convert.df)[colnames(Fimo.all.sum.sub2.convert.df) %in% as.character(linear.regression.output.sub.top.uniq$ID)])>=1){
          #     kept_motifs_MER11_F2.tmp$is_top = ifelse(kept_motifs_MER11_F2.tmp$Motif %in% colnames(Fimo.all.sum.sub2.convert.df.uniq),"yes",kept_motifs_MER11_F2.tmp$is_top)
          #   }
          #   
          #   if (length(colnames(Fimo.all.sum.sub2.convert.df)[colnames(Fimo.all.sum.sub2.convert.df) %in% as.character(linear.regression.output.sub.top.uniq$ID)])>1){
          #     summary.tmp = data.frame(apply(Fimo.all.sum.sub2.convert.df.uniq, 2, function(x) length(unique(x))))
          #     colnames(summary.tmp) = "uniq.count"
          #     summary.tmp$motif = rownames(summary.tmp)
          #     if (nrow(summary.tmp[summary.tmp$uniq.count>1,]) == 0){
          #       next
          #     }
          #     res.logP_10.uniq <- cor(as.matrix(Fimo.all.sum.sub2.convert.df.uniq),use = "complete.obs")
          #     pdf(paste("Step9_3_",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-corrplot.logP_10.uniq.pdf",sep=""),                      
          #         width = 12,
          #         height = 12,
          #         pointsize = 10)        # smaller font size
          #     p.corr2.uniq = corrplot(res.logP_10.uniq, 
          #                             addCoef.col = 'black', 
          #                             #tl.srt = 45,
          #                             #p.mat = testRes$p,
          #                             p.mat = NULL,
          #                             #number.cex = 1.5, 
          #                             method = 'square', ## or 'circle'
          #                             diag = TRUE,
          #                             order = 'hclust',
          #                             sig.level = 0.005,
          #                             #insig='blank',
          #                             type = 'upper',
          #                             bg="white",
          #                             hclust.method = "complete",
          #                             cl.pos = 'b', col = COL2('RdBu'))
          #     dev.off()  
          #   } 
          #   kept_motifs_MER11_F2 = rbind(kept_motifs_MER11_F2,kept_motifs_MER11_F2.tmp)
          # }
          
          #### plot the branch and cluster
          # a. all motifs logP>=5
          linear.regression.output.sub.plot1 = linear.regression.output.sub[as.character(linear.regression.output.sub$ID) %in% motifOrder.logP_5,]
          linear.regression.output.sub.plot1$ID = factor(linear.regression.output.sub.plot1$ID,levels=motifOrder.logP_5)
          
          p.corr1_branch<-ggplot(linear.regression.output.sub.plot1, aes(x=ID,y=1)) +
            geom_tile(aes(fill = branch_color),col = NA,size=0.1) + 
            scale_fill_manual(values =branch.color)+
            scale_x_discrete(drop = FALSE)+
            xlab("x")+
            ylab("brahch")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  #axis.line = element_line(colour = "black"),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(), 
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1)))
          
          p.corr1_cluster<-ggplot(linear.regression.output.sub.plot1, aes(x=ID,y=1)) +
            geom_tile(aes(fill = cluster_color),col = NA,size=0.1) + 
            scale_fill_manual(values =cluster.color)+
            scale_x_discrete(drop = FALSE)+
            xlab("x")+
            ylab("cluster")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  #axis.line = element_line(colour = "black"),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(), 
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1)))
          
          p.corr1_pvalue<-ggplot(linear.regression.output.sub.plot1, aes(x=ID,y=1)) +
            geom_tile(aes(fill = log10.P),col = NA,size=0.1) + 
            scale_fill_gradient2(low = "white", high = "#e41a1c",na.value = "grey50")+
            geom_text(aes(label = is_sig_P10)) +
            scale_x_discrete(drop = FALSE)+
            xlab("x")+
            ylab("pvalue")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  #axis.line = element_line(colour = "black"),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(), 
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1)))
          p.corr1_exp<-ggplot(linear.regression.output.sub.plot1, aes(x=ID,y=1)) +
            geom_tile(aes(fill = X0h),col = NA,size=0.1) + 
            scale_fill_gradient2(low="white", high = "#377eb8",na.value = "grey50")+
            scale_x_discrete(drop = FALSE)+
            xlab("x")+
            ylab("brahch")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  #axis.line = element_line(colour = "black"),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(), 
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1)))
          
          # b. all motifs logP>=10
          # linear.regression.output.sub.plot2 = linear.regression.output.sub[as.character(linear.regression.output.sub$ID) %in% motifOrder.logP_10,]
          # linear.regression.output.sub.plot2$ID = factor(linear.regression.output.sub.plot2$ID,levels=motifOrder.logP_10)
          # if (nrow(linear.regression.output.sub.plot2)>0){
          #   p.corr2_branch<-ggplot(linear.regression.output.sub.plot2, aes(x=ID,y=1)) +
          #     geom_tile(aes(fill = branch_color),col = NA,size=0.1) + 
          #     scale_fill_manual(values =branch.color)+
          #     scale_x_discrete(drop = FALSE)+
          #     xlab("x")+
          #     ylab("brahch")+
          #     theme(panel.grid.major = element_blank(),
          #           panel.grid.minor = element_blank(),
          #           panel.background = element_blank(), 
          #           #axis.line = element_line(colour = "black"),
          #           axis.line.x = element_blank(),
          #           axis.line.y = element_blank(),
          #           axis.text.y=element_blank(), 
          #           axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
          #           axis.title=element_text(colour="black",size=rel(1)),
          #           legend.position="right",
          #           legend.background = element_blank(),
          #           legend.text=element_text(size=rel(1)))
          #   
          #   p.corr2_cluster<-ggplot(linear.regression.output.sub.plot2, aes(x=ID,y=1)) +
          #     geom_tile(aes(fill = cluster_color),col = NA,size=0.1) + 
          #     scale_fill_manual(values =cluster.color)+
          #     scale_x_discrete(drop = FALSE)+
          #     xlab("x")+
          #     ylab("cluster")+
          #     theme(panel.grid.major = element_blank(),
          #           panel.grid.minor = element_blank(),
          #           panel.background = element_blank(), 
          #           #axis.line = element_line(colour = "black"),
          #           axis.line.x = element_blank(),
          #           axis.line.y = element_blank(),
          #           axis.text.y=element_blank(), 
          #           axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
          #           axis.title=element_text(colour="black",size=rel(1)),
          #           legend.position="right",
          #           legend.background = element_blank(),
          #           legend.text=element_text(size=rel(1)))
          #   p.corr2_pvalue<-ggplot(linear.regression.output.sub.plot2, aes(x=ID,y=1)) +
          #     geom_tile(aes(fill = log10.P),col = NA,size=0.1) + 
          #     scale_fill_gradient2(low = "white", high = "#e41a1c",na.value = "grey50")+
          #     scale_x_discrete(drop = FALSE)+
          #     geom_text(aes(label = is_sig_P10)) +
          #     xlab("x")+
          #     ylab("pvalue")+
          #     theme(panel.grid.major = element_blank(),
          #           panel.grid.minor = element_blank(),
          #           panel.background = element_blank(), 
          #           #axis.line = element_line(colour = "black"),
          #           axis.line.x = element_blank(),
          #           axis.line.y = element_blank(),
          #           axis.text.y=element_blank(), 
          #           axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
          #           axis.title=element_text(colour="black",size=rel(1)),
          #           legend.position="right",
          #           legend.background = element_blank(),
          #           legend.text=element_text(size=rel(1)))
          #   
          #   p.corr2_exp<-ggplot(linear.regression.output.sub.plot2, aes(x=ID,y=1)) +
          #     geom_tile(aes(fill = X0h),col = NA,size=0.1) + 
          #     scale_fill_gradient2(low="white", high = "#377eb8",na.value = "grey50")+
          #     scale_x_discrete(drop = FALSE)+
          #     xlab("x")+
          #     ylab("exp")+
          #     theme(panel.grid.major = element_blank(),
          #           panel.grid.minor = element_blank(),
          #           panel.background = element_blank(), 
          #           #axis.line = element_line(colour = "black"),
          #           axis.line.x = element_blank(),
          #           axis.line.y = element_blank(),
          #           axis.text.y=element_blank(), 
          #           axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.9),
          #           axis.title=element_text(colour="black",size=rel(1)),
          #           legend.position="right",
          #           legend.background = element_blank(),
          #           legend.text=element_text(size=rel(1)))
          #   
          #   gA2 <- ggplotGrob(p.corr2_branch)
          #   gB2 <- ggplotGrob(p.corr2_cluster)
          #   gC2 <- ggplotGrob(p.corr2_pvalue)
          #   gD2 <- ggplotGrob(p.corr2_exp)
            # 
            # gA1 <- ggplotGrob(p.corr1_branch)
            # gB1 <- ggplotGrob(p.corr1_cluster)
            # gC1 <- ggplotGrob(p.corr1_pvalue)
            # gD1 <- ggplotGrob(p.corr1_exp)
            # g = rbind(gA1,gB1,gC1,gD1,gA2,gB2,gC2,gD2, size = "last")
            # pdf(paste("Step9_3_",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-corrplot.pdf",sep=""),                      
            #     width = 12,        # 5 x 300 pixels
            #     height = 16,
            #     pointsize = 10)        # smaller font size
            # grid.draw(g)
            # dev.off()
        }
#          } else {
            gA1 <- ggplotGrob(p.corr1_branch)
            gB1 <- ggplotGrob(p.corr1_cluster)
            gC1 <- ggplotGrob(p.corr1_pvalue)
            gD1 <- ggplotGrob(p.corr1_exp)
            g = rbind(gA1,gB1,gC1,gD1, size = "last")
            pdf(paste("Figure_S7-",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-corrplot.pdf",sep=""),
                width = 12,        # 5 x 300 pixels
                height = 8,
                pointsize = 5)        # smaller font size
            grid.draw(g)
            dev.off()
#          }
          
          ### 5. plots
          # if (grepl("logFC",Type)){
          #   p1<-ggplot(Annotation_plot.group.sub, aes(uniqueID_new,logFC)) +
          #     geom_bar(position="dodge", stat="identity",aes(fill=iPSC.alpha.Zscore.isActive.group))+
          #     scale_fill_manual(values= iPSC.alpha.Zscore.isActive.group.color)+
          #     scale_x_discrete(drop = FALSE)+
          #     geom_hline(yintercept = c(-2,-1,1,2), linetype="dotted",color = "black", size=.5)+
          #     ylab("Activity (Z-scaled alpha value)")+
          #     xlab(paste("N=",length(unique(Annotation_plot.group.sub$uniqueID_new))))+
          #     theme(panel.grid.major = element_blank(),
          #           panel.grid.minor = element_blank(),
          #           panel.background = element_blank(),
          #           axis.line = element_line(colour = "black"),
          #           axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          #           axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          #           #axis.text.y=element_text(colour="black",size=rel(0.7)),
          #           axis.text.y=element_text(colour="black",size=rel(1)),
          #           axis.text.x=element_blank(),
          #           axis.title=element_text(colour="black",size=rel(1)),
          #           legend.position="right",
          #           legend.background = element_blank(),
          #           legend.text=element_text(size=rel(1))) 
          # } else {
          #   p1<-ggplot(Annotation_plot.group.sub, aes(uniqueID_new,iPSC.alpha.Zscore)) +
          #     geom_bar(position="dodge", stat="identity",aes(fill=iPSC.alpha.Zscore.isActive.group))+
          #     scale_fill_manual(values= iPSC.alpha.Zscore.isActive.group.color)+
          #     scale_x_discrete(drop = FALSE)+
          #     geom_hline(yintercept = c(2,4,6,8), linetype="dotted",color = "black", size=.5)+
          #     ylab("Activity (Z-scaled alpha value)")+
          #     xlab(paste("N=",length(unique(Annotation_plot.group.sub$uniqueID_new))))+
          #     theme(panel.grid.major = element_blank(),
          #           panel.grid.minor = element_blank(),
          #           panel.background = element_blank(),
          #           axis.line = element_line(colour = "black"),
          #           axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          #           axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          #           #axis.text.y=element_text(colour="black",size=rel(0.7)),
          #           axis.text.y=element_text(colour="black",size=rel(1)),
          #           axis.text.x=element_blank(),
          #           axis.title=element_text(colour="black",size=rel(1)),
          #           legend.position="right",
          #           legend.background = element_blank(),
          #           legend.text=element_text(size=rel(1))) 
          # # }
          # Fimo.all.sum.sub$motif_name = factor(Fimo.all.sum.sub$motif_name,levels=rev(motifOrder.logP_5))
          # p2 <- ggplot(Fimo.all.sum.sub, aes(uniqueID.new,motif_name)) +
          #   geom_tile(aes(fill = n.group),col = NA) + 
          #   scale_fill_manual(values= color_motif)+
          #   scale_x_discrete(drop = FALSE)+
          #   xlab(Group2)+
          #   theme(panel.grid.major = element_blank(),
          #         panel.grid.minor = element_blank(),
          #         #panel.background = element_rect(fill = "#f0f0f0"),
          #         panel.background = element_blank(),
          #         axis.line = element_line(colour = "black"),
          #         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          #         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          #         axis.text.y=element_text(colour="black",size=rel(0.7)),
          #         axis.text.x=element_blank(),
          #         axis.title=element_text(colour="black",size=rel(1)),
          #         legend.position="right",
          #         legend.background = element_blank(),
          #         legend.text=element_text(size=rel(1))) 
          # if(length(motifOrder.logP_5.uniq)>0) {
          #   Fimo.all.sum.sub.tmp = Fimo.all.sum.sub[Fimo.all.sum.sub$motif_name %in% motifOrder.logP_5.uniq,]
          #   Fimo.all.sum.sub.tmp$motif_name = factor(Fimo.all.sum.sub.tmp$motif_name,levels=rev(motifOrder.logP_5.uniq))
          # } else {
          #   Fimo.all.sum.sub.tmp = Fimo.all.sum.sub
          # }
          # p2.uniq <- ggplot(Fimo.all.sum.sub.tmp, aes(uniqueID.new,motif_name)) +
          #   geom_tile(aes(fill = n.group),col = NA) + 
          #   scale_fill_manual(values= color_motif)+
          #   scale_x_discrete(drop = FALSE)+
          #   xlab(paste(Group2,"unique"))+
          #   theme(panel.grid.major = element_blank(),
          #         panel.grid.minor = element_blank(),
          #         #panel.background = element_rect(fill = "#f0f0f0"),
          #         panel.background = element_blank(),
          #         axis.line = element_line(colour = "black"),
          #         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          #         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          #         axis.text.y=element_text(colour="black",size=rel(0.7)),
          #         axis.text.x=element_blank(),
          #         axis.title=element_text(colour="black",size=rel(1)),
          #         legend.position="right",
          #         legend.background = element_blank(),
          #         legend.text=element_text(size=rel(1)))
          # rm(Fimo.all.sum.sub.tmp)
          # Fimo.all.sum.sub2 = Fimo.all.sum.sub[Fimo.all.sum.sub$motif_name %in% motifOrder.logP_10,]
          # if (nrow(Fimo.all.sum.sub2) > 0){
          #   Fimo.all.sum.sub2$motif_name = factor(Fimo.all.sum.sub2$motif_name,levels=rev(motifOrder.logP_10))
          #   p2_10<-ggplot(Fimo.all.sum.sub2, aes(uniqueID.new,motif_name)) +
          #     geom_tile(aes(fill = n.group),col = NA) + 
          #     scale_fill_manual(values= color_motif)+
          #     scale_x_discrete(drop = FALSE)+
          #     xlab(paste(Group,"-log10 padj >=10"))+
          #     theme(panel.grid.major = element_blank(),
          #           panel.grid.minor = element_blank(),
          #           #panel.background = element_rect(fill = "#f0f0f0"),
          #           panel.background = element_blank(),
          #           axis.line = element_line(colour = "black"),
          #           axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          #           axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          #           axis.text.y=element_text(colour="black",size=rel(0.7)),
          #           axis.text.x=element_blank(),
          #           axis.title=element_text(colour="black",size=rel(1)),
          #           legend.position="right",
          #           legend.background = element_blank(),
          #           legend.text=element_text(size=rel(1))) 
          #   p2_10.uniq<-ggplot(Fimo.all.sum.sub2[Fimo.all.sum.sub2$motif_name %in% linear.regression.output.sub.top.uniq$ID,], aes(uniqueID.new,motif_name)) +
          #     geom_tile(aes(fill = n.group),col = NA) + 
          #     scale_fill_manual(values= color_motif)+
          #     scale_x_discrete(drop = FALSE)+
          #     xlab(paste(Group,"-log10 padj >=10 unique"))+
          #     theme(panel.grid.major = element_blank(),
          #           panel.grid.minor = element_blank(),
          #           #panel.background = element_rect(fill = "#f0f0f0"),
          #           panel.background = element_blank(),
          #           axis.line = element_line(colour = "black"),
          #           axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          #           axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          #           axis.text.y=element_text(colour="black",size=rel(0.7)),
          #           axis.text.x=element_blank(),
          #           axis.title=element_text(colour="black",size=rel(1)),
          #           legend.position="right",
          #           legend.background = element_blank(),
          #           legend.text=element_text(size=rel(1))) 
          # }
          # gA <- ggplotGrob(p1)
          # gB <- ggplotGrob(p2)
          # gB_10.uniq <- ggplotGrob(p2_10.uniq)
          # if (nrow(Fimo.all.sum.sub2)>0){
          #   gB.uniq <- ggplotGrob(p2.uniq)
          #   gB_10 <- ggplotGrob(p2_10)
          #   g = rbind(gA,gB,gB_10,gB.uniq,gB_10.uniq, size = "last")
          #   Height = 16
          # } else {
          #   g = rbind(gA,gB,gB.uniq, size = "last")
          #   Height = 10
          # }
          # pdf(paste("Step9_3_",Family,"_",Species,"_",Group,"_",Type,"_motif_",Date,"-distribution.pdf",sep=""),    # create PNG for the heat map
          #     width = 16,        # 5 x 300 pixels
          #     height = Height,
          #     pointsize = 10)        # smaller font size
          # grid.draw(g)
          # dev.off()
      # }
        #####
    #}
    } ## each combination
  }  ## each species
}

#kept_motifs_MER11_F2.kept = kept_motifs_MER11_F2
#write.csv(kept_motifs_MER11_F2.kept,file=paste("motifs_kept_",Date,".csv",sep=""))

