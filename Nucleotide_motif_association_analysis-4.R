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
##################################################
######################## step 5: motif plots
Family = "MER11_combined_F1"
Frame = "F2"
#Families = c("MER11","MER34","MER52")
topN = 10
topN2 = 5
Order.type = "subfamily"
motifType = "motif"
variantType = "variantAll"

##### motifs
MER11_F2.motifs = c("RARA::RXRG","INSM1","DMBX1","CRX","NR5A1","FOXF2","RFX6","ZNF331","ZNF136",
                    "ZIC3","TEAD4","TEAD1","TFCP2","GATA5","POU2F1::SOX2","SOX17",
                    "IKZF1","ZNF317")
# kept result
kept_motifs.final = read.csv(paste("input/motifs_kept_2023_7_6.csv",sep=""))
kept_motifs.final = kept_motifs.final[kept_motifs.final$Family == "All" & 
                                        kept_motifs.final$thredhold == "logP_10" &
                                        kept_motifs.final$Species == "All",]
kept_motifs.final$is_kept = NA
kept_motifs.final$is_kept = ifelse(kept_motifs.final$Frame == "MER11_combined_F1" & 
                                     kept_motifs.final$Motif %in% c("SOX9","MSGN1","PRDM1","ZNF354A","RXRG","POU5F1::SOX2","SP3","GLIS2",
                                                                    "ZFP335","NR5A1","ZIC2","MYF6"),"kept",NA)
kept_motifs.final$is_kept = ifelse(kept_motifs.final$Frame == "MER11_combined_F2" & 
                                     kept_motifs.final$Motif %in% c("RARA::RXRG","INSM1","DMBX1","CRX","NR5A1","FOXF2","RFX6","ZNF331","ZNF136",
                                                                    "ZIC3","TEAD4","TEAD1","TFCP2","GATA5","POU2F1::SOX2","SOX17",
                                                                    "IKZF1","ZNF317"),"kept",kept_motifs.final$is_kept)
kept_motifs.final$is_kept = ifelse(kept_motifs.final$Frame == "MER34_combined_F1" & 
                                     kept_motifs.final$Motif %in% c("POU3F3"),"kept",kept_motifs.final$is_kept)
kept_motifs.final$is_kept = ifelse(kept_motifs.final$Frame == "MER52_combined_F1" & 
                                     kept_motifs.final$Motif %in% c("SP2","E2F6","THRB","ZNF281"),"kept",kept_motifs.final$is_kept)
kept_motifs.final$is_kept = ifelse(kept_motifs.final$Frame == "MER52_combined_F2" & 
                                     kept_motifs.final$Motif %in% c("SP2"),"kept",kept_motifs.final$is_kept)
kept_motifs.final = kept_motifs.final[!is.na(kept_motifs.final$is_kept),]
# assign colors
kept_motifs.final$motifcolor = NA

kept_motifs.final$Motif = factor(kept_motifs.final$Motif)
motifColor = distinctColorPalette(length(unique(kept_motifs.final$Motif)))
names(motifColor) = levels(kept_motifs.final$Motif)
motifColor['gap'] = "white"
motifColor['nt'] = "black"

# motif enrichment analysis output
Family = "MER11_combined_F2"
Order.type = "subfamily"

list_MER11.macFas5.distance.tmp = list_MER11.macFas5.distance[,c("subfamily","Node1.","Order.rerooted.tree","FG")]
colnames(list_MER11.macFas5.distance.tmp) = c("subfamily","distance","order.tree","FG")
list_MER11_hg19_macFas5 = rbind(list_MER11.distance[,c("subfamily","distance","order.tree","FG")],list_MER11.macFas5.distance.tmp)

##
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

#
linear.regression.output.variant = linear.regression.output.kept[grepl("variant",linear.regression.output.kept$Type),]
linear.regression.output.variant$Shape = ifelse(linear.regression.output.variant$log10.P>=5,"sig","not_sig")
linear.regression.output.variant$Shape = ifelse(linear.regression.output.variant$log10.P>=10,"strong_sig",linear.regression.output.variant$Shape)
linear.regression.output.variant$Label = ifelse(linear.regression.output.variant$Shape %in% c("strong_sig"),linear.regression.output.variant$POS,NA)
linear.regression.output.variant.unique = linear.regression.output.variant[!duplicated(linear.regression.output.variant[,c(1,2,3,4)]),c(1,2,3,4)]
linear.regression.output.variant.unique$Group2 = ifelse(grepl("MER",linear.regression.output.variant.unique$Combination),"family","subfamily")
linear.regression.output.variant.unique$Group2 = ifelse(grepl("MER",linear.regression.output.variant.unique$Combination) & grepl("_g",linear.regression.output.variant.unique$Combination),"group",linear.regression.output.variant.unique$Group2)
linear.regression.output.variant.unique$Group2 = ifelse(linear.regression.output.variant.unique$Combination == "All","All",linear.regression.output.variant.unique$Group2)


## loops
for(Family in c("MER11_combined_F1","MER11_combined_F2")){
    Summary_Table2.plot = Summary_Table2.kept[!is.na(Summary_Table2.kept$Family.frame) & Summary_Table2.kept$Family.frame == Family,]
    Summary_Table2.plot = Summary_Table2.plot[!is.na(Summary_Table2.plot$iPSC.alpha.range),]   ## exclude positive and negative
    #  
    Summary_Table2.plot = Summary_Table2.plot[Summary_Table2.plot$Family.frame == Family & 
                                                !is.na(Summary_Table2.plot$Family.frame) & 
                                                !(Summary_Table2.plot$Family %in% "MER11D"),]    ## there is no MER11D frame was designed
    # exclude instances without subfamilies
    Summary_Table2.plot = Summary_Table2.plot[!is.na(Summary_Table2.plot$label.final.correctedName),]
    
      # order by the original annotation group
      if (grepl("MER11",Family)){
        Summary_Table2.plot.tmp = merge(Summary_Table2.plot,list_MER11_hg19_macFas5,by.x="label.final.correctedName",by.y="subfamily",all.x=T)
        for (nrow_tmp in nrow(list_MER11_hg19_macFas5):1){          #### the order of
          if (nrow_tmp == nrow(list_MER11_hg19_macFas5)){
            Summary_Table2.plot2 = Summary_Table2.plot.tmp[Summary_Table2.plot.tmp$label.final.correctedName == list_MER11_hg19_macFas5[nrow_tmp,]$subfamily,]
          } else {
            Summary_Table2.plot2 = rbind(Summary_Table2.plot2,Summary_Table2.plot.tmp[Summary_Table2.plot.tmp$label.final.correctedName == list_MER11_hg19_macFas5[nrow_tmp,]$subfamily,])
          }
        }
        Summary_Table2.plot = Summary_Table2.plot2
      } 
      # or order by the subfamily
      Summary_Table2.plot$uniqueID_new = factor(Summary_Table2.plot$uniqueID_new,levels=unique(Summary_Table2.plot$uniqueID_new))
    
    # each group
    #Groupings = data.frame("Group2"=c("All","Combined",unique(Summary_Table2.plot$Family)))
    Groupings = data.frame("Group2"=c("All"))
    Groupings$Group1 = ifelse(Groupings$Group2 == "All","All",NA)
    Groupings$Group1 = ifelse(Groupings$Group2 == "Combined","Combined",Groupings$Group1)
    Groupings$Group1 = ifelse(grepl("^MER",Groupings$Group2),"family",Groupings$Group1)
    Groupings$Group1 = ifelse(grepl("^[135][142]",Groupings$Group2),"subfamily",Groupings$Group1)
    Groupings = Groupings[Groupings$Group1 != "subfamily",]
    # load the frame tree and position of tips
    # load the motif info and position
    motif_file = read.delim(paste("input/",Family,".mafft.JASPAR.fimo.table.gz",sep=""),header=T,sep="")
    motif_file$Posi.gt99_0bp = as.numeric(as.character(motif_file$Posi.gt99_0bp))
    motif_file$Posi.gt99_1bp = motif_file$Posi.gt99_0bp+1
    motif_file$unique_posi = paste(motif_file$motifName,motif_file$uniqueID_new,motif_file$Posi.gt99_1bp)
    motif_file$nt = "m"
    motif_file$motifName = toupper(motif_file$motifName)
    
    # load the nt posi and gap info
    posi.info = read.delim(paste("input/",Family,".mafft.posi.gz",sep=""),sep="\t",header=T)
    posi.info = posi.info[!grepl("^0|^1|^2|^3|^4|^5|^6|^7|^8|^9",posi.info$uniqueID_new),]
    posi.info = posi.info[posi.info$Posi.gt99_0bp!="-",]
    posi.info$motifName = "nt"
    posi.info$Posi.gt99_1bp = as.numeric(as.character(posi.info$Posi.gt99_0bp))+1
    
    # variant summary table
    variant.sum = read.csv(file=paste("input/Association_",Family,"_2023_9_8_","All",".variant.sum.csv",sep=""))
    variant.sum.all = variant.sum[variant.sum$Group1 == "All",]
    variant.sum.all$pos.nt = paste(gsub("X","",variant.sum.all$nt.pos),variant.sum.all$majority)
    Group = "Combined"
    Species = "hg19"
      # 1. Annotation plot group
      for (Species in c("hg19","macFas5","panTro4")){
        #for (Species in c("hg19","combined")){
        # 1.1 achieve the summary table
        if (Species == "combined"){
          if (Group == "All" ){
            Summary_Table2.plot.sub = Summary_Table2.plot
          } else if (Group == "Combined") {
            Summary_Table2.plot.sub = Summary_Table2.plot[with(Summary_Table2.plot, order(Family, Order.final)),]
            if (grepl("MER11",Family)){        ### order by the MJN order
              Summary_Table2.plot.sub = Summary_Table2.plot
            }
          } else {
            Summary_Table2.plot.sub = Summary_Table2.plot[Summary_Table2.plot$Family == Group,]
          }
        } else {
          if (Group == "All"){
            head(Summary_Table2.plot)
            Summary_Table2.plot.sub = Summary_Table2.plot[Summary_Table2.plot$Species == Species,]
          } else if (Group == "Combined") {
            Summary_Table2.plot.sub = Summary_Table2.plot[Summary_Table2.plot$Species == Species,]
            Summary_Table2.plot.sub = Summary_Table2.plot.sub[with(Summary_Table2.plot.sub, order(Family, Order.final)),]
            if (grepl("MER11",Family)){        ### order by the MJN order
              Summary_Table2.plot.sub = Summary_Table2.plot[Summary_Table2.plot$Species == Species,]
            }
          } else {
            Summary_Table2.plot.sub = Summary_Table2.plot[Summary_Table2.plot$Family == Group & Summary_Table2.plot$Species == Species,]
          }
        }
        # exclude group2 column
        Summary_Table2.plot.sub = Summary_Table2.plot.sub[,!colnames(Summary_Table2.plot.sub) %in% "group2"]
        # only kept instances
        Summary_Table2.plot.sub = Summary_Table2.plot.sub[Summary_Table2.plot.sub$Group == "Instance",]    ## added on 2023/2/21, it does not impact the hg19 plots
        # order of instances
        Summary_Table2.plot.sub$uniqueID_new = factor(Summary_Table2.plot.sub$uniqueID_new,levels=unique(as.character(Summary_Table2.plot.sub$uniqueID_new)))
        
        # skip empty results
        if (nrow(Summary_Table2.plot.sub) == 0){
          next
        }
        # 1.2 select the candidate motifs
        if (Group == "Combined"){
          linear.regression.output.sub = linear.regression.output[linear.regression.output$Family == Family & 
                                                                    linear.regression.output$Combination %in% Groupings[Groupings$Group1=="family",]$Group2 &
                                                                    linear.regression.output$Type == motifType,]
        } else {
          linear.regression.output.sub = linear.regression.output[linear.regression.output$Family == Family & 
                                                                    linear.regression.output$Combination == Group &
                                                                    linear.regression.output$Type == motifType,]
        }
        
        ## plot candidate motifs updated 2023/7/12
        linear.regression.output.sub = linear.regression.output.sub[linear.regression.output.sub$ID %in% kept_motifs.final[kept_motifs.final$Frame == Family,]$Motif,]
        
        linear.regression.output.sub$ID = factor(linear.regression.output.sub$ID,levels=kept_motifs.final[kept_motifs.final$Frame == Family,]$Motif)
        linear.regression.output.sub.uniq = linear.regression.output.sub
        # kept all candidate motifs updated 2023/9/8
        #linear.regression.output.sub.uniq = linear.regression.output.sub[!duplicated(linear.regression.output.sub$branch_color),]
        
        # 1.3 linear regression output and candidate nucleotides
        if (Group == "Combined"){
          linear.regression.output.each = linear.regression.output.variant[linear.regression.output.variant$Family == Family &
                                                                             linear.regression.output.variant$Combination %in% Groupings[Groupings$Group1=="family",]$Group2 &
                                                                             linear.regression.output.variant$Type == variantType,]
        } else {
          linear.regression.output.each = linear.regression.output.variant[linear.regression.output.variant$Family == Family &
                                                                             linear.regression.output.variant$Combination == Group &
                                                                             linear.regression.output.variant$Type == variantType,]
        }
        
        linear.regression.output.each = linear.regression.output.each[linear.regression.output.each$Shape!="not_sig" & !is.na(linear.regression.output.each$Shape),]
        # associated nt position
        posi.info.nt = posi.info[posi.info$Posi.noGap_1bp %in% linear.regression.output.each$POS,c("motifName","uniqueID_new","Posi.gt99_1bp","nt")]  ## gt99 position
        posi.info.nt = posi.info.nt[!(paste(posi.info.nt$Posi.gt99_1bp,posi.info.nt$nt) %in% variant.sum.all$pos.nt),]
        
        # gap nt position
        posi.info.gap = posi.info[posi.info$Posi.noGap_1bp=="-",c("motifName","uniqueID_new","Posi.gt99_1bp","nt")]  ## gt99 position
        posi.info.gap$motifName = "gap"
        posi.info.both = rbind(posi.info.nt,posi.info.gap)
        enrichmentType = "association"
        
          motif_file.plot=motif_file[motif_file$motifName %in% as.character(linear.regression.output.sub.uniq$ID),]
          # 2. plot
          # select kept frames
          posi.info.plot = posi.info.both[posi.info.both$uniqueID_new %in% as.character(Summary_Table2.plot.sub$uniqueID_new),]
          motif_file.plot = motif_file.plot[motif_file.plot$uniqueID_new %in% as.character(Summary_Table2.plot.sub$uniqueID_new),]
          
          # combine with the nt table
          motif_file.plot = motif_file.plot[!duplicated(motif_file.plot$unique_posi),c("motifName","uniqueID_new","Posi.gt99_1bp","nt")]
          motif.nt.plot = rbind(motif_file.plot,posi.info.plot)
          
          # motif color
          motif.nt.plot$motifName = factor(motif.nt.plot$motifName)

          # frame re-order (only kept orthologous frames)
          motif.nt.plot$uniqueID_new = factor(motif.nt.plot$uniqueID_new,levels=unique(as.character(Summary_Table2.plot.sub$uniqueID_new)))
          
          ##################### plots
          my_pal <- colorRampPalette(rev(brewer.pal(n = 10, name = "Set3")))
          # motif plot 
          p1_motif <-ggplot(motif.nt.plot, aes(Posi.gt99_1bp,uniqueID_new)) +
            geom_tile(aes(fill = motifName),alpha=0.8,col = NA,size=0.1) + 
            scale_fill_manual(values= motifColor)+
            scale_y_discrete(drop = FALSE)+
            geom_vline(xintercept = c(0,(nrow(variant.sum.all)+1)), linetype="dotted",color = "black", size=.5)+
            xlab(paste("frame alignment",nrow(variant.sum.all),"bp"))+
            ylab(paste("N=",length(levels(motif.nt.plot$uniqueID_new))))+
            ggtitle(paste(Family,Group,enrichmentType))+
            theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "#f0f0f0"),
                  axis.line = element_line(colour = "black"),
                  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1))) 
          p2_motif <-ggplot(motif.nt.plot[motif.nt.plot$motifName!="nt",], aes(Posi.gt99_1bp,uniqueID_new)) +
            geom_tile(aes(fill = motifName),alpha=0.8,col = NA,size=0.1) + 
            scale_fill_manual(values= motifColor)+
            scale_y_discrete(drop = FALSE)+
            geom_vline(xintercept = c(0,(nrow(variant.sum.all)+1)), linetype="dotted",color = "black", size=.5)+
            xlab(paste("frame alignment",nrow(variant.sum.all),"bp"))+
            ylab(paste("N=",length(levels(motif.nt.plot$uniqueID_new))))+
            ggtitle(paste(Family,Group,enrichmentType))+
            theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "#f0f0f0"),
                  axis.line = element_line(colour = "black"),
                  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1))) 
          
          gA <- ggplotGrob(p1_motif)
          gB <- ggplotGrob(p2_motif)
          g = rbind(gA,gB, size = "last")
          
          pdf(paste("Figure_S8C-",Species,"_",Family,"_",Group,"_",Order.type,"_",enrichmentType,"_",Date,"_motif.pdf",sep=""),    # create PNG for the heat map
              width = 24,        # 5 x 300 pixels
              height = 24,
              pointsize = 10)        # smaller font size
          grid.draw(g)
          dev.off()
        # plot the alpha distribution
        p2_alpha<-ggplot(Summary_Table2.plot.sub, aes(log2.iPSC.alpha,uniqueID_new)) +
          geom_col(aes(fill = iPSC.alpha.Zscore.isActive.group)) + 
          scale_fill_manual(values =iPSC.alpha.Zscore.isActive.group.color)+
          scale_y_discrete(drop = FALSE)+
          #geom_vline(xintercept = c(1,2), linetype="dotted",color = "grey44", size=1.5)+
          xlab("x")+
          ylab(paste("N=",length(levels(Summary_Table2.plot.sub$uniqueID_new))))+
          theme_void()+
          ggtitle(paste(Family,Group,"iPSC.alpha"))+
          theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.text.y=element_blank(), 
                axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                axis.title=element_text(colour="black",size=rel(1)),
                axis.title.y = element_text(angle=90),
                legend.position="right",
                legend.background = element_blank(),
                legend.text=element_text(size=rel(1)))
        p2_alpha.Zcore<-ggplot(Summary_Table2.plot.sub, aes(iPSC.alpha.Zscore,uniqueID_new)) +
          geom_col(aes(fill = iPSC.alpha.Zscore.isActive.group)) + 
          scale_fill_manual(values =iPSC.alpha.Zscore.isActive.group.color)+
          scale_y_discrete(drop = FALSE)+
          #geom_vline(xintercept = c(1,2), linetype="dotted",color = "grey44", size=1.5)+
          xlab("x")+
          ylab(paste("N=",length(levels(Summary_Table2.plot.sub$uniqueID_new))))+
          theme_void()+
          ggtitle(paste(Family,Group,"iPSC.alpha.Zscore"))+
          theme(plot.title = element_text(hjust = 0.5, size = rel(1)),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                axis.text.y=element_blank(), 
                axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                axis.title=element_text(colour="black",size=rel(1)),
                axis.title.y = element_text(angle=90),
                legend.position="right",
                legend.background = element_blank(),
                legend.text=element_text(size=rel(1)))
        gA <- ggplotGrob(p2_alpha)
        gB <- ggplotGrob(p2_alpha.Zcore)
        g = cbind(gA,gB, size = "last")
        
        pdf(paste("Figure_S8C-",Species,"_",Family,"_",Group,"_",Order.type,"_",Date,"_alpha.pdf",sep=""),    # create PNG for the heat map
            width = 12,        # 5 x 300 pixels
            height = 12,
            pointsize = 10)        # smaller font size
        grid.draw( g)
        dev.off()
        
        # plot the grouping info
        Summary_Table2.plot.color = Summary_Table2.plot.sub[!duplicated(Summary_Table2.plot.sub$label.final.correctedName),]
        col_family = colorRampPalette(brewer.pal(8, "Set3"))(n = length(unique(Summary_Table2.plot.sub$Family)))
        # col_family = distinctColorPalette(length(unique(Summary_Table2.plot$Family)))
        names(col_family) = unique(Summary_Table2.plot.sub$Family)
        cluster.color = Summary_Table2.plot.color$cluster.color
        names(cluster.color) = Summary_Table2.plot.color$label.final.correctedName
        
        p3_family<-ggplot(Summary_Table2.plot.sub, aes(x=1,y=uniqueID_new)) +
          geom_tile(aes(fill = TEfamily.x),col = NA,size=0.1) + 
          #scale_fill_manual(values =col_family)+
          scale_y_discrete(drop = FALSE)+
          xlab("x")+
          ylab(paste("N=",length(levels(Summary_Table2.plot.sub$uniqueID_new))))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                #axis.line = element_line(colour = "black"),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                axis.text.y=element_blank(), 
                axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                axis.title=element_text(colour="black",size=rel(1)),
                legend.position="right",
                legend.background = element_blank(),
                legend.text=element_text(size=rel(1)))
        
        p4_subfamily<-ggplot(Summary_Table2.plot.sub, aes(x=1,y=uniqueID_new)) +
          geom_tile(aes(fill = label.final.correctedName),col = NA,size=0.1) + 
          scale_fill_manual(values = cluster.color)+
          scale_y_discrete(drop = FALSE)+
          xlab("x")+
          ylab(paste("N=",length(levels(Summary_Table2.plot.sub$uniqueID_new))))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                #axis.line = element_line(colour = "black"),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                axis.text.y=element_blank(), 
                axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                axis.title=element_text(colour="black",size=rel(1)),
                legend.position="right",
                legend.background = element_blank(),
                legend.text=element_text(size=rel(1))) 
        
        if (grepl("MER11",Family)){
          if (Order.type == "subfamily"){
            Summary_Table2.plot.sub.dist = Summary_Table2.plot.sub
          } else {
            Summary_Table2.plot.sub.dist = merge(Summary_Table2.plot.sub,list_MER11_hg19_macFas5,by.x="label.final.correctedName",by.y="subfamily",all.x=T)
            Summary_Table2.plot.sub.dist$uniqueID_new  = factor(Summary_Table2.plot.sub.dist$uniqueID_new,levels = unique(as.character(Summary_Table2.plot.sub$uniqueID_new)))
          }
          p7_distance<-ggplot(Summary_Table2.plot.sub.dist, aes(x=1,y=uniqueID_new)) +
            geom_tile(aes(fill = distance),col = NA,size=0.1) + 
            scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",
                                 midpoint = median(list_MER11_hg19_macFas5$distance,na.rm=TRUE),
                                 limits = c(min(list_MER11_hg19_macFas5$distance,na.rm=TRUE),max(list_MER11_hg19_macFas5$distance,na.rm=TRUE)))+
            scale_y_discrete(drop = FALSE)+
            xlab("x")+
            ylab(paste("N=",length(levels(Summary_Table2.plot.sub.dist$uniqueID_new))))+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  #axis.line = element_line(colour = "black"),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.y=element_blank(), 
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1))) 
        } else {
          p7_distance = p4_subfamily
        }
        
        p5_activity<-ggplot(Summary_Table2.plot.sub, aes(x=1,y=uniqueID_new)) +
          geom_tile(aes(fill = iPSC.alpha.range),col = NA,size=0.1) + 
          scale_fill_manual(values =color_alpha.group)+
          scale_y_discrete(drop = FALSE)+
          xlab("x")+
          ylab(paste("N=",length(levels(Summary_Table2.plot.sub$uniqueID_new))))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                #axis.line = element_line(colour = "black"),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                axis.text.y=element_blank(), 
                axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                axis.title=element_text(colour="black",size=rel(1)),
                legend.position="none",
                legend.background = element_blank(),
                legend.text=element_text(size=rel(1))) 
        
        p6_div<-ggplot(Summary_Table2.plot.sub, aes(x=1,y=uniqueID_new)) +
          geom_tile(aes(fill = div_rate),col = NA,size=0.1) + 
          scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",
                               midpoint = median(Summary_Table2.plot$div_rate,na.rm=TRUE),
                               limits = c(min(Summary_Table2.plot$div_rate,na.rm=TRUE),max(Summary_Table2.plot$div_rate,na.rm=TRUE)))+
          scale_y_discrete(drop = FALSE)+
          xlab("x")+
          ylab(paste("N=",length(levels(Summary_Table2.plot.sub$uniqueID_new))))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                #axis.line = element_line(colour = "black"),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                axis.text.y=element_blank(), 
                axis.text.x=element_text(colour="black",size=rel(1),angle = 0,vjust=0.5, hjust=0.5),
                axis.title=element_text(colour="black",size=rel(1)),
                legend.position="right",
                legend.background = element_blank(),
                legend.text=element_text(size=rel(1))) 
        
        g3 <- ggplotGrob(p3_family)
        g4 <- ggplotGrob(p4_subfamily)
        g5 <- ggplotGrob(p5_activity)
        g6 <- ggplotGrob(p6_div)
        g7 <- ggplotGrob(p7_distance)
        
        # Combine the plots
        g = cbind(g3,g4,g5,g6,g7, size = "first")
        pdf(paste("Figure_S8C-",Species,"_",Family,"_",Group,"_",Order.type,"_",Date,"_family.pdf",sep=""),    # create PNG for the heat map
            width = 16,        # 5 x 300 pixels
            height = 12,
            pointsize = 10)        # smaller font size
        grid.draw(g)
        dev.off()
        
        ## plot the line distribution of proportion of motifs per subfamily
        if (nrow(linear.regression.output.sub) == 0) {
          next
        }
        
        Fimo.all.sum.plot.sub = Fimo.all.sum[Fimo.all.sum$uniqueID.new %in% as.character(Summary_Table2.plot.sub$uniqueID_new),]
        Fimo.all.sum.plot.sub$motif_name = toupper(Fimo.all.sum.plot.sub$motif_name)
        Fimo.all.sum.plot.sub = Fimo.all.sum.plot.sub[!duplicated(Fimo.all.sum.plot.sub[,c("motif_name","uniqueID.new")]),]
        Fimo.all.sum.plot.sub = Fimo.all.sum.plot.sub[Fimo.all.sum.plot.sub$motif_name %in% as.character(linear.regression.output.sub$ID),]
        Fimo.all.sum.plot.sub = merge(Fimo.all.sum.plot.sub,Summary_Table2.plot.sub[,c("uniqueID_new","species","cluster.final","label.final.correctedName","div_rate","Instance_len")],by.x="uniqueID.new",by.y="uniqueID_new",all.x=T)
        Fimo.all.sum.plot.sub.sum = data.frame(Fimo.all.sum.plot.sub %>% 
                                                 group_by(species,Family.frame,DB,label.final.correctedName,motif_name) %>% 
                                                 dplyr::summarise(count.motif = n()))
        Fimo.all.sum.plot.sub.sum$label.final.correctedName.motif = paste(Fimo.all.sum.plot.sub.sum$species,Fimo.all.sum.plot.sub.sum$label.final.correctedName,Fimo.all.sum.plot.sub.sum$motif_name)
        Summary_Table2.plot.sub.sum = data.frame(Summary_Table2.plot.sub %>% 
                                                   group_by(species,label.final.correctedName) %>% 
                                                   dplyr::summarise(count.instance = n(),mean.div = mean(div_rate,na.rm = TRUE),median.div = median(div_rate,na.rm = TRUE)))
        Summary_Table2.plot.sub.sum.all = data.frame("tmp1"=NA,"tmp2"=NA)
        for (motifID in unique(as.character(linear.regression.output.sub$ID))){
          Summary_Table2.plot.sub.sum.tmp = Summary_Table2.plot.sub.sum
          Summary_Table2.plot.sub.sum.tmp$motif = motifID
          if (nrow(Summary_Table2.plot.sub.sum.all) == 1){
            Summary_Table2.plot.sub.sum.all = Summary_Table2.plot.sub.sum.tmp
          } else {
            Summary_Table2.plot.sub.sum.all = rbind(Summary_Table2.plot.sub.sum.all,Summary_Table2.plot.sub.sum.tmp)
          }
        }
        Summary_Table2.plot.sub.sum.all$label.final.correctedName.motif = paste(Summary_Table2.plot.sub.sum.all$species,Summary_Table2.plot.sub.sum.all$label.final.correctedName,Summary_Table2.plot.sub.sum.all$motif)
        Summary_Table2.plot.sub.sum.all = merge(Summary_Table2.plot.sub.sum.all,Fimo.all.sum.plot.sub.sum,by="label.final.correctedName.motif",all.x=T)
        
        # keep subfamilies >=10
        Summary_Table2.plot.sub.sum.all.10 = Summary_Table2.plot.sub.sum.all[Summary_Table2.plot.sub.sum.all$count.instance>=10,]
        Summary_Table2.plot.sub.sum.all.10$perC = ifelse(is.na(Summary_Table2.plot.sub.sum.all.10$count.motif),0,Summary_Table2.plot.sub.sum.all.10$count.motif/Summary_Table2.plot.sub.sum.all.10$count.instance)
        # motif color
        linear.regression.output.sub.color = linear.regression.output.sub[!duplicated(paste(linear.regression.output.sub$branch_color,linear.regression.output.sub$ID)),]
        motifColor2 = linear.regression.output.sub.color$branch_color
        names(motifColor2) = linear.regression.output.sub.color$ID
        Summary_Table2.plot.sub.sum.all.10 = merge(Summary_Table2.plot.sub.sum.all.10,linear.regression.output.sub.color[,c("ID","branch_color")],by.x="motif",by.y="ID",all.x=T)
        Summary_Table2.plot.sub.sum.all.10$Group = paste(Summary_Table2.plot.sub.sum.all.10$species.x,Summary_Table2.plot.sub.sum.all.10$motif)
        Summary_Table2.plot.sub.sum.all.10.convert = data.frame("tmp1"=NA,"tmp2"=NA)
        for(Motif in unique(Summary_Table2.plot.sub.sum.all.10$motif)){
          Summary_Table2.plot.sub.sum.all.10.convert.tmp = Summary_Table2.plot.sub.sum.all.10[Summary_Table2.plot.sub.sum.all.10$motif == Motif,c("label.final.correctedName.x","perC")]
          colnames(Summary_Table2.plot.sub.sum.all.10.convert.tmp)[2] = Motif
          if (colnames(Summary_Table2.plot.sub.sum.all.10.convert)[1] == "tmp1"){
            Summary_Table2.plot.sub.sum.all.10.convert = Summary_Table2.plot.sub.sum.all.10.convert.tmp
          } else {
            Summary_Table2.plot.sub.sum.all.10.convert = merge(Summary_Table2.plot.sub.sum.all.10.convert,Summary_Table2.plot.sub.sum.all.10.convert.tmp,by="label.final.correctedName.x",all=T)
          }
        }
        rownames(Summary_Table2.plot.sub.sum.all.10.convert) = Summary_Table2.plot.sub.sum.all.10.convert$label.final.correctedName.x
        Summary_Table2.plot.sub.sum.all.10.convert = as.matrix(t(Summary_Table2.plot.sub.sum.all.10.convert[,-1]))
        Summary_Table2.plot.sub.sum.all.10.convert[is.na(Summary_Table2.plot.sub.sum.all.10.convert)] <-0
        
        if (nrow(Summary_Table2.plot.sub.sum.all.10.convert) == 1) {
          motif.order.cluster = unique(as.character(Summary_Table2.plot.sub.sum.all.10$motif))
        } else {
          hclust.dist = hclust(dist(Summary_Table2.plot.sub.sum.all.10.convert))
          motif.order.cluster = hclust.dist$labels[hclust.dist$order]
        }
        if (Family == "MER11_combined_F2"){
          motif.order.cluster = rev(MER11_F2.motifs)
        }
        Summary_Table2.plot.sub.sum.all.10$motif = factor(Summary_Table2.plot.sub.sum.all.10$motif,levels=motif.order.cluster)
        
        ## plot
        if (grepl("MER11",Family) & Species == "macFas5"){
          Summary_Table2.plot.sub.sum.all.10$label.final.correctedName.x = factor(Summary_Table2.plot.sub.sum.all.10$label.final.correctedName.x,levels=list_MER11.macFas5.distance$subfamily)
        } else if (grepl("MER11",Family) & Species %in% c("hg19","panTro4")) {
          Summary_Table2.plot.sub.sum.all.10$label.final.correctedName.x = factor(Summary_Table2.plot.sub.sum.all.10$label.final.correctedName.x,levels=list_MER11$subfamily)
        }
        if (nrow(Summary_Table2.plot.sub.sum.all.10)>=1) {
          p1 = ggplot(Summary_Table2.plot.sub.sum.all.10, aes(x=label.final.correctedName.x,y=motif)) +
            geom_tile(aes(fill=perC*100)) + 
            scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
            xlab("divergent rate")+
            ylab("% instances containing the motif")+
            scale_x_discrete(drop = FALSE)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "#f0f0f0"),
                  axis.line = element_line(colour = "black"),
                  #axis.line.x = element_blank(),
                  #axis.line.y = element_blank(),
                  axis.text.y=element_text(colour="black",size=rel(1)),
                  axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=0.5, hjust=0.5),
                  axis.title=element_text(colour="black",size=rel(1)),
                  legend.position="right",
                  legend.background = element_blank(),
                  legend.text=element_text(size=rel(1))) 
          pdf(paste("Figure_3F-",Species,"_",Family,"_",Group,"_",Order.type,"_",Date,"_motif-1.pdf",sep=""),    # create PNG for the heat map
              width = 12,        # 5 x 300 pixels
              height = 6,
              pointsize = 10)        # smaller font size
          grid.draw(p1)
          dev.off()
        }
  }
}

