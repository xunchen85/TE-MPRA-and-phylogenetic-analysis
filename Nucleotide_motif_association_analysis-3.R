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

##
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

##################################################
######################## step 4: nucleotide association plots
linear.regression.output.variant = linear.regression.output.kept[grepl("variant",linear.regression.output.kept$Type),]
linear.regression.output.variant$Shape = ifelse(linear.regression.output.variant$log10.P>=5,"sig","not_sig")
linear.regression.output.variant$Shape = ifelse(linear.regression.output.variant$log10.P>=10,"strong_sig",linear.regression.output.variant$Shape)
linear.regression.output.variant$Label = ifelse(linear.regression.output.variant$Shape %in% c("strong_sig"),linear.regression.output.variant$POS,NA)
linear.regression.output.variant.unique = linear.regression.output.variant[!duplicated(linear.regression.output.variant[,c(1,2,3,4)]),c(1,2,3,4)]
linear.regression.output.variant.unique$Group2 = ifelse(grepl("MER",linear.regression.output.variant.unique$Combination),"family","subfamily")
linear.regression.output.variant.unique$Group2 = ifelse(grepl("MER",linear.regression.output.variant.unique$Combination) & grepl("_g",linear.regression.output.variant.unique$Combination),"group",linear.regression.output.variant.unique$Group2)
linear.regression.output.variant.unique$Group2 = ifelse(linear.regression.output.variant.unique$Combination == "All","All",linear.regression.output.variant.unique$Group2)
shape_genotype2 = c("indel"=23,"variant"=21)
shape_genotype3 = c("indel"=18,"variant"=19)

Family = "MER11_combined_F2"
Families = c("MER11_combined_F1","MER11_combined_F2","MER34_combined_F1","MER52_combined_F1","MER52_combined_F2")

for(Family in Families){
  glist1 = list()
  Order1 = 1
  if (grepl("^MER11",Family)){
    Groupings = data.frame("Group2"=c("All",
                                      unique(linear.regression.output.variant.unique[linear.regression.output.variant.unique$Family== Family & 
                                                                                       linear.regression.output.variant.unique$Group2 == "family",]$Combination)))
    Groupings$Group1 = ifelse(Groupings$Group2 == "All","All",NA)
    Groupings$Group1 = ifelse(grepl("^MER",Groupings$Group2) & !grepl("_FG",Groupings$Group2),"family",Groupings$Group1)
    Groupings$Group1 = ifelse(grepl("^MER",Groupings$Group2) & grepl("_FG",Groupings$Group2),"group",Groupings$Group1)
    Groupings = Groupings[Groupings$Group1 != "family",]
  } else {
    Groupings = data.frame("Group2"=c("All"))
    Groupings$Group1 = ifelse(Groupings$Group2 == "All","All",NA)
  }
  # load the variant summary table
  variant.sum = read.csv(paste("input/Association_",Family,"_2023_9_8","_all.variant.sum.csv",sep=""))
  variant.sum.all = variant.sum[variant.sum$Group1 == "All",]
  #
  Group = "Combined"
  Type = "variantAll.both"
  Species = "Comparison"
  for(Group in Groupings$Group2){
    for (Species in c(unique(linear.regression.output.variant.unique$Species),"Comparison")){
      for (Type in c("variantAll.both")){
        if (Species == "Comparison"){
          Species2 = c("hg19","panTro4","macFas5")
        } else {
          Species2 = Species
        }
        if (Group == "Combined") {
          Group2 = Groupings[Groupings$Group1 == "family",]$Group2
        } else {
          Group2 = Group
        }
        # 2. linear regression output and map file
        if (grepl("both$",Type) & grepl("variantAll",Type)) {
          Type2 = gsub(".both$","",Type)
          linear.regression.output.variant.plot = linear.regression.output.variant[linear.regression.output.variant$Species %in% Species2 &
                                                                                     linear.regression.output.variant$Family == Family &
                                                                                     linear.regression.output.variant$Combination %in% Group2 &
                                                                                     linear.regression.output.variant$Type %in% c(Type2,gsub("variantAll","variantAll.indel",Type2)),]
        } else if (grepl("both$",Type) & grepl("variantEach",Type)) {
          Type2 = gsub(".both$","",Type)
          linear.regression.output.variant.plot = linear.regression.output.variant[linear.regression.output.variant$Species %in% Species2 &
                                                                                     linear.regression.output.variant$Family == Family &
                                                                                     linear.regression.output.variant$Combination %in% Group2 &
                                                                                     linear.regression.output.variant$Type %in% c(Type2,gsub("variantEach","variantEach.indel",Type2)),]
        }
        linear.regression.output.variant.plot$is_indel = ifelse(grepl("indel",linear.regression.output.variant.plot$Type),"indel","variant")
        linear.regression.output.variant.plot$colorGroup = paste(linear.regression.output.variant.plot$Species,linear.regression.output.variant.plot$Combination)
        # 3. plot
        p1 = ggplot(linear.regression.output.variant.plot, aes(x=POS, y=log10.P)) +
          ylab("-log10(P-value)")+
          xlab(paste("Position along the frame alignment ",nrow(variant.sum.all)," (bp)",sep=""))+
          geom_point(aes(fill = abs(BETA),shape = is_indel),size=3,color = "black")+
          scale_shape_manual(values=shape_genotype2)+
          xlim(0,nrow(variant.sum.all))+
          #scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 0)+
          scale_fill_gradient2(low = "white", high = "#e41a1c")+
          geom_hline(yintercept = c(5,10), linetype="dotted",color = "black", size=0.5)+
          geom_vline(xintercept = c(1,nrow(variant.sum.all)), linetype="dotted", 
                     color = "black", size=0.5)+
          geom_text_repel(aes(label = Label,color = Combination),max.overlaps = 100) +
          #scale_color_manual(values = color_sig)+
          #scale_shape_manual(values = shape_variant)+
          ggtitle(paste(Family,Group,Type,Species))+
          theme(
            plot.title = element_text(hjust = 0.5, size = rel(1)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.text=element_text(colour="black",size=rel(1),angle = 0),
            axis.text.x=element_text(colour="black",vjust=.5,hjust=0.5,angle = 0),
            axis.title=element_text(colour="black",size=rel(1)),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.position ="right",
            #legend.position = "none",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1)))
        glist1[[Order1]] <- ggplotGrob(p1)
        Order1 = Order1 + 1
        p1 = ggplot(linear.regression.output.variant.plot, aes(x=POS, y=log10.P)) +
          ylab("-log10(P-value)")+
          xlab(paste("Position along the frame alignment ",nrow(variant.sum.all)," (bp)",sep=""))+
          geom_point(aes(color = colorGroup,shape = is_indel),size=3,alpha=0.5)+
          scale_shape_manual(values=shape_genotype3)+
          xlim(0,nrow(variant.sum.all))+
          #scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 0)+
          geom_hline(yintercept = c(5,10), linetype="dotted",color = "black", size=0.5)+
          geom_vline(xintercept = c(1,nrow(variant.sum.all)), linetype="dotted", 
                     color = "black", size=0.5)+
          #geom_text_repel(aes(label = Label,color = Combination),max.overlaps = 100) +
          #scale_color_manual(values = color_sig)+
          #scale_shape_manual(values = shape_variant)+
          ggtitle(paste(Family,Group,Type))+
          theme(
            plot.title = element_text(hjust = 0.5, size = rel(1)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.text=element_text(colour="black",size=rel(1),angle = 0),
            axis.text.x=element_text(colour="black",vjust=.5,hjust=0.5,angle = 0),
            axis.title=element_text(colour="black",size=rel(1)),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.position ="right",
            #legend.position = "none",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1)))
        glist1[[Order1]] <- ggplotGrob(p1)
        Order1 = Order1 + 1
      }
    }
  }
  pdf(paste("Figure_5D_S8AB-",Family,"_",Date,".pdf",sep=""),    # create PNG for the heat map        
      width = 48,        # 5 x 300 pixels
      height = 3 * (length(glist1)/4),
      pointsize = 10 )        # smaller font size
  do.call("grid.arrange",c(glist1,ncol=4))
  dev.off()
}
