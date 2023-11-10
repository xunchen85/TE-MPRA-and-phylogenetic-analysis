##### load pacakges
library(splitstackshape)
library(gplots)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid)
library(gridExtra)
library(pheatmap)
library(randomcoloR)
library(RColorBrewer)
library("rGREAT")

#Date = "2023_9_7"
############################ Enviroment
setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/Project_Neurogenesis/Final_edited_version_2022_11_25/Final_scripts/")

######################### read annotation
Summary_Table1 = read.csv("input/Summary_Table1_2022_8_9.csv")

Summary_Table2 = read.csv("input/Summary_Table2_2022_8_9.csv")

Summary_Table1$Instance_coordinate_species_1bp = paste(Summary_Table1$species,":",Summary_Table1$chr,":",Summary_Table1$start+1,"-",Summary_Table1$end,sep="")
Summary_Table1$Instance_coordinate_species_1bp = ifelse(Summary_Table1$species == "Ancient",Summary_Table1$Instance_coordinate_species,Summary_Table1$Instance_coordinate_species_1bp)

# MPRA activity 
Summary_Table2 = data.frame(Summary_Table2[,c("Instance_coordinate","uniqueID_new","Family","uniqueID_MPRA","Group","Frame","Species","iPSC.activity.mean","iPSC.alpha.Zscore","iPSC.alpha.Zscore.isActive","NPC.activity.mean","NPC.alpha.Zscore","NPC.alpha.Zscore.isActive","iPSC.alpha","NPC.alpha","logFC","is_kept")])
Summary_Table2$Instance_coordinate_1bp = ifelse(is.na(Summary_Table2$Instance_coordinate) & Summary_Table2$Species == "Ancient",Summary_Table2$Family,Summary_Table2$Instance_coordinate)
Summary_Table2$Instance_coordinate_species_1bp = paste(Summary_Table2$Species,Summary_Table2$Instance_coordinate_1bp,sep=":")
Summary_Table2$log2.iPSC.alpha = log2(Summary_Table2$iPSC.alpha)

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


# 
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha >=5 & !is.na(Summary_Table2$log2.iPSC.alpha),"5",NA)
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha >=4 & Summary_Table2$log2.iPSC.alpha < 5 &
                                                !is.na(Summary_Table2$log2.iPSC.alpha) &
                                                is.na(Summary_Table2$log2.iPSC.alpha.range),"4",Summary_Table2$log2.iPSC.alpha.range)
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha >=3 & Summary_Table2$log2.iPSC.alpha < 4 &
                                                !is.na(Summary_Table2$log2.iPSC.alpha) &
                                                is.na(Summary_Table2$log2.iPSC.alpha.range),"3",Summary_Table2$log2.iPSC.alpha.range)
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha >=2 & Summary_Table2$log2.iPSC.alpha < 3 &
                                                !is.na(Summary_Table2$log2.iPSC.alpha) &
                                                is.na(Summary_Table2$log2.iPSC.alpha.range),"2",Summary_Table2$log2.iPSC.alpha.range)
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha >=1 & Summary_Table2$log2.iPSC.alpha < 2 &
                                                !is.na(Summary_Table2$log2.iPSC.alpha) &
                                                is.na(Summary_Table2$log2.iPSC.alpha.range),"1",Summary_Table2$log2.iPSC.alpha.range)
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha >= 0 & Summary_Table2$log2.iPSC.alpha < 1 &
                                                !is.na(Summary_Table2$log2.iPSC.alpha) &
                                                is.na(Summary_Table2$log2.iPSC.alpha.range),"0",Summary_Table2$log2.iPSC.alpha.range)
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha >= -1 & Summary_Table2$log2.iPSC.alpha < 0 &
                                                !is.na(Summary_Table2$log2.iPSC.alpha) &
                                                is.na(Summary_Table2$log2.iPSC.alpha.range),"-1",Summary_Table2$log2.iPSC.alpha.range)
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha >= -2 & Summary_Table2$log2.iPSC.alpha < -1 &
                                                !is.na(Summary_Table2$log2.iPSC.alpha) &
                                                is.na(Summary_Table2$log2.iPSC.alpha.range),"-2",Summary_Table2$log2.iPSC.alpha.range)
Summary_Table2$log2.iPSC.alpha.range = ifelse(Summary_Table2$log2.iPSC.alpha < -2 &
                                                !is.na(Summary_Table2$log2.iPSC.alpha) &
                                                is.na(Summary_Table2$log2.iPSC.alpha.range),"-3",Summary_Table2$log2.iPSC.alpha.range)

Summary_Table2$iPSC.alpha.range = ifelse(Summary_Table2$iPSC.alpha >=4 & !is.na(Summary_Table2$log2.iPSC.alpha),"alpha>=4","None")
Summary_Table2$iPSC.alpha.range = ifelse(Summary_Table2$iPSC.alpha >=2 & Summary_Table2$iPSC.alpha < 4 &
                                           !is.na(Summary_Table2$iPSC.alpha) ,"alpha>=2",Summary_Table2$iPSC.alpha.range)
Summary_Table2$iPSC.alpha.range = ifelse(is.na(Summary_Table2$iPSC.alpha),"Missing",Summary_Table2$iPSC.alpha.range)

activity_range = c("5" = "#5e4fa2","4"="#3288bd","3"="#66c2a5","2"="#e6f598","1"="#ffffbf","0"="#f5f5f5","-1"="#f5f5f5","-2"="#f5f5f5","-3"="#f5f5f5")
color_alpha.group = c("alpha>=4"="#bd0026","alpha>=2"="#fd8d3c","None"="#f0f0f0","Missing"="#ffffff")

# each frame
Summary_Table2_F1 = Summary_Table2[Summary_Table2$Frame == "F1",c("Instance_coordinate_species_1bp","uniqueID_new","iPSC.alpha","iPSC.alpha.range","log2.iPSC.alpha","log2.iPSC.alpha.range","iPSC.alpha.Zscore","iPSC.alpha.Zscore.isActive.group","is_kept")]
Summary_Table2_F2 = Summary_Table2[Summary_Table2$Frame == "F2",c("Instance_coordinate_species_1bp","uniqueID_new","iPSC.alpha","iPSC.alpha.range","log2.iPSC.alpha","log2.iPSC.alpha.range","iPSC.alpha.Zscore","iPSC.alpha.Zscore.isActive.group","is_kept")]

colnames(Summary_Table2_F1) = c("Instance_coordinate_species_1bp","uniqueID_new.F1","iPSC.alpha.F1","iPSC.alpha.range.F1","log2.iPSC.alpha.F1","log2.iPSC.alpha.range.F1","iPSC.alpha.Zscore.F1","iPSC.alpha.Zscore.isActive.group.F1","is_kept.F1")
colnames(Summary_Table2_F2) = c("Instance_coordinate_species_1bp","uniqueID_new.F2","iPSC.alpha.F2","iPSC.alpha.range.F2","log2.iPSC.alpha.F2","log2.iPSC.alpha.range.F2","iPSC.alpha.Zscore.F2","iPSC.alpha.Zscore.isActive.group.F2","is_kept.F2")

# merge with annotation_plot
Summary_Table1 = merge(Summary_Table1,Summary_Table2_F1,by = "Instance_coordinate_species_1bp",all.x=T)
Summary_Table1 = merge(Summary_Table1,Summary_Table2_F2,by = "Instance_coordinate_species_1bp",all.x=T)

############################ step1.1 load nonshuffle ATAC and K27ac
# obtain H1 metadata
metadata.H1 = read.delim("input/metadata_H1.tsv",sep="\t",header=T)
metadata.H1 = metadata.H1[metadata.H1$File.type == "bed",c("File.accession","Experiment.target")]
metadata.H1$cell = "H1"
metadata.NPC = read.delim("input/metadata_NPC.tsv",sep="\t",header=T)
metadata.NPC = metadata.NPC[metadata.NPC$File.type == "bed",c("File.accession","Experiment.target")]
metadata.NPC$cell = "NPC"
metadata = rbind(metadata.H1,metadata.NPC)
rm(metadata.H1,metadata.NPC)
metadata$mark.type = ifelse(grepl("H[1234]K|H2[AB]K",metadata$Experiment.target),"histone","TF")
metadata$Experiment.target = gsub("-human","",metadata$Experiment.target)

# obtain the total number of peaks 
TF.totalPeaks = read.delim("input/hg19_Chipseq_total_peaks.sum_2022_12_28",sep="",header=F)
colnames(TF.totalPeaks) = c("fileID","total.peaks")
TF.totalPeaks$fileID = gsub("_hg19","",TF.totalPeaks$fileID)
TF.totalPeaks = merge(TF.totalPeaks,metadata,by.x="fileID",by.y="File.accession",all.x=T)
TF.totalPeaks$cell = ifelse(is.na(TF.totalPeaks$cell),"293T",TF.totalPeaks$cell)
TF.totalPeaks$Experiment.target = ifelse(is.na(TF.totalPeaks$Experiment.target),TF.totalPeaks$fileID,TF.totalPeaks$Experiment.target)
TF.totalPeaks$mark.type = ifelse(is.na(TF.totalPeaks$mark.type),"KZNF",TF.totalPeaks$mark.type)
TF.totalPeaks$Experiment.target = gsub("ZNF33A-rep2","ZNF33A.2",TF.totalPeaks$Experiment.target)
unique(TF.totalPeaks$Experiment.target)
# filter duplicated set with a lower counts 
TF.totalPeaks = TF.totalPeaks[order(-TF.totalPeaks$total.peaks),]
TF.totalPeaks = TF.totalPeaks[!duplicated(TF.totalPeaks[,c("Experiment.target","cell")]),]
# exclude NPC cells for now
TF.totalPeaks = TF.totalPeaks[TF.totalPeaks$cell != "NPC",]

##### number of TFs and chip-seq analyzed in H1 cells
nrow(TF.totalPeaks[TF.totalPeaks$mark.type !="KZNF" & TF.totalPeaks$mark.type == "TF",])
nrow(TF.totalPeaks[TF.totalPeaks$mark.type !="KZNF" & TF.totalPeaks$mark.type != "TF",])

# hg19 counts per family
TE_rmsk_0bp = read.delim("input/hg19_rmsk_TE_0bp_MER11_34_52.rename.bed",header=F,sep="")
TE_rmsk_0bp_200bp = TE_rmsk_0bp[TE_rmsk_0bp$V3-TE_rmsk_0bp$V2>=200,]
colnames(Summary_Table1)
TE_rmsk_0bp_200bp = merge(TE_rmsk_0bp_200bp,Summary_Table1[,c("instanceID.renamed","TEfamily")],by.x="V4",by.y="instanceID.renamed",all.x=T)
TE_rmsk_0bp_200bp.sum = data.frame(TE_rmsk_0bp_200bp %>% group_by(TEfamily) %>% dplyr::count())
colnames(TE_rmsk_0bp_200bp.sum) = c("TEfamily","total.instances.tree")

# # order of TE subfamilies from the tree; based on the circular plot and unrooted trees rather than the network
# list_MER11A_g1 = data.frame(group="MER11A_g1",subfamily = c("11A_h","11A_g","11A_f","11A_e","11A_d","11A_j","11A_c","11A_b"))
# list_MER11A_g2 = data.frame(group="MER11A_g2",subfamily = c("11A_l","11A_o","11A_a","11A_p","11A_r"))
# list_MER11A_g3 = data.frame(group="MER11A_g3",subfamily = c("11A_q","11A_k","11A_i","11A_m","11A_n"))
# 
# list_MER11B_g1 = data.frame(group="MER11B_g1",subfamily = c("11B_i","11B_h","11B_g","11B_k","11B_f"))
# list_MER11B_g2 = data.frame(group="MER11B_g2",subfamily = c("11B_q","11B_o","11B_r","11B_n"))
# list_MER11B_g3 = data.frame(group="MER11B_g3",subfamily = c("11B_e","11B_c","11B_b","11B_p","11B_m","11B_d","11B_l","11B_j","11B_a"))
# 
# list_MER11C_g1 = data.frame(group="MER11C_g1",subfamily = c("11C_r","11C_o","11C_w","11C_v"))
# list_MER11C_g2 = data.frame(group="MER11C_g2",subfamily = c("11C_x","11C_q","11C_a","11C_y"))
# list_MER11C_g3 = data.frame(group="MER11C_g3",subfamily = c("11C_z","11C_t","11C_a2","11C_s","11C_b2","11C_b","11C_c","11C_d","11C_e","11C_f","11C_g","11C_h","11C_p"))
# list_MER11C_g4 = data.frame(group="MER11C_g4",subfamily = c("11C_i","11C_j","11C_n","11C_m","11C_k","11C_l"))
# list_MER11D_g = data.frame(group="MER11D",subfamily = c("11D_b","11D_c","11D_a"))
# 
# subfamily.list_MER11 = rbind(list_MER11A_g1,list_MER11A_g2,list_MER11A_g3,
#                    list_MER11B_g1,list_MER11B_g2,list_MER11B_g3,
#                    list_MER11C_g1,list_MER11C_g2,list_MER11C_g3,list_MER11C_g4
#                    #,list_MER11D_g
#                    )


##### MER11 subfamily tree
#list_MER11.distance = read.csv(file="../Final_subfamilies_2023_2_17/Step14_1_hg19_MER11ABC_combined.consensus2.mafft.prank.best.fas.nonrev_dna_rerooted.csv")
#list_MER11.distance = list_MER11.distance[!is.na(list_MER11.distance$subfamily.name.final),c("Order.rerooted.tree","subfamily.name.final2","Node1.","Family")]
#colnames(list_MER11.distance) = c("tree.order","subfamily","distanceToRoot","family")

# updated on 2023/9/5
list_MER11.distance = read.csv(file="input/TEwide_group.list_2023_9_4.functional_group_edited.csv")
list_MER11.distance = list_MER11.distance[list_MER11.distance$subfamily2_1 %in% c("MER11A","MER11B","MER11C"),]
list_MER11.distance = list_MER11.distance[!is.na(list_MER11.distance$subfamily.species),c("Order","subfamily.species","branch.length.ToKeptNode","subfamily2_1","functional.group")]
colnames(list_MER11.distance) = c("tree.order","subfamily","distanceToRoot","family","functional.group")
list_MER11.distance$functional.group = gsub("PG","G",list_MER11.distance$functional.group)
head(list_MER11.distance)

# kept the order based on the combined tree
list_MER11.distance = list_MER11.distance[order(list_MER11.distance$tree.order),]
list_MER11.distance = list_MER11.distance[!is.na(list_MER11.distance$subfamily),]

##### obtain the list of TFs and histone marks of MER11
#"H3K9me3_H1" was removed
#TF.list_MER11 = c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3",
#            "BCL11A","RXRA","POU5F1","NANOG","EP300","YY1","TCF12","SP1","TEAD4","CEBPB","USF1","ZNF143",
#            "ZNF525","ZNF808","ZNF440","ZNF433","ZNF33A.2","ZNF468","ZNF611","KAP1","H3K9me3","H3K27me3")

# reordered within groups 2023/5/24
TF.list_MER11 = c("DNase","H3K27ac","H3K4me3","H3K4me2","H3K4me1",
                  "ZNF143","USF1","CEBPB","TEAD4","SP1","TCF12","EP300","POU5F1","YY1","BCL11A","RXRA","NANOG",
                  "ZNF440","ZNF433","ZNF468","ZNF611","ZNF33A.2","ZNF525","ZNF808","KAP1","H3K9me3","H3K27me3")
#                  "ZNF525","ZNF808","ZNF440","ZNF433","ZNF33A.2","ZNF468","ZNF611","KAP1","H3K9me3","H3K27me3")

TF.list_MER11.active = c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","BCL11A","RXRA","POU5F1","NANOG","EP300","YY1","TCF12","SP1","TEAD4","CEBPB","USF1","ZNF143")

TF.list_MER11.repressive = c("ZNF525","ZNF808","ZNF440","ZNF433","ZNF33A.2","ZNF468","ZNF611","KAP1","H3K9me3","H3K27me3")

TF.list_MER34 = c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","H2BK5ac","H3K9me3","H3K27me3",
                  "POU5F1","CHD7","H2AFZ",
                  "ZNF211","ZNF211.2")
TF.list_MER52 = c("DNase","H3K27ac","H3K4me1","H3K4me2","H3K4me3","H4K8ac","H2BK5ac","H4K91ac","H3K9me3","H3K27me3",
                  "REST","H2AFZ","RAD21","TBP","MAX","ZNF143",
                  "ZNF736","ZNF273","ZNF680","ZNF680.2","ZNF519","ZNF793","ZNF28","ZNF479")

# load the subfamily info from step 4
Subfamily_info1 = read.csv("input_trees/hg19_rmsk_TE_MER11A_0bp.mafft.prank.opt.gt99_group_final.csv")
Subfamily_info2 = read.csv("input_trees/hg19_rmsk_TE_MER11B_0bp.mafft.prank.opt.gt99_group_final.csv")
Subfamily_info3 = read.csv("input_trees/hg19_rmsk_TE_MER11C_0bp.mafft.prank.opt.gt99_group_final.csv")
#Subfamily_info4 = read.csv("input_trees/hg19_rmsk_TE_MER11D_0bp.mafft.prank.opt.gt99_group_final.csv")
Subfamily_info = rbind(Subfamily_info1,Subfamily_info2,Subfamily_info3)

# load the new subfamily names for MER11 for now
Subfamily_names = read.csv("input/MER11_subtree_info_2022_12_27.csv")
Subfamily_info = merge(Subfamily_info,Subfamily_names[,c("consensus.name","subfamily.name.final")],by="consensus.name",all.x=T)

# load the intersection result of ENCODE chip-seq data
TEs.intersected = read.delim("input/hg19_instance_final_2022_12_28.combined.bed",sep="",header=F)
colnames(TEs.intersected) = c("fileID","chr","start","end","TEname","anno1","direction")
TEs.intersected$fileID = gsub("_hg19","",TEs.intersected$fileID)

# load the intersection results of Dnase data E003
Roadmap.original.instance = read.delim("input/hg19_instance_final_2022_12_29.roadmap.combined.bed",header=F,sep="")
colnames(Roadmap.original.instance) = c("fileID","chr","start","end","TEname","anno1","direction")
Roadmap.original.instance = Roadmap.original.instance[Roadmap.original.instance$fileID == "DNase",]
unique(Roadmap.original.instance$fileID)

# add info to the total peaks of DNase
TF.totalPeaks[nrow(TF.totalPeaks)+1,] = c("DNase",339366,"DNase","H1","accessibility")
head(TF.totalPeaks)

# combined the DNase and Chip-seq data
TEs.intersected = rbind(TEs.intersected,Roadmap.original.instance)

# add the group info here
TEs.intersected = merge(TEs.intersected,Subfamily_info[,c("label","subfamily.name.final","TEfamily")],by.x="TEname",by.y="label",all.x=T)
TEs.intersected = TEs.intersected[!duplicated(TEs.intersected[,c("fileID","TEname")]),]
TEs.intersected = TEs.intersected[!is.na(TEs.intersected$subfamily.name.final),]

# summary
TEs.intersected.sum = data.frame(TEs.intersected %>% group_by(TEfamily,subfamily.name.final,fileID) %>% dplyr::count())
TEs.intersected.sum.family = data.frame(TEs.intersected %>% group_by(TEfamily,fileID) %>% dplyr::count())
colnames(TEs.intersected.sum.family) = c("Family","fileID","total.intersected.family")
TEs.intersected.sum.family$uniqueID = paste(TEs.intersected.sum.family$Family,TEs.intersected.sum.family$fileID)
#
TEs.intersected.sum = merge(TEs.intersected.sum,TF.totalPeaks,by="fileID",all.x=T)
TEs.intersected.sum = merge(TEs.intersected.sum,Subfamily_names[,c("subfamily.name.final","Count_with_consensus")],by="subfamily.name.final",all.x=T)   ## needed to be corrected
TEs.intersected.sum$uniqueID = paste(TEs.intersected.sum$TEfamily,TEs.intersected.sum$fileID)
TEs.intersected.sum = merge(TEs.intersected.sum,TEs.intersected.sum.family[,c("uniqueID","total.intersected.family")],by="uniqueID",all.x=T)
TEs.intersected.sum = merge(TEs.intersected.sum,TE_rmsk_0bp_200bp.sum,by.x="TEfamily",by.y="TEfamily",all.x=T)

# adjust the number of instances
TEs.intersected.sum$Count_without_consensus = TEs.intersected.sum$Count_with_consensus
TEs.intersected.sum$Count_without_consensus = ifelse(TEs.intersected.sum$subfamily.name.final %in% c("11A_q","11A_U","11B_i","11C_a","11C_i","11C_t","11C_w"),TEs.intersected.sum$Count_with_consensus - 1,TEs.intersected.sum$Count_without_consensus)
TEs.intersected.sum$Count_without_consensus = ifelse(TEs.intersected.sum$subfamily.name.final == "11A_h",TEs.intersected.sum$Count_with_consensus - 2,TEs.intersected.sum$Count_without_consensus)
TEs.intersected.sum$Count_without_consensus = ifelse(TEs.intersected.sum$subfamily.name.final == "11B_a",TEs.intersected.sum$Count_with_consensus - 3,TEs.intersected.sum$Count_without_consensus)

TEs.intersected.sum[TEs.intersected.sum$Count_with_consensus != TEs.intersected.sum$Count_without_consensus,]
# analyze the frame sequences
Subfamily_info.frame = merge(Subfamily_info,Summary_Table1,by.x="label",by.y="instanceID.renamed",all.x=T)
# exclude the consensus sequences
Subfamily_info.frame = Subfamily_info.frame[Subfamily_info.frame$is_consensus == "instance",]  # do not include the consensus sequences
Subfamily_info.frame$TEfamily = Subfamily_info.frame$TEfamily.x

subfamily = "11A_q"
for (subfamily in unique(Subfamily_info.frame$subfamily.name.final)){
  Subfamily_info.frame_sub = Subfamily_info.frame[Subfamily_info.frame$subfamily.name.final == subfamily,]
  colnames(Subfamily_info.frame_sub)
  Subfamily_info.frame_sub2 = Subfamily_info.frame_sub[!is.na(Subfamily_info.frame_sub$label),c("label","label","subfamily.name.final","consensus.name","cluster.final","count","is_kept.F1","uniqueID_new.F1","iPSC.alpha.F1","iPSC.alpha.Zscore.F1","is_kept.F2","uniqueID_new.F2","iPSC.alpha.F2","iPSC.alpha.Zscore.F2")]
  Subfamily_info.frame_sub2.F1 = Subfamily_info.frame_sub[!is.na(Subfamily_info.frame_sub$uniqueID_new.F1),c("uniqueID_new.F1","label","subfamily.name.final","consensus.name","cluster.final","count","is_kept.F1","uniqueID_new.F1","iPSC.alpha.F1","iPSC.alpha.Zscore.F1","is_kept.F2","uniqueID_new.F2","iPSC.alpha.F2","iPSC.alpha.Zscore.F2")]
  Subfamily_info.frame_sub2.F2 = Subfamily_info.frame_sub[!is.na(Subfamily_info.frame_sub$uniqueID_new.F2),c("uniqueID_new.F2","label","subfamily.name.final","consensus.name","cluster.final","count","is_kept.F1","uniqueID_new.F1","iPSC.alpha.F1","iPSC.alpha.Zscore.F1","is_kept.F2","uniqueID_new.F2","iPSC.alpha.F2","iPSC.alpha.Zscore.F2")]
  #write.table(Subfamily_info.frame_sub2, row.names = FALSE, col.names = FALSE, quote = FALSE,file=paste(subfamily,"_instance_",Date,".id",sep=""))
  #write.table(Subfamily_info.frame_sub2.F1, row.names = FALSE, col.names = FALSE, quote = FALSE,file=paste(subfamily,"_F1_",Date,".id",sep=""))
  #write.table(Subfamily_info.frame_sub2.F2, row.names = FALSE, col.names = FALSE, quote = FALSE,file=paste(subfamily,"_F2_",Date,".id",sep=""))
}

# write to the list files
Subfamily_info.frame.noUngrouped = Subfamily_info.frame[!grepl("_U$",Subfamily_info.frame$subfamily.name.final),]

# exclude low quality activities
Subfamily_info.frame.F1 = Subfamily_info.frame.noUngrouped[Subfamily_info.frame.noUngrouped$is_kept.F1 %in% c("highQuality.iPSC/NPC","highQuality.iPSC","highQuality.NPC"),]    # highQuality.NPC was also kept for the following analysis
Subfamily_info.frame.F2 = Subfamily_info.frame.noUngrouped[Subfamily_info.frame.noUngrouped$is_kept.F2 %in% c("highQuality.iPSC/NPC","highQuality.iPSC","highQuality.NPC"),]    # highQuality.NPC was also kept for the following analysis

# summary
Subfamily_info.sum.frame1 = data.frame(Subfamily_info.frame.F1 %>% 
                                         group_by(TEfamily,subfamily.name.final,log2.iPSC.alpha.range.F1) %>% 
                                         dplyr::summarise(n = n(),mean = mean(log2.iPSC.alpha.F1,na.rm = TRUE),median = median(log2.iPSC.alpha.F1,na.rm = TRUE)) %>% 
                                         mutate(total = sum(n)) %>% 
                                         mutate(freq = n/ sum(n)))
Subfamily_info.sum.frame2 = data.frame(Subfamily_info.frame.F2 %>% 
                                         group_by(TEfamily,subfamily.name.final,log2.iPSC.alpha.range.F2) %>% 
                                         dplyr::summarise(n = n(),mean = mean(log2.iPSC.alpha.F2,na.rm = TRUE),median = median(log2.iPSC.alpha.F2,na.rm = TRUE)) %>% 
                                         mutate(total = sum(n)) %>% 
                                         mutate(freq = n/ sum(n)))

# Zcore

Subfamily_info.sum.frame1 = data.frame(Subfamily_info.frame.F1[!is.na(Subfamily_info.frame.F1$iPSC.alpha.Zscore.isActive.group.F1),] %>% 
                                         group_by(TEfamily,subfamily.name.final,iPSC.alpha.Zscore.isActive.group.F1) %>% 
                                         dplyr::summarise(n = n(),mean = mean(iPSC.alpha.Zscore.F1,na.rm = TRUE),median = median(iPSC.alpha.Zscore.F1,na.rm = TRUE)) %>% 
                                         mutate(total = sum(n)) %>% 
                                         mutate(freq = n/ sum(n)))
Subfamily_info.sum.frame2 = data.frame(Subfamily_info.frame.F2[!is.na(Subfamily_info.frame.F2$iPSC.alpha.Zscore.isActive.group.F2),] %>% 
                                         group_by(TEfamily,subfamily.name.final,iPSC.alpha.Zscore.isActive.group.F2) %>% 
                                         dplyr::summarise(n = n(),mean = mean(iPSC.alpha.Zscore.F2,na.rm = TRUE),median = median(iPSC.alpha.Zscore.F2,na.rm = TRUE)) %>% 
                                         mutate(total = sum(n)) %>% 
                                         mutate(freq = n/ sum(n)))

Subfamily_info.sum.frame1_each = data.frame(Subfamily_info.frame.F1[Subfamily_info.frame.F1$iPSC.alpha.Zscore.isActive.group.F1!="No activity",] %>% 
                                              group_by(TEfamily,subfamily.name.final) %>% 
                                              dplyr::summarise(n = n(),mean = mean(iPSC.alpha.Zscore.F1,na.rm = TRUE),median = median(iPSC.alpha.Zscore.F1,na.rm = TRUE),sd = sd(iPSC.alpha.Zscore.F1,na.rm = TRUE)) %>% 
                                              mutate(total = sum(n)) %>% 
                                              mutate(freq = n/ sum(n)))

Subfamily_info.sum.frame2_each = data.frame(Subfamily_info.frame.F2[Subfamily_info.frame.F2$iPSC.alpha.Zscore.isActive.group.F2!="No activity",] %>% 
                                              group_by(TEfamily,subfamily.name.final) %>% 
                                              dplyr::summarise(n = n(),mean = mean(iPSC.alpha.Zscore.F2,na.rm = TRUE),median = median(iPSC.alpha.Zscore.F2,na.rm = TRUE),sd = sd(iPSC.alpha.Zscore.F2,na.rm = TRUE)) %>% 
                                              mutate(total = sum(n)) %>% 
                                              mutate(freq = n/ sum(n)))


################### Analyze each family
# selected family
TEs.intersected.sum$Experiment.target
TEs.intersected.sum.table = TEs.intersected.sum[,c("fileID","TEfamily","total.instances.tree","subfamily.name.final","subfamily.name.final","Count_without_consensus","n","n")]
TEs.intersected.sum.table$subfamily.color = NA
TEs.intersected.sum.table$node = NA
colnames(TEs.intersected.sum.table) = c("mark","family","family.count_200bp","subfamily","node","subfamily.count_200bp","intersect.count","intersect.count_200bp","subfamily.color")
#write.csv(TEs.intersected.sum.table,file = "MER11_TFs.subfamily.summary_2023_6_27.csv")
write.csv(TEs.intersected.sum.table,file = paste("MER11_TFs.subfamily.summary_",Date,".csv",sep=""))

TEs.intersected.sum_MER11 = TEs.intersected.sum[TEs.intersected.sum$Experiment.target %in% TF.list_MER11,]

TEs.intersected.sum_MER11$pValue = 1.0-phyper(TEs.intersected.sum_MER11$n-1, TEs.intersected.sum_MER11$total.intersected.family, TEs.intersected.sum_MER11$total.instances.tree-TEs.intersected.sum_MER11$total.intersected.family, TEs.intersected.sum_MER11$Count_without_consensus)
TEs.intersected.sum_MER11$padj = p.adjust(TEs.intersected.sum_MER11$pValue,method="BH",n=nrow(TEs.intersected.sum_MER11))
TEs.intersected.sum_MER11$padj.log10 = -log10(TEs.intersected.sum_MER11$padj)
TEs.intersected.sum_MER11$padj.log10 = ifelse(is.infinite(TEs.intersected.sum_MER11$padj.log10),min(TEs.intersected.sum_MER11$padj.log10),TEs.intersected.sum_MER11$padj.log10)
TEs.intersected.sum_MER11$perC = TEs.intersected.sum_MER11$n/TEs.intersected.sum_MER11$Count_without_consensus

TEs.intersected.sum_MER11[TEs.intersected.sum_MER11$Experiment.target=="H3K9me3" & TEs.intersected.sum_MER11$subfamily.name.final == "11A_b",]
TEs.intersected.sum_MER11$Experiment.target_axis = paste(TEs.intersected.sum_MER11$Experiment.target,TEs.intersected.sum_MER11$fileID)
TEs.intersected.sum_MER11$padj.log10 = ifelse(TEs.intersected.sum_MER11$padj.log10>=10,10,TEs.intersected.sum_MER11$padj.log10)

glist = list()
Order = 1
glist_frame = list()
Order_frame = 1
TEs.intersected.sum_MER11.plot = TEs.intersected.sum_MER11[!grepl("_U$",TEs.intersected.sum_MER11$subfamily.name.final) &
                                                             TEs.intersected.sum_MER11$TEfamily %in% c("MER11A","MER11B","MER11C") & 
                                                             TEs.intersected.sum_MER11$Experiment.target %in% TF.list_MER11,]
# order by the re-rooted tree
Subfamilies = list_MER11.distance$subfamily

Subfamily_info.sum.frame1.plot = Subfamily_info.sum.frame1[Subfamily_info.sum.frame1$TEfamily %in% c("MER11A","MER11B","MER11C"),]
Subfamily_info.sum.frame2.plot = Subfamily_info.sum.frame2[Subfamily_info.sum.frame2$TEfamily %in% c("MER11A","MER11B","MER11C"),]


# convert dataframe and obtain the clustered order

TEs.intersected.sum_MER11.convert = data.frame("Experiment.target" = TF.list_MER11)
  
Subfamilies = Subfamilies[!grepl("_U$",Subfamilies)]
for(Subfamily in Subfamilies){
  TEs.intersected.sum_MER11.tmp = TEs.intersected.sum_MER11.plot[TEs.intersected.sum_MER11.plot$subfamily.name.final == Subfamily,]
  TEs.intersected.sum_MER11.convert = merge(TEs.intersected.sum_MER11.convert,TEs.intersected.sum_MER11.tmp[,c("Experiment.target","padj.log10")],by="Experiment.target",all.x=T)
  colnames(TEs.intersected.sum_MER11.convert)[ncol(TEs.intersected.sum_MER11.convert)] = Subfamily
}
rownames(TEs.intersected.sum_MER11.convert) = TEs.intersected.sum_MER11.convert$Experiment.target
TEs.intersected.sum_MER11.convert = as.matrix(TEs.intersected.sum_MER11.convert[,-1])
TEs.intersected.sum_MER11.convert[is.na(TEs.intersected.sum_MER11.convert)] <-0
Order_hclust = hclust(dist(t(TEs.intersected.sum_MER11.convert)), method="median")
Subfamily_order = Order_hclust$labels[Order_hclust$order]

Subfamily_order.active = Subfamilies
Subfamily_order = Subfamilies
TF_order = TF.list_MER11
# 

### plot 1 activity
Subfamily_info.sum.frame1.plot$subfamily.name.final = factor(Subfamily_info.sum.frame1.plot$subfamily.name.final,levels=Subfamily_order.active)
Subfamily_info.sum.frame2.plot$subfamily.name.final = factor(Subfamily_info.sum.frame2.plot$subfamily.name.final,levels=Subfamily_order.active)
Subfamily_info.sum.frame1.plot = Subfamily_info.sum.frame1.plot[Subfamily_info.sum.frame1.plot$total>=10,]
Subfamily_info.sum.frame2.plot = Subfamily_info.sum.frame2.plot[Subfamily_info.sum.frame2.plot$total>=10,]
Subfamily_info.sum.frame1.plot$iPSC.alpha.Zscore.isActive.group.F1 = factor(Subfamily_info.sum.frame1.plot$iPSC.alpha.Zscore.isActive.group.F1,levels=rev(names(iPSC.alpha.Zscore.isActive.group.color)))
Subfamily_info.sum.frame2.plot$iPSC.alpha.Zscore.isActive.group.F2 = factor(Subfamily_info.sum.frame2.plot$iPSC.alpha.Zscore.isActive.group.F2,levels=rev(names(iPSC.alpha.Zscore.isActive.group.color)))
Subfamily_info.sum.frame1.plot$iPSC.alpha.Zscore = factor(Subfamily_info.sum.frame1.plot$iPSC.alpha.Zscore.isActive.group.F1)
Subfamily_info.sum.frame2.plot$iPSC.alpha.Zscore = factor(Subfamily_info.sum.frame2.plot$iPSC.alpha.Zscore.isActive.group.F2)

p_frame = ggplot(Subfamily_info.sum.frame1.plot, aes(x=subfamily.name.final,y=n,fill = iPSC.alpha.Zscore)) +
  geom_bar(stat='identity',position = "fill")+
  scale_x_discrete(drop=FALSE)+
  scale_fill_manual(values = iPSC.alpha.Zscore.isActive.group.color)+
  ylab("Proportion of frames")+
  xlab("Subfamily")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fff7ec"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
        axis.title=element_text(colour="black",size=rel(1)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        #      legend.title="(observed/random)",
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))
glist_frame[[Order_frame]] <- ggplotGrob(p_frame)
Order_frame = Order_frame + 1

p_frame = ggplot(Subfamily_info.sum.frame2.plot, aes(x=subfamily.name.final,y=n,fill = iPSC.alpha.Zscore)) +
  geom_bar(stat='identity',position = "fill")+
  scale_x_discrete(drop=FALSE)+
  scale_fill_manual(values = iPSC.alpha.Zscore.isActive.group.color)+
  ylab("Proportion of frames")+
  xlab("Subfamily")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fff7ec"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
        axis.title=element_text(colour="black",size=rel(1)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        #      legend.title="(observed/random)",
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))
glist_frame[[Order_frame]] <- ggplotGrob(p_frame)
Order_frame = Order_frame + 1 

### plot 2 TF and histone mark
TEs.intersected.sum_MER11.plot$subfamily.name.final = factor(TEs.intersected.sum_MER11.plot$subfamily.name.final,levels=Subfamily_order)
TEs.intersected.sum_MER11.plot$Experiment.target = factor(TEs.intersected.sum_MER11.plot$Experiment.target,levels=TF_order)
# p = ggplot(TEs.intersected.sum_MER11.plot, aes(subfamily.name.final,Experiment.target)) +
#   scale_x_discrete(drop=FALSE)+
#   scale_y_discrete(drop=FALSE)+
#   geom_point(aes(size=perC,fill = (padj.log10)),shape=21,color="black")+
#   scale_fill_gradient2(low = "white", high = "#b2182b",limits = c(0,max(TEs.intersected.sum_MER11$padj.log10)))+
#   ylab("Mark")+
#   xlab("Subfamily")+
#   scale_size(limits=c(0,1))+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.text.y=element_text(colour="black",size=rel(1)),
#         axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
#         axis.title=element_text(colour="black",size=rel(1)),
#         #        legend.position=c(0.8,0.8),
#         #        legend.position="bottom",
#         #      legend.title="(observed/random)",
#         legend.position="right",
#         legend.background = element_blank(),
#         legend.text=element_text(size=rel(1)))
# glist[[Order]] <- ggplotGrob(p)
# Order = Order + 1
# # 
# p = ggplot(TEs.intersected.sum_MER11.plot, aes(subfamily.name.final,Experiment.target)) +
#   geom_point(aes(size=perC,fill = (padj.log10)),shape=21,color="black")+
#   scale_x_discrete(drop=FALSE)+
#   scale_y_discrete(drop=FALSE)+
#   scale_fill_gradient2(low = "white", high = "#b2182b",limits = c(0,max(TEs.intersected.sum_MER11$padj.log10)))+
#   scale_size(limits=c(0,1))+
#   ylab("Mark")+
#   xlab("Subfamily")+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#         axis.text.y=element_text(colour="black",size=rel(1)),
#         axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
#         axis.title=element_text(colour="black",size=rel(1)),
#         #        legend.position=c(0.8,0.8),
#         #        legend.position="bottom",
#         #      legend.title="(observed/random)",
#         legend.position="right",
#         legend.background = element_blank(),
#         legend.text=element_text(size=rel(1)))
# glist[[Order]] <- ggplotGrob(p)
# Order = Order + 1
TEs.intersected.sum_MER11.plot$subfamily.name.final = factor(TEs.intersected.sum_MER11.plot$subfamily.name.final,levels=Subfamily_order.active)
p = ggplot(TEs.intersected.sum_MER11.plot, aes(subfamily.name.final,Experiment.target)) +
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)+
  geom_point(aes(size=perC,fill = (padj.log10)),shape=21,color="black")+
  scale_fill_gradient2(low = "white", high = "#b2182b",limits = c(0,max(TEs.intersected.sum_MER11$padj.log10))) +
  scale_size(limits=c(0,1))+
  ylab("Mark")+
  xlab("Subfamily")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
        axis.title=element_text(colour="black",size=rel(1)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        #      legend.title="(observed/random)",
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))
glist[[Order]] <- ggplotGrob(p)
Order = Order + 1
p = ggplot(TEs.intersected.sum_MER11.plot, aes(subfamily.name.final,Experiment.target)) +
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)+
  geom_point(aes(size=perC,fill = (padj.log10)),shape=21,color="black")+
  scale_fill_gradient2(low = "white", high = "#08306b",limits = c(0,max(TEs.intersected.sum_MER11$padj.log10))) +
  scale_size(limits=c(0,1))+
  ylab("Mark")+
  xlab("Subfamily")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
        axis.title=element_text(colour="black",size=rel(1)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        #      legend.title="(observed/random)",
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))
glist[[Order]] <- ggplotGrob(p)
Order = Order + 1
  Width = 27
  Width_frame = 12

Height = 7
pdf(paste("Figure_2B",".pdf",sep=""),    # create PNG for the heat map        
    height = Height,        # 5 x 300 pixels
    width = Width,
    pointsize = 10 )        # smaller font size
do.call("grid.arrange",c(glist,ncol=2))
dev.off()
pdf(paste("Figure_3E",".pdf",sep=""),    # create PNG for the heat map        
    height = 8,        # 5 x 300 pixels
    width = Width_frame,
    pointsize = 10 )        # smaller font size
do.call("grid.arrange",c(glist_frame,ncol=1))
dev.off()


## sum by phyletic group
TEs.intersected.sum_MER11.comparison = merge(TEs.intersected.sum_MER11.plot,list_MER11.distance[,c("subfamily","functional.group")],by.x="subfamily.name.final",by.y="subfamily",all.x=T)
# the total number of instances analyzed
TEs.intersected.sum_MER11.plot.total = TEs.intersected.sum_MER11.comparison[!duplicated(TEs.intersected.sum_MER11.comparison$subfamily.name.final),c("subfamily.name.final","functional.group","TEfamily","Count_without_consensus")]
TEs.intersected.sum_MER11.plot.total.group = data.frame(TEs.intersected.sum_MER11.plot.total %>% 
                                                          group_by(functional.group) %>% 
                                                          dplyr::summarise(total = sum(Count_without_consensus)))
TEs.intersected.sum_MER11.plot.total.family = data.frame(TEs.intersected.sum_MER11.plot.total %>% 
                                                           group_by(TEfamily) %>% 
                                                           dplyr::summarise(total = sum(Count_without_consensus)))

# group
TEs.intersected.sum_MER11.comparison.group =  data.frame(TEs.intersected.sum_MER11.comparison %>% 
                                                           group_by(Experiment.target,functional.group) %>% 
                                                           dplyr::summarise(intersect = sum(n)))
TEs.intersected.sum_MER11.comparison.group$x.axis = paste(TEs.intersected.sum_MER11.comparison.group$Experiment.target,TEs.intersected.sum_MER11.comparison.group$functional.group)
TEs.intersected.sum_MER11.comparison.group.full = data.frame("x.axis"=paste(rep(TF_order,each=4),c("G1","G2","G3","G4")),"Experiment.target"=rep(TF_order,each=4),"functional.group"=c("G1","G2","G3","G4"))
TEs.intersected.sum_MER11.comparison.group.full = merge(TEs.intersected.sum_MER11.comparison.group.full,TEs.intersected.sum_MER11.comparison.group[,c("x.axis","intersect")],by="x.axis",all.x=T)
TEs.intersected.sum_MER11.comparison.group.full = merge(TEs.intersected.sum_MER11.comparison.group.full,TEs.intersected.sum_MER11.plot.total.group,by="functional.group",all.x=T)
TEs.intersected.sum_MER11.comparison.group.full$x.axis = factor(TEs.intersected.sum_MER11.comparison.group.full$x.axis,
                                                                levels=paste(rep(TF_order,each=4),c("G1","G2","G3","G4")))
TEs.intersected.sum_MER11.comparison.group.full$intersect = ifelse(is.na(TEs.intersected.sum_MER11.comparison.group.full$intersect),0,TEs.intersected.sum_MER11.comparison.group.full$intersect)
TEs.intersected.sum_MER11.comparison.group.full$perC = TEs.intersected.sum_MER11.comparison.group.full$intersect/TEs.intersected.sum_MER11.comparison.group.full$total
TEs.intersected.sum_MER11.comparison.group.full$Type = "group"
colnames(TEs.intersected.sum_MER11.comparison.group.full)[1] = "Group"
TEs.intersected.sum_MER11.comparison.group.full2 = TEs.intersected.sum_MER11.comparison.group.full[TEs.intersected.sum_MER11.comparison.group.full$Experiment.target %in% TF_order[1:17],]
TEs.intersected.sum_MER11.comparison.group.full2$x.axis = factor(TEs.intersected.sum_MER11.comparison.group.full2$x.axis,
                                                                 levels=paste(rep(TF_order[1:17],each=4),c("G1","G2","G3","G4")))
TEs.intersected.sum_MER11.comparison.group.full2$Experiment.target = factor(TEs.intersected.sum_MER11.comparison.group.full2$Experiment.target,levels=TF_order[1:17])

# subfamily
TEs.intersected.sum_MER11.comparison.family =  data.frame(TEs.intersected.sum_MER11.comparison %>% 
                                                            group_by(Experiment.target,TEfamily) %>% 
                                                            dplyr::summarise(intersect = sum(n)))
TEs.intersected.sum_MER11.comparison.family$x.axis = paste(TEs.intersected.sum_MER11.comparison.family$Experiment.target,TEs.intersected.sum_MER11.comparison.family$TEfamily)
TEs.intersected.sum_MER11.comparison.family.full = data.frame("x.axis"=paste(rep(TF_order,each=3),c("MER11A","MER11B","MER11C")),"Experiment.target"=rep(TF_order,each=3),"TEfamily"=c("MER11A","MER11B","MER11C"))
TEs.intersected.sum_MER11.comparison.family.full = merge(TEs.intersected.sum_MER11.comparison.family.full,TEs.intersected.sum_MER11.comparison.family[,c("x.axis","intersect")],by="x.axis",all.x=T)
TEs.intersected.sum_MER11.comparison.family.full = merge(TEs.intersected.sum_MER11.comparison.family.full,TEs.intersected.sum_MER11.plot.total.family,by="TEfamily",all.x=T)
TEs.intersected.sum_MER11.comparison.family.full$x.axis = factor(TEs.intersected.sum_MER11.comparison.family.full$x.axis,
                                                                 levels=paste(rep(TF_order,each=3),c("MER11A","MER11B","MER11C")))
TEs.intersected.sum_MER11.comparison.family.full$intersect = ifelse(is.na(TEs.intersected.sum_MER11.comparison.family.full$intersect),0,TEs.intersected.sum_MER11.comparison.family.full$intersect)
TEs.intersected.sum_MER11.comparison.family.full$perC = TEs.intersected.sum_MER11.comparison.family.full$intersect/TEs.intersected.sum_MER11.comparison.family.full$total
TEs.intersected.sum_MER11.comparison.family.full$Type = "subfamily"
colnames(TEs.intersected.sum_MER11.comparison.family.full)[1] = "Group"
TEs.intersected.sum_MER11.comparison.family.full2 = TEs.intersected.sum_MER11.comparison.family.full[TEs.intersected.sum_MER11.comparison.family.full$Experiment.target %in% TF_order[1:17],]
TEs.intersected.sum_MER11.comparison.family.full2$x.axis = factor(TEs.intersected.sum_MER11.comparison.family.full2$x.axis,
                                                                  levels=paste(rep(TF_order[1:17],each=3),c("MER11A","MER11B","MER11C")))
TEs.intersected.sum_MER11.comparison.family.full2$Experiment.target = factor(TEs.intersected.sum_MER11.comparison.family.full2$Experiment.target,levels=TF_order[1:17])

TEs.intersected.sum_MER11.comparison.combined = rbind(TEs.intersected.sum_MER11.comparison.group.full,TEs.intersected.sum_MER11.comparison.family.full)
TEs.intersected.sum_MER11.comparison.combined$uniqueID = paste(TEs.intersected.sum_MER11.comparison.combined$Type,TEs.intersected.sum_MER11.comparison.combined$Experiment.target)
# combined
head(TEs.intersected.sum_MER11.comparison.combined)
head(TEs.intersected.sum_MER11.comparison.family)
head(TEs.intersected.sum_MER11.comparison.group)
TEs.intersected.sum_MER11.comparison.combined$perC

max_ones = data.frame(TEs.intersected.sum_MER11.comparison.combined %>%
                        group_by(uniqueID) %>%
                        dplyr::slice(which.max(perC)))
min_ones = data.frame(TEs.intersected.sum_MER11.comparison.combined %>%
                        group_by(uniqueID) %>%
                        dplyr::slice(which.min(perC)))
TEs.intersected.sum_MER11.comparison.combined2 = merge(max_ones[,c("uniqueID","Experiment.target","Type","perC")],min_ones[,c("uniqueID","perC")],by="uniqueID",all=T)
TEs.intersected.sum_MER11.comparison.combined2$perC.diff = TEs.intersected.sum_MER11.comparison.combined2$perC.x-TEs.intersected.sum_MER11.comparison.combined2$perC.y
TEs.intersected.sum_MER11.comparison.combined2$Experiment.target = factor(TEs.intersected.sum_MER11.comparison.combined2$Experiment.target,levels=TF_order)

# plot
my_pal <- colorRampPalette(brewer.pal(8, "Paired"))(n = length(TF_order))
my_pal2 <- colorRampPalette(brewer.pal(8, "Paired"))(n = length(TF_order[1:17]))
glist_com = list()
Order_com = 1
p1 = ggplot(TEs.intersected.sum_MER11.comparison.combined2, aes(x=Experiment.target,y=perC.diff*100,group=Type)) +
  geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(perC.diff*100,1)), vjust=-.2, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values =c("#016A70","#A2C579"))+
  geom_hline(yintercept=c(20,40), linetype="dashed", color = "black")+
  ylim(0,100)+
  #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
  scale_x_discrete(drop = FALSE)+
  ylab("specificity (max - min % peaks-associated instances)")+
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
glist_com[[Order_com]] <- ggplotGrob(p1)
Order_com = Order_com + 1

TEs.intersected.sum_MER11.comparison.combined3 = TEs.intersected.sum_MER11.comparison.combined2[TEs.intersected.sum_MER11.comparison.combined2$Experiment.target %in% TF_order[1:17],]
TEs.intersected.sum_MER11.comparison.combined3$Experiment.target = factor(TEs.intersected.sum_MER11.comparison.combined3$Experiment.target,levels=TF_order[1:17])
p1 = ggplot(TEs.intersected.sum_MER11.comparison.combined3, aes(x=Experiment.target,y=perC.diff*100,group=Type)) +
  geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(perC.diff*100,1)), vjust=-.2, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values =c("#016A70","#A2C579"))+
  geom_hline(yintercept=c(20,40), linetype="dashed", color = "black")+
  ylim(0,50)+
  #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
  scale_x_discrete(drop = FALSE)+
  ylab("specificity (max - min % peaks-associated instances)")+
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
glist_com[[Order_com]] <- ggplotGrob(p1)
Order_com = Order_com + 1
p1 = ggplot(TEs.intersected.sum_MER11.comparison.combined3, aes(x=Experiment.target,y=perC.x*100,group=Type)) +
  geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(perC.x*100,1)), vjust=-.2, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values =c("#016A70","#A2C579"))+
  geom_hline(yintercept=c(20,40), linetype="dashed", color = "black")+
  ylim(0,50)+
  #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
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
glist_com[[Order_com]] <- ggplotGrob(p1)
Order_com = Order_com + 1
p1 = ggplot(TEs.intersected.sum_MER11.comparison.combined3, aes(x=Experiment.target,y=perC.y*100,group=Type)) +
  geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
  geom_text(aes(label=round(perC.y*100,1)), vjust=-.2, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values =c("#016A70","#A2C579"))+
  geom_hline(yintercept=c(20,40), linetype="dashed", color = "black")+
  ylim(0,50)+
  #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
  scale_x_discrete(drop = FALSE)+
  ylab("min % peaks-associated instances")+
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
glist_com[[Order_com]] <- ggplotGrob(p1)
Order_com = Order_com + 1

p1 = ggplot(TEs.intersected.sum_MER11.comparison.group.full, aes(x=x.axis,y=perC)) +
  geom_col(aes(fill=Experiment.target)) + 
  scale_fill_manual(values =my_pal)+
  geom_hline(yintercept=c(0.2,0.4), linetype="dashed", color = "black")+
  ylim(0,1)+
  #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
  scale_x_discrete(drop = FALSE)+
  ylab("% peaks associated instances")+
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
glist_com[[Order_com]] <- ggplotGrob(p1)
Order_com = Order_com + 1

pdf(paste("Figure_2C","-comparison.pdf",sep=""),    # create PNG for the heat map        
    width = 12,        # 5 x 300 pixels
    height = 15,
    pointsize = 10)        # smaller font size
do.call("grid.arrange",c(glist_com,ncol=1))
dev.off()

## divergent rate
Subfamily_info.div = merge(Subfamily_info,Summary_Table1[,c("instanceID.renamed","Instance_len","div_rate")],by.x="label",by.y="instanceID.renamed",all.x=T)
Subfamily_info.div = Subfamily_info.div[Subfamily_info.div$subfamily.name.final %in% list_MER11.distance$subfamily &
                                          Subfamily_info.div$species == "hg19",]
Subfamily_info.div.sum = data.frame(Subfamily_info.div %>% 
                                      group_by(TEfamily,subfamily.name.final) %>% 
                                      dplyr::summarise(n = n(),mean.div = mean(div_rate,na.rm = TRUE),median.div = median(div_rate,na.rm = TRUE),
                                                       mean.len = mean(Instance_len,na.rm = TRUE),median.len = median(Instance_len,na.rm = TRUE)) %>% 
                                      mutate(total = sum(n)) %>% 
                                      mutate(freq = n/ sum(n)))
Subfamily_info.div.sum$subfamily.name.final = factor(Subfamily_info.div.sum$subfamily.name.final,levels=list_MER11.distance$subfamily)
# div
p1 = ggplot(Subfamily_info.div.sum, aes(subfamily.name.final,1)) +
  geom_tile(aes(fill = (median.div)),col = NA,size=0.1) + 
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 10,na.value = 'black')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
        axis.title=element_text(colour="black",size=rel(1)),
        #        legend.position=c(0.8,0.8),
        #        legend.position="bottom",
        #      legend.title="(observed/random)",
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))
p2 = ggplot(Subfamily_info.div.sum, aes(subfamily.name.final,1)) +
  geom_tile(aes(fill = (median.len)),col = NA,size=0.1) + 
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 1050,na.value = 'black')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))
list_MER11.distance$subfamily = factor(list_MER11.distance$subfamily,levels=list_MER11.distance$subfamily)
p3 = ggplot(list_MER11.distance, aes(subfamily,1)) +
  geom_tile(aes(fill = (distanceToRoot)),col = NA,size=0.1) + 
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 0.2,na.value = 'black')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.y=element_text(colour="black",size=rel(1)),
        axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
        axis.title=element_text(colour="black",size=rel(1)),
        legend.position="right",
        legend.background = element_blank(),
        legend.text=element_text(size=rel(1)))

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
gC <- ggplotGrob(p3)

g2 = rbind(gA,gB,gC, size = "first")
pdf(paste("Figure_3F","-annotation.pdf",sep=""),    # create PNG for the heat map        
    height = 5,        # 5 x 300 pixels
    width = Width,
    pointsize = 10 )        # smaller font size
grid.draw(g2)
dev.off()
