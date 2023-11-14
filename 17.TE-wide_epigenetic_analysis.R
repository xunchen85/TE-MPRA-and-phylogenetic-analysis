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

#########################
Date = "2023_9_13"
############################ step 1: prepare the group.list format 
### for MER11 (a non-used hg19 MER11C tree was included by mistake, which do not impact the results but should be excluded in the final figure)  ## corrected in the final input
hg19.target.families = read.csv("input/hg19_rmsk_TE_0bp_MER11_34_52_group_final_2022_11_9.csv")
hg19.target.families$is_consensus = ifelse(grepl("MER",hg19.target.families$node.info_2),hg19.target.families$consensus.name,hg19.target.families$is_consensus)
hg19.target.families$consensus.name = ifelse(grepl("MER",hg19.target.families$node.info_2),hg19.target.families$perC,hg19.target.families$consensus.name)
hg19.target.families$perC = ifelse(grepl("MER",hg19.target.families$node.info_2),hg19.target.families$count,hg19.target.families$perC)
hg19.target.families$count = ifelse(grepl("MER",hg19.target.families$node.info_2),hg19.target.families$TEfamily.representative,hg19.target.families$count)
hg19.target.families$TEfamily.representative = ifelse(grepl("MER",hg19.target.families$node.info_2),hg19.target.families$node.info_2,hg19.target.families$TEfamily.representative)
hg19.target.families$TEfamily.representative = hg19.target.families$consensus.name
hg19.target.families$is_consensus = hg19.target.families$TEfamily
hg19.target.families$X = hg19.target.families$consensus.name
hg19.target.families$edgeID = gsub(" ","\t",hg19.target.families$edgeID)
#write.table(hg19.target.families, sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE,file=paste("hg19_rmsk_TE_0bp_MER11_34_52_group_final_2022_11_9",".list",sep=""))

# 
# ### for other families with a small number of instances
# missed.families.hg19 = c("LTR5","HERV1_LTRb","HERV1_LTRd")
# missed.families.hg19.table  = hg19.target.families[1,]
# missed.families.hg19.table[1,] = NA
# missed.families.hg19.table.final = missed.families.hg19.table

# # LTR5
# for (Family in missed.families.hg19){
#   missed.families.hg19.fasta = read.delim(paste("input_trees/hg19_rmsk_",Family,"_0bp.mafft.prank.opt.gt99.fa",sep=""),header=F,sep="")
#   missed.families.hg19.fasta = gsub(">","",missed.families.hg19.fasta[grepl(">",missed.families.hg19.fasta$V1),])
#   missed.families.hg19.table.tmp = missed.families.hg19.table
#   missed.families.hg19.table.tmp[1:length(missed.families.hg19.fasta),]$label = missed.families.hg19.fasta
#   if (Family == "LTR5"){
#     Family.cluster = "LTR5"
#   } else {
#     Family.cluster = "LTR15"
#   }
#   missed.families.hg19.table.tmp[1:length(missed.families.hg19.fasta),]$is_consensus = Family
#   missed.families.hg19.table.tmp[1:length(missed.families.hg19.fasta),]$consensus.name = paste(Family,"_g0",sep="")
#   missed.families.hg19.table.tmp[1:length(missed.families.hg19.fasta),]$TEfamily.representative = paste(Family,"_g0",sep="")
#   missed.families.hg19.table.tmp[1:length(missed.families.hg19.fasta),]$X = paste(Family,"_g0",sep="")
#   missed.families.hg19.table.tmp[1:length(missed.families.hg19.fasta),]$count = length(missed.families.hg19.fasta)
#   missed.families.hg19.table.tmp[1:length(missed.families.hg19.fasta),]$perC = length(missed.families.hg19.fasta)
#   missed.families.hg19.table.final = rbind(missed.families.hg19.table.final,missed.families.hg19.table.tmp)
# }
# missed.families.hg19.table.final = missed.families.hg19.table.final[-1,]
# missed.families.hg19.table.final[is.na(missed.families.hg19.table.final)]<-"-"
# missed.families.hg19.table.final$branch.length = missed.families.hg19.table.final$label
# #write.table(missed.families.hg19.table.final, sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE,file=paste("hg19_rmsk_TE_0bp_missed_families_group_final",".list_2023_7_1",sep=""))

############################### load the metadata of chip-seq
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
TF.totalPeaks = read.delim("input/total.peaks_2023_6_29",sep="",header=F)
colnames(TF.totalPeaks) = c("fileID","total.peaks")
TF.totalPeaks$fileID = gsub("_hg19","",TF.totalPeaks$fileID)
TF.totalPeaks = merge(TF.totalPeaks,metadata,by.x="fileID",by.y="File.accession",all.x=T)

# Roadmap
TF.totalPeaks$Experiment.target = ifelse(grepl("^E00[34567]",TF.totalPeaks$fileID),gsub("^E00[34567]-","",TF.totalPeaks$fileID),TF.totalPeaks$Experiment.target)
TF.totalPeaks$cell = ifelse(grepl("^E003",TF.totalPeaks$fileID),"E003",TF.totalPeaks$cell)
TF.totalPeaks$cell = ifelse(grepl("^E004",TF.totalPeaks$fileID),"E004",TF.totalPeaks$cell)
TF.totalPeaks$cell = ifelse(grepl("^E005",TF.totalPeaks$fileID),"E005",TF.totalPeaks$cell)
TF.totalPeaks$cell = ifelse(grepl("^E006",TF.totalPeaks$fileID),"E006",TF.totalPeaks$cell)
TF.totalPeaks$cell = ifelse(grepl("^E007",TF.totalPeaks$fileID),"E007",TF.totalPeaks$cell)
TF.totalPeaks$mark.type = ifelse(grepl("^E00[34567]",TF.totalPeaks$fileID),"histone",TF.totalPeaks$mark.type)

# Taka's data
TF.totalPeaks$Experiment.target = ifelse(grepl("^Taka",TF.totalPeaks$fileID) & grepl("K27ac",TF.totalPeaks$fileID),
                                         paste("H3K27ac",gsub("Taka_7TP_|_K27ac|Taka_3TP_|Taka_","",TF.totalPeaks$fileID),sep="_"),TF.totalPeaks$Experiment.target)
TF.totalPeaks$Experiment.target = ifelse(grepl("^Taka",TF.totalPeaks$fileID) & grepl("ATAC",TF.totalPeaks$fileID),
                                         paste("ATAC",gsub("Taka_7TP_|_ATAC|Taka_3TP_|Taka_","",TF.totalPeaks$fileID),sep="_"),TF.totalPeaks$Experiment.target)
TF.totalPeaks$cell = ifelse(grepl("^Taka",TF.totalPeaks$fileID),"iPSC",TF.totalPeaks$cell)
TF.totalPeaks$mark.type = ifelse(grepl("^Taka",TF.totalPeaks$fileID),"histone",TF.totalPeaks$mark.type)

# KRAB-ZNFs
TF.totalPeaks$Experiment.target = ifelse(is.na(TF.totalPeaks$Experiment.target),TF.totalPeaks$fileID,TF.totalPeaks$Experiment.target)
TF.totalPeaks$cell = ifelse(is.na(TF.totalPeaks$cell),"HEK293T",TF.totalPeaks$cell)
TF.totalPeaks$mark.type = ifelse(TF.totalPeaks$cell == "HEK293T","KZNF",TF.totalPeaks$mark.type)
TF.totalPeaks$cell.mark = paste(TF.totalPeaks$cell,TF.totalPeaks$mark.type,TF.totalPeaks$Experiment.target,sep="|")
TF.totalPeaks$cell.mark = gsub("-rep",".rep",TF.totalPeaks$cell.mark)
TF.totalPeaks = TF.totalPeaks[order(-TF.totalPeaks$total.peaks),]
TF.totalPeaks = TF.totalPeaks[!duplicated(TF.totalPeaks$cell.mark),]

# exclude 
TF.totalPeaks.excluded = TF.totalPeaks[(TF.totalPeaks$cell == "E003" & !TF.totalPeaks$Experiment.target %in% c("DNase","H3K27ac")) |
                                         TF.totalPeaks$cell == "NPC",]

TF.totalPeaks = TF.totalPeaks[!TF.totalPeaks$cell.mark %in% TF.totalPeaks.excluded$cell.mark,]
#write.csv(TF.totalPeaks,file="TF.totalPeaks_2023_10_19.csv")

################################ step 2
### load the TE-wide chipseq summary table including MER11_34_52
TEwide_TFs_subfamily = read.delim("input/TEwide_TFs.subfamily.summary_2023_7_1.gz",header=F,sep="")
colnames(TEwide_TFs_subfamily) = c("mark","family","family.count_200bp","subfamily","node","subfamily.count_200bp","intersect.count","intersect.count_200bp","subfamily.color")
TEwide_TFs_subfamily$mark = gsub("MER11_combined.hg19.","",TEwide_TFs_subfamily$mark)
TEwide_TFs_subfamily$mark = gsub("_hg19","",TEwide_TFs_subfamily$mark)
TEwide_TFs_subfamily$uniqueID = paste(TEwide_TFs_subfamily$mark,TEwide_TFs_subfamily$family,TEwide_TFs_subfamily$subfamily)
TEwide_TFs_subfamily$subfamily.count_200bp = as.numeric(as.character(TEwide_TFs_subfamily$subfamily.count_200bp))

# corrected for MER11 families
TEwide_TFs_subfamily$subfamily.count_200bp = ifelse(TEwide_TFs_subfamily$subfamily %in% c("11A_q","11A_U","11B_i","11C_a","11C_i","11C_t","11C_w") | 
                                                      TEwide_TFs_subfamily$subfamily %in% c("11A|N8|303","11A|U|53","11B|N37|13","11C|N1|148","11C|N118|45","11C|N243|16","11C|N37|11"),TEwide_TFs_subfamily$subfamily.count_200bp - 1,TEwide_TFs_subfamily$subfamily.count_200bp)
TEwide_TFs_subfamily$subfamily.count_200bp = ifelse(TEwide_TFs_subfamily$subfamily == "11A_h" |
                                                      TEwide_TFs_subfamily$subfamily == "11A|N161|17",TEwide_TFs_subfamily$subfamily.count_200bp - 2,TEwide_TFs_subfamily$subfamily.count_200bp)
TEwide_TFs_subfamily$subfamily.count_200bp = ifelse(TEwide_TFs_subfamily$subfamily == "11B_a" |
                                                      TEwide_TFs_subfamily$subfamily == "11B|N1|30",TEwide_TFs_subfamily$subfamily.count_200bp - 3,TEwide_TFs_subfamily$subfamily.count_200bp)
TEwide_TFs_subfamily$family.count_200bp = ifelse(TEwide_TFs_subfamily$family == "MER11A",907,TEwide_TFs_subfamily$family.count_200bp)
TEwide_TFs_subfamily$family.count_200bp = ifelse(TEwide_TFs_subfamily$family == "MER11B",508,TEwide_TFs_subfamily$family.count_200bp)
TEwide_TFs_subfamily$family.count_200bp = ifelse(TEwide_TFs_subfamily$family == "MER11C",843,TEwide_TFs_subfamily$family.count_200bp)
TEwide_TFs_subfamily.unique = TEwide_TFs_subfamily[!duplicated(TEwide_TFs_subfamily$subfamily),]
head(TEwide_TFs_subfamily[TEwide_TFs_subfamily$mark == "Taka_72hr_K27ac_r1",])

### shuffle table
TEwide_TFs_subfamily.shuffle = read.delim("input/TEwide_TFs.subfamily.shuffle.summary_2023_6_30.gz",header=F,sep="")
colnames(TEwide_TFs_subfamily.shuffle) = c("mark","shuffle","family","family.count_200bp","subfamily","node","subfamily.count_200bp","intersect.count","intersect.count_200bp","subfamily.color")
TEwide_TFs_subfamily.shuffle$mark = gsub("MER11_combined.hg19.","",TEwide_TFs_subfamily.shuffle$mark)
TEwide_TFs_subfamily.shuffle$mark = gsub("_hg19","",TEwide_TFs_subfamily.shuffle$mark)
TEwide_TFs_subfamily.shuffle$uniqueID = paste(TEwide_TFs_subfamily.shuffle$mark,TEwide_TFs_subfamily.shuffle$family,TEwide_TFs_subfamily.shuffle$subfamily)
TEwide_TFs_subfamily.shuffle = TEwide_TFs_subfamily.shuffle[,c("uniqueID","shuffle","intersect.count_200bp")]

TEwide_TFs_subfamily.shuffle.sum = data.frame(TEwide_TFs_subfamily.shuffle %>% 
                                                group_by(uniqueID) %>% 
                                                dplyr::summarise(count.shuffle = n(),mean.value = mean(intersect.count_200bp,na.rm=T),sd.value = sd(intersect.count_200bp,na.rm=T)))
# 
TEwide_TFs_subfamily.shuffle.sum = TEwide_TFs_subfamily.shuffle.sum[TEwide_TFs_subfamily.shuffle.sum$count.shuffle == 100,]
head(TEwide_TFs_subfamily.shuffle.sum[grepl("Taka_72hr_K27ac_r1",TEwide_TFs_subfamily.shuffle.sum$uniqueID),])

## load the latest TE-wide table 
ordered_subfamilies.table = read.csv(file = "input/TEwide_group.list_2023_9_4.functional_group_edited.csv")

#################################################
################################## LTR7 subfamily level enrichment analysis
Roadmap.marks = c("E003|histone|DNase","E004|histone|DNase","E005|histone|DNase","E006|histone|DNase","E007|histone|DNase",
                  "E003|histone|H3K27ac","E004|histone|H3K27ac","E005|histone|H3K27ac","E006|histone|H3K27ac","E007|histone|H3K27ac")
iPSC.marks = c("iPSC|histone|ATAC_0hr_r1","iPSC|histone|ATAC_0hr_r2","iPSC|histone|ATAC_3hr_r1","iPSC|histone|ATAC_3hr_r2",
               "iPSC|histone|ATAC_6hr_r1","iPSC|histone|ATAC_6hr_r2","iPSC|histone|ATAC_12hr_r1","iPSC|histone|ATAC_12hr_r2",
               "iPSC|histone|ATAC_24hr_r1","iPSC|histone|ATAC_24hr_r2","iPSC|histone|ATAC_48hr_r1","iPSC|histone|ATAC_48hr_r2",
               "iPSC|histone|ATAC_72hr_r1","iPSC|histone|ATAC_72hr_r2",
               "iPSC|histone|H3K27ac_0hr_r1","iPSC|histone|H3K27ac_0hr_r2","iPSC|histone|H3K27ac_3hr_r1","iPSC|histone|H3K27ac_3hr_r2",
               "iPSC|histone|H3K27ac_6hr_r1","iPSC|histone|H3K27ac_6hr_r2","iPSC|histone|H3K27ac_12hr_r1","iPSC|histone|H3K27ac_12hr_r2",
               "iPSC|histone|H3K27ac_24hr_r1","iPSC|histone|H3K27ac_24hr_r2","iPSC|histone|H3K27ac_48hr_r1","iPSC|histone|H3K27ac_48hr_r2",
               "iPSC|histone|H3K27ac_72hr_r1","iPSC|histone|H3K27ac_72hr_r2")

TEwide_TFs_subfamily.final = TEwide_TFs_subfamily[TEwide_TFs_subfamily$family %in% c("LTR7C","LTR7B","LTR7_Rhe","LTR7","LTR7Y"),]

# more than 5 intersected peaks
TEwide_TFs_subfamily.final.marks = unique(TEwide_TFs_subfamily.final[TEwide_TFs_subfamily.final$intersect.count_200bp>=5,]$mark)
TEwide_TFs_subfamily.final = TEwide_TFs_subfamily.final[TEwide_TFs_subfamily.final$mark %in% c(TEwide_TFs_subfamily.final.marks,Roadmap.marks,iPSC.marks),]
TEwide_TFs_subfamily.final$pvalue = NA

# shuffle data
TEwide_TFs_subfamily.shuffle.tmp = TEwide_TFs_subfamily.shuffle[TEwide_TFs_subfamily.shuffle$uniqueID %in% TEwide_TFs_subfamily.final$uniqueID,]
# loop
for (i in 1:nrow(TEwide_TFs_subfamily.final)){
  TEwide_TFs_subfamily.shuffle.plot = TEwide_TFs_subfamily.shuffle.tmp[TEwide_TFs_subfamily.shuffle.tmp$uniqueID == TEwide_TFs_subfamily.final[i,]$uniqueID,]
  TEwide_TFs_subfamily.final[i,]$pvalue = 2*mean(as.numeric(as.character(TEwide_TFs_subfamily.shuffle.plot$intersect.count_200bp)) >= as.numeric(as.character(TEwide_TFs_subfamily.final[i,]$intersect.count_200bp)))
}
# merge with non-shuffled data
TEwide_TFs_subfamily.final = merge(TEwide_TFs_subfamily.final,TEwide_TFs_subfamily.shuffle.sum[,c("uniqueID","mean.value","sd.value")],by="uniqueID",all.x=T)
TEwide_TFs_subfamily.final$log2FC = log2((TEwide_TFs_subfamily.final$intersect.count_200bp+1)/(TEwide_TFs_subfamily.final$mean.value+1))
TEwide_TFs_subfamily.final$perC = TEwide_TFs_subfamily.final$intersect.count_200bp/as.numeric(as.character(TEwide_TFs_subfamily.final$subfamily.count_200bp))
#write.csv(TEwide_TFs_subfamily.final,file = "TEwide_TFs_subfamily.final_LTR7_2023_9_4.csv")

#################################################
################################## LTR5 subfamily level enrichment analysis
TEwide_TFs_subfamily.final = TEwide_TFs_subfamily[TEwide_TFs_subfamily$family %in% c("LTR5","LTR5_Hs","LTR5A","LTR5B"),]
# more than 5 intersected peaks
TEwide_TFs_subfamily.final.marks = unique(TEwide_TFs_subfamily.final[TEwide_TFs_subfamily.final$intersect.count_200bp>=5,]$mark)
TEwide_TFs_subfamily.final = TEwide_TFs_subfamily.final[TEwide_TFs_subfamily.final$mark %in% c(TEwide_TFs_subfamily.final.marks,Roadmap.marks,iPSC.marks),]
TEwide_TFs_subfamily.final$pvalue = NA

# shuffle data
TEwide_TFs_subfamily.shuffle.tmp = TEwide_TFs_subfamily.shuffle[TEwide_TFs_subfamily.shuffle$uniqueID %in% TEwide_TFs_subfamily.final$uniqueID,]
# loop
for (i in 1:nrow(TEwide_TFs_subfamily.final)){
  TEwide_TFs_subfamily.shuffle.plot = TEwide_TFs_subfamily.shuffle.tmp[TEwide_TFs_subfamily.shuffle.tmp$uniqueID == TEwide_TFs_subfamily.final[i,]$uniqueID,]
  TEwide_TFs_subfamily.final[i,]$pvalue = 2*mean(as.numeric(as.character(TEwide_TFs_subfamily.shuffle.plot$intersect.count_200bp)) >= as.numeric(as.character(TEwide_TFs_subfamily.final[i,]$intersect.count_200bp)))
}
# merge with non-shuffled data
TEwide_TFs_subfamily.final = merge(TEwide_TFs_subfamily.final,TEwide_TFs_subfamily.shuffle.sum[,c("uniqueID","mean.value","sd.value")],by="uniqueID",all.x=T)
TEwide_TFs_subfamily.final$log2FC = log2((TEwide_TFs_subfamily.final$intersect.count_200bp+1)/(TEwide_TFs_subfamily.final$mean.value+1))
TEwide_TFs_subfamily.final$perC = TEwide_TFs_subfamily.final$intersect.count_200bp/as.numeric(as.character(TEwide_TFs_subfamily.final$subfamily.count_200bp))
#write.csv(TEwide_TFs_subfamily.final,file = "TEwide_TFs_subfamily.final_LTR5_2023_9_4.csv")

################################################# 
##################################### TEwide and functional groups 
# select families
TEwide_TFs_subfamily.final = TEwide_TFs_subfamily[TEwide_TFs_subfamily$family %in% c(unique(ordered_subfamilies.table$subfamily2_1)),]

# combine with the family information
TEwide_TFs_subfamily.final = merge(TEwide_TFs_subfamily.final,ordered_subfamilies.table[,c("subfamily","functional.group","family.group","family.cluster","subfamily2_1")],by="subfamily",all.x=T)

# a new table for the family.cluster and functional groups
TEwide_TFs_subfamily.final$subfamily.count_200bp = as.numeric(as.character(TEwide_TFs_subfamily.final$subfamily.count_200bp))
TEwide_TFs_subfamily.final$intersect.count_200bp = as.numeric(as.character(TEwide_TFs_subfamily.final$intersect.count_200bp))
TEwide_TFs_subfamily.final$uniqueID.FG = paste(TEwide_TFs_subfamily.final$mark,TEwide_TFs_subfamily.final$family.group,TEwide_TFs_subfamily.final$functional.group,TEwide_TFs_subfamily.final$family.cluster)

# functional group table
TEwide_TFs_subfamily.final.FG = data.frame(TEwide_TFs_subfamily.final %>% 
                                             group_by(mark,family.group,functional.group,family.cluster) %>% 
                                             dplyr::summarise(no.subfamilies.FG = n(),
                                                              subfamily.count_200bp.FG = sum(subfamily.count_200bp,na.rm=T),
                                                              intersect.count_200bp.FG = sum(intersect.count_200bp,na.rm=T)))
TEwide_TFs_subfamily.final.FG$uniqueID.FG = paste(TEwide_TFs_subfamily.final.FG$mark,TEwide_TFs_subfamily.final.FG$family.group,TEwide_TFs_subfamily.final.FG$functional.group,TEwide_TFs_subfamily.final.FG$family.cluster)

# more than 5 intersected peaks and Taka's data
TEwide_TFs_subfamily.final.FG.marks = unique(TEwide_TFs_subfamily.final.FG[(TEwide_TFs_subfamily.final.FG$intersect.count_200bp>=5 & TEwide_TFs_subfamily.final.FG$subfamily.count_200bp.FG<100) |
                                                                             (TEwide_TFs_subfamily.final.FG$intersect.count_200bp.FG/TEwide_TFs_subfamily.final.FG$subfamily.count_200bp.FG>= 0.05 & 
                                                                                TEwide_TFs_subfamily.final.FG$subfamily.count_200bp.FG>=100) | grepl("^Taka_",TEwide_TFs_subfamily.final.FG$mark),]$mark)
TEwide_TFs_subfamily.final.FG = TEwide_TFs_subfamily.final.FG[!is.na(TEwide_TFs_subfamily.final.FG$mark),]
TEwide_TFs_subfamily.final.FG$pvalue = NA
TEwide_TFs_subfamily.final.FG$mean.value = NA
TEwide_TFs_subfamily.final.FG$sd.value = NA

# shuffle data
TEwide_TFs_subfamily.shuffle.tmp = TEwide_TFs_subfamily.shuffle[TEwide_TFs_subfamily.shuffle$uniqueID %in% TEwide_TFs_subfamily.final$uniqueID,]
i = 1
for (i in 1:nrow(TEwide_TFs_subfamily.final.FG)){
  TEwide_TFs_subfamily.final.tmp = TEwide_TFs_subfamily.final[TEwide_TFs_subfamily.final$uniqueID.FG == TEwide_TFs_subfamily.final.FG[i,]$uniqueID.FG,]
  TEwide_TFs_subfamily.shuffle.plot = TEwide_TFs_subfamily.shuffle.tmp[TEwide_TFs_subfamily.shuffle.tmp$uniqueID %in% TEwide_TFs_subfamily.final.tmp$uniqueID,]
  TEwide_TFs_subfamily.shuffle.plot2 = data.frame(TEwide_TFs_subfamily.shuffle.plot %>% 
                                                    group_by(shuffle) %>% 
                                                    dplyr::summarise(no.subfamilies.FG = n(),
                                                                     intersect.count_200bp.FG = sum(intersect.count_200bp,na.rm=T)))
  TEwide_TFs_subfamily.final.FG[i,]$pvalue = 2*mean(as.numeric(as.character(TEwide_TFs_subfamily.shuffle.plot2$intersect.count_200bp.FG)) >= as.numeric(as.character(TEwide_TFs_subfamily.final.FG[i,]$intersect.count_200bp.FG)))
  TEwide_TFs_subfamily.final.FG[i,]$mean.value = mean(TEwide_TFs_subfamily.shuffle.plot2$intersect.count_200bp.FG,na.rm=T)
  TEwide_TFs_subfamily.final.FG[i,]$sd.value = sd(TEwide_TFs_subfamily.shuffle.plot2$intersect.count_200bp.FG,na.rm=T)
}

TEwide_TFs_subfamily.final.FG$log2FC = log2((TEwide_TFs_subfamily.final.FG$intersect.count_200bp+1)/(TEwide_TFs_subfamily.final.FG$mean.value+1))
TEwide_TFs_subfamily.final.FG$perC = TEwide_TFs_subfamily.final.FG$intersect.count_200bp.FG/TEwide_TFs_subfamily.final.FG$subfamily.count_200bp.FG
#write.csv(TEwide_TFs_subfamily.final.FG,file = "TEwide_TFs_subfamily.final_all_2023_9_4.csv")

###################################################
################################# group summary
ordered_subfamilies.table_final = read.csv(file = "input/TEwide_group.list_2023_9_4.functional_group_edited.csv")

ordered_subfamilies.table_final.sum = data.frame(ordered_subfamilies.table_final %>% 
                                                   group_by(Y.axis,family.group,family.cluster,functional.group) %>% 
                                                   dplyr::summarise(no.subfamilies.FG = n(),
                                                                    subfamily.count_200bp.FG = sum(subfamily.count_200bp,na.rm=T)))
ordered_subfamilies.table_final.sum.tmp = data.frame(ordered_subfamilies.table_final.sum %>% 
                                                       group_by(family.group,family.cluster) %>% 
                                                       dplyr::summarise(count.FG.family.cluster = n(),
                                                                        count.200bp.family.cluster = sum(subfamily.count_200bp.FG,na.rm=T)))
ordered_subfamilies.table_final.sum.tmp = ordered_subfamilies.table_final.sum.tmp[order(-ordered_subfamilies.table_final.sum.tmp$count.FG.family.cluster),]
ordered_subfamilies.table_final.sum.tmp$family.cluster = factor(ordered_subfamilies.table_final.sum.tmp$family.cluster,levels=ordered_subfamilies.table_final.sum.tmp$family.cluster)
ordered_subfamilies.table_final.sum.tmp2 = data.frame(ordered_subfamilies.table_final[!duplicated(ordered_subfamilies.table_final$subfamily2_1),] %>% 
                                                        group_by(family.group,family.cluster) %>% 
                                                        dplyr::summarise(count.family = n()))
ordered_subfamilies.table_final.sum.tmp2 = ordered_subfamilies.table_final.sum.tmp2[order(-ordered_subfamilies.table_final.sum.tmp2$count.family),]

ordered_subfamilies.table_final.sum = merge(ordered_subfamilies.table_final.sum,ordered_subfamilies.table_final.sum.tmp[,c("family.group","count.FG.family.cluster","count.200bp.family.cluster")],by = "family.group",all.x=T)
ordered_subfamilies.table_final.sum = merge(ordered_subfamilies.table_final.sum,ordered_subfamilies.table_final.sum.tmp2[,c("family.group","count.family")],by = "family.group",all.x=T)
ordered_subfamilies.table_final.sum$family.cluster = factor(ordered_subfamilies.table_final.sum$family.cluster,levels = ordered_subfamilies.table_final.sum.tmp$family.cluster)
ordered_subfamilies.table_final.sum = ordered_subfamilies.table_final.sum[with(ordered_subfamilies.table_final.sum, order(-count.FG.family.cluster,family.group,-subfamily.count_200bp.FG)), ]
ordered_subfamilies.table_final.sum$Y.axis = factor(ordered_subfamilies.table_final.sum$Y.axis,levels = rev(ordered_subfamilies.table_final.sum$Y.axis))
ordered_subfamilies.table_final.sum.reordered = ordered_subfamilies.table_final.sum
ordered_subfamilies.table_final.sum.reordered$Group.order = as.numeric(as.character(gsub("PG","",ordered_subfamilies.table_final.sum.reordered$functional.group)))
ordered_subfamilies.table_final.sum.reordered = ordered_subfamilies.table_final.sum.reordered[with(ordered_subfamilies.table_final.sum.reordered, order(-count.FG.family.cluster,family.group,Group.order)), ]
ordered_subfamilies.table_final.sum.reordered$Y.axis = factor(ordered_subfamilies.table_final.sum.reordered$Y.axis,levels=ordered_subfamilies.table_final.sum.reordered$Y.axis)

sum(ordered_subfamilies.table_final.sum.tmp$count.FG.family.cluster)

########## barplot: how many functional groups
p1 = ggplot(ordered_subfamilies.table_final.sum.tmp, aes(x=family.cluster,y=count.FG.family.cluster)) + 
  geom_col()+
  geom_text(aes(label=count.FG.family.cluster))+
  ylab("# functional groups")+
  xlab("")+
  #scale_fill_manual(values=c("down"="#2166ac","up"="#b2182b"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.margin=unit(c(10,10,10,10),"mm"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1),angle = 0),
        axis.title=element_text(colour="black",size=rel(1.2)),
        legend.position="right",
        legend.text=element_text(size=rel(0.8)),
        axis.text.x=element_text(angle=90,hjust = 0.95,vjust = .5),
        legend.key = element_rect(fill = NA))

### plot for the proportion of families
ordered_subfamilies.table_final.plot = ordered_subfamilies.table_final
ordered_subfamilies.table_final.plot = ordered_subfamilies.table_final.plot[!is.na(ordered_subfamilies.table_final.plot$family.count_200bp),]
ordered_subfamilies.table_final.plot.sum = data.frame(ordered_subfamilies.table_final.plot %>% 
                                                        group_by(family.group,family.cluster,subfamily2_1,family.count_200bp,functional.group) %>% 
                                                        dplyr::summarise(count.FGs = n(),
                                                                         count.instances = sum(subfamily.count_200bp,na.rm=T)))
ordered_subfamilies.table_final.plot.sum = ordered_subfamilies.table_final.plot.sum[with(ordered_subfamilies.table_final.plot.sum, order(family.group,family.cluster,subfamily2_1,-count.instances)),]
ordered_subfamilies.table_final.plot.sum$Order.FG = NA
for(i in 1:nrow(ordered_subfamilies.table_final.plot.sum)){
  if (i == 1) {
    ordered_subfamilies.table_final.plot.sum[i,]$Order.FG = 1
  } else if (ordered_subfamilies.table_final.plot.sum[i,]$subfamily2_1 != ordered_subfamilies.table_final.plot.sum[i-1,]$subfamily2_1) {
    ordered_subfamilies.table_final.plot.sum[i,]$Order.FG = 1
  } else {
    ordered_subfamilies.table_final.plot.sum[i,]$Order.FG = ordered_subfamilies.table_final.plot.sum[i-1,]$Order.FG+ 1
  }
}
ordered_subfamilies.table_final.plot.sum$perC.FG = ordered_subfamilies.table_final.plot.sum$count.instances/ordered_subfamilies.table_final.plot.sum$family.count_200bp
ordered_subfamilies.table_final.plot.sum.order = ordered_subfamilies.table_final.plot.sum[ordered_subfamilies.table_final.plot.sum$Order.FG ==1,]
ordered_subfamilies.table_final.plot.sum.order = ordered_subfamilies.table_final.plot.sum.order[order(-ordered_subfamilies.table_final.plot.sum.order$perC.FG),]
ordered_subfamilies.table_final.plot.sum$subfamily2_1 = factor(ordered_subfamilies.table_final.plot.sum$subfamily2_1,levels = ordered_subfamilies.table_final.plot.sum.order$subfamily2_1)
ordered_subfamilies.table_final.plot.sum.uniq = ordered_subfamilies.table_final.plot.sum[!duplicated(ordered_subfamilies.table_final.plot.sum[,c("family.group","subfamily2_1")]),]
ordered_subfamilies.table_final.plot.sum.different_Groups = ordered_subfamilies.table_final.plot.sum[!duplicated(ordered_subfamilies.table_final.plot.sum[,c("family.group","subfamily2_1","Order.FG")]),]
ordered_subfamilies.table_final.plot.sum.different_Groups.families = ordered_subfamilies.table_final.plot.sum.different_Groups[ordered_subfamilies.table_final.plot.sum.different_Groups$Order.FG==2,]

# number of total instances that were grouped
totalCounts.FG = sum(ordered_subfamilies.table_final.plot.sum.different_Groups[ordered_subfamilies.table_final.plot.sum.different_Groups$subfamily2_1 %in% ordered_subfamilies.table_final.plot.sum.different_Groups.families$subfamily2_1,]$count.instances)

# number of instances that were not grouped into FG1
totalCounts.noFG1 = sum(ordered_subfamilies.table_final.plot.sum.different_Groups[ordered_subfamilies.table_final.plot.sum.different_Groups$Order.FG !=1 &
                                                                                    ordered_subfamilies.table_final.plot.sum.different_Groups$subfamily2_1 %in% ordered_subfamilies.table_final.plot.sum.different_Groups.families$subfamily2_1,]$count.instances)

p2.bar = ggplot(ordered_subfamilies.table_final.plot.sum, aes(x=subfamily2_1,y = perC.FG,fill = Order.FG)) + 
  geom_bar(position="stack", stat="identity")+
  ylab("# functional groups")+
  xlab("")+
  geom_text(data = ordered_subfamilies.table_final.plot.sum.uniq,aes(x=subfamily2_1,label=family.count_200bp,y=1),angle=90)+
  #scale_fill_manual(values=c("down"="#2166ac","up"="#b2182b"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_blank(), 
        plot.margin=unit(c(10,10,10,10),"mm"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1),angle = 0),
        axis.title=element_text(colour="black",size=rel(1.2)),
        legend.position="right",
        legend.text=element_text(size=rel(0.8)),
        axis.text.x=element_text(angle=90,hjust = 0.95,vjust = .5),
        legend.key = element_rect(fill = NA))

### plot for the phyletic groups
pdf(paste("Figure_6A_",Date,".pdf",sep=""),
    width = 4,
    height = 4,
    pointsize = 10)        # smaller font size
grid.draw(p1)
dev.off()  
pdf(paste("Figure_6B_",Date,".pdf",sep=""),
    width = 10,
    height = 4,
    pointsize = 10)        # smaller font size
grid.draw(p2.bar)
dev.off()  

################################# heatmap plot
combination = "LTR5"
#Versions = c("perC","count")
Versions = c("perC")
for (combination in c("LTR5","LTR7","all")){
  Version = "perC"
  for (Version in Versions) {
    TEwide_TFs_plot = read.csv(file = paste("input/TEwide_TFs_subfamily.final_",combination,"_2023_9_4.csv",sep=""))
    TEwide_TFs_plot = merge(TEwide_TFs_plot,TF.totalPeaks[,c("fileID","cell.mark")],by.x="mark",by.y="fileID",all.x=T)
    if (combination %in% c("LTR5","LTR7")){
      if (Version == "perC") {
        TEwide_TFs_plot$is_sig = ifelse(TEwide_TFs_plot$subfamily.count_200bp < 100 & 
                                          TEwide_TFs_plot$intersect.count_200bp >=5 & 
                                          TEwide_TFs_plot$pvalue<=0.05 & 
                                          TEwide_TFs_plot$log2FC >=1,"*",NA)
        TEwide_TFs_plot$is_sig = ifelse(TEwide_TFs_plot$subfamily.count_200bp >= 100 &
                                          TEwide_TFs_plot$intersect.count_200bp>=20 &
                                          TEwide_TFs_plot$pvalue<=0.05 &
                                          TEwide_TFs_plot$log2FC>=1,"*",TEwide_TFs_plot$is_sig)
      } else {
        TEwide_TFs_plot$is_sig = ifelse(TEwide_TFs_plot$intersect.count_200bp>=5 & 
                                          TEwide_TFs_plot$pvalue<=0.05 & 
                                          TEwide_TFs_plot$log2FC>=1,"*",NA)
      }
      if (combination == "LTR5"){ 
        groupID = "G2"
        Height = 8
        Width = 24
      } else {
        groupID = "G6"
        Height = 8
        Width = 24
      }
      TEwide_TFs_plot$subfamily = factor(TEwide_TFs_plot$subfamily,levels=ordered_subfamilies.table[ordered_subfamilies.table$family.group == groupID,]$subfamily)
      TEwide_TFs_plot$Y.axis = TEwide_TFs_plot$subfamily
      TEwide_TFs_plot = TEwide_TFs_plot[!is.na(TEwide_TFs_plot$Y.axis),]
      ### 
      TEwide_TFs_plot$perC = ifelse(TEwide_TFs_plot$perC == 0,NA,TEwide_TFs_plot$perC)
      unique(TEwide_TFs_plot$mark)
      TEwide_TFs_plot$mark
      # filter by the significance
      TEwide_TFs_plot.cell.mark = unique(TEwide_TFs_plot[!is.na(TEwide_TFs_plot$is_sig),]$cell.mark)
      TEwide_TFs_plot = TEwide_TFs_plot[TEwide_TFs_plot$cell.mark %in% c(TEwide_TFs_plot.cell.mark,Roadmap.marks,iPSC.marks),]
      # adjust the log2FC
      TEwide_TFs_plot$log2FC.adj = ifelse(TEwide_TFs_plot$log2FC>=0,TEwide_TFs_plot$log2FC,0)
    } else {
      if (Version == "perC") {
        TEwide_TFs_plot$is_sig = ifelse(TEwide_TFs_plot$subfamily.count_200bp.FG < 100 & 
                                          TEwide_TFs_plot$intersect.count_200bp>=5 & 
                                          TEwide_TFs_plot$pvalue<=0.05 & 
                                          TEwide_TFs_plot$log2FC>=1,"*",NA)
        TEwide_TFs_plot$is_sig = ifelse(TEwide_TFs_plot$subfamily.count_200bp.FG >= 100 &
                                          TEwide_TFs_plot$intersect.count_200bp.FG>=20 &
                                          TEwide_TFs_plot$pvalue<=0.05 &
                                          TEwide_TFs_plot$log2FC>=1,"*",TEwide_TFs_plot$is_sig)
      } else {
        TEwide_TFs_plot$is_sig = ifelse(TEwide_TFs_plot$intersect.count_200bp>=5 & 
                                          TEwide_TFs_plot$pvalue<=0.05 & 
                                          TEwide_TFs_plot$log2FC>=1,"*",NA)
      }
      TEwide_TFs_plot$Y.axis = paste(TEwide_TFs_plot$family.group,TEwide_TFs_plot$functional.group)
      # by the group
      TEwide_TFs_plot$Y.axis = factor(TEwide_TFs_plot$Y.axis,levels=rev(ordered_subfamilies.table_final.sum.reordered$Y.axis))
      # by the count
      #TEwide_TFs_plot$Y.axis = factor(TEwide_TFs_plot$Y.axis,levels = rev(ordered_subfamilies.table_final.sum$Y.axis))
      
      TEwide_TFs_plot = TEwide_TFs_plot[!is.na(TEwide_TFs_plot$Y.axis),]
      ### 
      TEwide_TFs_plot$perC = ifelse(TEwide_TFs_plot$perC == 0,NA,TEwide_TFs_plot$perC)
      
      # filter by the significance
      TEwide_TFs_plot.cell.mark = unique(TEwide_TFs_plot[!is.na(TEwide_TFs_plot$is_sig),]$cell.mark)
      TEwide_TFs_plot = TEwide_TFs_plot[TEwide_TFs_plot$cell.mark %in% c(TEwide_TFs_plot.cell.mark,Roadmap.marks,iPSC.marks),]
      # adjust the log2FC
      TEwide_TFs_plot$log2FC.adj = ifelse(TEwide_TFs_plot$log2FC>=0,TEwide_TFs_plot$log2FC,0)
      # width
      Width = 40
      Height = 14
    }
    TEwide_TFs_plot = TEwide_TFs_plot[!TEwide_TFs_plot$mark %in% TF.totalPeaks.excluded$fileID,]
    # convert table for the cluster analysis
    
    TEwide_TFs_plot.convert = data.frame(matrix(ncol=1,nrow=length(unique(TEwide_TFs_plot$Y.axis))))
    colnames(TEwide_TFs_plot.convert)[1] = "TE"
    TEwide_TFs_plot.convert$TE = unique(as.character(TEwide_TFs_plot$Y.axis))
    
    TEwide_TFs_plot = TEwide_TFs_plot[!is.na(TEwide_TFs_plot$cell.mark),]
    Marks = unique(TEwide_TFs_plot$cell.mark)
    for (Mark in Marks[!Marks %in% c(Roadmap.marks,iPSC.marks)]){
      tmp1 = TEwide_TFs_plot[TEwide_TFs_plot$cell.mark == Mark,c("Y.axis","log2FC")]
      TEwide_TFs_plot.convert = merge(TEwide_TFs_plot.convert,tmp1,by.x = "TE",by.y = "Y.axis",all.x=T)
      colnames(TEwide_TFs_plot.convert)[which(colnames(TEwide_TFs_plot.convert) == "log2FC")] = Mark
    }
    rownames(TEwide_TFs_plot.convert) = TEwide_TFs_plot.convert$TE
    TEwide_TFs_plot.convert = TEwide_TFs_plot.convert[,-1]
    TEwide_TFs_plot.convert = as.matrix(TEwide_TFs_plot.convert)
    hclust.mark = hclust(dist(t(TEwide_TFs_plot.convert)), method="complete")
    hclust.mark_order = hclust.mark$labels[hclust.mark$order]
    hclust.mark_reorder = c(iPSC.marks,Roadmap.marks,
                            hclust.mark_order[grepl("^H1\\|histone",hclust.mark_order)],
                            hclust.mark_order[grepl("^H1\\|TF",hclust.mark_order)],
                            hclust.mark_order[grepl("^HEK293T",hclust.mark_order)])
    ###
    TEwide_TFs_plot$cell.mark = factor(TEwide_TFs_plot$cell.mark,levels = hclust.mark_reorder)
    p1 = ggplot(TEwide_TFs_plot, aes(cell.mark,Y.axis)) +
      scale_x_discrete(drop=FALSE)+
      scale_y_discrete(drop=FALSE)+
      geom_point(aes(size=perC,fill = (log2FC.adj)),shape=21,color="black")+
      geom_text(aes(label=is_sig),color="blue")+
      scale_fill_gradient2(low = "white", high = "#b2182b") +
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
    
    p2 = ggplot(TEwide_TFs_plot, aes(cell.mark,Y.axis)) +
      scale_x_discrete(drop=FALSE)+
      scale_y_discrete(drop=FALSE)+
      geom_tile(aes(fill=log2FC.adj))+
      geom_text(aes(label=is_sig),color="blue")+
      scale_fill_gradient2(low = "white", high = "#b2182b") +
      #scale_size(limits=c(0,1))+
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
            legend.position="right",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1)))
    p3 = ggplot(ordered_subfamilies.table_final.sum.reordered, aes(subfamily.count_200bp.FG,Y.axis)) +
      scale_y_discrete(drop=FALSE)+
      geom_col()+
      geom_text(aes(label=subfamily.count_200bp.FG),col="blue")+
      xlab("#")+
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
    
    if (combination == "all"){
      gA1 <- ggplotGrob(p1)
      gA2 <- ggplotGrob(p2)
      gA3 <- ggplotGrob(p3)
      g = cbind(gA1,gA2,gA3, size = "last")
      pdf(paste("Figure_6C_",combination,"-",Version,"_",Date,"-1.pdf",sep=""),    # create PNG for the heat map
          width = (Width/2)*3,        # 5 x 300 pixels
          height = Height,
          pointsize = 10)        # smaller font size
      grid.draw(g)
      dev.off() 
    } else {
      gA1 <- ggplotGrob(p1)
      gA2 <- ggplotGrob(p2)
      g = cbind(gA1,gA2, size = "last")
      pdf(paste("Figure_S11C_S11D_",combination,"-",Version,"_",Date,"-1.pdf",sep=""),    # create PNG for the heat map
          width = Width,        # 5 x 300 pixels
          height = Height,
          pointsize = 10)        # smaller font size
      grid.draw(g)
      dev.off() 
    }
  }
}

#################################
######################
TEwide_TFs_subfamily.comparison = merge(TEwide_TFs_subfamily,ordered_subfamilies.table[,c("subfamily","tree","family.cluster","functional.group","instance.cluster.name","subfamily.cluster.order","treeOrder")],by="subfamily",all.x=T)
TEwide_TFs_subfamily.comparison = TEwide_TFs_subfamily.comparison[TEwide_TFs_subfamily.comparison$mark %in% TF.totalPeaks$fileID,]
TEwide_TFs_subfamily.comparison = merge(TEwide_TFs_subfamily.comparison,TF.totalPeaks[,c("fileID","Experiment.target","cell.mark")],by.x="mark",by.y="fileID",all.x=T)
TEwide_TFs_subfamily.comparison.candidate = TEwide_TFs_subfamily.comparison[!is.na(TEwide_TFs_subfamily.comparison$family.cluster),]

# summarized by each group
TEwide_TFs_subfamily.comparison.candidate.sum1 = data.frame(TEwide_TFs_subfamily.comparison.candidate %>% 
                                                              group_by(tree,family.cluster,functional.group,cell.mark,Experiment.target) %>% 
                                                              dplyr::summarise(count.cluster = n(),
                                                                               subfamily.sum = sum(subfamily.count_200bp,na.rm=T),
                                                                               subfamily.intersect.sum = sum(intersect.count_200bp,na.rm=T)))

TEwide_TFs_subfamily.comparison.candidate.sum1$perC = TEwide_TFs_subfamily.comparison.candidate.sum1$subfamily.intersect.sum/TEwide_TFs_subfamily.comparison.candidate.sum1$subfamily.sum
TEwide_TFs_subfamily.comparison.candidate.sum1$Type = "group"
TEwide_TFs_subfamily.comparison.candidate.sum1$uniqueID = paste(TEwide_TFs_subfamily.comparison.candidate.sum1$family.cluster,
                                                                TEwide_TFs_subfamily.comparison.candidate.sum1$cell.mark,
                                                                TEwide_TFs_subfamily.comparison.candidate.sum1$Type)

# summarized by each subfamily
TEwide_TFs_subfamily.comparison.candidate.sum2 = data.frame(TEwide_TFs_subfamily.comparison.candidate %>% 
                                                              group_by(tree,family.cluster,family,cell.mark,Experiment.target) %>% 
                                                              dplyr::summarise(count.cluster = n(),
                                                                               subfamily.sum = sum(subfamily.count_200bp,na.rm=T),
                                                                               subfamily.intersect.sum = sum(intersect.count_200bp,na.rm=T)))

TEwide_TFs_subfamily.comparison.candidate.sum2$perC = TEwide_TFs_subfamily.comparison.candidate.sum2$subfamily.intersect.sum/TEwide_TFs_subfamily.comparison.candidate.sum2$subfamily.sum
TEwide_TFs_subfamily.comparison.candidate.sum2$Type = "subfamily"
TEwide_TFs_subfamily.comparison.candidate.sum2$uniqueID = paste(TEwide_TFs_subfamily.comparison.candidate.sum2$family.cluster,
                                                                TEwide_TFs_subfamily.comparison.candidate.sum2$cell.mark,
                                                                TEwide_TFs_subfamily.comparison.candidate.sum2$Type)

TEwide_TFs_subfamily.comparison.candidate.sum = rbind(TEwide_TFs_subfamily.comparison.candidate.sum1[!is.na(TEwide_TFs_subfamily.comparison.candidate.sum1$family.cluster),c("tree","family.cluster","Experiment.target","cell.mark","perC","Type","uniqueID")],
                                                      TEwide_TFs_subfamily.comparison.candidate.sum2[!is.na(TEwide_TFs_subfamily.comparison.candidate.sum2$family.cluster),c("tree","family.cluster","Experiment.target","cell.mark","perC","Type","uniqueID")])
max_ones = data.frame(TEwide_TFs_subfamily.comparison.candidate.sum %>%
                        group_by(tree,uniqueID) %>%
                        dplyr::slice(which.max(perC)))
min_ones = data.frame(TEwide_TFs_subfamily.comparison.candidate.sum %>%
                        group_by(tree,uniqueID) %>%
                        dplyr::slice(which.min(perC)))

TEwide_TFs_subfamily.comparison.candidate.sum.plot = merge(max_ones[,c("uniqueID","family.cluster","Experiment.target","cell.mark","Type","perC")],min_ones[,c("uniqueID","perC")],by="uniqueID",all=T)
TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.diff = TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.x-TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.y
TEwide_TFs_subfamily.comparison.candidate.sum.plot$uniqueID2 = paste(TEwide_TFs_subfamily.comparison.candidate.sum.plot$family.cluster,
                                                                     TEwide_TFs_subfamily.comparison.candidate.sum.plot$Experiment.target)
TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.diff = ifelse(TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.diff == 0,NA,TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.diff)
TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.x = ifelse(TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.x == 0,NA,TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.x)
TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.y = ifelse(TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.y == 0,NA,TEwide_TFs_subfamily.comparison.candidate.sum.plot$perC.y)
TEwide_TFs_subfamily.comparison.candidate.sum.plot$Experiment.target = ifelse(grepl("^E00",TEwide_TFs_subfamily.comparison.candidate.sum.plot$cell.mark),gsub("\\|histone\\|","-",TEwide_TFs_subfamily.comparison.candidate.sum.plot$cell.mark),
                                                                              TEwide_TFs_subfamily.comparison.candidate.sum.plot$Experiment.target)

TF_order = c("YY1","USF1","NANOG","EP300","TCF12","SP1","TEAD4","TBP","POLR2AphosphoS5","POLR2A","TAF1","ASH2L","E2F6","MAX","CEBPB",
             "CHD7","POU5F1","RXRA","BCL11A","KDM1A","CBX5","RAD21","CTCF","H2AFZ","ZNF143")

TF_order2 = c("E003-DNase","E004-DNase","E005-DNase","E006-DNase","E007-Nase","E003-H3K27ac","E004-H3K27ac","E005-H3K27ac","E006-H3K27ac","E007-H3K27ac",
              "H3K4me3","H3K4me2","H4K5ac","H3K27ac","H2BK5ac","H2BK12ac","H2AK5ac","H3K18ac","H2BK20ac","H2BK120ac","H3K36me3","H4K20me1","H3K79me2","H4K8ac","H3K4me1","H4K91ac","H3K9me3")
for (Type in c("histone","TF")){
  if (Type == "TF"){
    TEwide_TFs_subfamily.comparison.candidate.sum.plot2 = TEwide_TFs_subfamily.comparison.candidate.sum.plot[TEwide_TFs_subfamily.comparison.candidate.sum.plot$Experiment.target %in%TF_order,]
    TEwide_TFs_subfamily.comparison.candidate.sum.plot2$Experiment.target = factor(TEwide_TFs_subfamily.comparison.candidate.sum.plot2$Experiment.target,
                                                                                   levels=TF_order)
  } else {
    TEwide_TFs_subfamily.comparison.candidate.sum.plot2 = TEwide_TFs_subfamily.comparison.candidate.sum.plot[TEwide_TFs_subfamily.comparison.candidate.sum.plot$Experiment.target %in%TF_order2,]
    TEwide_TFs_subfamily.comparison.candidate.sum.plot2$Experiment.target = factor(TEwide_TFs_subfamily.comparison.candidate.sum.plot2$Experiment.target,
                                                                                   levels=TF_order2)
  }
  ####### plot 
  Family = "MER11A"
  glist = list()
  Order = 1
  for(Family in unique(TEwide_TFs_subfamily.comparison.candidate.sum.plot2$family.cluster)){
    p1 = ggplot(TEwide_TFs_subfamily.comparison.candidate.sum.plot2[TEwide_TFs_subfamily.comparison.candidate.sum.plot2$family.cluster == Family,], aes(x=Experiment.target,y=perC.diff*100,group=Type)) +
      geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
      geom_text(aes(label=round(perC.diff*100,1)),  color="black",
                position = position_dodge(0.9), size=4,vjust=-0.2, hjust=-0.2,angle=45)+
      scale_fill_manual(values =c("#016A70","#A2C579"))+
      ggtitle(Family)+
      #geom_hline(yintercept=c(20,40), linetype="dashed", color = "black")+
      #ylim(0,100)+
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
            legend.position="right",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1))) 
    glist[[Order]] <- ggplotGrob(p1)
    Order = Order + 1
    p1 = ggplot(TEwide_TFs_subfamily.comparison.candidate.sum.plot2[TEwide_TFs_subfamily.comparison.candidate.sum.plot2$family.cluster == Family,], aes(x=Experiment.target,y=perC.x*100,group=Type)) +
      geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
      geom_text(aes(label=round(perC.x*100,1)), color="black",
                position = position_dodge(0.9), size=4,vjust=-0.2, hjust=-0.2,angle=45)+
      scale_fill_manual(values =c("#016A70","#A2C579"))+
      ggtitle(Family)+
      #geom_hline(yintercept=c(20,40), linetype="dashed", color = "black")+
      #ylim(0,100)+
      #scale_fill_gradient2(low="white",high = "#e41a1c",na.value = 'black')+
      scale_x_discrete(drop = FALSE)+
      ylab("max % peaks-associated instances)")+
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
            legend.position="right",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1))) 
    
    glist[[Order]] <- ggplotGrob(p1)
    Order = Order + 1
    p1 = ggplot(TEwide_TFs_subfamily.comparison.candidate.sum.plot2[TEwide_TFs_subfamily.comparison.candidate.sum.plot2$family.cluster == Family,], aes(x=Experiment.target,y=perC.y*100,group=Type)) +
      geom_bar(aes(fill=Type),stat="identity", position=position_dodge()) + 
      geom_text(aes(label=round(perC.y*100,1)), vjust=-.2, color="black",
                position = position_dodge(0.9), size=4,vjust=-0.2, hjust=-0.2,angle=45)+
      scale_fill_manual(values =c("#016A70","#A2C579"))+
      ggtitle(Family)+
      scale_x_discrete(drop = FALSE)+
      ylab("min % peaks-associated instances)")+
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
            legend.position="right",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1))) 
    glist[[Order]] <- ggplotGrob(p1)
    Order = Order + 1
  }
  
  pdf(paste("Figure_6D-",Type,"-",Date,"-1.pdf",sep=""),    # create PNG for the heat map        
      width = 36,        # 5 x 300 pixels
      height = 90,
      pointsize = 10)        # smaller font size
  do.call("grid.arrange",c(glist,ncol=3))
  dev.off()
}


