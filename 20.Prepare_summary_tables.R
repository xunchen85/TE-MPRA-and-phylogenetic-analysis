#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

library(splitstackshape)
library(dplyr)
library(ggplot2)
library(gplots)
library(grid)
library(gridExtra)
library(ggridges)

options(scipen=999)

########################### Table 1: instance info
div.hg19 = read.delim("input/hg19_rmsk_TE_1bp.div.bed.gz",header=F,sep="")
div.hg19$Instance_coordinate = paste("hg19:",div.hg19$V1,":",div.hg19$V2-1,"-",div.hg19$V3,sep="")
div.macFas5_20140131 = read.delim("input/macFas5_rmsk_TE_1bp.div.bed.gz",header=F,sep="")
div.macFas5_20140131$Instance_coordinate = paste("macFas5:",div.macFas5_20140131$V1,":",div.macFas5_20140131$V2-1,"-",div.macFas5_20140131$V3,sep="")
div.panTro4_20140131 = read.delim("input/panTro4_rmsk_TE_1bp.div.bed.gz",header=F,sep="")
div.panTro4_20140131$Instance_coordinate = paste("panTro4:",div.panTro4_20140131$V1,":",div.panTro4_20140131$V2-1,"-",div.panTro4_20140131$V3,sep="")

# load tables
Summary_Table1 = read.delim("input/TE_full_length.rename_2022_6_25",sep="",header=F)
colnames(Summary_Table1) = c("instanceID","instanceID.renamed")
Summary_Table1.tmp = read.delim("input/TE_instances_withFrames_2022_8_24.bed",sep="",header=F)
colnames(Summary_Table1.tmp) = c("fileName","chr","start","end","instanceID","anno","strand")

# is used for MPRA
Summary_Table1.tmp.F1 = Summary_Table1.tmp[grepl("_F1_",Summary_Table1.tmp$fileName),]
Summary_Table1.tmp.F2 = Summary_Table1.tmp[grepl("_F2_",Summary_Table1.tmp$fileName),]
colnames(Summary_Table1.tmp.F1)[1] = "uniqueID_F1"
colnames(Summary_Table1.tmp.F2)[1] = "uniqueID_F2"

Summary_Table1$frames_tested = ifelse(Summary_Table1$instanceID %in% Summary_Table1.tmp$instanceID,"Tested",NA)
Summary_Table1$F1_tested = ifelse(Summary_Table1$instanceID %in% Summary_Table1.tmp.F1$instanceID,"Tested",NA)
Summary_Table1$F2_tested = ifelse(Summary_Table1$instanceID %in% Summary_Table1.tmp.F2$instanceID,"Tested",NA)
Summary_Table1 = merge(Summary_Table1,Summary_Table1.tmp.F1[,c("instanceID","uniqueID_F1")],by="instanceID",all.x=T)
Summary_Table1 = merge(Summary_Table1,Summary_Table1.tmp.F2[,c("instanceID","uniqueID_F2")],by="instanceID",all.x=T)

rm(Summary_Table1.tmp,Summary_Table1.tmp.F1,Summary_Table1.tmp.F2)

# add species and family info
Summary_Table1$info = Summary_Table1$instanceID
Summary_Table1$info = gsub("-",":",Summary_Table1$info)
Summary_Table1 = data.frame(cSplit(Summary_Table1,"info",sep=":",type.convert = as.character))
colnames(Summary_Table1)[8:12] = c("species","TEfamily","chr","start","end")
Summary_Table1$Instance_coordinate = paste(Summary_Table1$chr,":",Summary_Table1$start,"-",Summary_Table1$end,sep="")

# add the div info
Summary_Table1$Instance_coordinate_species = paste(Summary_Table1$species,":",Summary_Table1$chr,":",Summary_Table1$start,"-",Summary_Table1$end,sep="")
div.hg19.candidate = div.hg19[div.hg19$Instance_coordinate %in% Summary_Table1$Instance_coordinate_species,]
div.macFas5_20140131.candidate = div.macFas5_20140131[div.macFas5_20140131$Instance_coordinate %in% Summary_Table1$Instance_coordinate_species,]
div.panTro4_20140131.candidate = div.panTro4_20140131[div.panTro4_20140131$Instance_coordinate %in% Summary_Table1$Instance_coordinate_species,]
div.combined = rbind(div.hg19.candidate,div.macFas5_20140131.candidate,div.panTro4_20140131.candidate)
Summary_Table1[!Summary_Table1$Instance_coordinate_species %in% div.combined$Instance_coordinate,]
rm(div.hg19.candidate,div.macFas5_20140131.candidate,div.panTro4_20140131.candidate)
div.combined = cSplit(div.combined,"V4",sep=":",type.convert = as.character)
div.combined = div.combined[,c("Instance_coordinate","V4_2","V7")]
colnames(div.combined) = c("Instance_coordinate_species","TEinstance","div_rate")
Summary_Table1 = merge(Summary_Table1,div.combined,by="Instance_coordinate_species",all=T)
rm(div.hg19,div.macFas5_20140131,div.panTro4_20140131,div.combined)
Summary_Table1$Instance_len = as.numeric(as.character(Summary_Table1$end)) - as.numeric(as.character(Summary_Table1$start))

# add consensus info
Consensus_families = c("MER11A","MER11B","MER11C","MER11D",
                       "MER34","MER34A","MER34A1","MER34B","MER34C","MER34C2","MER34C_","MER34D",
                       "MER52A","MER52C","MER52D")
rowID = nrow(Summary_Table1)
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),] = NA
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$instanceID = Consensus_families
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$instanceID.renamed = Consensus_families
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$Instance_coordinate_species = paste("Ancient:",Consensus_families,sep="")
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$frames_tested = "Tested"
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$F1_tested = "Tested"
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$F2_tested = "Tested"
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$species = "Ancient"
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$TEfamily = Consensus_families
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$div_rate = NA
Summary_Table1[(rowID+1):(rowID+length(Consensus_families)),]$Instance_len = c(1266,1236,1071,897,
                                                                               542,571,587,565,551,555,585,579,
                                                                               1755,1278,2123)
# write to table
write.csv(Summary_Table1,file = "Summary_Table1_2022_8_9.csv")

########################## Table 2: MPRA results
# 2.1 observed barcodes in association experiment
minDNA.Counts = 10
#
MPRA_experiment = read.delim("input/Combined.counts_2022_8_25.gz",sep="",header=T) 
MPRA_experiment$DNAcount = as.numeric(as.character(MPRA_experiment$DNAcount))
MPRA_experiment$RNAcount = as.numeric(as.character(MPRA_experiment$RNAcount))
MPRA_experiment$TotalBCperInsertOriginal = as.numeric(as.character(MPRA_experiment$TotalBCperInsertOriginal))
MPRA_experiment = MPRA_experiment[!is.na(MPRA_experiment$DNAcount),]
colnames(MPRA_experiment)[1] = "sampleID"
MPRA_experiment$sampleID = gsub("_insert_BCs.counts","",MPRA_experiment$sampleID)
#
MPRA_experiment.DNA.sum = data.frame(MPRA_experiment[MPRA_experiment$DNAcount>0,] %>% group_by(sampleID,Insert,TotalBCperInsertOriginal) 
                                     %>% dplyr::summarise(n = n()))
#
MPRA_experiment.DNA.sum.iPSC_1 = MPRA_experiment.DNA.sum[MPRA_experiment.DNA.sum$sampleID == "iPSC_1",c("Insert","n")]
colnames(MPRA_experiment.DNA.sum.iPSC_1) = c("Insert","bc_DNA.iPSC_1")
MPRA_experiment.DNA.sum.iPSC_2 = MPRA_experiment.DNA.sum[MPRA_experiment.DNA.sum$sampleID == "iPSC_2",c("Insert","n")]
colnames(MPRA_experiment.DNA.sum.iPSC_2) = c("Insert","bc_DNA.iPSC_2")
MPRA_experiment.DNA.sum.iPSC_3 = MPRA_experiment.DNA.sum[MPRA_experiment.DNA.sum$sampleID == "iPSC_3",c("Insert","n")]
colnames(MPRA_experiment.DNA.sum.iPSC_3) = c("Insert","bc_DNA.iPSC_3")

MPRA_experiment.DNA.sum.NPC_1 = MPRA_experiment.DNA.sum[MPRA_experiment.DNA.sum$sampleID == "NPC_1",c("Insert","n")]
colnames(MPRA_experiment.DNA.sum.NPC_1) = c("Insert","bc_DNA.NPC_1")
MPRA_experiment.DNA.sum.NPC_2 = MPRA_experiment.DNA.sum[MPRA_experiment.DNA.sum$sampleID == "NPC_2",c("Insert","n")]
colnames(MPRA_experiment.DNA.sum.NPC_2) = c("Insert","bc_DNA.NPC_2")
MPRA_experiment.DNA.sum.NPC_3 = MPRA_experiment.DNA.sum[MPRA_experiment.DNA.sum$sampleID == "NPC_3",c("Insert","n")]
colnames(MPRA_experiment.DNA.sum.NPC_3) = c("Insert","bc_DNA.NPC_3")

MPRA_experiment.DNA.sum.final  = MPRA_experiment[!duplicated(MPRA_experiment[,c(2,3)]),c("Insert","TotalBCperInsertOriginal")]
colnames(MPRA_experiment.DNA.sum.final) = c("Insert","bc_asso")
MPRA_experiment.DNA.sum.final = merge(MPRA_experiment.DNA.sum.final,MPRA_experiment.DNA.sum.iPSC_1,by="Insert",all=T)
MPRA_experiment.DNA.sum.final = merge(MPRA_experiment.DNA.sum.final,MPRA_experiment.DNA.sum.iPSC_2,by="Insert",all=T)
MPRA_experiment.DNA.sum.final = merge(MPRA_experiment.DNA.sum.final,MPRA_experiment.DNA.sum.iPSC_3,by="Insert",all=T)
MPRA_experiment.DNA.sum.final = merge(MPRA_experiment.DNA.sum.final,MPRA_experiment.DNA.sum.NPC_1,by="Insert",all=T)
MPRA_experiment.DNA.sum.final = merge(MPRA_experiment.DNA.sum.final,MPRA_experiment.DNA.sum.NPC_2,by="Insert",all=T)
MPRA_experiment.DNA.sum.final = merge(MPRA_experiment.DNA.sum.final,MPRA_experiment.DNA.sum.NPC_3,by="Insert",all=T)

# 2.2 combined ratio in DNA/RNA experiment
# 2.1 observed barcodes in association experiment
MPRA_ratio = read.delim("input/Combined.ratio_2022_8_25.gz",sep="",header=T) 
MPRA_ratio$dna_count = as.numeric(as.character(MPRA_ratio$dna_count))
MPRA_ratio$rna_count = as.numeric(as.character(MPRA_ratio$rna_count))
MPRA_ratio$ratio = as.numeric(as.character(MPRA_ratio$ratio))
MPRA_ratio$log2 = as.numeric(as.character(MPRA_ratio$log2))
MPRA_ratio$n_obs_bc = as.numeric(as.character(MPRA_ratio$n_obs_bc))
colnames(MPRA_ratio)[c(1,2)] = c("sampleID","Insert")
MPRA_ratio = MPRA_ratio[!is.na(MPRA_ratio$ratio),]
MPRA_ratio$uniqueName = paste(MPRA_ratio$sampleID,MPRA_ratio$Insert)
# only keep ratio with minimum DNA counts
MPRA_ratio.minCounts = MPRA_ratio
MPRA_ratio.iPSC_1 = MPRA_ratio.minCounts[MPRA_ratio.minCounts$sampleID == "iPSC_1",c("Insert","ratio")]
colnames(MPRA_ratio.iPSC_1) = c("Insert","ratio.iPSC_1")
MPRA_ratio.iPSC_2 = MPRA_ratio.minCounts[MPRA_ratio.minCounts$sampleID == "iPSC_2",c("Insert","ratio")]
colnames(MPRA_ratio.iPSC_2) = c("Insert","ratio.iPSC_2")
MPRA_ratio.iPSC_3 = MPRA_ratio.minCounts[MPRA_ratio.minCounts$sampleID == "iPSC_3",c("Insert","ratio")]
colnames(MPRA_ratio.iPSC_3) = c("Insert","ratio.iPSC_3")
MPRA_ratio.NPC_1 = MPRA_ratio.minCounts[MPRA_ratio.minCounts$sampleID == "NPC_1",c("Insert","ratio")]
colnames(MPRA_ratio.NPC_1) = c("Insert","ratio.NPC_1")
MPRA_ratio.NPC_2 = MPRA_ratio.minCounts[MPRA_ratio.minCounts$sampleID == "NPC_2",c("Insert","ratio")]
colnames(MPRA_ratio.NPC_2) = c("Insert","ratio.NPC_2")
MPRA_ratio.NPC_3 = MPRA_ratio.minCounts[MPRA_ratio.minCounts$sampleID == "NPC_3",c("Insert","ratio")]
colnames(MPRA_ratio.NPC_3) = c("Insert","ratio.NPC_3")

MPRA_ratio.final = merge(MPRA_ratio.iPSC_1,MPRA_ratio.iPSC_2,by="Insert",all=T)
MPRA_ratio.final = merge(MPRA_ratio.final,MPRA_ratio.iPSC_3,by="Insert",all=T)
MPRA_ratio.final = merge(MPRA_ratio.final,MPRA_ratio.NPC_1,by="Insert",all=T)
MPRA_ratio.final = merge(MPRA_ratio.final,MPRA_ratio.NPC_2,by="Insert",all=T)
MPRA_ratio.final = merge(MPRA_ratio.final,MPRA_ratio.NPC_3,by="Insert",all=T)

MPRA_ratio.final[!MPRA_ratio.final$Insert %in% MPRA_experiment.DNA.sum.final$Insert,]

### there are few inserts have 0 DNA counts in the MPRA table
MPRA_results = merge(MPRA_experiment.DNA.sum.final,MPRA_ratio.final,by="Insert",all=T)
MPRA_results = MPRA_results[MPRA_results$Insert != "no_BC",]

# count how many replicates have more than 10 barcodes
MPRA_results$minCount_10.DNA.iPSC = apply(MPRA_results[,grepl("^bc_DNA.iPSC",colnames(MPRA_results))],1,function(a){sum(a>=minDNA.Counts &!is.na(a))})
MPRA_results$minCount_10.DNA.NPC = apply(MPRA_results[,grepl("^bc_DNA.NPC",colnames(MPRA_results))],1,function(a){sum(a>=minDNA.Counts &!is.na(a))})

MPRA_results$average.ratio.iPSC = apply(MPRA_results[,grepl("^ratio.iPSC",colnames(MPRA_results))],1,mean,na.rm=T)
MPRA_results$average.ratio.NPC = apply(MPRA_results[,grepl("^ratio.NPC",colnames(MPRA_results))],1,mean,na.rm=T)

## exclude inaccurate ratio values
MPRA_results$ratio.iPSC_1 = ifelse(MPRA_results$bc_DNA.iPSC_1>=minDNA.Counts & !is.na(MPRA_results$bc_DNA.iPSC_1),MPRA_results$ratio.iPSC_1,NA)
MPRA_results$ratio.iPSC_2 = ifelse(MPRA_results$bc_DNA.iPSC_2>=minDNA.Counts & !is.na(MPRA_results$bc_DNA.iPSC_2),MPRA_results$ratio.iPSC_2,NA)
MPRA_results$ratio.iPSC_3 = ifelse(MPRA_results$bc_DNA.iPSC_3>=minDNA.Counts & !is.na(MPRA_results$bc_DNA.iPSC_3),MPRA_results$ratio.iPSC_3,NA)
MPRA_results$ratio.NPC_1 = ifelse(MPRA_results$bc_DNA.NPC_1>=minDNA.Counts & !is.na(MPRA_results$bc_DNA.NPC_1),MPRA_results$ratio.NPC_1,NA)
MPRA_results$ratio.NPC_2 = ifelse(MPRA_results$bc_DNA.NPC_2>=minDNA.Counts & !is.na(MPRA_results$bc_DNA.NPC_2),MPRA_results$ratio.NPC_2,NA)
MPRA_results$ratio.NPC_3 = ifelse(MPRA_results$bc_DNA.NPC_3>=minDNA.Counts & !is.na(MPRA_results$bc_DNA.NPC_3),MPRA_results$ratio.NPC_3,NA)
MPRA_results$iPSC.activity.mean = apply(MPRA_results[,grepl("^ratio.iPSC",colnames(MPRA_results))],1,mean,na.rm=T)
MPRA_results$NPC.activity.mean = apply(MPRA_results[,grepl("^ratio.NPC",colnames(MPRA_results))],1,mean,na.rm=T)

##
MPRA_results$is_kept = ifelse(MPRA_results$bc_asso>=minDNA.Counts & MPRA_results$minCount_10.DNA.iPSC>=2 & MPRA_results$minCount_10.DNA.NPC>=2,"highQuality.iPSC/NPC",NA)
MPRA_results$is_kept = ifelse(MPRA_results$bc_asso>=minDNA.Counts & (MPRA_results$minCount_10.DNA.iPSC>=2 & MPRA_results$minCount_10.DNA.NPC<2),"highQuality.iPSC",MPRA_results$is_kept)
MPRA_results$is_kept = ifelse(MPRA_results$bc_asso>=minDNA.Counts & (MPRA_results$minCount_10.DNA.iPSC<2 & MPRA_results$minCount_10.DNA.NPC>=2),"highQuality.NPC",MPRA_results$is_kept)
MPRA_results$is_kept = ifelse(MPRA_results$bc_asso<minDNA.Counts | (MPRA_results$minCount_10.DNA.iPSC<2 & MPRA_results$minCount_10.DNA.NPC<2),"lowQuality",MPRA_results$is_kept)

# alpha value
Alpha_value = read.csv("input/iPSC_NPC_combined_alphaValue_2021_12_15.tsv")
colnames(Alpha_value) = c("uniqueID","iPSC.alpha","NPC.alpha")

# log2FC
log2FC = read.csv("input/iPSC_NPC_combined_logFC_2022_3_20.tsv")

# 
Summary_Table2.tmp = merge(Alpha_value,MPRA_results[,c("Insert","is_kept","iPSC.activity.mean","NPC.activity.mean")],by.x="uniqueID",by.y="Insert",all=T)
Summary_Table2.tmp = merge(Summary_Table2.tmp,log2FC[,2:9],by="uniqueID",by.y="Insert",all=T)
Summary_Table2.tmp.plot = merge(Alpha_value,MPRA_results,by.x="uniqueID",by.y="Insert",all=T)

# combined into a table
Summary_Table2 = read.delim("input/MPRA_combined_final_2022_8_25.tsv",header=T,sep="\t")
Summary_Table2.final = merge(Summary_Table2,Summary_Table2.tmp,by.x="uniqueID_MPRA",by.y="uniqueID",all.x=T)

# normalize the activity score using MAD 2023/1/4
Summary_Table2.final.negative = Summary_Table2.final[Summary_Table2.final$Family == "Negative",]
Summary_Table2.final.negative.iPSC.alpha = Summary_Table2.final.negative[Summary_Table2.final.negative$is_kept %in% c("highQuality.iPSC/NPC","highQuality.iPSC"),]$iPSC.alpha
Summary_Table2.final.negative.iPSC.alpha.mad = mad(Summary_Table2.final.negative.iPSC.alpha,center = median(Summary_Table2.final.negative.iPSC.alpha), constant = 1.4826, na.rm = TRUE)
Summary_Table2.final$iPSC.alpha.Zscore = (Summary_Table2.final$iPSC.alpha - median(Summary_Table2.final.negative[Summary_Table2.final.negative$is_kept %in% c("highQuality.iPSC/NPC","highQuality.iPSC"),]$iPSC.alpha))/Summary_Table2.final.negative.iPSC.alpha.mad
Summary_Table2.final$iPSC.alpha.Zscore.pvalue = pnorm(Summary_Table2.final$iPSC.alpha.Zscore,
                                                      mean=mean(Summary_Table2.final.negative.iPSC.alpha),
                                                      sd = sd(Summary_Table2.final.negative.iPSC.alpha),lower.tail = FALSE)
Summary_Table2.final.iPSC.padj = Summary_Table2.final[Summary_Table2.final$Family != "Negative" & !is.na(Summary_Table2.final$iPSC.alpha.Zscore.pvalue) & 
                                                        (Summary_Table2.final$is_kept %in% c("highQuality.iPSC/NPC","highQuality.iPSC")),]
Summary_Table2.final.iPSC.padj$iPSC.alpha.Zscore.padj = p.adjust(Summary_Table2.final.iPSC.padj$iPSC.alpha.Zscore.pvalue,
                                                                 method = "fdr",
                                                                 n=nrow(Summary_Table2.final.iPSC.padj))
Summary_Table2.final.iPSC.padj$iPSC.alpha.Zscore.isActive = ifelse(Summary_Table2.final.iPSC.padj$iPSC.alpha.Zscore.padj <= 0.05,"Active","None")

# NPC
Summary_Table2.final.negative.NPC.alpha = Summary_Table2.final.negative[Summary_Table2.final.negative$is_kept %in% c("highQuality.iPSC/NPC","highQuality.NPC"),]$NPC.alpha
Summary_Table2.final.negative.NPC.alpha.mad = mad(Summary_Table2.final.negative.NPC.alpha,center = median(Summary_Table2.final.negative.NPC.alpha), constant = 1.4826, na.rm = TRUE)
Summary_Table2.final$NPC.alpha.Zscore = (Summary_Table2.final$NPC.alpha - median(Summary_Table2.final.negative[Summary_Table2.final.negative$is_kept %in% c("highQuality.iPSC/NPC","highQuality.NPC"),]$NPC.alpha))/Summary_Table2.final.negative.NPC.alpha.mad
Summary_Table2.final$NPC.alpha.Zscore.pvalue = pnorm(Summary_Table2.final$NPC.alpha.Zscore,
                                                     mean=mean(Summary_Table2.final.negative.NPC.alpha),
                                                     sd = sd(Summary_Table2.final.negative.NPC.alpha),lower.tail = FALSE)

Summary_Table2.final.NPC.padj = Summary_Table2.final[Summary_Table2.final$Family != "Negative" & !is.na(Summary_Table2.final$NPC.alpha.Zscore.pvalue) &
                                                       (Summary_Table2.final$is_kept %in% c("highQuality.iPSC/NPC","highQuality.iPSC")),]
Summary_Table2.final.NPC.padj$NPC.alpha.Zscore.padj = p.adjust(Summary_Table2.final.NPC.padj$NPC.alpha.Zscore.pvalue,
                                                               method = "fdr",
                                                               n=nrow(Summary_Table2.final.NPC.padj))
Summary_Table2.final.NPC.padj$NPC.alpha.Zscore.isActive = ifelse(Summary_Table2.final.NPC.padj$NPC.alpha.Zscore.padj <= 0.05,"Active","None")

# Combined with the table
Summary_Table2.final = merge(Summary_Table2.final,Summary_Table2.final.iPSC.padj[!duplicated(Summary_Table2.final.iPSC.padj$uniqueID_MPRA),c("uniqueID_MPRA","iPSC.alpha.Zscore.padj","iPSC.alpha.Zscore.isActive")],by.x="uniqueID_MPRA",all.x=T)
Summary_Table2.final = merge(Summary_Table2.final,Summary_Table2.final.NPC.padj[!duplicated(Summary_Table2.final.NPC.padj$uniqueID_MPRA),c("uniqueID_MPRA","NPC.alpha.Zscore.padj","NPC.alpha.Zscore.isActive")],by.x="uniqueID_MPRA",all.x=T)
# write to table
#write.csv(Summary_Table2.final,file = "Summary_Table2_2022_8_9.csv")

