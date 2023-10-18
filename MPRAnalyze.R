## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load and examine data, include=TRUE--------------------------------------
#BiocManager::install("MPRAnalyze")

#memory.limit(size=16000)

#data("ChrEpi")
#summary(ce.colAnnot)
#head(ce.rnaCounts)
library(MPRAnalyze)
library(splitstackshape)
library(gplots)
library(ggplot2)
library(grid)
library(gridExtra)

###
setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/Project_Neurogenesis/Final_edited_version_2022_11_25/Final_scripts/")

################################################# step 1 prepare our own data
sampleListDNA = c("iPSC_1_insert_BCs.DNAcounts","iPSC_2_insert_BCs.DNAcounts","iPSC_3_insert_BCs.DNAcounts",
                  "NPC_1_insert_BCs.DNAcounts","NPC_2_insert_BCs.DNAcounts","NPC_3_insert_BCs.DNAcounts")
sampleListRNA = c("iPSC_1_insert_BCs.RNAcounts","iPSC_2_insert_BCs.RNAcounts","iPSC_3_insert_BCs.RNAcounts",
                  "NPC_1_insert_BCs.RNAcounts","NPC_2_insert_BCs.RNAcounts","NPC_3_insert_BCs.RNAcounts")
### DNA count table 
ce.dnaCounts_TE = data.frame("tmp1"=NA,"tmp2"=NA)
for (sampleID in sampleListDNA) {
  inputTmp = read.delim(paste("input/",sampleID,sep=""),sep="\t",header=T)
  sampleID_tmp = gsub("_insert_BCs.DNAcounts|_insert_BCs.RNAcounts","",sampleID)
  sampleID_tmp = gsub("_",".",sampleID_tmp)
  colnames(inputTmp) = gsub("^X",paste(sampleID_tmp,".",sep=""),colnames(inputTmp))
  if (nrow(ce.dnaCounts_TE)==1){
    ce.dnaCounts_TE = inputTmp
  } else {
    inputTmp = inputTmp[,2:ncol(inputTmp)]
    ce.dnaCounts_TE = cbind(ce.dnaCounts_TE,inputTmp)
  }
}
rownames(ce.dnaCounts_TE) = ce.dnaCounts_TE$insert
ce.dnaCounts_TE = ce.dnaCounts_TE[,-1]

### RNA count table
ce.rnaCounts_TE = data.frame("tmp1"=NA,"tmp2"=NA)
for (sampleID in sampleListRNA) {
  inputTmp = read.delim(paste("input/",sampleID,sep=""),sep="\t",header=T)
  sampleID_tmp = gsub("_insert_BCs.DNAcounts|_insert_BCs.RNAcounts","",sampleID)
  sampleID_tmp = gsub("_",".",sampleID_tmp)
  colnames(inputTmp) = gsub("^X",paste(sampleID_tmp,".",sep=""),colnames(inputTmp))
  if (nrow(ce.rnaCounts_TE)==1){
    ce.rnaCounts_TE = inputTmp
  } else {
    inputTmp = inputTmp[,2:ncol(inputTmp)]
    ce.rnaCounts_TE = cbind(ce.rnaCounts_TE,inputTmp)
  }
}
rownames(ce.rnaCounts_TE) = ce.rnaCounts_TE$insert
ce.rnaCounts_TE = ce.rnaCounts_TE[,-1]

#### annotation table 
ce.colAnnot_TE = data.frame(colnames(ce.dnaCounts_TE))
colnames(ce.colAnnot_TE)[1] = "Info"
ce.colAnnot_TE$rowName = ce.colAnnot_TE$Info
rownames(ce.colAnnot_TE) = ce.colAnnot_TE$Info
ce.colAnnot_TE = data.frame(cSplit(ce.colAnnot_TE,"Info",sep=".",type.convert = as.character))
ce.colAnnot_TE = ce.colAnnot_TE[ce.colAnnot_TE$Info_1 !="insert",]
rownames(ce.colAnnot_TE) = ce.colAnnot_TE$rowName
head(ce.colAnnot_TE)
colnames(ce.colAnnot_TE) = c("rowName","condition","batch","barcode") 
ce.colAnnot_TE = ce.colAnnot_TE[,c("batch","condition","barcode")]
ce.colAnnot_TE$batch = factor(ce.colAnnot_TE$batch)
ce.colAnnot_TE$condition = factor(ce.colAnnot_TE$condition)
ce.colAnnot_TE$barcode = factor(ce.colAnnot_TE$barcode)

#### list of control
ce.control_TE = ifelse(grepl("ctrl",rownames(ce.dnaCounts_TE)),TRUE,FALSE)

##################################
################################## step 2, normalization final version
## step 2.1 ----init object, include=TRUE------------------------------------------------
obj_TE <- MpraObject(dnaCounts = as.matrix(ce.dnaCounts_TE), rnaCounts = as.matrix(ce.rnaCounts_TE), 
                 dnaAnnot = ce.colAnnot_TE, rnaAnnot = ce.colAnnot_TE, 
                 controls = ce.control_TE)
# save(obj_TE, file = "input/MPRAnalyze.obj_TE_2022_2_2")
# load(file = "input/MPRAnalyze.obj_TE_2022_2_2")

## step 2.2 normalization same factor
## In this case, the factors are the same - each combination of batch and 
## condition is a single library, and we'll use the default estimation
obj_TE <- estimateDepthFactors(obj_TE, lib.factor = c("batch", "condition"),
                            which.lib = "both",
                            depth.estimator = "uq")
# save(obj_TE, file = "input/MPRAnalyze.obj_TE.depthfactors_2022_2_2")
# load(file = "input/MPRAnalyze.obj_TE.depthfactors_2022_2_2")

## step 2.3 fit the quant model 
obj_TE <- analyzeQuantification(obj = obj_TE,
                            dnaDesign = ~ batch + condition,
                            rnaDesign = ~ batch + condition)
# save(obj_TE, file = "input/MPRAnalyze.obj_TE.quantmodel_2022_2_2")
# load(file = "input/MPRAnalyze.obj_TE.quantmodel_2022_2_2")

## ----quant extract and viz
##extract alpha values from the fitted model
alpha <- getAlpha(obj_TE, by.factor = "condition")
write.csv(alpha,file="input/iPSC_NPC_combined_alphaValue_2021_12_15.tsv")
