library(gplots)
library(ggplot2)
library(grid)
library(dplyr)
library(splitstackshape)
library(gridExtra)
library(fastcluster)
library(RColorBrewer)

# Z-transformation
ztran <- function(x, na.rm = TRUE) {
  mns <- colMeans(x, na.rm = na.rm)
  sds <- apply(x, 2, sd, na.rm = na.rm)
  x <- sweep(x, 2, mns, "-")
  x <- sweep(x, 2, sds, "/")
  x
}

### path
setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/Project_Neurogenesis/Final_edited_version_2022_11_25/Final_scripts/")

### file names

TEfamily = "MER11A"                               ### family name
File = "hg19_MER11A_1bp.bed.group"                ### combined .group matrix file
infoFile="metadata.sampleList"            ### no header, but the first two columns should be sampleID followed by group
accessibleTEFile = "Neuro_TE_enrichment_2022_7_12.TE2.bed"          ### can be the entire TE2.bed file
ConsensusLenFile = "hg19.fa.align.seq.len"           ### length of consensus file
is_accessible = "Yes"

TEfamilies = c("MER11A","MER11B","MER11C","MER11D",
               "MER52A","MER52C","MER52D",
               "MER34","MER34A","MER34A1","MER34B","MER34C","MER34C2","MER34C_","MER34D")
for (TEfamily in TEfamilies){
  File = paste("hg19_",TEfamily,"_1bp.bed.group",sep="")                 ### combined .group matrix file
  
  ### Input files (group)
  Input_peak= read.delim(paste("input/",File,sep=""),header=T,sep="")
  colnames(Input_peak) = gsub(".bigWig","",colnames(Input_peak))
  
  ### sample info table with group info
  sampleInfo = read.delim(paste("input/",infoFile,sep=""),header=F,sep="")
  colnames(sampleInfo)[1:2] = c("sampleID","group")
  sampleInfo$group = factor(sampleInfo$group,levels=unique(sampleInfo$group))
  sampleInfo_sum = data.frame(sampleInfo %>% group_by(group) %>% dplyr::count(sort = TRUE))
  
  ### accessible instances per sample
  accessibleTEs = read.delim(paste("input/",accessibleTEFile,sep=""),header=F,sep="")
  accessibleTEs = data.frame(cSplit(accessibleTEs[,1:10],"V5",sep=":",type.convert = as.character))
  accessibleTEs$V1 = gsub("_hg19.overlap_summit.|_mm10.overlap_summit.|.TE2.bed|_summit."," ",accessibleTEs$V1)
  accessibleTEs = data.frame(cSplit(accessibleTEs,"V1",sep=" ",type.convert = as.character))
  colnames(accessibleTEs) = c("chr","start","end","score","strand","chr_summit","start_summit","end_summit","TEfamily","TEinstance","TEsuperfamily","TEsubclass","sampleID","anno")
  
  # add ref posi
  accessibleTEs$posi_ref = paste(accessibleTEs$chr_summit,accessibleTEs$start_summit,sep=":")
  # add consensus posi
  Input_peak_tmp = Input_peak[Input_peak$posi_ref %in% accessibleTEs$posi_ref,] 
  head(Input_peak_tmp)
  accessibleTEs = merge(accessibleTEs,Input_peak_tmp[,c("posi_ref","TEconsensus","posi_consensus")],by = "posi_ref",all.x=T)
  accessibleTEs[is.na(accessibleTEs$posi_consensus),]
  # add group info
  accessibleTEs_plot = merge(accessibleTEs,sampleInfo,by="sampleID",all.x=T)
  # select candidate family
  accessibleTEs_plot = accessibleTEs_plot[accessibleTEs_plot$TEfamily == TEfamily,]
  # select candidate samples
  accessibleTEs_plot = accessibleTEs_plot[accessibleTEs_plot$sampleID %in% sampleInfo$sampleID,]
  # remove samples without group info
  accessibleTEs_plot = accessibleTEs_plot[!is.na(accessibleTEs_plot$group),]
  ### obtain consensus positions for peak summits
  accessibleTEs_plot$posi_consensus = gsub("_NA","",accessibleTEs_plot$posi_consensus)
  accessibleTEs_plot$posi_consensus = as.numeric(as.character(accessibleTEs_plot$posi_consensus))
  accessibleTEs_plot$group_posi = paste(accessibleTEs_plot$group,accessibleTEs_plot$posi_consensus)
  accessibleTEs_plot$group_instance = paste(accessibleTEs_plot$group,accessibleTEs_plot$TEinstance)
  
  ### 
  Consensus_len = read.delim(paste("input/",ConsensusLenFile,sep=""),header=F,sep="")
  colnames(Consensus_len) = c("TEfamily","Len")
  Consensus_len = Consensus_len[order(-Consensus_len$Len),]
  Consensus_len = Consensus_len[!duplicated(Consensus_len$TEfamily),]
  
  #args <- commandArgs(trailingOnly = TRUE)
  #File <- args[1]
  #species <- args[2]
  #sampleInfo <- args[3]
  
  #############################          ##### load data have some issue
  ##### step 1: read file
  #color_stages_hg19 = c("iPSC_ATAC" = "#3288bd","iMeLC_ATAC" = "#92c5de","d6hPGCLC_ATAC" = "#41ab5d","ag77_ATAC" = "#a50f15")
  #color_stages_mm10 = c("ESC_ATAC" = "#3288bd","EpiLC_ATAC" = "#92c5de","d4PGCLC_ATAC" = "#41ab5d","d4c7RAB2_ATAC" = "#a50f15")
  #color_stages_mm10_full = c("ESC_ATAC" = "#3288bd","EpiLC_ATAC" = "#92c5de","d2PGCLC_ATAC" = "#a8ddb5", "d4PGCLC_ATAC" = "#41ab5d","d4c7RAB2_ATAC" = "#a50f15","d4c7PGCLC_ATAC" = "#8c6bb1","GSC_ATAC" = "#d9d9d9")
  
  ### select accessible instances and TE consensus
  Input_peak_plot = Input_peak[Input_peak$TEinstance %in% unique(accessibleTEs_plot$TEinstance),]
  # keep unique posi_consensus (some TE annotation may be overlapped)
  Input_peak_plot = Input_peak_plot[!duplicated(Input_peak_plot[,c("TEconsensus","TEinstance","posi_consensus")]),]
  
  ### sum to keep top consensus
  Input_peak_plot_sum = data.frame(Input_peak[!duplicated(Input_peak[,c("TEconsensus","TEinstance")]),] %>% group_by(TEconsensus) %>% dplyr::count(sort = TRUE))
  Input_peak_plot_sum_sub = data.frame(Input_peak_plot[!duplicated(Input_peak_plot[,c("TEconsensus","TEinstance")]),] %>% group_by(TEconsensus) %>% dplyr::count(sort = TRUE))
  Input_peak_plot_sum = merge(Input_peak_plot_sum,Input_peak_plot_sum_sub,by="TEconsensus",all.x=T)
  colnames(Input_peak_plot_sum) = c("TEconsensus","accessibleTEs.total","accessibleTEs.candidate")
  Input_peak_plot_sum = Input_peak_plot_sum[order(-Input_peak_plot_sum$accessibleTEs.candidate),]
  # write.csv(Input_peak_plot_sum,file = paste(File,".summary.csv",sep=""))
  Input_peak_plot_sum_sub = Input_peak_plot_sum_sub[!is.na(Input_peak_plot_sum_sub$TEconsensus),]
  Input_peak_plot_sum_sub = Input_peak_plot_sum_sub[1,]
  Input_peak_plot_sum_sub$Len = Consensus_len[Consensus_len$TEfamily == Input_peak_plot_sum_sub$TEconsensus,]$Len
  rm(Input_peak_plot_sum)
  
  # exclude insertions relative to the consensus
  Input_peak_plot = Input_peak_plot[!grepl("_NA",Input_peak_plot$posi_consensus) & !is.na(Input_peak_plot$TEconsensus) & Input_peak_plot$TEconsensus == Input_peak_plot_sum_sub$TEconsensus,]
  
  ##### step 2: data reformat
  ## combine if there are multiple replicates (average)
  GroupID = "iPSC_ATAC"
  colnames(Input_peak)
  Input_peak_plot[1:10,]
  
  ### convert the data format
  Input_peak_convert = data.frame("tmp1"=NA,"tmp2"=NA)
  GroupID = "T_0hr"
  colnames(Input_peak_plot)
  for (GroupID in levels(sampleInfo$group)) {
    if (!(GroupID %in% colnames(Input_peak_plot))){
      next
    }
    Input_peak_tmp = Input_peak_plot[,c("TEinstance","posi_consensus",GroupID)]
    Input_peak_tmp$group = GroupID
    colnames(Input_peak_tmp)[3] = "rpm"
    colnames(Input_peak_tmp)[2] = "posi"
    if (nrow(Input_peak_convert) == 1){
      Input_peak_convert = Input_peak_tmp
    } else {
      Input_peak_convert = rbind(Input_peak_convert,Input_peak_tmp)
    }
  }
  Input_peak_convert$group_posi = paste(Input_peak_convert$group,Input_peak_convert$posi)
  
  ## prepare the dataset for the clustering analysis add the deletion
  Input_peak_convert_plot = data.frame("group" = rep(levels(sampleInfo$group),each = Input_peak_plot_sum_sub$Len),"posi" = (1:Input_peak_plot_sum_sub$Len))
  Input_peak_convert_plot$order = 1:nrow(Input_peak_convert_plot)
  Input_peak_convert_plot$group_posi = paste(Input_peak_convert_plot$group,Input_peak_convert_plot$posi)
  Instance = "THE1B-int_dup9"
  rowID = 1
  
  Instances_summits = data.frame("TEinstance"=NA,"group"=NA,"samplesWithPeak" = NA,"posi"=NA,"group_posi"=NA)
  for (Instance in unique(Input_peak_convert$TEinstance)){
    ## convert the table one by one
    Input_peak_convert_tmp = Input_peak_convert[Input_peak_convert$TEinstance == Instance,]
    Input_peak_convert_plot = merge(Input_peak_convert_plot,Input_peak_convert_tmp[,c("group_posi","rpm")],all.x=T,by="group_posi")
    colnames(Input_peak_convert_plot)[ncol(Input_peak_convert_plot)] = Instance
    
    ## achieve Summits per instance
    Instances_summits_tmp = data.frame("TEinstance"=NA,"group"=NA,"samplesWithPeak" = NA,"posi"=NA,"group_posi"=NA)
    Instances_summits_tmp[1:length(levels(sampleInfo$group)),]$TEinstance = Instance
    groupOrder = 1
    for (groupOrder in 1:length(levels(sampleInfo$group))){
      Instances_summits_tmp[groupOrder,]$group = levels(sampleInfo$group)[groupOrder]
      Instances_summits_tmp[groupOrder,]$TEinstance = Instance
      # obtain peak summits per group
      accessibleTEs_plot_tmp = accessibleTEs_plot[accessibleTEs_plot$TEinstance == Instance & accessibleTEs_plot$group == levels(sampleInfo$group)[groupOrder],]
      Instances_summits_tmp[groupOrder,]$samplesWithPeak = length(unique(accessibleTEs_plot_tmp$sampleID))
      if (nrow(accessibleTEs_plot_tmp) > 0) {
        accessibleTEs_plot_tmp = accessibleTEs_plot_tmp[order(accessibleTEs_plot_tmp$posi_consensus),]
        Instances_summits_tmp[groupOrder,]$posi = accessibleTEs_plot_tmp[round(nrow(accessibleTEs_plot_tmp)/2,0)+1,]$posi_consensus
        Instances_summits_tmp[groupOrder,]$group_posi = paste(Instances_summits_tmp[groupOrder,]$group,Instances_summits_tmp[groupOrder,]$posi)
      }
    }
    Instances_summits = rbind(Instances_summits,Instances_summits_tmp)
  }
  Instances_summits$group_instance = paste(Instances_summits$group,Instances_summits$TEinstance)
  Instances_summits = Instances_summits[Instances_summits$group_instance %in% accessibleTEs_plot$group_instance,]
  
  ## 
  Input_peak_convert_plot = Input_peak_convert_plot[order(Input_peak_convert_plot$order),]
  rm(Input_peak_convert_tmp,Input_peak_tmp)
  Input_peak_convert_plot$rpm_aggregated = apply(Input_peak_convert_plot[,5:ncol(Input_peak_convert_plot)],1,sum,na.rm=T)
  Input_peak_convert_plot$rpm_mean = apply(Input_peak_convert_plot[,5:ncol(Input_peak_convert_plot)],1,mean,na.rm=T)
  colnames(Input_peak_convert_plot) = gsub("-|:",".",colnames(Input_peak_convert_plot))
  Input_peak_convert_plot$group = factor(Input_peak_convert_plot$group,levels=levels(sampleInfo$group))
  
  ###### step 3: plot the aggregated plot
  ## barplot
  p1 = ggplot(Input_peak_convert_plot, aes(x=posi,y=rpm_aggregated,group=group))+
    geom_bar(aes(fill = group),alpha = 0.5,stat="identity",position="identity") + 
    ylab("Aggregated RPM value")+
    xlab(paste(Input_peak_plot_sum_sub$TEconsensus,Input_peak_plot_sum_sub$Len,"(bp)"))+
    scale_x_continuous(breaks = seq(0,Input_peak_plot_sum_sub$Len, by = 250))+
    ggtitle(paste(TEfamily," (n=",Input_peak_plot_sum_sub$n,")",sep=""))+
    expand_limits(x = 0, y = 0)+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          legend.key = element_rect(fill = NA),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title=element_text(colour="black",size=rel(1)),
          legend.title=element_blank(),
          legend.position="right",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1)))
  pdf(paste("Figure_3A-",File,"aggregated.plot.pdf",sep="."),    # create PNG for the heat map
      width = 6,        # 5 x 300 pixels
      height = 4,
      pointsize = 12)        # smaller font size
  grid.draw(p1)
  dev.off()    
  
  #write.csv(Input_peak_convert_plot,file = paste(File,".aggregated.plot.csv",sep=""))
  
  ####### step 4: plot the heatmap
  ## 4.1
  rownames(Input_peak_convert_plot) = Input_peak_convert_plot$group_posi
  Input_peak_convert_heatmap = Input_peak_convert_plot[,5:(ncol(Input_peak_convert_plot)-1)]
  Input_peak_convert_heatmap_Z = ztran((Input_peak_convert_heatmap))
  Input_peak_convert_heatmap_Z_plot = as.matrix(t(Input_peak_convert_heatmap_Z))
  Input_peak_convert_heatmap_Z_plot[Input_peak_convert_heatmap_Z_plot<0 | is.na(Input_peak_convert_heatmap_Z_plot)] <- 0
  
  ## 4.2
  if (max(Input_peak_convert_heatmap_Z_plot) > 0 | min(Input_peak_convert_heatmap_Z_plot) < 0){
    ## clust
    Input_peak_convert_dist<-dist(Input_peak_convert_heatmap_Z_plot, method="euclidean")
    Input_peak_convert_cluster<-hclust(Input_peak_convert_dist, method='ward.D')
    Input_peak_convert_cluster_order = Input_peak_convert_cluster$labels[c(Input_peak_convert_cluster$order)]
    
    ## cut trees into 2:6 groups
    ClusterGroups_tmp = data.frame(cutree(Input_peak_convert_cluster,k = 2:6))
    colnames(ClusterGroups_tmp) = paste("clusters_",2:6,sep="")
    ClusterGroups_tmp = ClusterGroups_tmp[Input_peak_convert_cluster_order,]
    ClusterGroups_tmp = ClusterGroups_tmp[Input_peak_convert_cluster_order,]
  }
  
  ## 4.3 make the heatmap using geom_tile
  Col1 = 1
  Input_peak_convert_heatmap_Z.t = data.frame((Input_peak_convert_heatmap_Z))
  for (Col1 in 1:ncol(Input_peak_convert_heatmap_Z.t)){
    dis_plot_heat_final_Z_plot_tmp = data.frame(Input_peak_convert_heatmap_Z.t[,Col1])
    colnames(dis_plot_heat_final_Z_plot_tmp)[1] = "RPM_Z"
    dis_plot_heat_final_Z_plot_tmp$group_posi = rownames(Input_peak_convert_heatmap_Z.t)
    dis_plot_heat_final_Z_plot_tmp$TEinstance = colnames(Input_peak_convert_heatmap_Z.t)[Col1]
    if (Col1 == 1){
      Input_peak_convert_heatmap_Z_plot2 = dis_plot_heat_final_Z_plot_tmp
    } else {
      Input_peak_convert_heatmap_Z_plot2 = rbind(Input_peak_convert_heatmap_Z_plot2,dis_plot_heat_final_Z_plot_tmp)
    }
  }
  
  #Input_peak_convert_heatmap_Z_plot2$TEinstance = gsub("\\.","-",Input_peak_convert_heatmap_Z_plot2$TEinstance)
  Input_peak_convert_heatmap_Z_plot2$TEinstance = factor(Input_peak_convert_heatmap_Z_plot2$TEinstance,levels = rev(Input_peak_convert_cluster_order))
  Input_peak_convert_heatmap_Z_plot2$group_posi = factor(Input_peak_convert_heatmap_Z_plot2$group_posi,levels = rownames(Input_peak_convert_heatmap_Z))
  Input_peak_convert_heatmap_Z_plot2$Order = Input_peak_convert_heatmap_Z_plot2$group_posi
  levels(Input_peak_convert_heatmap_Z_plot2$Order) = 1:length(levels(Input_peak_convert_heatmap_Z_plot2$Order))
  Input_peak_convert_heatmap_Z_plot2$Order = as.numeric(as.character(Input_peak_convert_heatmap_Z_plot2$Order))
  #write.csv(Input_peak_convert_heatmap_Z_plot2,file = paste(File,".heatmap.rpm.csv",sep=""))
  
  ### visualize all potential minor peak regions
  Input_peak_convert_heatmap_tmp = Input_peak_convert_heatmap[,c(1,2)]
  Input_peak_convert_heatmap_tmp[,1] = rownames(Input_peak_convert_heatmap_tmp)
  Input_peak_convert_heatmap_tmp[,2] = 1:nrow(Input_peak_convert_heatmap_tmp)
  colnames(Input_peak_convert_heatmap_tmp) = c("group_posi","Order")
  Instances_summits_plot = merge(Instances_summits,Input_peak_convert_heatmap_tmp,by="group_posi",all.x=T)
  #write.csv(Instances_summits,file = paste(File,".heatmap.summit.csv",sep=""))
  Instances_summits_plot$TEinstance = gsub("-",".",Instances_summits_plot$TEinstance)
  
  ##
  p1 = ggplot(Input_peak_convert_heatmap_Z_plot2, aes(x=Order,y=TEinstance)) +
    geom_tile(aes(fill = RPM_Z,colour = RPM_Z),size=0.1) + 
    scale_fill_gradient2(low="white", high = "#08306b",na.value = "#f0f0f0")+
    scale_color_gradient2(low="white", high = "#08306b",na.value = "#f0f0f0")+
    xlab(paste(TEfamily,Input_peak_plot_sum_sub[1,]$TEconsensus,Input_peak_plot_sum_sub[1,]$Len,"(bp)"))+
    ylab("")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          #axis.text=element_blank(),
          #axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title=element_text(colour="black",size=rel(1.5)),
          #        legend.position=c(0.8,0.8),
          #        legend.position="bottom",
          legend.title=element_blank(),
          legend.position="right",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1.2)))
  ## 6.5.3 put the summit positions
  p1 = p1 + geom_point(data=Instances_summits_plot,aes(x=Order,y=TEinstance),shape=17,size=2,color="red")
  png(paste("Figure_S5A-",File,".heatmap.png",sep=""),    # create PNG for the heat map
      width = 40*800,        # 5 x 300 pixels
      height = 6*800,
      res = 800,
      pointsize = 10)
  grid.draw(p1)
  dev.off()  
}




