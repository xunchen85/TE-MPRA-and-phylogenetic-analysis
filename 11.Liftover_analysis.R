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
# load FASTA
library("Biostrings")
library(randomcoloR)

myTheme = theme(
  plot.title = element_text(hjust = 0.5, size = rel(1)),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(), 
  axis.line = element_line(colour = "black"),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.text=element_text(colour="black",size=rel(1),angle = 0),
  axis.text.x=element_text(colour="black",hjust=0.95,vjust=0.2,angle = 90),
  axis.title=element_text(colour="black",size=rel(1)),
  legend.key = element_rect(colour = "transparent", fill = "white"),
  #legend.position="right",
  legend.position = "right",
  legend.background = element_blank(),
  legend.text=element_text(size=rel(1)))

#########################

Family="MER11A_F1"
color_species = c("Ancient" = "#6a3d9a","Consensus" = "#6a3d9a","hg19"="#e31a1c","macFas5"="#1f78b4","panTro4"="#33a02c")
shape_species = c("Ancient" = 18,"Consensus" = 18,"hg19"=NA,"macFas5"=NA,"panTro4"=NA)

######################### 1. summary table
Summary_Table1 = read.csv("input/Summary_Table1_2022_8_9.csv")
Summary_Table1$Instance_coordinate_species_0bp = paste(Summary_Table1$species,":",Summary_Table1$chr,":",Summary_Table1$start+1,"-",Summary_Table1$end,sep="")
Summary_Table1$Instance_coordinate_species_0bp = ifelse(Summary_Table1$species == "Ancient" | Summary_Table1$species == "Consensus",paste("Ancient",Summary_Table1$TEfamily,sep=":"),Summary_Table1$Instance_coordinate_species_0bp)
Summary_Table2 = read.csv("input/Summary_Table2_2022_8_9.csv")
Summary_Table2$Instance_coordinate_species_0bp = paste(Summary_Table2$Species,":",Summary_Table2$Instance_coordinate,sep="")
Summary_Table2$Instance_coordinate_species_0bp = ifelse(Summary_Table2$Group == "Consensus",paste("Ancient",Summary_Table2$Family,sep=":"),Summary_Table2$Instance_coordinate_species_0bp)
Summary_Table2 = merge(Summary_Table2,Summary_Table1[,c("Instance_coordinate_species_0bp","uniqueID_F1","uniqueID_F2","TEinstance","instanceID.renamed")],by="Instance_coordinate_species_0bp",all.x=T)

########################## 2. load each information file
# 2.1 final subfamily info of human 
subfamily_info = read.csv("input/MER11_subtree_info_2022_12_27.csv")

# 2.2 TEs each species
hg19_rmsk_0bp = read.delim("input/hg19_rmsk_TE_0bp_combined.rename.bed.gz",header=F,sep="")
macFas5_rmsk_0bp = read.delim("input/macFas5_rmsk_TE_0bp_combined.rename.bed.gz",header=F,sep="")
panTro4_rmsk_0bp = read.delim("input/panTro4_rmsk_TE_0bp_combined.rename.bed.gz",header=F,sep="")

# 2.3 hg19 liftover and intersect
hg19_liftover_intersect.raw = read.delim("input/hg19_rmsk_TE_0bp_combined.rename.liftover.intersect.bed.gz",sep="",header=F)

# missing
# Proboscis monkey; Golden snub-nosed monkey no liftover back:
# in total we retrive 39 pairs between hg19 to other species

# 2.4 corss-species between hg19, macFas5 and panTro4
crossspecies_intersect.raw = read.delim("input/hg19_macFas5_panTro4_rmsk_TE_0bp_combined.rename.pair.intersect_2023_2_18.bed.gz",sep="",header=F)

########################### 3. liftover hg19 to each species 
N.leaf = 10
M.dist = 0.02
# 
hg19_liftover_intersect = hg19_liftover_intersect.raw[hg19_liftover_intersect.raw$V9 != -1,]
hg19_liftover_intersect = data.frame(cSplit(hg19_liftover_intersect,"V11",sep=":",type.convert = as.character))
colnames(hg19_liftover_intersect) = c("species.liftover","chr_hg19.original","start_hg19.original","end_hg19.original",
                                      "name_hg19.original","anno.original","direction.original",
                                      "chr_hg19.intersect","start_hg19.intersect","end_hg19.intersect","len_hg19.intersect","name.liftover",
                                      "chr_species.liftover","start_species.liftover","end_species.liftover")
# kept unique row
hg19_liftover_intersect = hg19_liftover_intersect[hg19_liftover_intersect$name_hg19.original == hg19_liftover_intersect$name.liftover,]
hg19_liftover_intersect = hg19_liftover_intersect[order(-hg19_liftover_intersect$len_hg19.intersect),]
hg19_liftover_intersect = hg19_liftover_intersect[!duplicated(hg19_liftover_intersect[,c("species.liftover","name_hg19.original")]),]

# Ordered species
# 1. human:hg19ToHg38; 2. Chimp:hg19ToPanTro6 (5Myr); 3. Bonobo:hg19ToPanPan2 (5Myr); 4. Gorilla:hg19ToGorGor3 (10Myr); 5. Orangutan:hg19ToPonAbe2 (15Myr); 
# 6. Gibbon:hg19ToNomLeu3 (20 Myr); 7. Crab-eating macaque:hg19ToMacFas5 (25 Myr); 8. Rhesus:hg19ToRheMac10 (25 Myr); 9. Baboon (anubis):hg19ToPapAnu2 (30 Myr); 10. Marmoset:hg19ToCalJac3 (35Myr); 
# 11. Squirrel monkey:hg19ToSaiBol1 (35 Myr); 12. Mouse lemur:hg19ToMicMur1 (65Myr); 13. Mouse:hg19ToMm10 (75Myr); 14. Hedgehog:hg19ToEriEur2; 15. Cow:hg19ToSorAra2;
# 16. shrew:hg19ToSusScr11; 17. pig:hg19ToBosTau9 (80 Myr); 18. horse:hg19ToEquCab2, 19. dog:hg19ToCanFam; 20. cat:hg19ToFelCat5

# reordering
liftover_species = c("hg19ToHg38","hg19ToPanTro2","hg19ToPanTro3","hg19ToPanTro4","hg19ToPanTro5","hg19ToPanTro6","hg19ToPanPan2","hg19ToGorGor3","hg19ToPonAbe2",
                     "hg19ToNomLeu1","hg19ToNomLeu3","hg19ToMacFas5","hg19ToRheMac10","hg19ToRheMac3","hg19ToRheMac2","hg19ToRheMac8","hg19ToPapAnu2",
                     "hg19ToCalJac1","hg19ToCalJac3","hg19ToSaiBol1","hg19ToMicMur1","hg19ToMm9","hg19ToMm10",
                     "hg19ToEriEur2","hg19ToSorAra2","hg19ToSusScr11","hg19ToSusScr2","hg19ToSusScr3",
                     "hg19ToBosTau4","hg19ToBosTau6","hg19ToBosTau9","hg19ToBosTau8","hg19ToBosTau7",
                     "hg19ToEquCab2","hg19ToCanFam2","hg19ToCanFam3","hg19ToFelCat5","hg19ToFelCat3","hg19ToFelCat4")
### plot
Kept_species = liftover_species

Kept_species = Kept_species[c(grep("hg19ToHg38",Kept_species),   ## human 0 Myr hg38
                              grep("hg19ToPanTro6",Kept_species),  ## chimpanzee (hominoids) 5 Myr
                              grep("hg19ToGorGor3",Kept_species),  ## Gorilla (hominoids) 10 Myr
                              grep("hg19ToPonAbe2",Kept_species),  ## Orangutan (hominoids) 15 Myr
                              grep("hg19ToNomLeu3",Kept_species),  ## Gibbon (hominoids) 20 Myr
                              grep("hg19ToMacFas5",Kept_species),  ## Macaque (old world) 25 Myr
                              grep("hg19ToPapAnu2",Kept_species),  ## Baboon (old world) 30 Myr
                              grep("hg19ToCalJac3",Kept_species),  ## Marmoset (new world) 35 Myr
                              grep("hg19ToMicMur1",Kept_species), ## Mouse lemur 65 Myr
                              grep("hg19ToMm10",Kept_species))]   ## mouse 75 Myr

#Kept_species = liftover_species
########################### 4. liftover between hg19, macFas5 and panTro4
crossspecies_intersect.raw2 = data.frame(cSplit(crossspecies_intersect.raw,"V6",sep=":"))
crossspecies_intersect.raw2 = crossspecies_intersect.raw2[,c("V1","V2","V6_1","V10","V13")]
colnames(crossspecies_intersect.raw2) = c("ref","liftover","original.instance","liftover.instance","intersect.length")
# order by length
crossspecies_intersect.raw2 = crossspecies_intersect.raw2[order(-crossspecies_intersect.raw2$intersect.length),]
crossspecies_intersect.raw2$uniqueID = paste(crossspecies_intersect.raw2$original.instance,crossspecies_intersect.raw2$liftover.instance)
# exclude unique original and liftover paired
crossspecies_intersect.raw2 = crossspecies_intersect.raw2[!duplicated(crossspecies_intersect.raw2$uniqueID),]
# excluded non-intersected fragments
crossspecies_intersect.raw2.intersect = crossspecies_intersect.raw2[crossspecies_intersect.raw2$liftover.instance!=".",]
crossspecies_intersect.raw2.nonintersect = crossspecies_intersect.raw2[crossspecies_intersect.raw2$liftover.instance=="." & 
                                                                         !(crossspecies_intersect.raw2$original.instance %in% crossspecies_intersect.raw2.intersect$original.instance),]
crossspecies_intersect.raw2 = rbind(crossspecies_intersect.raw2.intersect,crossspecies_intersect.raw2.nonintersect)
crossspecies_intersect.raw2 = crossspecies_intersect.raw2[order(-crossspecies_intersect.raw2$intersect.length),]
# kept unique original instance with longer length
crossspecies_intersect.raw2.unique = crossspecies_intersect.raw2[!duplicated(crossspecies_intersect.raw2[,c("liftover","original.instance")]),]
crossspecies_intersect.raw2.unique$uniqueID2 = paste(crossspecies_intersect.raw2.unique$liftover,crossspecies_intersect.raw2.unique$original.instance)

# achieve 
crossspecies_intersect.raw2$is_liftover_intersect = ifelse(crossspecies_intersect.raw2$liftover.instance == ".","liftover","intersect")
crossspecies_intersect.raw2.count = data.frame(crossspecies_intersect.raw2 %>%
                                                 group_by(liftover,is_liftover_intersect,original.instance) %>%
                                                 dplyr::summarise(intersect.count = n()))
crossspecies_intersect.raw2.count$uniqueID2 = paste(crossspecies_intersect.raw2.count$liftover,crossspecies_intersect.raw2.count$original.instance)
crossspecies_intersect.raw2.count = merge(crossspecies_intersect.raw2.count,crossspecies_intersect.raw2.unique[,c("uniqueID2","liftover.instance","intersect.length")],by="uniqueID2",all.x=T)
crossspecies_intersect.raw2.count = merge(crossspecies_intersect.raw2.count,Summary_Table1[,c("instanceID.renamed","TEfamily","species")],by.x="liftover.instance",by.y="instanceID.renamed",all.x=T)
crossspecies_intersect.raw2.count$intersect.count = ifelse(crossspecies_intersect.raw2.count$intersect.length == 0,0,crossspecies_intersect.raw2.count$intersect.count)

# Summary 1
Summary_Table1.plot = Summary_Table1[Summary_Table1$species %in% c("hg19","panTro4","macFas5"),]
Summary_Table1.plot$liftover = ifelse(Summary_Table1.plot$species == "hg19","hg19ToMacFas5",NA)
Summary_Table1.plot$liftover = ifelse(Summary_Table1.plot$species == "panTro4","panTro4ToHg19",Summary_Table1.plot$liftover)
Summary_Table1.plot$liftover = ifelse(Summary_Table1.plot$species == "macFas5","macFas5ToHg19",Summary_Table1.plot$liftover)
Summary_Table1.plot.tmp = Summary_Table1[Summary_Table1$species %in% c("hg19"),]
Summary_Table1.plot.tmp$liftover = "hg19ToPanTro4"
Summary_Table1.plot = rbind(Summary_Table1.plot,Summary_Table1.plot.tmp)
Summary_Table1.plot$uniqueID2 = paste(Summary_Table1.plot$liftover,Summary_Table1.plot$instanceID.renamed)
Summary_Table1.plot = merge(Summary_Table1.plot,crossspecies_intersect.raw2.count,by="uniqueID2",all.x=T)

Summary_Table1.plot$is_liftover_intersect = ifelse(is.na(Summary_Table1.plot$is_liftover_intersect),"no_ortholog",Summary_Table1.plot$is_liftover_intersect)
Summary_Table1.plot$intersect.count = ifelse(!is.na(Summary_Table1.plot$intersect.count) & Summary_Table1.plot$intersect.count >=2,2,Summary_Table1.plot$intersect.count)
head(Summary_Table1.plot)
Summary_Table1.plot.sum = data.frame(Summary_Table1.plot %>%
                                       group_by(liftover.x,species.x,TEfamily.x,TEfamily.y,is_liftover_intersect,intersect.count) %>%
                                       dplyr::summarise(count = n()))
Summary_Table1.plot.sum$uniqueID3 = paste(Summary_Table1.plot.sum$liftover.x,Summary_Table1.plot.sum$species.x,Summary_Table1.plot.sum$TEfamily.x)
Summary_Table1.plot.sum.total = data.frame(Summary_Table1.plot %>%
                                             group_by(liftover.x,species.x,TEfamily.x) %>%
                                             dplyr::summarise(total = n()))
Summary_Table1.plot.sum.total$uniqueID3 = paste(Summary_Table1.plot.sum.total$liftover.x,Summary_Table1.plot.sum.total$species.x,Summary_Table1.plot.sum.total$TEfamily.x)

Summary_Table1.plot.sum = merge(Summary_Table1.plot.sum,Summary_Table1.plot.sum.total[,c("uniqueID3","total")],by="uniqueID3",all.x=T)
Summary_Table1.plot.sum$perC = Summary_Table1.plot.sum$count/Summary_Table1.plot.sum$total
Summary_Table1.plot.sum$is_sameFamily = ifelse(!is.na(Summary_Table1.plot.sum$TEfamily.y) & Summary_Table1.plot.sum$TEfamily.x == Summary_Table1.plot.sum$TEfamily.y,"Consistent","Different")
Summary_Table1.plot.sum$group = paste(Summary_Table1.plot.sum$is_liftover_intersect,Summary_Table1.plot.sum$intersect.count,Summary_Table1.plot.sum$is_sameFamily)
Summary_Table1.plot.sum$group = factor(Summary_Table1.plot.sum$group,levels = rev(c("no_ortholog NA Different","liftover 0 Different","intersect 1 Different","intersect 2 Different","intersect 2 Consistent","intersect 1 Consistent")))
Summary_Table1.plot.sum = Summary_Table1.plot.sum[Summary_Table1.plot.sum$TEfamily.x %in%c("MER11A","MER11B","MER11C","MER11D","MER34A1","MER34C_","MER52A","MER52C"),]
Summary_Table1.plot.sum$TEfamily.x = factor(Summary_Table1.plot.sum$TEfamily.x,levels=c("MER11A","MER11B","MER11C","MER11D","MER34A1","MER34C_","MER52A","MER52C"))
## proportion of pantro4 species-specific
Summary_Table1.plot.sum[Summary_Table1.plot.sum$species.x == "panTro4" & Summary_Table1.plot.sum$TEfamily.x %in% c("MER11A","MER11B","MER11C","MER11D") & !Summary_Table1.plot.sum$is_liftover_intersect %in% "intersect",]
Summary_Table1.plot.sum[Summary_Table1.plot.sum$species.x == "hg19" & Summary_Table1.plot.sum$TEfamily.x %in% c("MER11A","MER11B","MER11C","MER11D") & !Summary_Table1.plot.sum$is_liftover_intersect %in% "intersect",]
Summary_Table1.plot.sum[Summary_Table1.plot.sum$species.x == "macFas5" & Summary_Table1.plot.sum$TEfamily.x %in% c("MER11A","MER11B","MER11C","MER11D") & !Summary_Table1.plot.sum$is_liftover_intersect %in% "intersect",]

### original color
#color_group = c("no_ortholog NA Different" = "#f7f7f7","liftover 0 Different" = "#bababa",
#                "intersect 2 Different" = "#bababa","intersect 1 Consistent" = "#1f78b4","intersect 1 Different"="#bababa","intersect 2 Consistent"="#a6cee3")
color_group = c("no_ortholog NA Different" = "#1f78b4","liftover 0 Different" = "#a6cee3",
                "intersect 2 Different" = "#a6cee3","intersect 1 Consistent" = "#f7f7f7","intersect 1 Different"="#a6cee3","intersect 2 Consistent"="#bababa")

p1 = ggplot(Summary_Table1.plot.sum[Summary_Table1.plot.sum$species.x == "hg19" &Summary_Table1.plot.sum$liftover.x == "hg19ToPanTro4",], aes(fill=group, y=perC, x=TEfamily.x)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=color_group)+
  xlab("")+ylab("Proportion")+
  scale_x_discrete(drop = FALSE)+
  ggtitle("hg19 vs panTro4")+
  myTheme
p2 = ggplot(Summary_Table1.plot.sum[Summary_Table1.plot.sum$species.x == "hg19" &Summary_Table1.plot.sum$liftover.x == "hg19ToMacFas5" ,], aes(fill=group, y=perC, x=TEfamily.x)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=color_group)+
  xlab("")+ylab("Proportion")+
  scale_x_discrete(drop = FALSE)+
  ggtitle("hg19 vs macFas5")+
  myTheme
p3 = ggplot(Summary_Table1.plot.sum[Summary_Table1.plot.sum$species.x == "macFas5",], aes(fill=group, y=perC, x=TEfamily.x)) + 
  geom_bar(position="stack", stat="identity")+ 
  scale_fill_manual(values=color_group)+
  xlab("")+ylab("Proportion")+
  ggtitle("macFas5 vs hg19")+
  scale_x_discrete(drop = FALSE)+
  myTheme
p4 = ggplot(Summary_Table1.plot.sum[Summary_Table1.plot.sum$species.x == "panTro4",], aes(fill=group, y=perC, x=TEfamily.x)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=color_group)+
  xlab("")+ylab("Proportion")+
  ggtitle("panTro4 vs hg19")+
  scale_x_discrete(drop = FALSE)+
  myTheme
gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
gC <- ggplotGrob(p3)
gD <- ggplotGrob(p4)

g2 = rbind(gA,gB,gC,gD, size = "first")

pdf(paste("Figure_4A_S9A-liftover","_2023_2_19-1.pdf",sep=""),    # create PNG for the heat map
    width = 6,        # 5 x 300 pixels
    height = 12,
    pointsize = 10)        # smaller font size
grid.draw(g2)
dev.off()

#################### step 5: kept the top candidates; consistent liftover
Summary_Table1.plot$uniqueID4 = paste(Summary_Table1.plot$instanceID.renamed,Summary_Table1.plot$liftover.instance)
Summary_Table1.plot$uniqueID4_reciprocal = paste(Summary_Table1.plot$liftover.instance,Summary_Table1.plot$instanceID.renamed)
Summary_Table1.plot$is_intersect_reciprocal = ifelse(!is.na(Summary_Table1.plot$intersect.count) &
                                                       Summary_Table1.plot$intersect.count >0 &
                                                       Summary_Table1.plot$uniqueID4 %in% Summary_Table1.plot$uniqueID4_reciprocal,paste(Summary_Table1.plot$is_liftover_intersect,"reciprocal",sep="_"),Summary_Table1.plot$is_liftover_intersect)

#write.csv(Summary_Table1.plot,file="Summary_Table1.intersect_reciprocal_hu_ch_ma_2023_2_20.csv")







