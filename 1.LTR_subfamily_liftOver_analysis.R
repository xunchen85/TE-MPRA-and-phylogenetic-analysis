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
library("Biostrings")
library(randomcoloR)
library("see")

#########################
hg19.div.sum = read.csv(file = "input/hg19_TE_rmsk_UCSC.age_2023_6_12.csv")

##
hg19.div.sum = hg19.div.sum[order(hg19.div.sum$mean.div),]
hg19.div.sum$TE_family2 = hg19.div.sum$TE_family
hg19.div.sum = data.frame(cSplit(hg19.div.sum,"TE_family2",sep=":"))

# load the liftover intersect results
liftover.summary = read.delim("input/TEwide_liftover.subfamily.summary_2023_8_29",header=F,sep="")

head(liftover.summary)

colnames(liftover.summary) = c("species","TE_family","total.count","intersected.count","intersected.count.200bp")
liftover.summary = merge(liftover.summary,hg19.div.sum[,c("TE_family","count.200")],by="TE_family",all.x=T)
liftover.summary$perC.200bp = liftover.summary$intersected.count.200bp/liftover.summary$count.200
liftover.summary$perC = liftover.summary$intersected.count/liftover.summary$total.count

# liftover species
liftover_species = c("hg19ToHg38","hg19ToPanTro2","hg19ToPanTro3","hg19ToPanTro4","hg19ToPanTro5","hg19ToPanTro6","hg19ToPanPan2","hg19ToGorGor3","hg19ToPonAbe2",
                     "hg19ToNomLeu1","hg19ToNomLeu3","hg19ToMacFas5","hg19ToRheMac10","hg19ToRheMac3","hg19ToRheMac2","hg19ToRheMac8","hg19ToPapAnu2",
                     "hg19ToCalJac1","hg19ToCalJac3","hg19ToSaiBol1","hg19ToMicMur1","hg19ToMm9","hg19ToMm10",
                     "hg19ToEriEur2","hg19ToSorAra2","hg19ToSusScr11","hg19ToSusScr2","hg19ToSusScr3",
                     "hg19ToBosTau4","hg19ToBosTau6","hg19ToBosTau9","hg19ToBosTau8","hg19ToBosTau7",
                     "hg19ToEquCab2","hg19ToCanFam2","hg19ToCanFam3","hg19ToFelCat5","hg19ToFelCat3","hg19ToFelCat4")
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

################ combined
liftover.summary.plot = merge(liftover.summary,hg19.div.sum,by="TE_family",all=T)

# filtering
liftover.summary.plot.sub = liftover.summary.plot[liftover.summary.plot$species %in% Kept_species,]
liftover.summary.plot.sub$species = factor(liftover.summary.plot.sub$species,levels=Kept_species)

# order by divergent rate
liftover.summary.plot.sub$TE_family = factor(liftover.summary.plot.sub$TE_family,levels=rev(hg19.div.sum$TE_family))
liftover.summary.plot.sub$TE_family2_1 = factor(liftover.summary.plot.sub$TE_family2_1,levels=rev(hg19.div.sum$TE_family2_1))

# kept family >= 100; LTRs; no int sequences
liftover.summary.plot.sub = liftover.summary.plot.sub[grepl(":LTR$",liftover.summary.plot.sub$TE_family) & 
                                                        liftover.summary.plot.sub$total.count>=100 &
                                                        !grepl("-int:",liftover.summary.plot.sub$TE_family),]

# kept recently evolved ones, <=20% in MicMur1
primateSpecific.families.hg38 = as.character(liftover.summary.plot.sub[liftover.summary.plot.sub$species == "hg19ToHg38" & liftover.summary.plot.sub$perC>=0.8,]$TE_family)
primateSpecific.families.hg38.200 = as.character(liftover.summary.plot.sub[liftover.summary.plot.sub$species == "hg19ToHg38" & liftover.summary.plot.sub$perC.200bp>=0.8,]$TE_family)

primateSpecific.families = as.character(liftover.summary.plot.sub[liftover.summary.plot.sub$species == "hg19ToMicMur1" & 
                                                                    liftover.summary.plot.sub$perC<=0.2 & 
                                                                    liftover.summary.plot.sub$TE_family %in% primateSpecific.families.hg38,]$TE_family)
primateSpecific.families.200bp = as.character(liftover.summary.plot.sub[liftover.summary.plot.sub$species == "hg19ToMicMur1" & 
                                                                          liftover.summary.plot.sub$perC.200bp<=0.2 & 
                                                                          liftover.summary.plot.sub$TE_family %in% primateSpecific.families.hg38.200,]$TE_family)

write.table(primateSpecific.families,file="primateSpecific.families_2023_4_10.txt",sep="",row.names = FALSE,quote = FALSE,col.names = FALSE)

liftover.summary.plot.sub.candidate = liftover.summary.plot.sub[liftover.summary.plot.sub$TE_family %in% primateSpecific.families,]
liftover.summary.plot.sub$is_kept = ifelse(liftover.summary.plot.sub$TE_family %in% primateSpecific.families,"kept","not")

### plot1: 351 family >= 100; LTRs; no int sequences
p2 <- ggplot(liftover.summary.plot.sub, aes(x=species, y=TE_family)) + 
  geom_tile(aes(fill=perC)) + 
  #geom_text(aes(label=round(perC.200bp,1)),color="white") + 
  scale_x_discrete(drop = FALSE)+
  #scale_y_discrete(drop = FALSE)+
  scale_fill_viridis_c() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",hjust=0.95,vjust=0.2,angle = 90),
    axis.text.y=element_blank(),
    axis.title=element_text(colour="black",size=rel(1.2)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position="right",
    #legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(.9)))

p2.family <- ggplot(liftover.summary.plot.sub[!duplicated(liftover.summary.plot.sub$TE_family),], aes(x=1, y=TE_family)) + 
  geom_tile(aes(fill=TE_family2_2)) + 
  #geom_text(aes(label=round(perC.200bp,1)),color="white") + 
  scale_x_discrete(drop = FALSE)+
  #scale_y_discrete(drop = FALSE)+
  scale_fill_brewer(palette="Set3")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",hjust=0.95,vjust=0.2,angle = 90),
    axis.text.y=element_blank(),
    axis.title=element_text(colour="black",size=rel(1.2)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position="right",
    #legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(.9)))
p2.kept <- ggplot(liftover.summary.plot.sub[!duplicated(liftover.summary.plot.sub$TE_family),], aes(x=1, y=TE_family)) + 
  geom_tile(aes(fill=is_kept)) + 
  #geom_text(aes(label=round(perC.200bp,1)),color="white") + 
  scale_x_discrete(drop = FALSE)+
  #scale_y_discrete(drop = FALSE)+
  scale_fill_brewer(palette="Set3")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",hjust=0.95,vjust=0.2,angle = 90),
    axis.text.y=element_blank(),
    axis.title=element_text(colour="black",size=rel(1.2)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position="right",
    #legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(.9)))
p2.div <- ggplot(liftover.summary.plot.sub[!duplicated(liftover.summary.plot.sub$TE_family),], aes(x=1, y=TE_family)) + 
  geom_tile(aes(fill=mean.age)) + 
  #geom_text(aes(label=round(perC.200bp,1)),color="white") + 
  scale_x_discrete(drop = FALSE)+
  #scale_y_discrete(drop = FALSE)+
  #scale_fill_viridis_c() +
  scale_fill_gradient2(low="#377eb8",mid = "white", high = "#e41a1c",midpoint = 90)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",hjust=0.95,vjust=0.2,angle = 90),
    axis.text.y=element_blank(),
    axis.title=element_text(colour="black",size=rel(1.2)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position="right",
    #legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(.9)))


gB <- ggplotGrob(p2.kept)
gC <- ggplotGrob(p2.family)
gA <- ggplotGrob(p2)
gD <- ggplotGrob(p2.div)

g2 = cbind(gA,gD,gB,gC, size = "first")

pdf("Figure_1A.pdf",    # create PNG for the heat map
    width = 16,        # 5 x 300 pixels
    height = 6,
    pointsize = 10)        # smaller font size
grid.draw(g2)
dev.off()   



#######
liftover.summary.plot.sub.macaque = liftover.summary.plot.sub[liftover.summary.plot.sub$species == "hg19ToMacFas5" & 
                                                                liftover.summary.plot.sub$is_kept == "kept",]
liftover.summary.plot.sub.macaque$TE_family2_1 = as.character(liftover.summary.plot.sub.macaque$TE_family2_1)
liftover.summary.plot.sub.macaque$Label2 = ifelse(liftover.summary.plot.sub.macaque$perC.200<=0.6,as.character(liftover.summary.plot.sub.macaque$TE_family2_1),NA)
liftover.summary.plot.sub.macaque$family = ifelse(!is.na(liftover.summary.plot.sub.macaque$Label2),liftover.summary.plot.sub.macaque$TE_family2_2,NA)
liftover.summary.plot.sub.macaque$count.final = ifelse(liftover.summary.plot.sub.macaque$count.200.x>=3000,3000,liftover.summary.plot.sub.macaque$count.200.x)

color_family = c("ERVL"="#F8766D","ERV1"="#00BFC4","ERVK"="#C77CFF","ERVL-MaLR"="#7CAE00")

p1 = ggplot(liftover.summary.plot.sub.macaque, aes(x=mean.age, y=perC.200bp*100)) +
  geom_point(aes(size=(count.final)),shape=21,fill="#3182bd",color="black")+
  geom_hline(yintercept=60, linetype="dashed", color = "grey", size=1)+
  geom_text_repel(aes(label = Label2,color = TE_family2_2),max.overlaps = 100)+
  scale_color_manual(values=color_family)+
  ylab("% human instances shared with macaque")+
  xlab("TE evolutionary age (Myrs)")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
    axis.text.x=element_text(colour="black",hjust=0.5,vjust=0.5,angle = 0),
    axis.title=element_text(colour="black",size=rel(1.2)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position="right",
    #legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(.9)))
pdf("Figure_S1A.pdf",    # create PNG for the heat map        
    width = 8,        # 5 x 300 pixels
    height = 6,
    pointsize = 10 )        # smaller font size
grid.draw(p1)
dev.off()


##################
cyto_group = read.delim("input/TEwide_subfamily_2023_6_7.cyto.group",header=T,sep="")
cyto_group = cyto_group[!is.na(cyto_group$score_combined),]
cyto_group.sum = data.frame(cyto_group %>% group_by(score_combined) %>% dplyr::count())
cyto_group = merge(cyto_group,cyto_group.sum,by="score_combined",all.x=T)

########## divergent rates
## load the div rate
hg.div = read.delim("input/hg19_rmsk_TE_1bp.div.bed.gz",header=F,sep="")
hg.div = cSplit(hg.div,"V4",sep=":",type.convert=as.character)
hg.div$TE_family = paste(hg.div$V4_1,hg.div$V4_3,hg.div$V4_4,sep=":")
hg.div$Len = hg.div$V3-hg.div$V2+1
hg.div$Age = ((hg.div$V7/100)/0.0000000022)/1000000
hg.div.kept = hg.div[hg.div$V4_1 %in% cyto_group$family,]
hg.div.kept = merge(hg.div.kept,cyto_group,by.x="V4_1",by.y="family",all.x=T)
hg.div.kept$V4_1 = factor(hg.div.kept$V4_1,levels = cyto_group$family)

### summary
hg19.div.sum = data.frame(hg.div.kept %>%
                            group_by(TE_family) %>%
                            dplyr::summarise(mean.div = mean(V7,na.rm=T),sd.div = sd(V7,na.rm=T),md.div=sum(abs(V7-mean(V7,na.rm=T))),
                                             quan_25=quantile(V7,na.rm = TRUE,probs=0.25),quan_50=quantile(V7,na.rm = TRUE,probs=0.5),quan_75=quantile(V7,na.rm = TRUE,probs=0.75),
                                             mean.age = mean(Age,na.rm=T),sd.age = mean(Age,na.rm=T),
                                             mean.Len = mean(Len,na.rm=T),sd.Len = sd(Len,na.rm=T),sum.Len = sum(Len,na.rm=T),
                                             count=n()))
hg19.div.sum$md.div = hg19.div.sum$md.div/hg19.div.sum$count
hg19.div.sum$cv.div = hg19.div.sum$sd.div/hg19.div.sum$mean.div
hg19.div.sum$cv.Len = hg19.div.sum$sd.Len/hg19.div.sum$mean.Len
hg19.div.sum$md.mean.div = hg19.div.sum$md.div/hg19.div.sum$mean.div
hg19.div.sum$p.shapiro = NA
hg19.div.sum$p.chi = NA
i = 1
for(i in 1:nrow(hg19.div.sum)){
  hg19.div.sum[i,]$p.shapiro = shapiro.test(hg.div.kept[hg.div.kept$TE_family == hg19.div.sum[i,]$TE_family,]$V7)[[2]]
  hg19.div.sum[i,]$p.chi = chisq.test(hg.div.kept[hg.div.kept$TE_family == hg19.div.sum[i,]$TE_family,]$V7)[[3]]
  hg19.div.sum[i,]$p.ks = ks.test(hg.div.kept[hg.div.kept$TE_family == hg19.div.sum[i,]$TE_family,]$V7,'pnorm')[[2]]
}

hg19.div.sum$p.shapiro.adj = p.adjust(hg19.div.sum$p.shapiro,n=nrow(hg19.div.sum))
hg19.div.sum$p.chi.adj = p.adjust(hg19.div.sum$p.chi,n=nrow(hg19.div.sum))
hg19.div.sum$family = hg19.div.sum$TE_family
hg19.div.sum = data.frame(cSplit(hg19.div.sum,"family",sep=":",type.convert=as.character))
hg19.div.sum = merge(hg19.div.sum,cyto_group,by.x="family_1",by.y="family",all=T)
hg19.div.sum = hg19.div.sum[!is.na(hg19.div.sum$TE_family),]
hg19.div.sum.tmp = data.frame(hg19.div.sum %>% group_by(score_combined) %>% dplyr::count())
colnames(hg19.div.sum.tmp)[2] = "n.actual"
hg19.div.sum = merge(hg19.div.sum,hg19.div.sum.tmp,by="score_combined",all.x=T)
hg19.div.sum = hg19.div.sum[with(hg19.div.sum, order(-n.actual,score_combined,-mean.div)),]


hg19.div.sum$dis.color = ifelse(hg19.div.sum$p.chi.adj<=0.001,"#d95f0e","#662506")
col_group = hg19.div.sum$dis.color
names(col_group) = hg19.div.sum$family_1
hg19.div.sum$family_1 = factor(hg19.div.sum$family_1,levels=hg19.div.sum$family_1)
hg.div.kept$V4_1 = factor(hg.div.kept$V4_1,levels=rev(hg19.div.sum$family_1))
hg19.div.sum$p.chi.adj.formated = ifelse(hg19.div.sum$p.chi.adj<0.001,formatC(hg19.div.sum$p.chi.adj, format = "e", digits = 2),round(hg19.div.sum$p.chi.adj,3))

#####
hg19.div.sum.separate = hg19.div.sum[order(hg19.div.sum$p.chi.adj),]

hg.div.kept.separate1 = hg.div.kept[hg.div.kept$V4_1 %in% hg19.div.sum.separate[hg19.div.sum.separate$dis.color == "#d95f0e",]$family_1,]
hg.div.kept.separate1$V4_1 = factor(hg.div.kept.separate1$V4_1,levels=as.character(hg19.div.sum.separate[hg19.div.sum.separate$dis.color == "#d95f0e",]$family_1))
hg.div.kept.separate2 = hg.div.kept[hg.div.kept$V4_1 %in% hg19.div.sum.separate[hg19.div.sum.separate$dis.color == "#662506",]$family_1,]
hg.div.kept.separate2$V4_1 = factor(hg.div.kept.separate2$V4_1,levels=as.character(hg19.div.sum.separate[hg19.div.sum.separate$dis.color == "#662506",]$family_1))

p1 = ggplot(hg.div.kept.separate1, aes(x=V4_1, y=V7)) +
  geom_violinhalf(width = 3.5,aes(fill=V4_1),linewidth=0.3,alpha=0.5) +
  scale_fill_manual(values = col_group)+
  #geom_text(data=hg19.div.sum,aes(x=family_1,y=40,label = p.chi.adj.formated))+
  coord_flip()+
  ylab("divergent rate")+
  xlab("putative simian-specific LTR subfamily group")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
    axis.text.x=element_text(colour="black",hjust=0.5,vjust=0.95,angle = 0),
    axis.title=element_text(colour="black",size=rel(1.2)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    #legend.position="right",
    legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1.2)))
pdf(paste("Figure_1B-divergent_final-1.pdf",sep=""),    # create PNG for the heat map        
    width = 7,        # 5 x 300 pixels
    height = 7,
    pointsize = 10 )        # smaller font size
grid.draw(p1)
dev.off()

p1 = ggplot(hg.div.kept.separate2, aes(x=V4_1, y=V7)) +
  geom_violinhalf(width = 3.5,aes(fill=V4_1),linewidth=0.3,alpha=0.5) +
  scale_fill_manual(values = col_group)+
  #geom_text(data=hg19.div.sum,aes(x=family_1,y=40,label = p.chi.adj.formated))+
  coord_flip()+
  ylab("divergent rate")+
  xlab("putative simian-specific LTR subfamily group")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
    axis.text.x=element_text(colour="black",hjust=0.5,vjust=0.95,angle = 0),
    axis.title=element_text(colour="black",size=rel(1.2)),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    #legend.position="right",
    legend.position = "none",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1.2)))
pdf(paste("Figure_S3A-divergent_final-2.pdf",sep=""),    # create PNG for the heat map        
    width = 7,        # 5 x 300 pixels
    height = 7,
    pointsize = 10 )        # smaller font size
grid.draw(p1)
dev.off()
