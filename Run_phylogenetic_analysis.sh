#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

########## step 1: subfamily groups
# Subfamily consensus sequences could were downloaded from the repbase or DFAM databases. 
# Here we privided the blastn output file "TEwide_subfamily_2023_6_7.outfmt6.tsv".

#### 1.1 transposable element subfamily consensus sequences similarity analysis using blastn (contributed by Dr. Zicong Zhang)
makeblastdb -in TEwide_subfamily_2023_6_7.fa -dbtype nucl

blastn -task dc-megablast \
	-outfmt 6 \
	-num_threads 4 \
	-query TEwide_subfamily_2023_6_7.fa \
	-db TEwide_subfamily_2023_6_7.fa \
	-out TEwide_subfamily_2023_6_7.outfmt6.tsv

echo -e "from\tto\tscore" >TEwide_subfamily_2023_6_7.outfmt6.cyto.tsv
awk '{if($1!=$2)print$12"\t"$1"|"$2}' TEwide_subfamily_2023_6_7.outfmt6.tsv | sort -k2 -k 1nr | uniq -f 1 | sed 's/\|/\t/' | awk '{print$2"\t"$3"\t"$1}' >>TEwide_subfamily_2023_6_7.outfmt6.cyto.tsv

#### 1.2 use edge scores 200 to detect subfamily groups containing the 35 simina-specific LTR subfamilies.
# using all edges to retrieve subfamily groups containing a single candidate LTR subfamily
# Here we got the list of subfamily group indicated in the "score_combined" column in the file "TEwide_subfamily_2023_6_7.cyto.group"



########## step 2: subdivide each subfamily into instance clusters
# first downloaded the human reference genome "hg19.fa.gz" and in TE annotation file "hg19.fa.out.gz" from UCSC (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/).

#### 2.1 convert hg19.fa.out.gz file to BED format
# reformat the downloaded out file
zless hg19.fa.out.gz | awk '{print$1,$5,$6,$7,$9,$10,$11}' | grep chr |sed 's/ALR\/Alpha/ALR_Alpha/' | sed 's/BSR\/Beta/BSR_Beta/'| sed 's/\// /'| awk '{if($8=="") print $1,$2,$3,$4,$5,$6,"Unknown",$7; else print$0}' | sed 's/ /\t/g' | sed 's/\tC\t/\t-\t/' | sed 's/ALR_Alpha/ALR\/Alpha/' | sed 's/BSR_Beta/BSR\/Beta/' >hg19.fa.out_2

# convert from the out file to GTF format; makeTEgtf.pl perl script was downloaded from https://github.com/mhammell-laboratory/TEtranscripts/issues/83
perl makeTEgtf.pl -c 2 -s 3 -e 4 -o 5 -t 6 -n hg19_rmsk -f 8 -C 7 -S 1 -1 hg19.fa.out_2 >hg19_rmsk_TE.gtf

# convert from GTF to BED format
awk '{print$1"\t"$4-1"\t"$5"\t"$10":"$12":"$14":"$16"\t"$6"\t"$7}' hg19_rmsk_TE.gtf | sed 's/\"//g' | sed 's/\;//g' >hg19_rmsk_TE_0bp.bed

# achieve the divergence rate
zless hg19.fa.out.gz | awk '{if($9=="C"||$9=="+")print$5":"$6"-"$7"\t"$9"\t"$10"\t"$11"\t"$2}' > hg19.fa.out2
awk '{print$1"\t"$4"\t"$5"\t"$10":"$12":"$14":"$16"\t"$6"\t"$7}' hg19_rmsk_TE.gtf | sed 's/\"//g' | sed 's/\;//g' >hg19_rmsk_TE_1bp.bed
python Extract_div_by_coordinate.py -i hg19_rmsk_TE_1bp.bed -l hg19.fa.out2 >hg19_rmsk_TE_1bp.div.bed

rm hg19.fa.out_2
rm hg19.fa.out2
rm hg19_rmsk_TE.gtf

#### 2.2 extract candidate LTR subfamilies, such as MER11A
# also downloaded the "hg19.fa.gz" human reference genome from UCSC and then unzipped (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)
python Extract_by_TEfamilyname_each.py -i hg19_rmsk_TE_0bp.bed -l MER11A >hg19_rmsk_MER11A_0bp.bed

bedtools getfasta -nameOnly \
	-bed hg19_rmsk_MER11A_0bp.bed \
	-fi hg19.fa -fo hg19_rmsk_MER11A_0bp.fa -s

# we also extracted the sequences for other MER11 subfamilies, MER11A, MER11B, MER11C

#### 2.3 unrooted trees construction, such as MER11A
# run mafft v7.505-without-extensions
mafft --localpair --maxiterate 1000 --thread 4 hg19_rmsk_MER11A_0bp.fa >hg19_rmsk_MER11A_0bp.mafft.fa
# run prank v.170427
prank -showanc -njtree -uselogs -prunetree -F -showevents -d=hg19_rmsk_MER11A_0bp.mafft.fa -o=hg19_rmsk_MER11A_0bp.mafft.prank
sed 's/_+_$//' hg19_rmsk_MER11A_0bp.mafft.prank.best.fas | sed 's/_-_$//' >hg19_rmsk_MER11A_0bp.mafft.prank.best.fa
# run trimal v1.4.1
trimal -in hg19_rmsk_MER11A_0bp.mafft.prank.best.fa -gt 0.01 -out hg19_rmsk_MER11A_0bp.mafft.prank.gt99.fa -fasta
# remove sequences containing only gaps and "N" after trimming
perl Exclude_emptyFasta.pl hg19_rmsk_MER11A_0bp.mafft.prank.gt99.fa >hg19_rmsk_MER11A_0bp.mafft.prank.opt.gt99.fa
cp hg19_rmsk_MER11A_0bp.mafft.prank.opt.gt99.fa hg19_rmsk_MER11A_0bp.mafft.prank.gt99.fa
# run iqtree2
iqtree2 -s hg19_rmsk_MER11A_0bp.mafft.prank.opt.gt99.fa -nt AUTO -m MFP -bb 6000 -asr -minsup .95 -T 4

# we performed similar analysis for other subfamilies

#### 2.4 subdivide MER11A into instance clusters using "Tree_subdivision.R" R script and achieve the MER11A cluster consensus sequences



######## step 3: determine phyletic groups among clusters per subfamily group
#### 3.1 after renaming the consensus sequences, we then gather all of these consensus sequences from a subfamily group for the following rooted tree analysis.
# e.g., hg19_MER11ABC_v1.fa which containg MER11A, MER11B, MER11C consensus sequences in human

# mafft
mafft --localpair \
	--maxiterate 1000 \
	--thread 4 hg19_MER11ABC_v1.fa >hg19_MER11ABC_v1.mafft.fa

# prank without gaps
prank -showanc \
	-njtree \
	-uselogs \
	-prunetree \
	-F \
	-showevents \
	-d=hg19_MER11ABC_v1.mafft.fa \
	-o=hg19_MER11ABC_v1.mafft.prank

# iqtree2
iqtree2 -s hg19_MER11ABC_v1.mafft.prank.best.fas \
        --model-joint 12.12 \
        -B 1000 \
        -T AUTO \
        --root-test -zb 1000 -au \
        --prefix ${inFile}.mafft.prank.best.fas.nonrev_dna

# we did the same for other subfamily groups.

#### 3.2 divergence rate calculation between cluster consensus sequences from each subfamily group (contributed by Dr. Zicong Zhang)
cd input_trees_consensus/

for SEQ in `ls *i.v2.fa`
do
  NAME=`basename -s .fa ${SEQ}`
  # multiple sequence alignment
  ginsi --thread 4 ${NAME}.fa > ${NAME}.ginsi.fa
  # divergence rate calculation
  raxmlHPC-PTHREADS-AVX -f x \
	-p 12345 \
	-m GTRGAMMA \
	-T 4 \
	-s ${NAME}.ginsi.fa \
	-n ${NAME}
done

# we also computed the divergence rates between cluster consensus sequences we identified and the ones in the DFAM databases.
# we also computed the divergence rates between human and macaque MER11 cluster consensus sequences.

cd ../

#### 3.3 liftOver analysis
# lift over human TEs to other species, e.g., macaque genome; the chain files, such as hg19ToMacFas5.over.chain was downloaded from UCSC.
# hg19 versus other species could be found here (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/).
# macaque versus hg19 could be found here (https://hgdownload.soe.ucsc.edu/goldenPath/macFas5/liftOver/).

species1="hg19"
species1_u="Hg19"
species2="Macfas5"
species2_l="macfas5"

bnMapper.py -k \
	-t 0.5 \
	input/hg19_rmsk_TE_0bp.bed \
	input_chain_files/${species1}To${species2}.over.chain.pkl >hg19_rmsk_TE_0bp.${species1}To${species2}.bed

# format conversion
awk '{print$0":"$1":"$2":"$3}' hg19_rmsk_TE_0bp.${species1}To${species2}.bed >hg19_rmsk_TE_0bp.${species1}To${species2}.bed2

# liftover reciprocally against human
bnMapper.py -k \
	-t 0.5 \
	hg19_rmsk_TE_0bp.${species1}To${species2}.bed2 \
	chain_files/${species2_l}To${species1_u}.over.chain >hg19_rmsk_TE_0bp.${species1}To${species2}.back.bed

# intersect
bedtools intersect \
	-wao \
	-a input/hg19_rmsk_TE_0bp.rename2.bed \
	-b hg19_rmsk_TE_0bp.${species1}To${species2}.back.bed >hg19_rmsk_TE_0bp.${species1}To${species2}.back.intersect.bed

# summarize the results
python Summarize_liftoverIntersect_by_TEFamily.py \
	-l input/hg19_rmsk_TE.bed.TEfamilyCounts \
	-i hg19_rmsk_TE_0bp.${species1}To${species2}.back.intersect.bed >hg19_rmsk_TE_0bp.hg19ToMacFas5.back.intersect.out

python Summarize_liftoverIntersect_by_TEFamily_Group.py \
	-l hg19_rmsk_TE.bed.TEfamilyCounts \
	-i hg19_rmsk_TE_0bp.${species1}To${$species2}.back.intersect.bed \
	-n TEwide_instance.list_2023_8_29.hg19 \
	-g TEwide_group.list_final_2023_8_29.hg19 >hg19_rmsk_TE_0bp.${species1}To${species2}.back.intersect.group.out

# we then combine the liftOver results into a summary table at the subfamily and a table at the phyletic group level separate.
cat hg19_rmsk_TE_0bp.hg19To*.back.intersect.out >hg19_rmsk_TE_0bp.liftover.summary_2023_4_10


#### 3.4 phyletic groups determination

# We determined the phyletic groups based on the top-selected rooted tree. These are the rules:
#	1) we first examined the internal branch lengths of the tree and grouped adjacent clusters based on the top branch lengths (bootstrap values) manually. 
#	2) We also examined the heatmap of divergence rates to look at extreme values between every adjacent clusters to validate the phyletic groups. 
#	3) We kept the phyletic groups after we confirmed that the clusters from a phyletic group were evolutionary close to each other.
