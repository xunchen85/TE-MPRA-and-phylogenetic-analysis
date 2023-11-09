##########
## Authors: Xun Chen
## Email: xunchen85@gmail.com or xchen@outlook.com
## Date: 2023/10/17
##########

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
cd input_tree_consensus/

for SEQ in `ls *.fa`
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

#### 3.3 liftOver analysis

