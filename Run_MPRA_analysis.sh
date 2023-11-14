#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2023/10/17
###
#######################

####################### MPRAflow analysis
# optimize the "nf_ori_map_barcodes.py" included in MPRAflow to keep reads fully matched to the designed insert sequences


########### step 1: association analysis
# 1.1
nextflow run association.nf -w TE_MPRA_i1_work_opt \
	--mapq 5 \
	--min-frac 0.5 \
	--cigar 140M \
	--fastq-insert "/home/xchen/software/MPRAflow/Assoc_Basic/data/TE_MPRA_i1_S1_R1_001.fastq.gz" \
	--fastq-insertPE "/home/xchen/software/MPRAflow/Assoc_Basic/data/TE_MPRA_i1_S1_R3_001.fastq.gz" \
	--fastq-bc "/home/xchen/software/MPRAflow/Assoc_Basic/data/TE_MPRA_i1_S1_R2_001.fastq.gz" \
	--design "/home/xchen/software/MPRAflow/Assoc_Basic/data/MPRA_combined_final_2021_6_10_unique_withoutAdapter.fa" \
	--name TE_MPRA_i1_assoc_basic_opt_identical5

# 1.2 pickle file conversion
python ConvertPickleToCsv.py -i TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle \
	-f pickle

########### step 2: counts analysis
# 2.1 MPRAflow count analysis for all samples, including iPSC_1, iPSC_2, iPSC_3, iPSC_4, NPC_1, NPC_2, NPC_3
nextflow run count.nf --dir "/home/xchen/software/MPRAflow/Count_Basic/TE_MPRA" \
        --experiment-file "/home/xchen/software/MPRAflow/Count_Basic/TE_MPRA/experiment.csv" \
        --design "/home/xchen/software/MPRAflow/Assoc_Basic/data/MPRA_combined_final_2021_6_10_unique_withoutAdapter.fa" \
        --association "/home/xchen/software/MPRAflow/outs/TE_MPRA_i1_assoc_basic_opt_identical5/TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle" \
        --outdir "/home/xchen/software/MPRAflow/outs/TE_MPRA_i1_counts_basic_opt_identical5_F" \
        -w "/home/xchen/software/MPRAflow/work_tmp_F" \
        --bc-length 15 \
        --umi-length 15 \
        --thresh 5 \
        --merge_intersect FALSE

########### step 3: outputs format conversion
# 3.1
python Insert_barcodes_inspection.py -a TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle.out \
	-d NPC_1_DNA_counts.tsv \
	-r NPC_1_RNA_counts.tsv >NPC_1_insert_BCs.counts
python Insert_barcodes_inspection.py -a TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle.out \
	-d NPC_2_DNA_counts.tsv \
	-r NPC_2_RNA_counts.tsv >NPC_2_insert_BCs.counts
python Insert_barcodes_inspection.py -a TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle.out \
	-d NPC_3_DNA_counts.tsv \
	-r NPC_3_RNA_counts.tsv >NPC_3_insert_BCs.counts
python Insert_barcodes_inspection.py -a TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle.out \
	-d iPSC_1_DNA_counts.tsv \
	-r iPSC_1_RNA_counts.tsv >iPSC_1_insert_BCs.counts
python Insert_barcodes_inspection.py -a TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle.out \
	-d iPSC_2_DNA_counts.tsv \
	-r iPSC_2_RNA_counts.tsv >iPSC_2_insert_BCs.counts
python Insert_barcodes_inspection.py -a TE_MPRA_i1_assoc_basic_opt_identical5_filtered_coords_to_barcodes.pickle.out \
	-d iPSC_3_DNA_counts.tsv \
	-r iPSC_3_RNA_counts.tsv >iPSC_3_insert_BCs.counts

# 3.2 obtain unique insert and barcodes
cat *_BCs.counts | grep -v '^Insert' | awk '{print$1,$3}' | sort | uniq >insert_BCs_2022_2_1.list
python ../Insert_barcodes_ids.py -i insert_BCs_2022_2_1.list >insert_BCs_2022_2_1.list.ids
awk '{print$1}' insert_BCs_2022_2_1.list.ids | grep -v '^Insert' | sort | uniq >insert_BCs_2022_2_1.list.uniqInserts
awk '{print$3}' insert_BCs_2022_2_1.list.ids |grep -v '^BCid' | sort -n | uniq >insert_BCs_2022_2_1.list.uniqBCids

# 3.3.1 NPC
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c NPC_1_insert_BCs.counts \
	-t 3 >NPC_1_insert_BCs.DNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c NPC_1_insert_BCs.counts \
	-t 4 >NPC_1_insert_BCs.RNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c NPC_2_insert_BCs.counts \
	-t 3 >NPC_2_insert_BCs.DNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c NPC_2_insert_BCs.counts \
	-t 4 >NPC_2_insert_BCs.RNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c NPC_3_insert_BCs.counts \
	-t 3 >NPC_3_insert_BCs.DNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c NPC_3_insert_BCs.counts \
	-t 4 >NPC_3_insert_BCs.RNAcounts

# 3.3.2 iPSC
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c iPSC_1_insert_BCs.counts \
	-t 3 >iPSC_1_insert_BCs.DNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts 
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c iPSC_1_insert_BCs.counts \
	-t 4 >iPSC_1_insert_BCs.RNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c iPSC_2_insert_BCs.counts \
	-t 3 >iPSC_2_insert_BCs.DNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c iPSC_2_insert_BCs.counts \
	-t 4 >iPSC_2_insert_BCs.RNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c iPSC_3_insert_BCs.counts \
	-t 3 >iPSC_3_insert_BCs.DNAcounts
python Reformat_MPRAnalyze.py -i insert_BCs_2022_2_1.list.uniqInserts \
	-b insert_BCs_2022_2_1.list.ids \
	-m 2161 \
	-c iPSC_3_insert_BCs.counts \
	-t 4 >iPSC_3_insert_BCs.RNAcounts

# 3.4
grep '[A-Z0-9a-z]' *_insert_BCs.counts | sed 's/:/\t/' >Combined.counts_2022_8_25
grep '[0-9A-Za-z]' *_[123]_counts.tsv | sed 's/_counts.tsv:/\t/' >Combined.ratio_2022_8_25

############### step 4: ran the MPRAnalyze.R script with the these converted count tables as the inputs to compute the MPRA activity (alpha values)


############### step 5: run the Prepare_summary_tables.R to prepare a summary table containing the LTR instance information and a table containing the sequence and MPRA activity
# MPRA activity (alpha values) was also normalized by the negative controls using the MAD.
# results were organized into these csv tables "Summary_Table1_2022_8_9.csv" and Summary_Table2_2022_8_9.csv"
