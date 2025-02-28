# TE-MPRA-and-phylogenetic-analysis

## Introduction
Scripts included here are used for the study of "Cryptic endogenous retrovirus subfamilies in the primate lineage". 

## Input data
All input files are under the "input/", "input_trees" and "input_trees_consensus" folders (Zenodo database at DOI: 10.5281/zenodo.10016500). We could not upload all the input files under "input/" due to the file size limitation. However, all inputs are available upon request (xchen@siii.cas.cn) and will be submitted to Zenodo database at DOI:10.5281/zenodo.10016500.

## Shell scripts
### MPRA analysis
1) Run_MPRA_analysis.sh: it includes the command lines and parameters to compute the MPRA activity from the raw FASTQ files using MPRAflow, MPRAnalyze, and R scripts.

### Phylogenetic analysis
2) Run_phylogenetic_analysis.sh: it contains the command lines and parameters to perform the TE phylogenetic analysis based on multiple sequence alignments. This script contains the full steps on how to determine the TE phyletic groups:

#### 2.1 obtain the subfamily groups;
   
#### 2.2 subdivide all instances (copies) into clusters per subfamily;

#### 2.3 determine phyletic groups among clusters per subfamily group based on constructed rooted trees and liftOver rate:
   
•	2.3.1 Select the best root based on statistical tests of every trees and liftOver rate;
   
•	2.3.2 Examine the internal branch lengths of the tree and group adjacent clusters based on the top branch lengths (bootstrap values) manually. 
   
•	2.3.3 Examine the heatmap of divergence rates to look at extreme values between every adjacent clusters to validate the phyletic groups. 
   
•	2.3.4 Keep the phyletic groups after we confirmed that the clusters from a phyletic group were evolutionary close to each other.

## Python and perl scripts
These python and perl scripts are used in the shell scripts to prepare the inputs for the following R scripts

   •	Exclude_emptyFasta.pl: exclude empty FASTA sequences or with ambiguous nucleotides only in the multiple sequence alignments.

   •	makeTEgtf.pl: this script was obtained from https://github.com/mhammell-laboratory/TEtranscripts/issues/83.

   •	Extract_by_TEfamilyname_each.py: extract TE instances in BED format using the names of subfamilies.

   •	Extract_div_by_coordinate.py: extract the divergence rate of every TE instance using TE coordinates.

   •	Organize_seqFile_to_consensus.py: reformat the downloaded repeatmasker alignment file relative to each TE consensus sequence.

   •	Reformat_MPRAnalyze.py: convert the format of MPRAflow output files to the inputs for the MPRAnalyze.

   •	Summarize_liftoverIntersect_by_TEFamily.py and Summarize_liftoverIntersect_by_TEFamily_Group.py: summarize the liftOver ratio per subfamily and per phyletic group.

   •	nf_ori_map_barcodes_cx.py: extract the mappable reads with a variable length rather than a fixed CIGAR string. it need to replace the original nf_ori_map_barcodes.py included in the MPRAflow pipeline.
   

## R scripts (scripts for making supplementary figures are also availbale upon request at xchen@siii.cas.cn)
### Figure 1:
   •	1.LTR_subfamily_liftOver_analysis.R: Figure 1A; Figure 1B; Figure S1A; Figure S3A
   
   •	2.Tree_subdivision.R: Figure 1C; Figure S3C
   
   •	3.TEwide_phyletic_groups_visualization.R: Figure 1E; Figure 1F; Figure S11C
   
   •	4.TEwide_phyletic_groups_summary.R: Figure 1G; Figure S11A

### Figure 2:
   •	5.Phylo_regulatory_analysis-1.R: Figure 2A
   
   •	6.Phylo_regulatory_analysis-2.R: Figure 2B; Figure 2C; Figure 3F-annotation; Figure S5G
     
### Figure 3:
   •	7.TE_RPM_distribution_visualization.R: Figure 3A; Figure S5A
   
   •	8.Prepare_summary_tables-plots.R: Figure 3D; Figure S5C; Figure S5D
   
   •	6.Phylo_regulatory_analysis-2.R: Figure 3E; Figure 3F-annotation; Figure S5G
   
   •	9.Nucleotide_motif_association_analysis-4.R: Figure 3F; Figure S8C
   
   •	10.Motifs_and_activity_phyleticGroups.R: Figure 3G; Figure S5H; Figure S10A

### Figure 4:
   •	11.Liftover_analysis.R: Figure 4A; Figure S9A
   
   •	12.Macaque_tree_analysis-1.R: Figure 4B
   
   •	13.Macaque_tree_analysis-2.R: Figure 4C; Figure S9B; Figure S9C
   
   •	14.Macaque_tree_analysis-3.R: Figure 4D; Figure 4E; Figure 4F

### Figure 5:
   •	15.Nucleotide_motif_association_analysis-2.R: Figure 5A; Figure 5B; Figure 5C; Figure S6; Figure S7
   
   •	16.Nucleotide_motif_association_analysis-3.R: Figure 5D; Figure S8A; Figure S8B
   
   •	10.Motifs_and_activity_phyleticGroups.R: Figure 5F

### Figure 6:
   •	17.TE-wide_epigenetic_analysis.R: Figure 6A; Figure 6B; Figure 6C; Figure 6D; Figure S11C; Figure S11D

### Others:
   •	18.MPRA_activity_summary-final.R: Figure S5E; Figure S5F
   
   •	19.MPRAnalyze.R: run MPRAnalyze
   
   •	20.Prepare_summary_tables.R: prepare summary tables containing MPRA activity
   
   •	21.Nucleotide_motif_association_analysis-1.R: prepare the inputs for the motif and association analysis using plink2
   
   •	Figure 1D: manually revised based on the pdf file generated by the popart tool. The "hg19_rmsk_TE_0bp.MER11ABCD.MJN.aln.fa.phylip.nex" file under the input folder was used as the input
   
   •	Figure 3B: initially prepared by Fumitaka manually
   
   •	Figure 3C: manually prepared based on the length of TE consensus sequences and the locations of designed sequence frames
   
   •	Figure 5E: manually prepared based on the multiple sequence alignments of human "hg19_MER11ABC_v2.mafft.prank.best.fas.edited" and macaque "macFas5_MER11ABC_v2.mafft.prank.best.fas.edited" MER11 consensus sequences. Original MER11A/B/C consensus sequences were excluded


