# TE-MPRA-and-phylogenetic-analysis

## Introduction
Scripts included here are used for the study of "Cryptic endogenous retrovirus subfamilies in the primate lineage". We could not upload the input files due to the file size limitation. However, the inputs are available upon request (chen.xun.3r@kyoto-u.ac.jp) and the final version will be submitted to Zenodo database at DOI:10.5281/zenodo.10016500.

## Input data
All input files are kept under the "inputs/" folder (Zenodo database at DOI: 10.5281/zenodo.10016500).

## Shell scripts
1. Run_MPRA_analysis.sh: it includes the command lines and parameters to compute the MPRA activity from the raw FASTQ files using MPRAflow, MPRAnalyze, and R scripts.

2. Run_phylogenetic_analysis.sh: it contains the command lines and parameters to perform the TE phylogenetic analysis based on multiple sequence alignments. This script contains the full steps on how to determine the TE phyletic groups:

   2.1 obtain the subfamily groups;
   
   2.2 subdivide all instances (copies) into clusters per subfamily;

   2.3 determine phyletic groups among clusters per subfamily group based on constructed rooted trees and liftOver rate:
   
      2.3.1 Select the best root based on statistical tests of every trees and liftOver rate;
   
      2.3.2 Examine the internal branch lengths of the tree and group adjacent clusters based on the top branch lengths (bootstrap values) manually. 
   
      2.3.3 Examine the heatmap of divergence rates to look at extreme values between every adjacent clusters to validate the phyletic groups. 
   
      2.3.4 Keep the phyletic groups after we confirmed that the clusters from a phyletic group were evolutionary close to each other.

## Python and perl scripts
These python and perl scripts are used in the shell scripts to prepare the inputs for the following R scripts

   1. Exclude_emptyFasta.pl: exclude empty FASTA sequences or with ambiguous nucleotides only in the multiple sequence alignments.

   2. makeTEgtf.pl: this script was obtained from https://github.com/mhammell-laboratory/TEtranscripts/issues/83.

   3. Extract_by_TEfamilyname_each.py: extract TE instances in BED format using the names of subfamilies.

   4. Extract_div_by_coordinate.py: extract the divergence rate of every TE instance using TE coordinates.

   5. Organize_seqFile_to_consensus.py: reformat the downloaded repeatmasker alignment file relative to each TE consensus sequence.

   6. Reformat_MPRAnalyze.py: convert the format of MPRAflow output files to the inputs for the MPRAnalyze.

   7. Summarize_liftoverIntersect_by_TEFamily.py and Summarize_liftoverIntersect_by_TEFamily_Group.py: summarize the liftOver ratio per subfamily and per phyletic group.

## R scripts
   1. LTR_subfamily_liftOver_analysis.R: Figure 1A
   2. Tree_subdivision.R: Figure 1C
   3. TEwide_phyletic_groups_visualization.R: Figure 1E; Figure 1F
   4. TEwide_phyletic_groups_summary.R: Figure 1G
   5. Phylo_regulatory_analysis-1.R: Figure 2A
   6. Phylo_regulatory_analysis-2.R: Figure 2B; Figure 2C
   7. TE_RPM_distribution_visualization.R: Figure 3A
   8. Prepare_summary_tables-plots.R: Figure 3D
   9. Phylo_regulatory_analysis-2.R: Figure 3E; Figure 3F; Figure S5G
   10. Motifs_and_activity_phyleticGroups.R: Figure 3G; Figure S5H
   11. Liftover_analysis.R: Figure 4A
   12. Macaque_tree_analysis-1.R: Figure 4B
   13. Macaque_tree_analysis-2.R: Figure 4C
   14. Macaque_tree_analysis-3.R: Figure 4D; Figure 4E; Figure 4F
   15. Nucleotide_motif_association_analysis-2.R: Figure 5A; Figure 5B; Figure 5C; Figure S6
   16. Nucleotide_motif_association_analysis-3.R: Figure 5D; Figure S8A; Figure S8B
   17. 



