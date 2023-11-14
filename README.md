# TE-MPRA-and-phylogenetic-analysis

## Introduction
Scripts included here are used for the study of "Cryptic endogenous retrovirus subfamilies in the primate lineage". We could not upload the input files due to the file size limitation. However, the inputs are available upon request (chen.xun.3r@kyoto-u.ac.jp) and the final version will be submitted to Zenodo database at DOI:10.5281/zenodo.10016500.

## Input data
All input files are kept under the "inputs/" folder (Zenodo database at DOI: 10.5281/zenodo.10016500).

## Shell scripts
1. Run_MPRA_analysis.sh: it includes the command lines and parameters to compute the MPRA activity from the raw FASTQ files using MPRAflow, MPRAnalyze, and R scripts.

2. Run_phylogenetic_analysis.sh: it contains the command lines and parameters to perform the TE phylogenetic analysis based on multiple sequence alignments. This script contains the full steps on how to determine the TE phyletic groups:

>•1) obtain the subfamily groups;
   
>•2) subdivide all instances (copies) into clusters per subfamily;

>•3) determine phyletic groups among clusters per subfamily group based on constructed rooted trees and liftOver rate:
   
>•3.1) Select the best root based on statistical tests of every trees and liftOver rate.
>•3.2) Examine the internal branch lengths of the tree and group adjacent clusters based on the top branch lengths (bootstrap values) manually. 
>•3.3) Examine the heatmap of divergence rates to look at extreme values between every adjacent clusters to validate the phyletic groups. 
>•3.4) Keep the phyletic groups after we confirmed that the clusters from a phyletic group were evolutionary close to each other.

## Python and perl scripts
These python and perl scripts are used in the shell scripts to prepare the inputs for the following R scripts
>•Exclude_emptyFasta.pl: exclude empty FASTA sequences or with ambiguous nucleotides only in the multiple sequence alignments.
>•makeTEgtf.pl: this script was obtained from https://github.com/mhammell-laboratory/TEtranscripts/issues/83.
>•Extract_by_TEfamilyname_each.py: extract TE instances in BED format using the names of subfamilies.
>•Extract_div_by_coordinate.py: extract the divergence rate of every TE instance using TE coordinates.
>•Organize_seqFile_to_consensus.py: reformat the downloaded repeatmasker alignment file relative to each TE consensus sequence.
>•Reformat_MPRAnalyze.py: convert the format of MPRAflow output files to the inputs for the MPRAnalyze.
>•Summarize_liftoverIntersect_by_TEFamily.py and Summarize_liftoverIntersect_by_TEFamily_Group.py: summarize the liftOver ratio per subfamily and per phyletic group.

## R scripts




