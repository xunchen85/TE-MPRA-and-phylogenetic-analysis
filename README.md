# TE-MPRA-and-phylogenetic-analysis

The scripts here are used for the study of "Cryptic endogenous retrovirus subfamilies in the primate lineage". We could not upload the input files due to the file size limitation. However, the inputs are available upon request (chen.xun.3r@kyoto-u.ac.jp) and the final version will be submitted to Zenodo database at DOI:10.5281/zenodo.10016500.


All Python, R, and shell scripts will be finalized and uploaded asap...

1. Run_MPRA_analysis.sh script was used to compute the MPRA activity using MPRAflow, MPRAnalyze, and several R scripts.

2. Run_phylogenetic_analysis.sh script was used to perform the TE phylogenetic analysis based on multiple sequence alignment. This script contains the full steps how to determine the TE phyletic groups:


   •	1) obtain the subfamily groups;
   
   •	2) subdivide all instances (copies) into clusters per subfamily;

   •	3) determine phyletic groups among clusters per subfamily group based on constructed rooted trees and liftOver rate:
          3.1 we first examined the internal branch lengths of the tree and grouped adjacent clusters based on the top branch lengths (bootstrap values) manually. 
          3.2 We also examined the heatmap of divergence rates to look at extreme values between every adjacent clusters to validate the phyletic groups. 
          3.3 We kept the phyletic groups after we confirmed that the clusters from a phyletic group were evolutionary close to each other.

