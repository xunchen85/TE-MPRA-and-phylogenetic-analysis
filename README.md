# TE-MPRA-and-phylogenetic-analysis

The scripts here are used for the study of "Cryptic endogenous retrovirus subfamilies in the primate lineage". We could not upload the input files due to the file size limitation. However, the inputs are available upon request (chen.xun.3r@kyoto-u.ac.jp) and the final version will be submitted to Zenodo database at DOI:10.5281/zenodo.10016500.


All Python, R, and shell scripts will be finalized and uploaded asap...

1. Run_MPRAflow.sh script was used to compute the MPRA activity using MPRAflow, MPRAnalyze, and several R scripts.

2. Run_phylogenetic_analysis.sh script was used to perform the TE phylogenetic analysis based on multiple sequence alignment. It contains the full steps on how to determine the TE phyletic groups:


   •	1) obtain the subfamily groups;
   •	2) subdivide all instances (copies) into clusters per subfamily;
   •	3) determine phyletic groups among clusters per subfamily group based on the top selected rooted tree and liftOver rate;
