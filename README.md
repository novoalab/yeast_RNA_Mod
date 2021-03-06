
# Analysis of Direct RNA Sequencing data for :
-  Saccharomyces Cerevisiae Ribosomal RNAs
-  Saccharomyces Cerevisiae sn/sno RNAs
-  Saccharomyces Cerevisiae mRNAs

Files related to yeast RNA analysis can be found [here](https://public-docs.crg.es/enovoa/public/begik_lucas_2021_NatBiotech_Yeast/) 



![alt text](./images/readme_image.png "init_fig")

## What's included:

### Tools 
-  [Epinano](https://github.com/novoalab/yeast_RNA_Mod/tree/master/Softwares) base frequency version in order to analyse direct RNA sequencing bam files
-  [Bam2Stats](https://github.com/novoalab/yeast_RNA_Mod/tree/master/Softwares) tool in order to analyse NanoCMC-Seq bam files]

### Bash scripts 
-  [Pre-process the files for R analysis](https://github.com/novoalab/yeast_RNA_Mod/tree/master/Analysis/Epinano)

### R scripts
-  [Process ribosomal RNA runs](https://github.com/novoalab/yeast_RNA_Mod/tree/master/Analysis/rRNA)
-  [Process sn/snoRNA runs](https://github.com/novoalab/yeast_RNA_Mod/tree/master/Analysis/ncRNA)
-  [Process mRNA runs](https://github.com/novoalab/yeast_RNA_Mod/tree/master/Analysis/mRNA)
-  [Process NanoCMC-Seq runs](https://github.com/novoalab/yeast_RNA_Mod/tree/master/Analysis/NanoCMCSeq)


## Dependencies/requirements: 
Following R Packages are needed to run the scripts: 
|R packages|
|----|
|plyr|
|stringr|
|reshape2|
|dplyr|
|ggplot2|
|ggbeeswarm|
|ggpubr|
|data.table|


## Citation:
Begik, O., Lucas, M.C., Pryszcz, L.P. et al. Quantitative profiling of pseudouridylation dynamics in native RNAs with nanopore sequencing. Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-00915-6


## Citation:
Please open an issue in the GitHub repo if you have any questions/doubts/suggestions about how to use this software. Thanks!
