
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
Begik O*, Lucas MC*, Ramirez JM, Milenkovic I, Cruciani S, Vieira HGS, Medina R, Liu H, Sas-Chen A, Mattick JS, Schwartz S and Novoa EM. Decoding ribosomal RNA modification dynamics at single molecule resolution. bioRxiv 2020. doi: https://doi.org/10.1101/2020.07.06.189969


## Citation:
Please open an issue in the GitHub repo if you have any questions/doubts/suggestions about how to use this software. Thanks!
