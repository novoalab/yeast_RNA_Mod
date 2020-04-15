# YEAST RIBOSOMAL RNA MODIFICATION ANALYSIS
Analysis of Direct RNA Sequencing of Yeast Ribosomal RNA 
These scripts are used for the generation of results in the study below:
"Decoding ribosomal RNA modifications in translating and non-translating ribosomes at single molecule resolution"


## DATASETS

### WT and Pseudouridylation KO Dataset

Run ID : RNA814001

|Barcodes| Strain        | Source           | Modification  | KO Positions |
|----| ------------- |:-------------:| -----:| -------: |
|BC1|    snR3  | Weizmann | Pseudouridylation | LSU:2129, 2133, 2264 |
|BC2|   snR34  | Weizmann | Pseudouridylation | LSU:2880, 2826 |
|BC3|    snR36  | Weizmann | Pseudouridylation | SSU: 1187 |
|BC4|   BY4741  | Weizmann/CRG | NA | NA |

### Normal Condition Sucrose Gradient Fractions Datasets

Run ID :  RNA442567, RNA92741

|Barcodes| Strain        | Condition           | rRNA population | Fraction number |
|----| ------------- |:-------------:| -----:| -------: |
|BC1|    BY4741  | Normal | Unassembled | 1,2 |
|BC2|   BY4741  | Normal | Small subunit + large subunit| 3,4,5,6 |
|BC3|    BY4741  | Normal | Lowly translating | 7,8,9,10 |
|BC4|   BY4741  | Normal | Highly translating | 12,13,14,15,16,17 |

### Normal and Stress (H2O2) Conditions Sucrose Gradient Fractions Datasets

Run ID :  RNA639991, RNA563572

|Barcodes| Strain        | Condition           | rRNA population | Fraction number |
|----| ------------- |:------------:| -----:| -------: |
|BC1|    BY4741  | Normal | Input | NA |
|BC2|   BY4741  | H2O2 (Stress) | Input | NA |
|BC3|    BY4741  | Normal | Lowly and Highly translating | 7-17 |
|BC4|   BY4741  | H2O2 (Stress)| Lowly and Highly translating | 7-17 |



[Information on the fractions can be found here](https://github.com/novoalab/yeast_rRNA_Mod/blob/master/scripts/yeast_rRNA_fractions.md)


## TOOLS NEEDED

Following R Packages are needed : 
plyr
stringr
reshape2
dplyr
ggplot2
ggbeeswarm
ggpubr



