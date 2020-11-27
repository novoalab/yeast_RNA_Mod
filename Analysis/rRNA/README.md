# Analysis of Ribosomal RNA sequences

![alt text](../../images/rrna/rrna_igv_image.png "rrna_igv")

## 1. Base-called features of rRNA modifications

### Dot-plots of each distinct type of modification (Figure 2B)
```
Rscript mods_dotplot.R <epinano_5mer.csv> all_rrna_mod_status.tsv
```
Example using test data:

```
Rscript mods_dotplot.R test_data/wt_epinano_5mer.csv all_rrna_mod_status.tsv
```
<img src="../../images/rrna/dotplot_example.png " width="600">

### Ternary plots (base-frequency) of each distinct type of modification (Figure 2C)
```
Rscript mods_ternary.R <epinano.csv> all_rrna_mod_status.tsv
```
Example using test data:

```
Rscript mods_ternary.R test_data/wt_epinano.csv all_rrna_mod_status.tsv
```

<img src="../../images/rrna/ternary_example.png " width="400"> 


### 5-mer dotplots for each mods
```
Rscript dotplot_5mer <epinano_5mer.csv> rrna_mod_5mer.tsv
```
Example using test data:

```
Rscript dotplot_5mer wt.bam.tsv.per.site.var.per_site_var.5mer.csv rrna_mod_5mer.tsv
```

<img src="../../images/rrna/5mer_dotplot_example.png " width="1000"> 