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

![alt text](../../images/rrna/dotplot_example.png "dotplot_example")


### Ternary plots (base-frequency) of each distinct type of modification (Figure 2C)
```
Rscript mods_ternary.R <epinano.csv> all_rrna_mod_status.tsv
```
Example using test data:

```
Rscript mods_ternary.R test_data/wt_epinano.csv all_rrna_mod_status.tsv
```

![alt text](../../images/rrna/ternary_example.png "ternary _example")

<img src="../../images/rrna/ternary_example.png " width="200">