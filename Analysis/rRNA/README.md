# Analysis of Ribosomal RNA sequences

![alt text](../images/rrna/rrna_igv_image.png "rrna_igv")





#### 2.2. Mean current intensity plots centered in the modified sites
```
Rscript --vanilla nanopolish_meanlineplot.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```
Example using test data:

```
Rscript --vanilla nanopolish_meanlineplot.R test_data/sn34_window_file.tsv test_data/wt_window_file.tsv
```


![alt text](./img/mean_current.png "Mean_current")
