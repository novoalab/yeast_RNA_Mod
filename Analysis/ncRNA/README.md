# Analysis of ncRNA runs

# Dot-plots Y sites based on normal and stress conditions
```
Rscript ncRNA_stress_dotplot.R <normal_rep1_epinano> <stress_rep1_epinano> <normal_rep2_epinano> <stress_rep1_epinano> 
```
Example using test data:

```
Rscript ncRNA_stress_dotplot.R test_data/normal_rep1_epinano.csv test_data/heat_rep1_epinano.csv test_data/normal_rep2_epinano.csv test_data/heat_rep1_epinano.csv
```