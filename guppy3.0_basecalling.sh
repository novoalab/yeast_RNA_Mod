# Basecalling of the yeast ribosomal RNA runs Using Guppy 3.0.3
## Performed this on CLUSTER 

#Location to the old guppy 3.0.3
guppy=/users/enovoa/joramirez/software/ont-guppy-cpu_3/bin. 



input=/no_backup_isis/enovoa/data/ont/yeast_rna/path/to/fast5
output=~/nanopore_analysis/ribosomal_rna_mods/path/to/output


guppy_basecaller --device cuda:0 -c rna_r9.4.1_70bps_hac.cfg --fast5_out --num_callers 4 --compress_fastq -ri $input -s $output

Paralel basecalling script is located in the scripts folder 

