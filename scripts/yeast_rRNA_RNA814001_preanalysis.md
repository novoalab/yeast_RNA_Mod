# Yeast Ribosomal RNA 

## Pre-Analysis of RNA814001 Run (Pseudouridylation KO strains)

Here we only explain the following : 
	- Basecalling with Guppy (version 3.0.3)
	- Demultiplexing with Deeplexicon Final model, using MoP PreProcessSimple (0.95 Specificity threshold)
	- Mapping with different settings


### Demultiplexing filtering

```bash
	#Extracting the filtered read tags from the demultiplexing output
		
		module load Python/3.7.2-GCCcore-8.2.0 #Source the python3


   #Extract difference barcodes from the default demultiplexing
   		grep bc_1 RNA814001.all_demux.txt > bc1.reads.tsv
   		grep bc_2 RNA814001.all_demux.txt > bc2.reads.tsv
   		grep bc_3 RNA814001.all_demux.txt > bc3.reads.tsv
   		grep bc_4 RNA814001.all_demux.txt > bc4.reads.tsv

	#Exract the second column which contains read IDs
		cut -f2 bc1.reads.tsv > bc1.reads
		cut -f2 bc2.reads.tsv > bc2.reads
		cut -f2 bc3.reads.tsv> bc3.reads
		cut -f2 bc4.reads.tsv > bc4.reads


	#extract_sequence_from_fastq.py is used (placed in scripts folder)
		python3 ~/scripts/extract_sequence_from_fastq.py bc1.reads RNA814001.fastq > bc1.fastq
		python3 ~/scripts/extract_sequence_from_fastq.py bc2.reads RNA814001.fastq > bc2.fastq
		python3 ~/scripts/extract_sequence_from_fastq.py bc3.reads RNA814001.fastq > bc3.fastq
		python3 ~/scripts/extract_sequence_from_fastq.py bc4.reads RNA814001.fastq > bc4.fastq



   #Filtering the demultiplexing, using sort_probability_deviance.py (placed in scripts folder)
		python3.7 ~/scripts/sort_probability_deviance.py RNA927416.all_demux.txt 0.85

   #Extract difference barcodes from the default demultiplexing
   		grep bc_1 RNA814001.all_demux.txt.filt > bc1.reads.tsv.filt
   		grep bc_2 RNA814001.all_demux.txt.filt > bc2.reads.tsv.filt
   		grep bc_3 RNA814001.all_demux.txt.filt > bc3.reads.tsv.filt
   		grep bc_4 RNA814001.all_demux.txt.filt > bc4.reads.tsv.filt


	#Exract the second column which contains read IDs
		cut -f2 bc1.reads.tsv.filt > bc1.reads.filt
		cut -f2 bc2.reads.tsv.filt > bc2.reads.filt
		cut -f2 bc3.reads.tsv.filt > bc3.reads.filt
		cut -f2 bc4.reads.tsv.filt > bc4.reads.filt


	#extract_sequence_from_fastq.py is used (placed in scripts folder)
		python3 ~/scripts/extract_sequence_from_fastq.py bc1.reads.filt RNA814001.fastq > bc1.fastq.filt
		python3 ~/scripts/extract_sequence_from_fastq.py bc2.reads.filt RNA814001.fastq > bc2.fastq.filt
		python3 ~/scripts/extract_sequence_from_fastq.py bc3.reads.filt RNA814001.fastq > bc3.fastq.filt
		python3 ~/scripts/extract_sequence_from_fastq.py bc4.reads.filt RNA814001.fastq > bc4.fastq.filt


	#In order to count how many reads exist, divide the output of the script below by 4
		wc -l *fastq | awk '{x=$1/4; print x}'
```



### Mapping


```bash

#Minimap2 Default mapping
ref=/users/enovoa/boguzhan/references/yeast_rRNA_ref.fa

minimap2 -ax map-ont $ref bc1.fastq  -o bc1.mdef.bam
minimap2 -ax map-ont $ref bc2.fastq  -o bc2.mdef.bam
minimap2 -ax map-ont $ref bc3.fastq -o bc3.mdef.bam
minimap2 -ax map-ont $ref bc4.fastq  -o bc4.mdef.bam

#Sorting and indexing

samtools view  -Sb bc1.mdef.bam | samtools sort - bc1.mdef.sorted && samtools index bc1.mdef.sorted.bam
samtools view  -Sb bc2.mdef.bam | samtools sort - bc2.mdef.sorted && samtools index bc2.mdef.sorted.bam
samtools view  -Sb bc3.mdef.bam | samtools sort - bc3.mdef.sorted && samtools index bc3.mdef.sorted.bam
samtools view  -Sb bc4.mdef.bam | samtools sort - bc4.mdef.sorted && samtools index bc4.mdef.sorted.bam



#Minimap2 Sensitive mapping

minimap2 -ax map-ont -k9 -w4 -m20 -A3 -B1  $ref bc1.fastq  -o bc1.msens.bam
minimap2 -ax map-ont -k9 -w4 -m20 -A3 -B1  $ref bc2.fastq -o bc2.msens.bam
minimap2 -ax map-ont -k9 -w4 -m20 -A3 -B1  $ref bc3.fastq -o bc3.msens.bam
minimap2 -ax map-ont -k9 -w4 -m20 -A3 -B1  $ref bc4.fastq -o bc4.msens.bam

#Sorting and indexing

samtools view  -Sb bc1.msens.bam| samtools sort - bc1.msens.sorted && samtools index bc1.msens.sorted.bam
samtools view  -Sb bc2.msens.bam | samtools sort - bc2.msens.sorted && samtools index bc2.msens.sorted.bam
samtools view  -Sb bc3.msens.bam | samtools sort - bc3.msens.sorted && samtools index bc3.msens.sorted.bam
samtools view  -Sb bc4.msens.bam | samtools sort - bc4.msens.sorted && samtools index bc4.msens.sorted.bam




#U to T transformation

awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' bc1.fastq > bc1.reads_U2T.fastq
awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' bc2.fastq > bc2.reads_U2T.fastq
awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' bc3.fastq > bc3.reads_U2T.fastq
awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' bc4.fastq > bc4.reads_U2T.fastq


#Graphmap Default mapping

ref=/users/enovoa/boguzhan/references/yeast_rRNA_ref.fa

graphmap align -r $ref -d bc1.reads_U2T.fastq -o bc1.gdef.bam -v 1 -K fastq
graphmap align -r $ref -d bc2.reads_U2T.fastq -o bc2.gdef.bam -v 1 -K fastq
graphmap align -r $ref -d bc3.reads_U2T.fastq -o bc3.gdef.bam -v 1 -K fastq
graphmap align -r $ref -d bc4.reads_U2T.fastq -o bc4.gdef.bam -v 1 -K fastq

samtools view  -Sb bc1.gdef.bam | samtools sort - bc1.gdef.sorted && samtools index bc1.gdef.sorted.bam
samtools view  -Sb bc2.gdef.bam | samtools sort - bc2.gdef.sorted && samtools index bc2.gdef.sorted.bam
samtools view  -Sb bc3.gdef.bam | samtools sort - bc3.gdef.sorted && samtools index bc3.gdef.sorted.bam
samtools view  -Sb bc4.gdef.bam | samtools sort - bc4.gdef.sorted && samtools index bc4.gdef.sorted.bam

#Graphmap Sensitive mapping


ref=/users/enovoa/boguzhan/references/yeast_rRNA_ref.fa

graphmap align -r $ref -d bc1.reads_U2T.fastq -o bc1.gsen.bam -v 1 -K fastq --rebuild-index --double-index --mapq -1 -x sensitive -z -1 --min-read-len 0 -A 7 -k 5
graphmap align -r $ref -d bc2.reads_U2T.fastq -o bc2.gsen.bam -v 1 -K fastq --rebuild-index --double-index --mapq -1 -x sensitive -z -1 --min-read-len 0 -A 7 -k 5
graphmap align -r $ref -d bc3.reads_U2T.fastq -o bc3.gsen.bam -v 1 -K fastq --rebuild-index --double-index --mapq -1 -x sensitive -z -1 --min-read-len 0 -A 7 -k 5
graphmap align -r $ref -d bc4.reads_U2T.fastq -o bc4.gsen.bam -v 1 -K fastq --rebuild-index --double-index --mapq -1 -x sensitive -z -1 --min-read-len 0 -A 7 -k 5


samtools view  -Sb bc1.gsen.bam | samtools sort - bc1.gsen.sorted && samtools index bc1.gsen.sorted.bam
samtools view  -Sb bc2.gsen.bam | samtools sort - bc2.gsen.sorted && samtools index bc2.gsen.sorted.bam
samtools view  -Sb bc3.gsen.bam | samtools sort - bc3.gsen.sorted && samtools index bc3.gsen.sorted.bam
samtools view  -Sb bc4.gsen.bam | samtools sort - bc4.gsen.sorted && samtools index bc4.gsen.sorted.bam



	#Count the mapped reads

		for i in *.sorted.bam;do samtools view -F 4 $i | cut -f1 | sort | uniq | wc -l;echo $i; done
```



### Epinano Anaysis

```bash
sam2tsv=/users/enovoa/boguzhan/Software/jvarkit/dist/sam2tsv.jar
tsv_to_var=/users/enovoa/boguzhan/Software/EpiNano_NEW/scripts/TSV_to_Variants_Freq.py3
ref=/users/enovoa/boguzhan/references/yeast_rRNA_ref.fa

module load Python/3.7.2-GCCcore-8.2.0

#Calling variants for each single read-to-reference alignment
samtools view -h bc1.gdef.sorted.bam| java -jar $sam2tsv -r $ref  > sn3.bam.tsv 
samtools view -h bc2.gdef.sorted.bam| java -jar $sam2tsv -r $ref  > sn34.bam.tsv 
samtools view -h bc3.gdef.sorted.bam| java -jar $sam2tsv -r $ref  > sn36.bam.tsv 
samtools view -h bc4.gdef.sorted.bam| java -jar $sam2tsv -r $ref  > wt.bam.tsv


python3.7 $tsv_to_var -f sn3.bam.tsv -t 10 
python3.7 $tsv_to_var -f sn34.bam.tsv -t 10 
python3.7 $tsv_to_var -f sn36.bam.tsv -t 10 
python3.7 $tsv_to_var -f wt.bam.tsv -t 10
```

There are two output files :

  * per.site.var.csv : Contains features per site

 |Ref|pos|base|cov|q_mean|q_median|q_std|mis|ins|del|
 |---|---|---|---|---|---|---|---|---|---|
|18s|308|C|4460.0|15.11733|15.00000|6.55744|0.024|0.02|0.01|


  * per.site.var.per_site_var.5mer.csv : Contains features per site in a 5 mer context


|Kmer|Window|Ref|Coverage|q1|q2|q3|q4|q5|mis1|mis2|mis3|mis4|mis5|ins1|ins2|ins3|ins4|ins5|del1|del2|del3|del4|del5
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|GGTTG|2585:2586:2587:2588:2589|25s|6802.0:6806.0:6807.0:6806.0:6809.0|25.7|26.4|21.9|20.8|22.4|0.005|0.006|0.012|0.024|0.005|0.023|0.020|0.025|0.024|0.024|0.006|0.142|0.004|0.077|0.0145|