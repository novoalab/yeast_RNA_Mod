# Yeast Ribosomal RNA 

## Nanopolish Analysis of RNA814001 Run (Pseudouridylation KO strains)

### Processing using Nanopolish Tool

### Pre-processing the Nanopolish Output (Event Align) for per-read/pos 

#### Process the raw eventalign output in order to clean it a bit (Take mean of the multiple observations per-read/pos)
```bash
#Load the module of Python
module load Python/3.7.2-GCCcore-8.2.0
```

```Python
import pandas as pd

bc1=pd.read_csv("bc1.reads-ref.eventalign.txt",sep='\t')
bc2=pd.read_csv("bc2.reads-ref.eventalign.txt",sep='\t')
bc3=pd.read_csv("bc3.reads-ref.eventalign.txt",sep='\t')
bc4=pd.read_csv("bc4.reads-ref.eventalign.txt",sep='\t')

#For the simple figure of (per mean/read-pos) 
grouped_multiple_mean_bc1 = bc1.groupby(['contig', 'position','reference_kmer','read_index']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc1 = grouped_multiple_mean_bc1.reset_index()
grouped_multiple_mean_bc1.columns =  grouped_multiple_mean_bc1.columns.droplevel(-1)
grouped_multiple_mean_bc1.to_csv(r'bc1_processed_perread_mean.tsv', sep='\t', index = False)

grouped_multiple_mean_bc2 = bc2.groupby(['contig', 'position','reference_kmer','read_index']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc2 = grouped_multiple_mean_bc2.reset_index()
grouped_multiple_mean_bc2.columns =  grouped_multiple_mean_bc2.columns.droplevel(-1)
grouped_multiple_mean_bc2.to_csv(r'bc2_processed_perread_mean.tsv', sep='\t', index = False)

grouped_multiple_mean_bc3 = bc3.groupby(['contig', 'position','reference_kmer','read_index']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc3 = grouped_multiple_mean_bc3.reset_index()
grouped_multiple_mean_bc3.columns =  grouped_multiple_mean_bc3.columns.droplevel(-1)
grouped_multiple_mean_bc3.to_csv(r'bc3_processed_perread_mean.tsv', sep='\t', index = False)

grouped_multiple_mean_bc4 = bc4.groupby(['contig', 'position','reference_kmer','read_index']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc4 = grouped_multiple_mean_bc4.reset_index()
grouped_multiple_mean_bc4.columns =  grouped_multiple_mean_bc4.columns.droplevel(-1)
grouped_multiple_mean_bc4.to_csv(r'bc4_processed_perread_mean.tsv', sep='\t', index = False)
```

#### Process the raw eventalign output in order to collapse all the reads per position

```bash
#Load the module of Python
module load Python/3.7.2-GCCcore-8.2.0
```
```Python
import pandas as pd

bc1=pd.read_csv("bc1.reads-ref.eventalign.txt",sep='\t')
bc2=pd.read_csv("bc2.reads-ref.eventalign.txt",sep='\t')
bc3=pd.read_csv("bc3.reads-ref.eventalign.txt",sep='\t')
bc4=pd.read_csv("bc4.reads-ref.eventalign.txt",sep='\t')

#For the simple figure of (per median/position) 
grouped_multiple_mean_bc1 = bc1.groupby(['contig', 'position','reference_kmer']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc1 = grouped_multiple_mean_bc1.reset_index()
grouped_multiple_mean_bc1.columns =  grouped_multiple_mean_bc1.columns.droplevel(-1)
grouped_multiple_mean_bc1.to_csv(r'bc1_processed_perpos_mean.tsv', sep='\t', index = False)

grouped_multiple_mean_bc2 = bc2.groupby(['contig', 'position','reference_kmer']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc2 = grouped_multiple_mean_bc2.reset_index()
grouped_multiple_mean_bc2.columns =  grouped_multiple_mean_bc2.columns.droplevel(-1)
grouped_multiple_mean_bc2.to_csv(r'bc2_processed_perpos_mean.tsv', sep='\t', index = False)

grouped_multiple_mean_bc3 = bc3.groupby(['contig', 'position','reference_kmer']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc3 = grouped_multiple_mean_bc3.reset_index()
grouped_multiple_mean_bc3.columns =  grouped_multiple_mean_bc3.columns.droplevel(-1)
grouped_multiple_mean_bc3.to_csv(r'bc3_processed_perpos_mean.tsv', sep='\t', index = False)

grouped_multiple_mean_bc4 = bc4.groupby(['contig', 'position','reference_kmer']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc4 = grouped_multiple_mean_bc4.reset_index()
grouped_multiple_mean_bc4.columns =  grouped_multiple_mean_bc4.columns.droplevel(-1)
grouped_multiple_mean_bc4.to_csv(r'bc4_processed_perpos_mean.tsv', sep='\t', index = False)
```


#### Modify the perpos data, so that we can use it to create sliding window
```R
bc1<-read.delim("bc1_processed_perpos_mean.tsv")
bc1$readidx<- rep("1", nrow(bc1))
bc1_2<- bc1[,c(1,2,3,5,4)]
write.table(bc1_2, file="bc1_processed_perpos_modified.tsv", quote=FALSE, row.names=FALSE)

bc2<-read.delim("bc2_processed_perpos_mean.tsv")
bc2$readidx<- rep("1", nrow(bc2))
bc2_2<- bc2[,c(1,2,3,5,4)]
write.table(bc2_2, file="bc2_processed_perpos_modified.tsv", quote=FALSE, row.names=FALSE)

bc3<-read.delim("bc3_processed_perpos_mean.tsv")
bc3$readidx<- rep("1", nrow(bc3))
bc3_2<- bc3[,c(1,2,3,5,4)]
write.table(bc3_2, file="bc3_processed_perpos_modified.tsv", quote=FALSE, row.names=FALSE)

bc4<-read.delim("bc4_processed_perpos_mean.tsv")
bc4$readidx<- rep("1", nrow(bc4))
bc4_2<- bc4[,c(1,2,3,5,4)]
write.table(bc4_2, file="bc4_processed_perpos_modified.tsv", quote=FALSE, row.names=FALSE)

```

#### Modify the perpos data, so that we can use it to create sliding window


```bash
module load Python/3.7.2-GCCcore-8.2.0
python /users/enovoa/hliu/SHARE/huanle_useful_snippets/src/oz_sliding_win/5.Sliding_event_aln_tabl_on_ref_with_reads_info.py -f bc1_processed_perpos_modified.tsv -w 15 > bc1_perpos_sliding15.tsv
python /users/enovoa/hliu/SHARE/huanle_useful_snippets/src/oz_sliding_win/5.Sliding_event_aln_tabl_on_ref_with_reads_info.py -f bc2_processed_perpos_modified.tsv -w 15 > bc2_perpos_sliding15.tsv
python /users/enovoa/hliu/SHARE/huanle_useful_snippets/src/oz_sliding_win/5.Sliding_event_aln_tabl_on_ref_with_reads_info.py -f bc3_processed_perpos_modified.tsv -w 15 > bc3_perpos_sliding15.tsv
python /users/enovoa/hliu/SHARE/huanle_useful_snippets/src/oz_sliding_win/5.Sliding_event_aln_tabl_on_ref_with_reads_info.py -f bc4_processed_perpos_modified.tsv -w 15 > bc4_perpos_sliding15.tsv
```



### Plotting the processed outputs

#### Density plots for the single positions


```R
#Libraries needed
library(dplyr)
library(ggplot2)

#Import the processed event align files
sn3<- read.csv("bc1_processed_perread_var.tsv", sep="\t")
sn34<- read.csv("bc2_processed_perread_var.tsv", sep="\t")
sn36<- read.csv("bc3_processed_perread_var.tsv", sep="\t")
wt<- read.csv("bc4_processed_perread_var.tsv", sep="\t"

#Add 3 to the positions because event align results are shifted by 2 nts and its zero based
sn3$position<- sn3$position+3
sn34$position<- sn34$position+3
sn36$position<- sn36$position+3
wt$position<- wt$position+3

#Create a column for unique positions
sn3$Pos<- paste(sn3$contig, sn3$position, sep="_")
sn34$Pos<- paste(sn34$contig, sn34$position, sep="_")
sn36$Pos<- paste(sn36$contig, sn36$position, sep="_")
wt$Pos<- paste(wt$contig, wt$position, sep="_")

#Add a column for the sample information
sn3$sample<- rep("sn3", nrow(sn3))
sn34$sample<- rep("sn34", nrow(sn34))
sn36$sample<- rep("sn36", nrow(sn36))
wt$sample<- rep("wt", nrow(wt))

#Merge all the tables
all<- rbind(wt,sn3,sn34,sn36)

#Import the table to use for 5mer extraction
mod<- read.delim("pU_5mer.tsv")

#For only Y positions in ALL strains
for (pos in unique(mod$Pos)) {
	subs<- subset(all, Pos==pos)
	pdf(file=paste(pos, "_Y_pos_all_strains.pdf", sep="."),height=4,width=9,onefile=FALSE)
		print(ggplot(subs, aes(x= event_level_mean, fill=sample,color=sample)) +
 		geom_density(alpha=0.3,adjust = 2)+
  		theme_bw()+
		ggtitle(paste(pos, "Y"))+
		xlab("Current Intensity (pA) ")+
      	ylab("Density") +
		theme(axis.text=element_text(size=15),strip.text = element_text(size=13),
    		axis.title=element_text(size=17,face="bold"),
    		legend.title = element_text(size = 20),
    		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
    		legend.text = element_text(color = "black", size=15)))
	dev.off()
}

#Create a table for only Y positions
#Join the table with mod positions/5mers by Pos
pu_all<- join(all,mod, by="Pos")
#Remove na rows
pu_all<- na.omit(pu_all)
#Remove redundant column
pu_all_2<-pu_all[,-12]
#Create a column that will give the base information
pu_all_2$base<- rep("Y",nrow(pu_all_2))


#Create a column in mod file in order to use for the join function
mod$reference_kmer<- mod$Kmer
#Join the table with mod positions/5mers by reference 5mer
u_all<- join(all,mod, by="reference_kmer")
#Remove all the NAs
u_all_2<- na.omit(u_all)
#Now join the table with mod file by Pos this time,
u_all_3<- join(u_all_2,mod, by="Pos")
#Keep the columns that ARE NA (Which is unmodified 5mer)
u_all_4 <- u_all_3[rowSums(is.na(u_all_3)) > 0,]
#Remove redundant columns
u_all_5<- u_all_4[,-c(10,13,14,15,16,17,18,19)]
#Create a column that will give the base information
u_all_5$base<- rep("U",nrow(u_all_5))

#Create a clumn in mod file in order to use for join function
mod2<- mod
mod2$reference_kmer<- mod$C.Kmer


c_all<- join(all,mod2, by="reference_kmer")
c_all_2<- na.omit(c_all)
c_all_3<-c_all_2[,-10]
c_all_3<-c_all_3[,-12]
c_all_3$base<- rep("C",nrow(c_all_3))


all_all<- rbind(pu_all_2,u_all_5, c_all_3)

sn3_all<- subset(all_all, sample=="wt" | sample=="sn3")
sn34_all<- subset(all_all, sample=="wt" | sample=="sn34")
sn36_all<- subset(all_all, sample=="wt" | sample=="sn36")


for (pos in unique(sn3_all$Posn.) ) {
	subs<- subset(sn3_all, Posn.==pos)
	subs2<-subset(sn3_all, Kmer ==as.character(unique(subs$Kmer)))
	subs3<-subset(sn3_all, Posn.==pos)
	library(ggplot2)
		pdf(file=paste(pos, "sn3_wt_with_U_Cmer.pdf", sep="."),height=3,width=10,onefile=FALSE)
			print(ggplot(subs3, aes(x= event_level_mean, fill=base,color=sample)) +
 			geom_density(alpha=0.3,adjust = 2)+
  			theme_bw()+
  			ggtitle(paste(pos, "Sn3KO"))+
			xlab("Current Intensity (pA) ")+
      		ylab("Density") +
  			theme(axis.text=element_text(size=15),strip.text = element_text(size=13),
    		axis.title=element_text(size=17,face="bold"),
    		legend.title = element_text(size = 20),
    		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
    		legend.text = element_text(color = "black", size=15)))
  		dev.off()
  	}




for (pos in unique(sn34_all$Posn.) ) {
	subs<- subset(sn34_all, Posn.==pos)
	subs2<-subset(sn34_all, Kmer ==as.character(unique(subs$Kmer)))
	subs3<-subset(sn34_all, Posn.==pos)
	library(ggplot2)
		pdf(file=paste(pos, "sn34_wt_with_U_Cmer.pdf", sep="."),height=3,width=10,onefile=FALSE)
			print(ggplot(subs3, aes(x= event_level_mean, fill=base,color=sample)) +
 			geom_density(alpha=0.3,adjust = 2)+
  			theme_bw()+
  			ggtitle(paste(pos, "Sn34KO"))+
			xlab("Current Intensity (pA) ")+
      		ylab("Density") +
  			theme(axis.text=element_text(size=15),strip.text = element_text(size=13),
    		axis.title=element_text(size=17,face="bold"),
    		legend.title = element_text(size = 20),
    		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
    		legend.text = element_text(color = "black", size=15)))
  		dev.off()
  	}



for (pos in unique(sn36_all$Posn.) ) {
	subs<- subset(sn36_all, Posn.==pos)
	subs2<-subset(sn36_all, Kmer ==as.character(unique(subs$Kmer)))
	subs3<-subset(sn36_all, Posn.==pos)
	library(ggplot2)
		pdf(file=paste(pos, "sn36_wt_with_U_Cmer.pdf", sep="."),height=3,width=10,onefile=FALSE)
			print(ggplot(subs3, aes(x= event_level_mean, fill=base,color=sample)) +
 			geom_density(alpha=0.3,adjust = 2)+
  			theme_bw()+
  			ggtitle(paste(pos, "Sn36KO"))+
			xlab("Current Intensity (pA) ")+
      		ylab("Density") +
  			theme(axis.text=element_text(size=15),strip.text = element_text(size=13),
    		axis.title=element_text(size=17,face="bold"),
    		legend.title = element_text(size = 20),
    		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
    		legend.text = element_text(color = "black", size=15)))
  		dev.off()
  	}
```

#### Scatter plots for the single position

```R
  #Scatter Plots
    library(plyr)
    library(ggplot2)
    library(ggrepel)


    #Import the data
      ##Import the table pre-made for the position labelling
      ##To label mod position and their neighbours
      status<- read.delim("neighbour_for_current.tsv")
      #Create a column for unique positions
      status$pos<- paste(status$Chr,status$Position)

      #Import the Per position table for Sn3
      sn3ko <- read.delim("bc1_processed_perpos_mean.tsv")
      #Add 3 nt to each position since Nanopolish output is 0 based and falls behind 2 nts
      sn3ko$position<- sn3ko$position+3 
      #Create a column for unique positions
      sn3ko$pos<- paste(sn3ko$contig, sn3ko$position)
      #Add a column replicating sn3ko
      sn3ko$sample <- rep("sn3ko",nrow(sn3ko)) 
      #use join function to add count number to the data
      sn3ko_st<- join(sn3ko, status, by="pos") 
      #Rename the columns
      colnames(sn3ko_st)<- c("contig","position", "kmer", "sn3ko_event_level_mean_mean", "pos", "sample", "Chr", "Position", "ModStatus","Status", "Neighbour")

      #Import the Per position table for Sn3
      sn34ko <- read.delim("bc2_processed_perpos_mean.tsv")
      #Add 3 nt to each position since Nanopolish output is 0 based and falls behind 2 nts
      sn34ko$position<- sn34ko$position+3 
      #Create a column for unique positions
      sn34ko$pos<- paste(sn34ko$contig, sn34ko$position)
      #Add a column replicating sn34ko
      sn34ko$sample <- rep("sn34ko",nrow(sn34ko)) 
      #use join function to add count number to the data
      sn34ko_st<- join(sn34ko, status, by="pos") 
      #Rename the columns
      colnames(sn34ko_st)<- c("contig","position", "kmer", "sn34ko_event_level_mean_mean", "pos", "sample", "Chr", "Position", "ModStatus","Status", "Neighbour")

      #Import the Per position table for Sn3
      sn36ko <- read.delim("bc3_processed_perpos_mean.tsv")
      #Add 3 nt to each position since Nanopolish output is 0 based and falls behind 2 nts
      sn36ko$position<- sn36ko$position+3 
      #Create a column for unique positions
      sn36ko$pos<- paste(sn36ko$contig, sn36ko$position)
      #Add a column replicating sn36ko
      sn36ko$sample <- rep("sn36ko",nrow(sn36ko)) 
      #use join function to add count number to the data
      sn36ko_st<- join(sn36ko, status, by="pos") 
      #Rename the columns
      colnames(sn36ko_st)<- c("contig","position", "kmer", "sn36ko_event_level_mean_mean", "pos", "sample", "Chr", "Position", "ModStatus","Status", "Neighbour")


      # Use join function to merge the strains
      wt_sn3<- join(wt_st, sn3ko_st, by="pos") 
      wt_sn34<- join(wt_st, sn34ko_st, by="pos")
      wt_sn36<- join(wt_st, sn36ko_st, by="pos")


      #Remove the columns that do not match
      wt_sn3<- subset(wt_sn3, sn3ko_event_level_mean_mean!="NA")
      wt_sn34<- subset(wt_sn34, sn34ko_event_level_mean_mean!="NA")
      wt_sn36<- subset(wt_sn36, sn36ko_event_level_mean_mean!="NA")

      #Remove the 20 nt from 5'UTR
      wt_sn3<- subset(wt_sn3, position>20)
      wt_sn34<- subset(wt_sn34, position>20)
      wt_sn36<- subset(wt_sn36, position>20)

      #Calculate the difference 
      wt_sn3$diff<- abs(wt_sn3$wt_event_level_mean_mean - wt_sn3$sn3ko_event_level_mean_mean)
      wt_sn34$diff<- abs(wt_sn34$wt_event_level_mean_mean -  wt_sn34$sn34ko_event_level_mean_mean)
      wt_sn36$diff<- abs(wt_sn36$wt_event_level_mean_mean - wt_sn36$sn36ko_event_level_mean_mean)


      #Plot for sn3KO
      for (rna in unique(wt_sn3$contig)) {
      subs<-subset(wt_sn3, contig==rna)  
      pdf(file=paste(rna, "wt_sn3_event_mean_xyplot.pdf", sep="_"),height=5,width=5,onefile=FALSE)
      print(ggplot(subs, aes(x=sn3ko_event_level_mean_mean, y=wt_event_level_mean_mean)) +
        geom_point(size=1, color="grey")+
        geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black")+
        geom_point(data=subset(subs, diff>0 & Neighbour == "Mod") , size=2, color="red")+
        geom_point(data=subset(subs, diff>0 & Neighbour == "Yes") , size=2, color="blue")+
        geom_text_repel(data=subset(subs, diff>0), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
        ggtitle(paste(rna,"snR3-KO"))+
        xlab("Mean Current Intensity (WT)")+
        ylab("Mean Current Intensity (KO)") +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold", color="black",size=11),
          axis.text.y = element_text(face="bold", color="black", size=11),
          plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
          axis.title.x = element_text(color="black", size=15, face="bold"),
          axis.title.y = element_text(color="black", size=15, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title = element_text(color = "black", size = 20,face="bold"),
          legend.text = element_text(color = "black", size=20)))
      dev.off()
      }


      #Plot for sn34KO
      for (rna in unique(wt_sn34$contig)) {
      subs<-subset(wt_sn34, contig==rna)   
      pdf(file=paste(rna, "wt_sn34_event_mean_xyplot.pdf", sep="_"),height=5,width=5,onefile=FALSE)
      print(ggplot(subs, aes(x=sn34ko_event_level_mean_mean, y=wt_event_level_mean_mean)) +
        geom_point(size=1, color="grey")+
        geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black")+
        geom_point(data=subset(subs, diff>0 & Neighbour == "Mod") , size=2, color="red")+
        geom_point(data=subset(subs, diff>0 & Neighbour == "Yes") , size=2, color="blue")+
        geom_text_repel(data=subset(subs, diff>0), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
        ggtitle(paste(rna,"snR34-KO"))+
        xlab("Mean Current Intensity (WT)")+
        ylab("Mean Current Intensity (KO)") +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold", color="black",size=11),
          axis.text.y = element_text(face="bold", color="black", size=11),
          plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
          axis.title.x = element_text(color="black", size=15, face="bold"),
          axis.title.y = element_text(color="black", size=15, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title = element_text(color = "black", size = 20,face="bold"),
          legend.text = element_text(color = "black", size=20)))
      dev.off()
      }


      #Plot for sn36KO
      for (rna in unique(wt_sn36$contig)) {
      subs<-subset(wt_sn36, contig==rna)   
      pdf(file=paste(rna, "wt_sn36_event_mean_xyplot.pdf", sep="_"),height=5,width=5,onefile=FALSE)
      print(ggplot(subs, aes(x=sn36ko_event_level_mean_mean, y=wt_event_level_mean_mean)) +
        geom_point(size=1, color="grey")+
        geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black")+
        geom_point(data=subset(subs, diff>0 & Neighbour == "Mod") , size=2, color="red")+
        geom_point(data=subset(subs, diff>0 & Neighbour == "Yes") , size=2, color="blue")+
        geom_text_repel(data=subset(subs, diff>0), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
        ggtitle(paste(rna,"snR36-KO"))+
        xlab("Mean Current Intensity (WT)")+
        ylab("Mean Current Intensity (KO)") +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold", color="black",size=11),
          axis.text.y = element_text(face="bold", color="black", size=11),
          plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
          axis.title.x = element_text(color="black", size=15, face="bold"),
          axis.title.y = element_text(color="black", size=15, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title = element_text(color = "black", size = 20,face="bold"),
          legend.text = element_text(color = "black", size=20)))
      dev.off()
      }
    ```

### Scatter plots for the single position

```R
library(plyr)
library(ggrepel)
library(dplyr)
library(stringr)
#Importing the tables processed in python (median event value per position)
sn3<- read.delim("bc1_processed_perpos_mean.tsv")
sn34<- read.delim("bc2_processed_perpos_mean.tsv")
sn36<- read.delim("bc3_processed_perpos_mean.tsv")
wt<- read.delim("bc4_processed_perpos_mean.tsv")


#Add 3 to the position, because thats how we built the modification tables
sn3$position<- sn3$position+3
sn34$position<- sn34$position+3
sn36$position<- sn36$position+3
wt$position<- wt$position+3


#Include a column that contains chr-position informatio
sn3$ref<- paste(sn3$contig, sn3$position)
sn34$ref<- paste(sn34$contig, sn34$position)
sn36$ref<- paste(sn36$contig, sn36$position)
wt$ref<- paste(wt$contig, wt$position)

#Rename the columns
colnames(wt)<- c("contig", "position","reference_kmer","wt_event_level_mean", "ref")
colnames(sn3)<- c("contig", "position","reference_kmer","sn3_event_level_mean", "ref")
colnames(sn34)<- c("contig", "position","reference_kmer","sn34_event_level_mean", "ref")
colnames(sn36)<- c("contig", "position","reference_kmer","sn36_event_level_mean", "ref")

#Create a column with base information
wt$base<- substring(wt$reference_kmer, 3,3)
sn3$base<- substring(sn3$reference_kmer, 3,3)
sn34$base<- substring(sn34$reference_kmer, 3,3)
sn36$base<- substring(sn36$reference_kmer, 3,3)

#Join the files by the chr-position information
wt_sn3<-join(wt,sn3, by="ref")
wt_sn34<-join(wt,sn34, by="ref")
wt_sn36<-join(wt,sn36, by="ref")

#remove the NAs
wt_sn3<- na.omit(wt_sn3)
wt_sn34<- na.omit(wt_sn34)
wt_sn36<- na.omit(wt_sn36)

#calculate delta current intensity
wt_sn3$diff<-abs(wt_sn3$sn3_event_level_mean- wt_sn3$wt_event_level_mean)
wt_sn34$diff<-abs(wt_sn34$sn34_event_level_mean- wt_sn34$wt_event_level_mean)
wt_sn36$diff<-abs(wt_sn36$sn36_event_level_mean- wt_sn36$wt_event_level_mean)


###### 25s 
## sn3 
wt_sn3_25s<-subset(wt_sn3, contig=="25s")
    pdf(file="wt_sn3_25s_delta.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn3_25s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_text_repel(data=subset(wt_sn3_25s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 0.4,segment.color = "grey50")+
          geom_vline(xintercept = 2129, color = "#f64b3c", alpha=0.4, size=0.3)+
          geom_vline(xintercept = 2133, color = "#f64b3c", alpha=0.4, size=0.3)+
          geom_vline(xintercept = 2264, color = "#f64b3c", alpha=0.4, size=0.3)+
          scale_x_continuous(limits=c(50, 3300), breaks=c(50, 400, 800, 1200, 1600,2000,2400,2800,3200))+
        ggtitle("25s rRNA snR3-KO")+
        xlab("Positions")+
        ylab("Delta Current Intensity (WT-KO)") +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()



    pdf(file="wt_sn3_25s_delta.zoomed.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn3_25s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_text_repel(data=subset(wt_sn3_25s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 1,segment.color = "black")+
          geom_vline(xintercept = 2129, color = "#f64b3c", alpha=0.4, size=2)+
          geom_vline(xintercept = 2133, color = "#f64b3c", alpha=0.4, size=2)+
          geom_vline(xintercept = 2264, color = "#f64b3c", alpha=0.4, size=2)+
          scale_x_continuous(limits=c(2100, 2300), breaks=c(2100, 2150, 2200, 2250, 2300))+
          ggtitle("25s rRNA snR3-KO")+
          xlab("Positions")+
          ylab("Delta Current Intensity (WT-KO)") +
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()

## sn34
wt_sn34_25s<-subset(wt_sn34, contig=="25s")
    
    pdf(file="wt_sn34_25s_delta.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn34_25s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_text_repel(data=subset(wt_sn34_25s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 1,segment.color = "black")+
          geom_vline(xintercept = 2826, color = "#f64b3c", alpha=0.4, size=0.3)+
          geom_vline(xintercept = 2880, color = "#f64b3c", alpha=0.4, size=0.3)+
          scale_x_continuous(limits=c(50, 3300), breaks=c(50, 400, 800, 1200, 1600,2000,2400,2800,3200))+
          ggtitle("25s rRNA snR34-KO")+
          xlab("Positions")+
          ylab("Delta Current Intensity (WT-KO)") +
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()


    pdf(file="wt_sn34_25s_delta.zoomed.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn34_25s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_text_repel(data=subset(wt_sn34_25s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 1,segment.color = "black")+
          geom_vline(xintercept = 2826, color = "#f64b3c", alpha=0.4, size=3)+
          geom_vline(xintercept = 2880, color = "#f64b3c", alpha=0.4, size=3)+
          scale_x_continuous(limits=c(2800, 2900), breaks=c(2800, 2850, 2900))+
          ggtitle("25s rRNA snR34-KO")+
          xlab("Positions")+
          ylab("Delta Current Intensity (WT-KO)") +
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()


## sn36
wt_sn36_25s<-subset(wt_sn36, contig=="25s")
    pdf(file="wt_sn36_25s_delta.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn36_25s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_text_repel(data=subset(wt_sn36_25s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 1,segment.color = "black")+
          scale_x_continuous(limits=c(50, 3300), breaks=c(50, 400, 800, 1200, 1600,2000,2400,2800,3200))+
           ggtitle("25s rRNA snR36-KO")+
          xlab("Positions")+
          ylab("Delta Current Intensity (WT-KO)") +
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()


###### 18s 
wt_sn3_18s <-subset(wt_sn3, contig=="18s")
    pdf(file="wt_sn3_18s_delta.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn3_18s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_text_repel(data=subset(wt_sn3_18s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 0.4,segment.color = "grey50")+
          scale_x_continuous(limits=c(50, 1750), breaks=c(50, 400, 800, 1200, 1600,1750))+
          ggtitle("18s rRNA snR3-KO")+
          xlab("Positions")+
          ylab("Delta Current Intensity (WT-KO)") +
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()


wt_sn34_18s <-subset(wt_sn34, contig=="18s")
    pdf(file="wt_sn34_18s_delta.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn34_18s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_text_repel(data=subset(wt_sn34_18s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 0.4,segment.color = "grey50")+
          scale_x_continuous(limits=c(50, 1750), breaks=c(50, 400, 800, 1200, 1600,1750))+
          ggtitle("18s rRNA snR34-KO")+
          xlab("Positions")+
          ylab("Delta Current Intensity (WT-KO)") +
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()



wt_sn36_18s<-subset(wt_sn36, contig=="18s")
    pdf(file="wt_sn36_18s_delta.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn36_18s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886")+
          geom_vline(xintercept = 1187, color = "#f64b3c", alpha=0.4, size=0.6)+
          geom_text_repel(data=subset(wt_sn36_18s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 1,segment.color = "black")+
          scale_x_continuous(limits=c(50, 1750), breaks=c(50, 400, 800, 1200, 1600,1750))+
          ggtitle("18s rRNA snR36-KO")+
          xlab("Positions")+
          ylab("Delta Current Intensity (WT-KO)") +
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()


wt_sn36_18s<-subset(wt_sn36, contig=="18s")
    pdf(file="wt_sn36_18s_delta.zoomed.pdf",height=5,width=20,onefile=FALSE)
      print(ggplot(wt_sn36_18s, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886")+
          geom_vline(xintercept = 1187, color = "#f64b3c", alpha=0.4, size=5)+
          geom_text_repel(data=subset(wt_sn36_18s, diff > 2), aes(position, diff, label=position,size=3, color="red"), segment.size  = 1,segment.color = "black")+
          scale_x_continuous(limits=c(1100, 1200), breaks=c(1100,1155,1200))+
          ggtitle("18s rRNA snR36-KO")+
          xlab("Positions")+
          ylab("Delta Current Intensity (WT-KO)") +
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()
      ```

### Scatter plots for the single position

```R
#We just need to ignore the 15mer sequence information because it is shifted 2 nucleotides
#Its because Nanopolish event means are based on 5 mers (base in the middle) but the position is reported for the starting base (position1)
library(stringr)
library(ggplot2)
library(reshape2)
library(data.table)


sn3<- read.delim("bc1_perpos_sliding15.tsv")
sn3_2<- sn3[!duplicated(sn3[,1:2]),]
ref <- str_split_fixed(sn3_2$windown, ":" , 15)
ref_2<- as.numeric(ref[,11])
sn3_2$windown<- paste(ref_2)
sn3_2$base<- substring(sn3_2$kmer, 10,10)
event_level_mean<- str_split_fixed(sn3_2$mean_current, ":" ,15)
sn3_3<- cbind(sn3_2, event_level_mean)
sn3_3$mean_current<-NULL
sn3_3$read<- NULL
sn3_3[,-c(1,2,3,4)] <- data.frame(sapply(sn3_3[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
colnames(sn3_3)<- c("ref", "windown", "kmer", "base", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1","2","3","4","5","6" , "7" )

sn3_final<- melt(sn3_3)
sn3_final$reference<-paste(sn3_final$ref, sn3_final$windown)
sn3_final$Strain<- rep("snR3-KO", nrow(sn3_final))


sn34<- read.delim("bc2_perpos_sliding15.tsv")
sn34_2<- sn34[!duplicated(sn34[,1:2]),]
ref <- str_split_fixed(sn34_2$windown, ":" , 15)
ref_2<- as.numeric(ref[,11])
sn34_2$windown<- paste(ref_2)
sn34_2$base<- substring(sn34_2$kmer, 10,10)
event_level_mean<- str_split_fixed(sn34_2$mean_current, ":" ,15)
sn34_3<- cbind(sn34_2, event_level_mean)
sn34_3$mean_current<-NULL
sn34_3$read<- NULL
sn34_3[,-c(1,2,3,4)] <- data.frame(sapply(sn34_3[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
colnames(sn34_3)<- c("ref", "windown", "kmer", "base", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1","2","3","4","5","6" , "7" )
sn34_final<- melt(sn34_3)
sn34_final$reference<-paste(sn34_final$ref, sn34_final$windown)
sn34_final$Strain<- rep("snR34-KO", nrow(sn34_final))



sn36<- read.delim("bc3_perpos_sliding15.tsv")
sn36_2<- sn36[!duplicated(sn36[,1:2]),]
ref <- str_split_fixed(sn36_2$windown, ":" , 15)
ref_2<- as.numeric(ref[,11])
sn36_2$windown<- paste(ref_2)
sn36_2$base<- substring(sn36_2$kmer, 10,10)
event_level_mean<- str_split_fixed(sn36_2$mean_current, ":" ,15)
sn36_3<- cbind(sn36_2, event_level_mean)
sn36_3$mean_current<-NULL
sn36_3$read<- NULL
sn36_3[,-c(1,2,3,4)] <- data.frame(sapply(sn36_3[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
colnames(sn36_3)<- c("ref", "windown", "kmer", "base", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1","2","3","4","5","6" , "7" )
sn36_final<- melt(sn36_3)
sn36_final$reference<-paste(sn36_final$ref, sn36_final$windown)
sn36_final$Strain<- rep("snR36-KO", nrow(sn36_final))



wt<- read.delim("bc4_perpos_sliding15.tsv")
wt_2<- wt[!duplicated(wt[,1:2]),]
ref <- str_split_fixed(wt_2$windown, ":" , 15)
ref_2<- as.numeric(ref[,11])
wt_2$windown<- paste(ref_2)
wt_2$base<- substring(wt_2$kmer, 10,10)
event_level_mean<- str_split_fixed(wt_2$mean_current, ":" ,15)
wt_3<- cbind(wt_2, event_level_mean)
wt_3$mean_current<-NULL
wt_3$read<- NULL
wt_3[,-c(1,2,3,4)] <- data.frame(sapply(wt_3[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
colnames(wt_3)<- c("ref", "windown", "kmer", "base", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1","2","3","4","5","6" , "7" )
wt_final<- melt(wt_3)
wt_final$reference<-paste(wt_final$ref, wt_final$windown)
wt_final$Strain<- rep("WT", nrow(wt_final))


all_final<- rbind(wt_final, sn3_final, sn34_final, sn36_final)


for (i in seq_along(unique(all_final$reference))) { 
  subs<- subset(all_final, all_final$reference == unique(all_final$reference)[i])
  pdf(file=paste(unique(all_final$reference)[i],"15nt_window.pdf",sep="."),height=10,width=25,onefile=FALSE)
  print(ggplot(subs, aes(x=variable, y=value, group=Strain)) +
    geom_line(aes(color=Strain), size=2)+
    geom_point(aes(color=Strain))+
    ggtitle(paste(unique(all_final$reference)[i]))+
    xlab("Relative position") +
    ylab("Mean Event Level")+
    theme(axis.text.x = element_text(face="bold", color="black",size=40),
      axis.text.y = element_text(face="bold", color="black", size=40),
      plot.title = element_text(color="black", size=40, face="bold.italic",hjust = 0.5),
      axis.title.x = element_text(color="black", size=40, face="bold"),
      axis.title.y = element_text(color="black", size=40, face="bold"),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size=1),
      panel.grid.major = element_line(colour = "black", size=0.01, linetype="dashed"),
      legend.key.size = unit(1.5, "cm"),
      legend.title = element_text(color = "black", size = 30,face="bold"),
      legend.text = element_text(color = "black", size=25)))
    dev.off()
}



for (i in seq_along(unique(all_final$reference))) { 
  subs<- subset(all_final, all_final$reference == unique(all_final$reference)[i])
  png(file=paste(unique(all_final$reference)[i],"15nt_window.png",sep="."),height=800,width=2400)
  print(ggplot(subs, aes(x=variable, y=value, group=Strain)) +
    geom_line(aes(color=Strain), size=2)+
    geom_point(aes(color=Strain))+
    ggtitle(paste(unique(all_final$reference)[i]))+
    xlab("Relative position") +
    ylab("Mean Event Level")+
    theme(axis.text.x = element_text(face="bold", color="black",size=40),
          axis.text.y = element_text(face="bold", color="black", size=40),
      plot.title = element_text(color="black", size=40, face="bold.italic",hjust = 0.5),
      axis.title.x = element_text(color="black", size=40, face="bold"),
      axis.title.y = element_text(color="black", size=40, face="bold"),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size=1),
      panel.grid.major = element_line(colour = "black", size=0.1),
      legend.key.size = unit(1.5, "cm"),
      legend.title = element_text(color = "black", size = 30,face="bold"),
            legend.text = element_text(color = "black", size=25)))
    dev.off()
}
```




