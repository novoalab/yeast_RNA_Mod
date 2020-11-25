library(plyr)
library(ggrepel)
library(dplyr)
library(stringr)


 args <- commandArgs(trailingOnly = TRUE)
cond1 <- args[1] #1st variable
cond2 <- args[2] #2nd variable
threshold<- args[3] 

#Importing the tables processed in python (median event value per position)
condition1<- read.delim(cond1,sep=" ")
condition2<- read.delim(cond2,sep=" ")


#Add 3 to the position, because thats how we built the modification tables
condition1$position<- condition1$position+3
condition2$position<- condition2$position+3


#Include a column that contains chr-position informatio
condition1$ref<- paste(condition1$contig, condition1$position)
condition2$ref<- paste(condition2$contig, condition2$position)

#Rename the columns
colnames(condition2)<- c("contig", "position","reference_kmer","read_idx", "condition2_event_level_mean", "ref")
colnames(condition1)<- c("contig", "position","reference_kmer","read_idx","condition1_event_level_mean", "ref")

#Create a column with base information
condition2$base<- substring(condition2$reference_kmer, 3,3)
condition1$base<- substring(condition1$reference_kmer, 3,3)

#Join the files by the chr-position information
condition2_condition1<-join(condition2,condition1, by="ref")

#remove the NAs
condition2_condition1<- na.omit(condition2_condition1)

#calculate delta current intensity
condition2_condition1$diff<-abs(condition2_condition1$condition1_event_level_mean- condition2_condition1$condition2_event_level_mean)


###### 25s 
## condition1 
for (chr in unique(condition2_condition1$contig)) {
condition2_condition1_subs<-subset(condition2_condition1, contig==chr)
    pdf(file=paste(chr,"condition2_condition1_delta.pdf",sep="."),height=5,width=20,onefile=FALSE)
      print(ggplot(condition2_condition1_subs, aes(x=position, y=diff)) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_text_repel(data=subset(condition2_condition1_subs, diff > threshold), aes(position, diff, label=position,size=3, color="red"), segment.size  = 0.4,segment.color = "grey50")+
        ggtitle(paste(chr))+
        xlab("Positions")+
        ylab("Delta Current Intensity (Condition 1 vs 2)") +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold", color="black",size=11),
            axis.text.y = element_text(face="bold", color="black", size=11),
            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size=0.5)))
      dev.off()
}
