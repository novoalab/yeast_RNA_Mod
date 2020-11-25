#Libraries needed
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
cond1 <- args[1] #1st variable
cond2 <- args[2] #2nd variable



#Import the processed event align files
condition1<- read.csv(cond1, sep="\t")
condition2<- read.csv(cond2, sep="\t")

#Add 3 to the positions because event align results are shifted by 2 nts and its zero based
condition1$position<- condition1$position+3
condition2$position<- condition2$position+3


#Create a column for unique positions
condition1$Pos<- paste(condition1$contig, condition1$position, sep="_")
condition2$Pos<- paste(condition2$contig, condition2$position, sep="_")


#Add a column for the sample information
condition1$sample<- rep("condition1", nrow(condition1))
condition2$sample<- rep("condition2", nrow(condition2))


#Merge all the tables
all<- rbind(condition1,condition2)

#Import the table to use for 5mer extraction
mod<- read.delim("pU_5mer.tsv")

#For only Y positions in ALL strains
for (pos in unique(mod$Pos)) {
	subs<- subset(all, Pos==pos)
	pdf(file=paste(pos, "_Y_pos_two_strains.pdf", sep="."),height=4,width=9,onefile=FALSE)
		print(ggplot(subs, aes(x= event_level_mean, fill=sample)) +
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