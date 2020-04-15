#Libraries needed
library(dplyr)
library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)
cond1 <- args[1] #1st variable
cond2 <- args[2] #2nd variable

#Import the processed event align files
condition1<- read.csv(cond1, sep="\t")
condition2<- read.csv(cond1, sep="\t")

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
all<- rbind(condition2,condition1)

#Import the table to use for 5mer extraction
mod<- read.delim("pU_5mer.tsv")

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

merged<- subset(all_all, sample=="condition2" | sample=="condition1")

for (pos in unique(merged$Posn.) ) {
	subs<- subset(merged, Posn.==pos)
	subs2<-subset(merged, Kmer ==as.character(unique(subs$Kmer)))
	subs3<-subset(merged, Posn.==pos)
	library(ggplot2)
		pdf(file=paste(pos, "condition1_condition2_with_U_Cmer.pdf", sep="."),height=3,width=10,onefile=FALSE)
			print(ggplot(subs3, aes(x= event_level_mean, fill=base,color=sample)) +
 			geom_density(alpha=0.3,adjust = 2)+
  			theme_bw()+
  			ggtitle(paste(pos, "Condition 1 vs 2"))+
			xlab("Current Intensity (pA) ")+
      		ylab("Density") +
  			theme(axis.text=element_text(size=15),strip.text = element_text(size=13),
    		axis.title=element_text(size=17,face="bold"),
    		legend.title = element_text(size = 20),
    		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
    		legend.text = element_text(color = "black", size=15)))
  		dev.off()
  	}

