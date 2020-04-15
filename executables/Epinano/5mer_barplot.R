#Load the libraries needed for this script
library(plyr)
library(ggplot2)
library(stringr)
library(reshape2)
#Rscript 5mer_barplot.R wt.bam.tsv.per.site.var.per_site_var.5mer.csv sn3.bam.tsv.per.site.var.per_site_var.5mer.csv


#Input
args <- commandArgs(trailingOnly = TRUE)
cond1 <- args[1] #1st variable
cond2 <- args[2] #2nd variable

#Importing and manipulating the rna mod 5mer file
	#Read the table for modification positions		
	mod_positions<- read.delim("rrna_mod_5mer.tsv")
	#Rename the columns
	colnames(mod_positions)<- c("rRNA","Ref_Pos","Mod", "X.Kmer")
	#Create a column for unique positions (18s 145)
	mod_positions$chr_pos<- paste(mod_positions$rRNA, mod_positions$Ref_Pos)

#Import the Epinano 5mer window outputs 
	#For cond1
		#Importing the epinano output
		cond1<- read.delim(cond1, sep=",")
		#Creating a vector 
		Ref_Pos<- str_split_fixed(cond1$Window, ":", 5)
		#Adding it as a column
		cond1$Ref_Pos<- Ref_Pos[,3]
		#Adding another column that is chromosome and position
		cond1$chr_pos<-paste(cond1$Ref,cond1$Ref_Pos)
		#Coverage of the mid position
		Cov<- str_split_fixed(cond1$Coverage, ":", 5)
		#Adding the coverage as column
		cond1$Cov<- as.numeric(Cov[,3])
		#Take only the onces that have higher coverage than 30
		cond1_2<-subset(cond1, Cov>30)
		#Join the two tables
		cond1_mod<- join(mod_positions,cond1_2, by="chr_pos") 
		#Keep only the complete cases
		cond1_mod<- cond1_mod[complete.cases(cond1_mod), ]
		#Create a new column with position and modification
		cond1_mod$site<- paste(cond1_mod$chr_pos,cond1_mod$Mod)
		#Remove unncessary columns
		cond1_mod$Ref_Pos<- NULL
		cond1_mod$Cov<- NULL
		cond1_mod$Ref_Pos<- NULL
		cond1_mod$Coverage<- NULL
		#Melt the table into a format we could use
		cond1_mod2<- melt(cond1_mod)
		#Create a column that contains sample information
		cond1_mod2$Sample<- rep("cond1",nrow(cond1_mod2))
		#Change some characters in the feature column
		cond1_mod2$feature<- gsub('.{1}$', '', cond1_mod2$variable)
		cond1_mod2$feature <- gsub("q", "Quality Score",cond1_mod2$feature)
		cond1_mod2$feature <- gsub("del", "Deletion Frequency",cond1_mod2$feature)
		cond1_mod2$feature <- gsub("ins", "Insertion Frequency",cond1_mod2$feature)
		cond1_mod2$feature <- gsub("mis", "Mismatch Frequency",cond1_mod2$feature)
		#Create the column that will be unique and be used for the title of the plots
		cond1_mod2$title <- paste(cond1_mod2$site, cond1_mod2$X.Kmer)

	#For cond2
		#Importing the epinano output
		cond2<- read.delim(cond2, sep=",")
		#Creating a vector 
		Ref_Pos<- str_split_fixed(cond2$Window, ":", 5)
		#Adding it as a column
		cond2$Ref_Pos<- Ref_Pos[,3]
		#Adding another column that is chromosome and position
		cond2$chr_pos<-paste(cond2$Ref,cond2$Ref_Pos)
		#Coverage of the mid position
		Cov<- str_split_fixed(cond2$Coverage, ":", 5)
		#Adding the coverage as column
		cond2$Cov<- as.numeric(Cov[,3])
		#Take only the onces that have higher coverage than 30
		cond2_2<-subset(cond2, Cov>30)
		#Join the two tables
		cond2_mod<- join(mod_positions,cond2_2, by="chr_pos") 
		#Keep only the complete cases
		cond2_mod<- cond2_mod[complete.cases(cond2_mod), ]
		#Create a new column with position and modification
		cond2_mod$site<- paste(cond2_mod$chr_pos,cond2_mod$Mod)
		#Remove unncessary columns
		cond2_mod$Ref_Pos<- NULL
		cond2_mod$Cov<- NULL
		cond2_mod$Ref_Pos<- NULL
		cond2_mod$Coverage<- NULL
		#Melt the table into a format we could use
		cond2_mod2<- melt(cond2_mod)
		#Create a column that contains sample information
		cond2_mod2$Sample<- rep("cond2",nrow(cond2_mod2))
		#Change some characters in the feature column
		cond2_mod2$feature<- gsub('.{1}$', '', cond2_mod2$variable)
		cond2_mod2$feature <- gsub("q", "Quality Score",cond2_mod2$feature)
		cond2_mod2$feature <- gsub("del", "Deletion Frequency",cond2_mod2$feature)
		cond2_mod2$feature <- gsub("ins", "Insertion Frequency",cond2_mod2$feature)
		cond2_mod2$feature <- gsub("mis", "Mismatch Frequency",cond2_mod2$feature)
		#Create the column that will be unique and be used for the title of the plots
		cond2_mod2$title <- paste(cond2_mod2$site, cond2_mod2$X.Kmer)

	



#Plotting barplots using ggplot2
	#For Condition 1 ONLY
	#Create a for loop
	for (mod in unique(cond1_mod2$title)) {
		subset_cond1_mod2<- subset(cond1_mod2, title==mod)
		pdf(file=paste(mod, "barplot_cond1_ONLY.pdf", sep="_"),height=4,width=10,onefile=FALSE)
			print(ggplot(subset_cond1_mod2, aes(x=variable, y=value, fill=Sample)) +
  				geom_bar(position=position_dodge(0.8), stat = "summary", fun.y = "mean", width=0.7) +
  				ggtitle(paste(mod))+
          		theme_bw()+
  				geom_jitter( position = position_dodge(1),
              		color = "black",size=0.02) + 
              	xlab("5-mer positions")+
              	ylab("Features") +
  				facet_wrap(~feature,scales="free")+
         		theme(axis.text.x = element_text(face="bold", color="black",size=11),
                	axis.text.y = element_text(face="bold", color="black", size=11),
              		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
              		axis.title.x = element_text(color="black", size=15, face="bold"),
              		axis.title.y = element_text(color="black", size=15, face="bold"),
              		panel.background = element_blank(),
              		strip.text = element_text(size=15),
              		axis.line = element_line(colour = "black", size=0.5),
              		legend.title = element_text(color = "black", size = 20,face="bold"),
                    legend.text = element_text(color = "black", size=20)))
  		dev.off()
	}


	#For Condition 1 ONLY
	#Create a for loop
	for (mod in unique(cond1_mod2$title)) {
		subset_cond2_mod2<- subset(cond2_mod2, title==mod)
		pdf(file=paste(mod, "barplot_cond2_ONLY.pdf", sep="_"),height=4,width=10,onefile=FALSE)
			print(ggplot(subset_cond2_mod2, aes(x=variable, y=value, fill=Sample)) +
  				geom_bar(position=position_dodge(0.8), stat = "summary", fun.y = "mean", width=0.7) +
  				ggtitle(paste(mod))+
          		theme_bw()+
  				geom_jitter( position = position_dodge(1),
              		color = "black",size=0.02) + 
              	xlab("5-mer positions")+
              	ylab("Features") +
  				facet_wrap(~feature,scales="free")+
         		theme(axis.text.x = element_text(face="bold", color="black",size=11),
                	axis.text.y = element_text(face="bold", color="black", size=11),
              		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
              		axis.title.x = element_text(color="black", size=15, face="bold"),
              		axis.title.y = element_text(color="black", size=15, face="bold"),
              		panel.background = element_blank(),
              		strip.text = element_text(size=15),
              		axis.line = element_line(colour = "black", size=0.5),
              		legend.title = element_text(color = "black", size = 20,face="bold"),
                    legend.text = element_text(color = "black", size=20)))
  		dev.off()
	}




	#For WT and Sn3KO
	#Create a for loop
	cond1_cond2<- rbind(cond1_mod2,cond2_mod2)
	for (mod in unique(cond1_cond2$title)) {
		subset_cond1_cond2 <- subset(cond1_cond2, title==mod)
		pdf(file=paste(mod, "barplot_cond1_cond2.pdf", sep="_"),height=4,width=10,onefile=FALSE)
			print(ggplot(subset_cond1_cond2, aes(x=variable, y=value, fill=Sample)) +
  				geom_bar(position=position_dodge(0.8), stat = "summary", fun.y = "mean", width=0.7) +
  				ggtitle(paste(mod))+
         		theme_bw()+
  				geom_jitter( position = position_dodge(1),
              		color = "black",size=0.02) + 
              	xlab("5-mer positions")+
              	ylab("Features") +
  				facet_wrap(~feature,scales="free")+
        	  	theme(axis.text.x = element_text(face="bold", color="black",size=11),
                	axis.text.y = element_text(face="bold", color="black", size=11),
              		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
              		axis.title.x = element_text(color="black", size=15, face="bold"),
              		axis.title.y = element_text(color="black", size=15, face="bold"),
              		panel.background = element_blank(),
              		strip.text = element_text(size=15),
              		axis.line = element_line(colour = "black", size=0.5),
              		legend.title = element_text(color = "black", size = 20,face="bold"),
                    legend.text = element_text(color = "black", size=20)))
  		dev.off()
	}
