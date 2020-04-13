# Executable post-processing

```R
#This is for scatter plot between WT and KO
#Rscript xy_plot.R wt.bam.tsv.per.site.var.csv ko.bam.tsv.per.site.var.csv 50
#How to use the script
#Load the libraries needed for this script
library(plyr)
library(ggplot2)
library(ggrepel) 

#Importing and manipulating the Epinano outputs
		#For wt (Barcode 4)
		#Read the Epinano per site table
		args <- commandArgs(trailingOnly = TRUE)
		input1 <- args[1] #1st variable
		input2 <- args[2] #1st variable
		input3 <- args[3]

		#Coverage filter
		coverage<- as.numeric(input3)
	
		#For WT 
		wt <- read.delim(input1,sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		wt<-subset(wt, cov>coverage)
		#Add a column replicating wt
		wt$sample <- rep("wt",nrow(wt)) 
		#Create a column for unique positions (18s 145)
		wt$position<- paste(wt$X.Ref,wt$pos)
		#Remove excess columns		
		wt<- wt[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position")]
		#Rename columns
		colnames(wt)<- c("base","wt_q_mean", "wt_q_median", "wt_q_std", "wt_mis", "wt_del", "wt_ins", "position") #Rename the columns



		#For KO
		ko <- read.delim(input1,sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		ko<-subset(ko, cov>coverage)
		#Add a column replicating ko
		ko$sample <- rep("ko",nrow(ko)) 
		#Create a column for unique positions (18s 145)
		ko$position<- paste(ko$X.Ref,ko$pos)
		#Remove excess columns		
		ko<- ko[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position")]
		#Rename columns
		colnames(ko)<- c("base","ko_q_mean", "ko_q_median", "ko_q_std", "ko_mis", "ko_del", "ko_ins", "position") #Rename the columns




#Creating comparison files
	#Use join function to merge WT and KO
	wt_ko<- join(wt, ko, by="position")
	#Remove the columns that contain NA
	wt_ko<- subset(wt_ko, ko_q_mean!="NA")



	# Scatter plots of Mismatch frequencies
		#Calculate the thresholds for the significant positions (to be labeled)
			#Calculate threshold for WT and KO
				#Calculate residuals
					res_mis<- lm(wt_ko$wt_mis ~ wt_ko$ko_mis)
				#Add the residuals as columns
					wt_ko$res_mis<- res_mis$residuals	
				#Calculate difference between a point and threshold		
					wt_ko$diff_mis<- abs(wt_ko$res_mis)-0.1

		#Plotting the scatter plots using ggplot2
			#For WT and SN3KO
				#Create a subset of the significantly different positions
				diff_ko_mis<-subset(wt_ko, diff_mis>0)
				#Export as PDF
				pdf(file="wt_ko_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_ko, aes(x=wt_mis, y=ko_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_ko, diff_mis>0), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_ko, diff_mis>0), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("WT-KO Mismatch Frequency ")+
			  			xlab("Mismatch Frequency (WT)")+
			  			ylab("Mismatch Frequency (KO)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


	# Scatter plots of Insertion frequencies
		#Calculate the thresholds for the significant positions (to be labeled)
			#Calculate threshold for WT and koKO
				#Calculate residuals
					res_ins<- lm(wt_ko$wt_ins ~ wt_ko$ko_ins)
				#Add the residuals as columns
					wt_ko$res_ins<- res_ins$residuals
				#Calculate difference between a point and threshold	
					wt_ko$diff_ins<- abs(wt_ko$res_ins)-0.1

		#Plotting the scatter plots using ggplot2
			#For WT and KO
				#Create a subset of the significantly different positions
				diff_ko_ins<-subset(wt_ko, diff_ins>0)
				#Export as PDF
				pdf(file="wt_ko_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_ko, aes(x=wt_ins, y=ko_ins)) +
			 			geom_point(size=1, color="grey")+
			  			geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_ko, diff_ins>0), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_ko, diff_ins>0), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("WT-KO Insertion Frequency")+
			  			xlab("Insertion Frequency (WT)")+
			  			ylab("Insertion Frequency (KO)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


	# Scatter plots of Deletion frequencies
		#Calculate the thresholds for the significant positions (to be labeled)
			#Calculate threshold for WT and SN3KO
				#Calculate residuals
					res_del<- lm(wt_ko$wt_del ~ wt_ko$ko_del)
				#Add the residuals as columns
					wt_ko$res_del<- res_del$residuals
				#Calculate difference between a point and threshold	
					wt_ko$diff_del<- abs(wt_ko$res_del)-0.1
		#Plotting the scatter plots using ggplot2
			#For WT and KO
				#Create a subset of the significantly different positions
				diff_ko_del<-subset(wt_ko, diff_del>0)
				#Export as PDF
				pdf(file="wt_ko_scatter_deletion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_ko, aes(x=wt_del, y=ko_del)) +
			 			geom_point(size=1, color="grey")+
			  			geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_ko, diff_del>0), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_ko, diff_del>0), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("WT-KO Deletion Frequency")+
			  			xlab("Deletion Frequency (WT)")+
			  			ylab("Deletion Frequency (KO)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()
```





### 5-Mer Barplots for all modifications

```R
#Load the libraries needed for this script
library(plyr)
library(ggplot2)
library(stringr)
library(reshape2)


#Importing and manipulating the rna mod 5mer file
	#Read the table for modification positions		
	mod_positions<- read.delim("rrna_mod_5mer.tsv")
	#Rename the columns
	colnames(mod_positions)<- c("rRNA","Ref_Pos","Mod", "X.Kmer")
	#Create a column for unique positions (18s 145)
	mod_positions$chr_pos<- paste(mod_positions$rRNA, mod_positions$Ref_Pos)

#Import the Epinano 5mer window outputs 
	#For sn3KO
		#Importing the epinano output
		sn3<- read.delim("sn3.bam.tsv.per.site.var.per_site_var.5mer.csv", sep=",")
		#Creating a vector 
		Ref_Pos<- str_split_fixed(sn3$Window, ":", 5)
		#Adding it as a column
		sn3$Ref_Pos<- Ref_Pos[,3]
		#Adding another column that is chromosome and position
		sn3$chr_pos<-paste(sn3$Ref,sn3$Ref_Pos)
		#Coverage of the mid position
		Cov<- str_split_fixed(sn3$Coverage, ":", 5)
		#Adding the coverage as column
		sn3$Cov<- as.numeric(Cov[,3])
		#Take only the onces that have higher coverage than 30
		sn3_2<-subset(sn3, Cov>30)
		#Join the two tables
		sn3_mod<- join(mod_positions,sn3_2, by="chr_pos") 
		#Keep only the complete cases
		sn3_mod<- sn3_mod[complete.cases(sn3_mod), ]
		#Create a new column with position and modification
		sn3_mod$site<- paste(sn3_mod$chr_pos,sn3_mod$Mod)
		#Remove unncessary columns
		sn3_mod$Ref_Pos<- NULL
		sn3_mod$Cov<- NULL
		sn3_mod$Ref_Pos<- NULL
		sn3_mod$Coverage<- NULL
		#Melt the table into a format we could use
		sn3_mod2<- melt(sn3_mod)
		#Create a column that contains sample information
		sn3_mod2$Sample<- rep("sn3",nrow(sn3_mod2))
		#Change some characters in the feature column
		sn3_mod2$feature<- gsub('.{1}$', '', sn3_mod2$variable)
		sn3_mod2$feature <- gsub("q", "Quality Score",sn3_mod2$feature)
		sn3_mod2$feature <- gsub("del", "Deletion Frequency",sn3_mod2$feature)
		sn3_mod2$feature <- gsub("ins", "Insertion Frequency",sn3_mod2$feature)
		sn3_mod2$feature <- gsub("mis", "Mismatch Frequency",sn3_mod2$feature)
		#Create the column that will be unique and be used for the title of the plots
		sn3_mod2$title <- paste(sn3_mod2$site, sn3_mod2$X.Kmer)


	#For sn34KO
		#Importing the epinano output
		sn34<- read.delim("sn34.bam.tsv.per.site.var.per_site_var.5mer.csv", sep=",")
		#Creating a vector 
		Ref_Pos<- str_split_fixed(sn34$Window, ":", 5)
		#Adding it as a column
		sn34$Ref_Pos<- Ref_Pos[,3]
		#Adding another column that is chromosome and position
		sn34$chr_pos<-paste(sn34$Ref,sn34$Ref_Pos)
		#Coverage of the mid position
		Cov<- str_split_fixed(sn34$Coverage, ":", 5)
		#Adding the coverage as column
		sn34$Cov<- as.numeric(Cov[,3])
		#Take only the onces that have higher coverage than 30
		sn34_2<-subset(sn34, Cov>30)
		#Join the two tables
		sn34_mod<- join(mod_positions,sn34_2, by="chr_pos") 
		#Keep only the complete cases
		sn34_mod<- sn34_mod[complete.cases(sn34_mod), ]
		#Create a new column with position and modification
		sn34_mod$site<- paste(sn34_mod$chr_pos,sn34_mod$Mod)
		#Remove unncessary columns
		sn34_mod$Ref_Pos<- NULL
		sn34_mod$Cov<- NULL
		sn34_mod$Ref_Pos<- NULL
		sn34_mod$Coverage<- NULL
		#Melt the table into a format we could use
		sn34_mod2<- melt(sn34_mod)
		#Create a column that contains sample information
		sn34_mod2$Sample<- rep("sn34",nrow(sn34_mod2))
		#Change some characters in the feature column
		sn34_mod2$feature<- gsub('.{1}$', '', sn34_mod2$variable)
		sn34_mod2$feature <- gsub("q", "Quality Score",sn34_mod2$feature)
		sn34_mod2$feature <- gsub("del", "Deletion Frequency",sn34_mod2$feature)
		sn34_mod2$feature <- gsub("ins", "Insertion Frequency",sn34_mod2$feature)
		sn34_mod2$feature <- gsub("mis", "Mismatch Frequency",sn34_mod2$feature)
		#Create the column that will be unique and be used for the title of the plots
		sn34_mod2$title <- paste(sn34_mod2$site, sn34_mod2$X.Kmer)

	#For sn36KO
		#Importing the epinano output
		sn36<- read.delim("sn36.bam.tsv.per.site.var.per_site_var.5mer.csv", sep=",")
		#Creating a vector 
		Ref_Pos<- str_split_fixed(sn36$Window, ":", 5)
		#Adding it as a column
		sn36$Ref_Pos<- Ref_Pos[,3]
		#Adding another column that is chromosome and position
		sn36$chr_pos<-paste(sn36$Ref,sn36$Ref_Pos)
		#Coverage of the mid position
		Cov<- str_split_fixed(sn36$Coverage, ":", 5)
		#Adding the coverage as column
		sn36$Cov<- as.numeric(Cov[,3])
		#Take only the onces that have higher coverage than 30
		sn36_2<-subset(sn36, Cov>30)
		#Join the two tables
		sn36_mod<- join(mod_positions,sn36_2, by="chr_pos") 
		#Keep only the complete cases
		sn36_mod<- sn36_mod[complete.cases(sn36_mod), ]
		#Create a new column with position and modification
		sn36_mod$site<- paste(sn36_mod$chr_pos,sn36_mod$Mod)
		#Remove unncessary columns
		sn36_mod$Ref_Pos<- NULL
		sn36_mod$Cov<- NULL
		sn36_mod$Ref_Pos<- NULL
		sn36_mod$Coverage<- NULL
		#Melt the table into a format we could use
		sn36_mod2<- melt(sn36_mod)
		#Create a column that contains sample information
		sn36_mod2$Sample<- rep("sn36",nrow(sn36_mod2))
		#Change some characters in the feature column
		sn36_mod2$feature<- gsub('.{1}$', '', sn36_mod2$variable)
		sn36_mod2$feature <- gsub("q", "Quality Score",sn36_mod2$feature)
		sn36_mod2$feature <- gsub("del", "Deletion Frequency",sn36_mod2$feature)
		sn36_mod2$feature <- gsub("ins", "Insertion Frequency",sn36_mod2$feature)
		sn36_mod2$feature <- gsub("mis", "Mismatch Frequency",sn36_mod2$feature)
		#Create the column that will be unique and be used for the title of the plots
		sn36_mod2$title <- paste(sn36_mod2$site, sn36_mod2$X.Kmer)

	#For WT
		#Importing the epinano output
		wt<- read.delim("wt.bam.tsv.per.site.var.per_site_var.5mer.csv", sep=",")
		#Creating a vector 
		Ref_Pos<- str_split_fixed(wt$Window, ":", 5)
		#Adding it as a column
		wt$Ref_Pos<- Ref_Pos[,3]
		#Adding another column that is chromosome and position
		wt$chr_pos<-paste(wt$Ref,wt$Ref_Pos)
		#Coverage of the mid position
		Cov<- str_split_fixed(wt$Coverage, ":", 5)
		#Adding the coverage as column
		wt$Cov<- as.numeric(Cov[,3])
		#Take only the onces that have higher coverage than 30
		wt_2<-subset(wt, Cov>30)
		#Join the two tables
		wt_mod<- join(mod_positions,wt_2, by="chr_pos") 
		#Keep only the complete cases
		wt_mod<- wt_mod[complete.cases(wt_mod), ]
		#Create a new column with position and modification
		wt_mod$site<- paste(wt_mod$chr_pos,wt_mod$Mod)
		#Remove unncessary columns
		wt_mod$Ref_Pos<- NULL
		wt_mod$Cov<- NULL
		wt_mod$Ref_Pos<- NULL
		wt_mod$Coverage<- NULL
		#Melt the table into a format we could use
		wt_mod2<- melt(wt_mod)
		#Create a column that contains sample information
		wt_mod2$Sample<- rep("wt",nrow(wt_mod2))
		#Change some characters in the feature column
		wt_mod2$feature<- gsub('.{1}$', '', wt_mod2$variable)
		wt_mod2$feature <- gsub("q", "Quality Score",wt_mod2$feature)
		wt_mod2$feature <- gsub("del", "Deletion Frequency",wt_mod2$feature)
		wt_mod2$feature <- gsub("ins", "Insertion Frequency",wt_mod2$feature)
		wt_mod2$feature <- gsub("mis", "Mismatch Frequency",wt_mod2$feature)
		#Create the column that will be unique and be used for the title of the plots
		wt_mod2$title <- paste(wt_mod2$site, wt_mod2$X.Kmer)




#Plotting barplots using ggplot2
	#For WT ONLY
	#Create a for loop
	for (mod in unique(wt_mod2$title)) {
		subset_wt_mod2<- subset(wt_mod2, title==mod)
		pdf(file=paste(mod, "barplot_wt_ONLY.pdf", sep="_"),height=4,width=10,onefile=FALSE)
			print(ggplot(subset_wt_mod2, aes(x=variable, y=value, fill=Sample)) +
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
	wt_sn3<- rbind(wt_mod2,sn3_mod2)
	for (mod in unique(wt_sn3$title)) {
		subset_wt_sn3 <- subset(wt_sn3, title==mod)
		pdf(file=paste(mod, "barplot_wt_sn3.pdf", sep="_"),height=4,width=10,onefile=FALSE)
			print(ggplot(subset_wt_sn3, aes(x=variable, y=value, fill=Sample)) +
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


	#For WT and Sn34KO
	#Create a for loop
	wt_sn34<- rbind(wt_mod2,sn34_mod2)
	for (mod in unique(wt_sn34$title)) {
		subset_wt_sn34<- subset(wt_sn34, title==mod)
		pdf(file=paste(mod, "barplot_wt_sn34.pdf", sep="_"),height=4,width=10,onefile=FALSE)
			print(ggplot(subset_wt_sn34, aes(x=variable, y=value, fill=Sample)) +
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




	#For WT and Sn34KO
	#Create a for loop
	wt_sn36<- rbind(wt_mod2,sn36_mod2)
	for (mod in unique(wt_sn36$title)) {
		subset_wt_sn36<- subset(wt_sn36, title==mod)
		pdf(file=paste(mod, "barplot_wt_sn36.pdf", sep="_"),height=4,width=10,onefile=FALSE)
			print(ggplot(subset_wt_sn36, aes(x=variable, y=value, fill=Sample)) +
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
```





### Dot-Plots of Modified vs Unmodified bases
  * We are using a table called rrna_mod_status.tsv, which contains ModStatus and Status values for each position. 
  * ModStatus only indicates whether a position is modified or not
  * Status indicates whether a position is AFFECTED by a neighboruing position (i.e. they are within the 5-mer of another modification)
  * So we will plot ALL positions and only Unaffected positions
```R
#Libraries needed
	library(plyr)
	library(stringr)
	library(reshape2)
	library(dplyr)
	library(ggplot2)
	library(ggbeeswarm)

#Plotting UNAFFECTED Positions
	#Import status file
	status<-read.delim("rrna_mod_status.tsv")
	status$chr_pos<- paste(status$Chr, status$Position)
	#Import yeast sequencing data (WT)
	yeast<-read.delim("wt.bam.tsv.per.site.var.per_site_var.5mer.csv", sep=",")
	#Create a vector for positions
	Ref_Pos<- str_split_fixed(yeast$Window, ":", 5) 
	#Create a vector for coverage
	Coverage<- str_split_fixed(yeast$Coverage, ":", 5)
	#Add it to the table
	yeast$CoverageX<- as.numeric(Coverage[,3])
	#Add it to the table
	yeast$Ref_Pos<- Ref_Pos[,3] 
	#Create a new column with two columns pasted
	yeast$chr_pos<-paste(yeast$Ref,yeast$Ref_Pos) 
	#Select the positions that has coverage higher than 30
	yeast<- subset(yeast, CoverageX>50)
	#Remove the CoverageX column
	yeast$CoverageX<- NULL

	#Add the status to the table
	#We need to use the chr_pos to join yeast table and status table
	yeast_mod<- join(yeast,status, by="chr_pos")
	#Order by position
	yeast_mod<-yeast_mod[order(yeast_mod$Position),]
	#Order by Chromosome
	yeast_mod<-yeast_mod[order(yeast_mod$Ref),]
	#Extract the base information
	yeast_mod$Nuc<- substring(yeast_mod$X.Kmer, 3,3)

	#Extract mod positions which are UNAFFECTED by other mods
	mod_only<- subset(yeast_mod, Status != "Unm") 
	#Remove some columns
	mod_only_stats<- mod_only[,c("Nuc","chr_pos", "Status", "q3", "mis3", "del3", "ins3")]
	#Create a new column as identifier
	mod_only_stats$identifier<-paste(mod_only_stats$chr_pos,mod_only_stats$Status)
	#Place it to the beginning
	mod_only_stats <- mod_only_stats %>% select(identifier, everything())
	#Transform all the columns except 1:4 into numeric
    mod_only_stats[,-c(1:4)] <- data.frame(sapply(mod_only_stats[,-c(1:4)] , function(x) as.numeric(as.character(x))))
    #Log transform the mismatch frequency
	mod_only_stats$mis3<-log(mod_only_stats$mis3)
	#Log transform the deletion frequency
	mod_only_stats$del3<-log(mod_only_stats$del3)
	#Log transform the insertion frequency
	mod_only_stats$ins3<-log(mod_only_stats$ins3)
	#Rename the columns
	colnames(mod_only_stats)<- c("identifier", "Nuc", "chr_pos", "Status", "Quality Score", "log (Mismatch Frequency)", "log (Deletion Frequency)", "log (Insertion Frequency)")

	#Extract unmod positions which are UNAFFECTED by other mods
	unmod_stats<-subset(yeast_mod, Status == "Unm")
	#Remove some columns
	unmod_stats2<- unmod_stats[,c("Nuc", "chr_pos", "Status", "q3", "mis3", "del3","ins3")]
	#Create a column as identified
	unmod_stats2$identifier<-paste(unmod_stats$Nuc)
	#Place it to the beginning
	unmod_stats2 <- unmod_stats2 %>% select(identifier, everything())
	#Transform all the columns except 1:4 into numeric
    unmod_stats2[,-c(1:4)] <- data.frame(sapply(unmod_stats2[,-c(1:4)] , function(x) as.numeric(as.character(x))))
	#Log transform the deletion frequency
	unmod_stats2$mis3<-log(unmod_stats2$mis3)
	#Log transform the deletion frequency
	unmod_stats2$del3<-log(unmod_stats2$del3)
	#Log transform the deletion frequency
	unmod_stats2$ins3<-log(unmod_stats2$ins3)
	colnames(unmod_stats2)<- c("identifier", "Nuc", "chr_pos", "Status", "Quality Score", "log (Mismatch Frequency)", "log (Deletion Frequency)", "log (Insertion Frequency)")



	#Density Plots 
	for (mod in unique(mod_only_stats$Status)) {
		subset_mod <- subset(mod_only_stats, Status==mod)
		subset_unmod<-subset(unmod_stats2, Nuc==unique(subset_mod$Nuc))
		binded<- rbind(subset_mod,subset_unmod)
		mbinded<- melt(binded)
		mbinded$Base<- gsub("Unm", paste(unique(mbinded$Nuc)), mbinded$Status)
		mbinded$Base<- gsub("T", "U", mbinded$Base)
		mbinded$Base <- factor(mbinded$Base, levels = unique(mbinded$Base))
		pdf(file=paste(mod, "density.unaffected.pdf", sep="_"),height=3,width=10,onefile=FALSE)
			print(ggplot(mbinded, aes(x= value, fill=Base,color=Base)) +
					geom_density(alpha=0.3,adjust = 2)+
					facet_wrap(~variable,scales = "free", nrow=1)+
					theme_bw())
		dev.off()
	}

	#Dot Plots
	for (mod in unique(mod_only_stats$Status)) {
		subset_mod <- subset(mod_only_stats, Status==mod)
		subset_unmod<-subset(unmod_stats2, Nuc==unique(subset_mod$Nuc))
		binded<- rbind(subset_mod,subset_unmod)
		mbinded<- melt(binded)
		mbinded$Base<- gsub("Unm", paste(unique(mbinded$Nuc)), mbinded$Status)
		mbinded$Base<- gsub("T", "U", mbinded$Base)
		mbinded$Base <- factor(mbinded$Base, levels = unique(mbinded$Base))
		pdf(file=paste(mod, "dotplot.ALL.pdf", sep="_"),height=4,width=12,onefile=FALSE)
			print(ggplot(mbinded, aes(x=Base, y=value)) + 
				geom_quasirandom(varwidth = TRUE, aes(color=Base))+
				scale_color_manual(values=c("#ffa41b", "#79bac1"))+
				geom_boxplot(aes(alpha=0), outlier.size=NA)+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
				facet_wrap(~variable,scales = "free", nrow=1)+
				theme_bw()+
				xlab("Positions")+
              	ylab("Features") +
				theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            		axis.title=element_text(size=17,face="bold"),
            		legend.title = element_text(size = 20),
            		legend.text = element_text(color = "black", size=15)))
		dev.off()
	}




#Plotting UNAFFECTED Positions
	#Import status file
	status<-read.delim("rrna_mod_status.tsv")
	status$chr_pos<- paste(status$Chr, status$Position)
	#Import yeast sequencing data (WT)
	yeast<-read.delim("wt.bam.tsv.per.site.var.per_site_var.5mer.csv", sep=",")
	#Create a vector for positions
	Ref_Pos<- str_split_fixed(yeast$Window, ":", 5) 
	#Create a vector for coverage
	Coverage<- str_split_fixed(yeast$Coverage, ":", 5)
	#Add it to the table
	yeast$CoverageX<- as.numeric(Coverage[,3])
	#Add it to the table
	yeast$Ref_Pos<- Ref_Pos[,3] 
	#Create a new column with two columns pasted
	yeast$chr_pos<-paste(yeast$Ref,yeast$Ref_Pos) 
	#Select the positions that has coverage higher than 30
	yeast<- subset(yeast, CoverageX>50)
	#Remove the CoverageX column
	yeast$CoverageX<- NULL

	#Add the status to the table
	#We need to use the chr_pos to join yeast table and status table
	yeast_mod<- join(yeast,status, by="chr_pos")
	#Order by position
	yeast_mod<-yeast_mod[order(yeast_mod$Position),]
	#Order by Chromosome
	yeast_mod<-yeast_mod[order(yeast_mod$Ref),]
	#Extract the base information
	yeast_mod$Nuc<- substring(yeast_mod$X.Kmer, 3,3)

	#Extract mod positions which are UNAFFECTED by other mods
	mod_only<- subset(yeast_mod, ModStatus != "Unm") 
	#Remove some columns
	mod_only_stats<- mod_only[,c("Nuc","chr_pos", "ModStatus", "q3", "mis3", "del3", "ins3")]
	#Create a new column as identifier
	mod_only_stats$identifier<-paste(mod_only_stats$chr_pos,mod_only_stats$ModStatus)
	#Place it to the beginning
	mod_only_stats <- mod_only_stats %>% select(identifier, everything())
	#Transform all the columns except 1:4 into numeric
    mod_only_stats[,-c(1:4)] <- data.frame(sapply(mod_only_stats[,-c(1:4)] , function(x) as.numeric(as.character(x))))
    #Log transform the mismatch frequency
	mod_only_stats$mis3<-log(mod_only_stats$mis3)
	#Log transform the deletion frequency
	mod_only_stats$del3<-log(mod_only_stats$del3)
	#Log transform the insertion frequency
	mod_only_stats$ins3<-log(mod_only_stats$ins3)
	#Rename the columns
	colnames(mod_only_stats)<- c("identifier", "Nuc", "chr_pos", "ModStatus", "Quality Score", "log (Mismatch Frequency)", "log (Deletion Frequency)", "log (Insertion Frequency)")

	#Extract unmod positions which are UNAFFECTED by other mods
	unmod_stats<-subset(yeast_mod, ModStatus == "Unm")
	#Remove some columns
	unmod_stats2<- unmod_stats[,c("Nuc", "chr_pos", "ModStatus", "q3", "mis3", "del3","ins3")]
	#Create a column as identified
	unmod_stats2$identifier<-paste(unmod_stats$Nuc)
	#Place it to the beginning
	unmod_stats2 <- unmod_stats2 %>% select(identifier, everything())
	#Transform all the columns except 1:4 into numeric
    unmod_stats2[,-c(1:4)] <- data.frame(sapply(unmod_stats2[,-c(1:4)] , function(x) as.numeric(as.character(x))))
	#Log transform the deletion frequency
	unmod_stats2$mis3<-log(unmod_stats2$mis3)
	#Log transform the deletion frequency
	unmod_stats2$del3<-log(unmod_stats2$del3)
	#Log transform the deletion frequency
	unmod_stats2$ins3<-log(unmod_stats2$ins3)
	colnames(unmod_stats2)<- c("identifier", "Nuc", "chr_pos", "ModStatus", "Quality Score", "log (Mismatch Frequency)", "log (Deletion Frequency)", "log (Insertion Frequency)")



	#Density Plots 
	for (mod in unique(mod_only_stats$ModStatus)) {
		subset_mod <- subset(mod_only_stats, ModStatus==mod)
		subset_unmod<-subset(unmod_stats2, Nuc==unique(subset_mod$Nuc))
		binded<- rbind(subset_mod,subset_unmod)
		mbinded<- melt(binded)
		mbinded$Base<- gsub("Unm", paste(unique(mbinded$Nuc)), mbinded$ModStatus)
		mbinded$Base<- gsub("T", "U", mbinded$Base)
		mbinded$Base <- factor(mbinded$Base, levels = unique(mbinded$Base))
		pdf(file=paste(mod, "density.ALL.pdf", sep="_"),height=3,width=10,onefile=FALSE)
			print(ggplot(mbinded, aes(x= value, fill=Base,color=Base)) +
					geom_density(alpha=0.3,adjust = 2)+
					facet_wrap(~variable,scales = "free", nrow=1)+
					theme_bw())
		dev.off()
	}

	#Dot Plots
	for (mod in unique(mod_only_stats$ModStatus)) {
		subset_mod <- subset(mod_only_stats, ModStatus==mod)
		subset_unmod<-subset(unmod_stats2, Nuc==unique(subset_mod$Nuc))
		binded<- rbind(subset_mod,subset_unmod)
		mbinded<- melt(binded)
		mbinded$Base<- gsub("Unm", paste(unique(mbinded$Nuc)), mbinded$ModStatus)
		mbinded$Base<- gsub("T", "U", mbinded$Base)
		mbinded$Base <- factor(mbinded$Base, levels = unique(mbinded$Base))
		pdf(file=paste(mod, "dotplot.ALL.pdf", sep="_"),height=4,width=12,onefile=FALSE)
			print(ggplot(mbinded, aes(x=Base, y=value)) + 
				geom_quasirandom(varwidth = TRUE, aes(color=Base))+
				scale_color_manual(values=c("#ffa41b", "#79bac1"))+
				geom_boxplot(aes(alpha=0), outlier.size=NA)+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
				facet_wrap(~variable,scales = "free", nrow=1)+
				theme_bw()+
				xlab("Positions")+
              	ylab("Features") +
				theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            		axis.title=element_text(size=17,face="bold"),
            		legend.title = element_text(size = 20),
            		legend.text = element_text(color = "black", size=15)))
		dev.off()
	}
```


