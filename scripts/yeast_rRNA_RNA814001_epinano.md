# Yeast Ribosomal RNA 

## Epinano Analysis of RNA814001 Run (Pseudouridylation KO strains)

Here we only explain the analysis of Epinano output of the Graphmap Default mapped reads

For this, we created a table that contains positions and modifications, called rrna_mod_status.tsv

### Scatter plots (WT vs KO)

```R
#Load the libraries needed for this script
library(plyr)
library(ggplot2)
library(ggrepel) 

#Importing and manipulating the rna mod status file
	#Read the table for modification positions		
	status<- read.delim("rrna_mod_status.tsv")
	#Create a column for unique positions (18s 145)
	status$position<- paste(status$Chr,status$Position)

#Importing and manipulating the Epinano outputs
	#For sn3 (Barcode 1)
		#Read the Epinano per site table
		sn3 <- read.delim("sn3.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		sn3<-subset(sn3, cov>50)
		#Add a column replicating sn3
		sn3$sample <- rep("sn3",nrow(sn3)) 
		#Create a column for unique positions (18s 145)
		sn3$position<- paste(sn3$X.Ref,sn3$pos)
		#Use join function to add the modification information to the Epinano outputs
		sn3_st<- join(sn3, status, by="position")
		#Remove redundant columns
		sn3_st<- sn3_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(sn3_st)<- c("base","sn3_q_mean", "sn3_q_median", "sn3_q_std", "sn3_mis", "sn3_del", "sn3_ins", "position", "ModStatus","Status") #Rename the columns

	#For sn34 (Barcode 2)
		#Read the Epinano per site table
		sn34 <- read.delim("sn34.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		sn34<-subset(sn34, cov>50)
		#Add a column replicating sn34
		sn34$sample <- rep("sn34",nrow(sn34)) 
		#Create a column for unique positions (18s 145)
		sn34$position<- paste(sn34$X.Ref,sn34$pos)
		#Use join function to add the modification information to the Epinano outputs
		sn34_st<- join(sn34, status, by="position")
		#Remove redundant columns
		sn34_st<- sn34_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(sn34_st)<- c("base","sn34_q_mean", "sn34_q_median", "sn34_q_std", "sn34_mis", "sn34_del", "sn34_ins", "position", "ModStatus","Status") #Rename the columns

	#For sn36 (Barcode 3)
		#Read the Epinano per site table
		sn36 <- read.delim("sn36.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		sn36<-subset(sn36, cov>50)
		#Add a column replicating sn36
		sn36$sample <- rep("sn36",nrow(sn36)) 
		#Create a column for unique positions (18s 145)
		sn36$position<- paste(sn36$X.Ref,sn36$pos)
		#Use join function to add the modification information to the Epinano outputs
		sn36_st<- join(sn36, status, by="position")
		#Remove redundant columns
		sn36_st<- sn36_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(sn36_st)<- c("base","sn36_q_mean", "sn36_q_median", "sn36_q_std", "sn36_mis", "sn36_del", "sn36_ins", "position", "ModStatus","Status") #Rename the columns

	#For wt (Barcode 4)
		#Read the Epinano per site table
		wt <- read.delim("wt.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		wt<-subset(wt, cov>50)
		#Add a column replicating wt
		wt$sample <- rep("wt",nrow(wt)) 
		#Create a column for unique positions (18s 145)
		wt$position<- paste(wt$X.Ref,wt$pos)
		#Use join function to add the modification information to the Epinano outputs
		wt_st<- join(wt, status, by="position")
		#Remove redundant columns
		wt_st<- wt_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(wt_st)<- c("base","wt_q_mean", "wt_q_median", "wt_q_std", "wt_mis", "wt_del", "wt_ins", "position", "ModStatus","Status") #Rename the columns



#Creating comparison files
	#Use join function to merge WT and SN3KO
	wt_sn3<- join(wt_st, sn3_st, by="position")
	#Remove the columns that contain NA
	wt_sn3<- subset(wt_sn3, sn3_q_mean!="NA")

	#Use join function to merge WT and SN34KO
	wt_sn34<- join(wt_st, sn34_st, by="position") 
	#Remove the columns that contain NA
	wt_sn34<- subset(wt_sn34, sn34_q_mean!="NA")


	#Use join function to merge WT and SN36KO
	wt_sn36<- join(wt_st, sn36_st, by="position")
	#Remove the columns that contain NA
	wt_sn36<- subset(wt_sn36, sn36_q_mean!="NA")



	# Scatter plots of Mismatch frequencies
		#Calculate the thresholds for the significant positions (to be labeled)

					wt_sn3$diff_mis<- abs(wt_sn3$wt_mis - wt_sn3$sn3_mis)
					wt_sn34$diff_mis<- abs(wt_sn34$wt_mis - wt_sn34$sn34_mis)
					wt_sn36$diff_mis<- abs(wt_sn36$wt_mis - wt_sn36$sn36_mis)

					wt_sn3$diff_ins<- abs(wt_sn3$wt_ins - wt_sn3$sn3_ins)
					wt_sn34$diff_ins<- abs(wt_sn34$wt_ins - wt_sn34$sn34_ins)
					wt_sn36$diff_ins<- abs(wt_sn36$wt_ins - wt_sn36$sn36_ins)

					wt_sn3$diff_del<- abs(wt_sn3$wt_del - wt_sn3$sn3_del)
					wt_sn34$diff_del<- abs(wt_sn34$wt_del - wt_sn34$sn34_del)
					wt_sn36$diff_del<- abs(wt_sn36$wt_del - wt_sn36$sn36_del)

		#Plotting the scatter plots using ggplot2
			#For WT and SN3KO
				#Export as PDF
				pdf(file="wt_sn3_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn3, aes(x=wt_mis, y=sn3_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn3, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn3, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR3-KO")+
			  			xlim(0,1)+
			  			ylim(0,1)+
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

			#For WT and SN34KO
				#Export as PDF
				pdf(file="wt_sn34_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn34, aes(x=wt_mis, y=sn34_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn34, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn34, diff_mis>0.1 & ModStatus == "Y"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR34-KO")+
			  			xlab("Mismatch Frequency (WT)")+
			  			ylab("Mismatch Frequency (KO)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
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

	 	#For WT and SN36KO
				#Export as PDF
				pdf(file="wt_sn36_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn36, aes(x=wt_mis, y=sn36_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn36, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn36, diff_mis>0.1 & ModStatus == "Y"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR36-KO")+
			  			xlab("Mismatch Frequency (WT)")+
			  			ylab("Mismatch Frequency (KO)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
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



		#Plotting the scatter plots using ggplot2
			#For WT and SN3KO
				#Export as PDF
				pdf(file="wt_sn3_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn3, aes(x=wt_ins, y=sn3_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn3, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn3, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR3-KO")+
			  			xlab("Insertion Frequency (WT)")+
			  			ylab("Insertion Frequency (KO)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
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

			#For WT and SN34KO
				#Export as PDF
				pdf(file="wt_sn34_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn34, aes(x=wt_ins, y=sn34_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn34, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn34, diff_ins>0.1 & ModStatus == "Y"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR34-KO")+
			  			xlab("Insertion Frequency (WT)")+
			  			ylab("Insertion Frequency (KO)")+
			  			xlim(0,1)+
			  			ylim(0,1)+
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

	 	#For WT and SN36KO
				#Export as PDF
				pdf(file="wt_sn36_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn36, aes(x=wt_ins, y=sn36_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn36, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn36, diff_ins>0.1 & ModStatus == "Y"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR36-KO")+
			  			xlab("Insertion Frequency (WT)")+
			  			ylab("Insertion Frequency (KO)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
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




		#Plotting the scatter plots using ggplot2
			#For WT and SN3KO
				#Export as PDF
				pdf(file="wt_sn3_scatter_deletion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn3, aes(x=wt_del, y=sn3_del)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn3, diff_del>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn3, diff_del>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR3-KO")+
			  			xlab("Deletion Frequency (WT)")+
			  			ylab("Deletion Frequency (KO)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
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

			#For WT and SN34KO
				#Export as PDF
				pdf(file="wt_sn34_scatter_deletion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn34, aes(x=wt_del, y=sn34_del)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn34, diff_del>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn34, diff_del>0.1 & ModStatus == "Y"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR34-KO")+
			  			xlab("Deletion Frequency (WT)")+
			  			ylab("Deletion Frequency (KO)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
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

	 	#For WT and SN36KO
				#Export as PDF
				pdf(file="wt_sn36_scatter_deletion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_sn36, aes(x=wt_del, y=sn36_del)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			 			geom_point(data=subset(wt_sn36, diff_del>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_sn36, diff_del>0.1 & ModStatus == "Y"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("snR36-KO")+
			  			xlab("Deletion Frequency (WT)")+
			  			ylab("Deletion Frequency (KO)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
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
		pdf(file=paste(mod, "dotplot.unaffected.pdf", sep="_"),height=4,width=12,onefile=FALSE)
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




#Plotting ALL Positions
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


