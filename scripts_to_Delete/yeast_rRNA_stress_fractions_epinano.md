# Yeast Ribosomal RNA 

## Epinano Analysis of  Run (Fractions in Normal and Stress Condition)

Here we only explain the analysis of Epinano output of the Graphmap Default mapped reads

For this, we created a table that contains positions and modifications, called rrna_mod_status.tsv



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
	#For wtinputrep1 (Rep 1)
		#Read the Epinano per site table
		wtinputrep1 <- read.delim("wtinput_rep1.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		wtinputrep1<-subset(wtinputrep1, cov>30)
		#Add a column replicating wtinputrep1
		wtinputrep1$sample <- rep("wtinputrep1",nrow(wtinputrep1)) 
		#Create a column for unique positions (18s 145)
		wtinputrep1$position<- paste(wtinputrep1$X.Ref,wtinputrep1$pos)
		#Use join function to add the modification information to the Epinano outputs
		wtinputrep1_st<- join(wtinputrep1, status, by="position")
		#Remove redundant columns
		wtinputrep1_st<- wtinputrep1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(wtinputrep1_st)<- c("base","wtinputrep1_q_mean", "wtinputrep1_q_median", "wtinputrep1_q_std", "wtinputrep1_mis", "wtinputrep1_del", "wtinputrep1_ins", "position", "ModStatus","Status") #Rename the columns

	#For h2o2inputrep1 (Rep 2)
		#Read the Epinano per site table
		h2o2inputrep1 <- read.delim("h2o2input_rep1.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		h2o2inputrep1<-subset(h2o2inputrep1, cov>30)
		#Add a column replicating h2o2inputrep1
		h2o2inputrep1$sample <- rep("h2o2inputrep1",nrow(h2o2inputrep1)) 
		#Create a column for unique positions (18s 145)
		h2o2inputrep1$position<- paste(h2o2inputrep1$X.Ref,h2o2inputrep1$pos)
		#Use join function to add the modification information to the Epinano outputs
		h2o2inputrep1_st<- join(h2o2inputrep1, status, by="position")
		#Remove redundant columns
		h2o2inputrep1_st<- h2o2inputrep1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(h2o2inputrep1_st)<- c("base","h2o2inputrep1_q_mean", "h2o2inputrep1_q_median", "h2o2inputrep1_q_std", "h2o2inputrep1_mis", "h2o2inputrep1_del", "h2o2inputrep1_ins", "position", "ModStatus","Status") #Rename the columns

	#For wt_f34rep1 (Rep 3)
		#Read the Epinano per site table
		wt_f34rep1 <- read.delim("wt_f34_rep1.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		wt_f34rep1<-subset(wt_f34rep1, cov>30)
		#Add a column replicating wt_f34rep1
		wt_f34rep1$sample <- rep("wt_f34rep1",nrow(wt_f34rep1)) 
		#Create a column for unique positions (18s 145)
		wt_f34rep1$position<- paste(wt_f34rep1$X.Ref,wt_f34rep1$pos)
		#Use join function to add the modification information to the Epinano outputs
		wt_f34rep1_st<- join(wt_f34rep1, status, by="position")
		#Remove redundant columns
		wt_f34rep1_st<- wt_f34rep1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(wt_f34rep1_st)<- c("base","wt_f34rep1_q_mean", "wt_f34rep1_q_median", "wt_f34rep1_q_std", "wt_f34rep1_mis", "wt_f34rep1_del", "wt_f34rep1_ins", "position", "ModStatus","Status") #Rename the columns

	#For h2o2_f34rep1 (Rep 4)
		#Read the Epinano per site table
		h2o2_f34rep1 <- read.delim("h2o2_f34_rep1.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		h2o2_f34rep1<-subset(h2o2_f34rep1, cov>30)
		#Add a column replicating h2o2_f34rep1
		h2o2_f34rep1$sample <- rep("h2o2_f34rep1",nrow(h2o2_f34rep1)) 
		#Create a column for unique positions (18s 145)
		h2o2_f34rep1$position<- paste(h2o2_f34rep1$X.Ref,h2o2_f34rep1$pos)
		#Use join function to add the modification information to the Epinano outputs
		h2o2_f34rep1_st<- join(h2o2_f34rep1, status, by="position")
		#Remove redundant columns
		h2o2_f34rep1_st<- h2o2_f34rep1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(h2o2_f34rep1_st)<- c("base","h2o2_f34rep1_q_mean", "h2o2_f34rep1_q_median", "h2o2_f34rep1_q_std", "h2o2_f34rep1_mis", "h2o2_f34rep1_del", "h2o2_f34rep1_ins", "position", "ModStatus","Status") #Rename the columns




#Importing and manipulating the Epinano outputs
	#For wtinputrep2 (Rep 2)
		#Read the Epinano per site table
		wtinputrep2 <- read.delim("wtinput_rep2.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		wtinputrep2<-subset(wtinputrep2, cov>30)
		#Add a column replicating wtinputrep2
		wtinputrep2$sample <- rep("wtinputrep2",nrow(wtinputrep2)) 
		#Create a column for unique positions (18s 145)
		wtinputrep2$position<- paste(wtinputrep2$X.Ref,wtinputrep2$pos)
		#Use join function to add the modification information to the Epinano outputs
		wtinputrep2_st<- join(wtinputrep2, status, by="position")
		#Remove redundant columns
		wtinputrep2_st<- wtinputrep2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(wtinputrep2_st)<- c("base","wtinputrep2_q_mean", "wtinputrep2_q_median", "wtinputrep2_q_std", "wtinputrep2_mis", "wtinputrep2_del", "wtinputrep2_ins", "position", "ModStatus","Status") #Rename the columns

	#For h2o2inputrep2 (Rep 2)
		#Read the Epinano per site table
		h2o2inputrep2 <- read.delim("h2o2input_rep2.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		h2o2inputrep2<-subset(h2o2inputrep2, cov>30)
		#Add a column replicating h2o2inputrep2
		h2o2inputrep2$sample <- rep("h2o2inputrep2",nrow(h2o2inputrep2)) 
		#Create a column for unique positions (18s 145)
		h2o2inputrep2$position<- paste(h2o2inputrep2$X.Ref,h2o2inputrep2$pos)
		#Use join function to add the modification information to the Epinano outputs
		h2o2inputrep2_st<- join(h2o2inputrep2, status, by="position")
		#Remove redundant columns
		h2o2inputrep2_st<- h2o2inputrep2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(h2o2inputrep2_st)<- c("base","h2o2inputrep2_q_mean", "h2o2inputrep2_q_median", "h2o2inputrep2_q_std", "h2o2inputrep2_mis", "h2o2inputrep2_del", "h2o2inputrep2_ins", "position", "ModStatus","Status") #Rename the columns

	#For wt_f34rep2 (Rep 3)
		#Read the Epinano per site table
		wt_f34rep2 <- read.delim("wt_f34_rep2.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		wt_f34rep2<-subset(wt_f34rep2, cov>30)
		#Add a column replicating wt_f34rep2
		wt_f34rep2$sample <- rep("wt_f34rep2",nrow(wt_f34rep2)) 
		#Create a column for unique positions (18s 145)
		wt_f34rep2$position<- paste(wt_f34rep2$X.Ref,wt_f34rep2$pos)
		#Use join function to add the modification information to the Epinano outputs
		wt_f34rep2_st<- join(wt_f34rep2, status, by="position")
		#Remove redundant columns
		wt_f34rep2_st<- wt_f34rep2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(wt_f34rep2_st)<- c("base","wt_f34rep2_q_mean", "wt_f34rep2_q_median", "wt_f34rep2_q_std", "wt_f34rep2_mis", "wt_f34rep2_del", "wt_f34rep2_ins", "position", "ModStatus","Status") #Rename the columns

	#For h2o2_f34rep2 (Rep 4)
		#Read the Epinano per site table
		h2o2_f34rep2 <- read.delim("h2o2_f34_rep2.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		h2o2_f34rep2<-subset(h2o2_f34rep2, cov>30)
		#Add a column replicating h2o2_f34rep2
		h2o2_f34rep2$sample <- rep("h2o2_f34rep2",nrow(h2o2_f34rep2)) 
		#Create a column for unique positions (18s 145)
		h2o2_f34rep2$position<- paste(h2o2_f34rep2$X.Ref,h2o2_f34rep2$pos)
		#Use join function to add the modification information to the Epinano outputs
		h2o2_f34rep2_st<- join(h2o2_f34rep2, status, by="position")
		#Remove redundant columns
		h2o2_f34rep2_st<- h2o2_f34rep2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(h2o2_f34rep2_st)<- c("base","h2o2_f34rep2_q_mean", "h2o2_f34rep2_q_median", "h2o2_f34rep2_q_std", "h2o2_f34rep2_mis", "h2o2_f34rep2_del", "h2o2_f34rep2_ins", "position", "ModStatus","Status")  

#Creating comparison files
	#Rep1 Only
	#Use join function to merge wtinputrep1 and h2o2_f34rep1
	wtinputrep1_h2o2_f34rep1<- join(wtinputrep1_st, h2o2_f34rep1_st, by="position")
	#Remove the columns that contain NA
	wtinputrep1_h2o2_f34rep1<- subset(wtinputrep1_h2o2_f34rep1, h2o2_f34rep1_q_mean!="NA")

	#Use join function to merge h2o2inputrep1 and h2o2_f34rep1
	h2o2inputrep1_h2o2_f34rep1<- join(h2o2inputrep1_st, h2o2_f34rep1_st, by="position")
	#Remove the columns that contain NA
	h2o2inputrep1_h2o2_f34rep1<- subset(h2o2inputrep1_h2o2_f34rep1, h2o2_f34rep1_q_mean!="NA")

	#Use join function to merge wt_f34rep1 and h2o2_f34rep1
	wt_f34rep1_h2o2_f34rep1<- join(wt_f34rep1_st, h2o2_f34rep1_st, by="position")
	#Remove the columns that contain NA
	wt_f34rep1_h2o2_f34rep1<- subset(wt_f34rep1_h2o2_f34rep1, h2o2_f34rep1_q_mean!="NA")

	#Use join function to merge h2o2inputrep1 and h2o2_f34rep1
	h2o2inputrep1_wt_f34rep1<- join(h2o2inputrep1_st, wt_f34rep1_st, by="position")
	#Remove the columns that contain NA
	h2o2inputrep1_wt_f34rep1<- subset(h2o2inputrep1_wt_f34rep1, wt_f34rep1_q_mean!="NA")
	
	#Use join function to merge wtinputrep1 and h2o2_f34rep1
	wtinputrep1_wt_f34rep1<- join(wtinputrep1_st, wt_f34rep1_st, by="position")
	#Remove the columns that contain NA
	wtinputrep1_wt_f34rep1<- subset(wtinputrep1_wt_f34rep1, wt_f34rep1_q_mean!="NA")
	
	#Use join function to merge wtinputrep1 and h2o2_f34rep1
	wtinputrep1_h2o2inputrep1<- join(wtinputrep1_st, h2o2inputrep1_st, by="position")
	#Remove the columns that contain NA
	wtinputrep1_h2o2inputrep1<- subset(wtinputrep1_h2o2inputrep1, h2o2inputrep1_q_mean!="NA")

	#Rep2 Only
	#Use join function to merge wtinputrep2 and h2o2_f34rep2
	wtinputrep2_h2o2_f34rep2<- join(wtinputrep2_st, h2o2_f34rep2_st, by="position")
	#Remove the columns that contain NA
	wtinputrep2_h2o2_f34rep2<- subset(wtinputrep2_h2o2_f34rep2, h2o2_f34rep2_q_mean!="NA")

	#Use join function to merge h2o2inputrep2 and h2o2_f34rep2
	h2o2inputrep2_h2o2_f34rep2<- join(h2o2inputrep2_st, h2o2_f34rep2_st, by="position")
	#Remove the columns that contain NA
	h2o2inputrep2_h2o2_f34rep2<- subset(h2o2inputrep2_h2o2_f34rep2, h2o2_f34rep2_q_mean!="NA")

	#Use join function to merge wt_f34rep2 and h2o2_f34rep2
	wt_f34rep2_h2o2_f34rep2<- join(wt_f34rep2_st, h2o2_f34rep2_st, by="position")
	#Remove the columns that contain NA
	wt_f34rep2_h2o2_f34rep2<- subset(wt_f34rep2_h2o2_f34rep2, h2o2_f34rep2_q_mean!="NA")

	#Use join function to merge h2o2inputrep2 and h2o2_f34rep2
	h2o2inputrep2_wt_f34rep2<- join(h2o2inputrep2_st, wt_f34rep2_st, by="position")
	#Remove the columns that contain NA
	h2o2inputrep2_wt_f34rep2<- subset(h2o2inputrep2_wt_f34rep2, wt_f34rep2_q_mean!="NA")
	
	#Use join function to merge wtinputrep2 and h2o2_f34rep2
	wtinputrep2_wt_f34rep2<- join(wtinputrep2_st, wt_f34rep2_st, by="position")
	#Remove the columns that contain NA
	wtinputrep2_wt_f34rep2<- subset(wtinputrep2_wt_f34rep2, wt_f34rep2_q_mean!="NA")
	
	#Use join function to merge wtinputrep2 and h2o2_f34rep2
	wtinputrep2_h2o2inputrep2<- join(wtinputrep2_st, h2o2inputrep2_st, by="position")
	#Remove the columns that contain NA
	wtinputrep2_h2o2inputrep2<- subset(wtinputrep2_h2o2inputrep2, h2o2inputrep2_q_mean!="NA")


	#Rep1 VS Rep2
	#Use join function to merge wtinputrep1 and wtinputrep2
	wtinputrep1_wtinputrep2<- join(wtinputrep1_st, wtinputrep2_st, by="position")
	#Remove the columns that contain NA
	wtinputrep1_wtinputrep2<- subset(wtinputrep1_wtinputrep2, wtinputrep2_q_mean!="NA")

	#Use join function to merge h2o2inputrep1 and h2o2inputrep2
	h2o2inputrep1_h2o2inputrep2<- join(h2o2inputrep1_st, h2o2inputrep2_st, by="position")
	#Remove the columns that contain NA
	h2o2inputrep1_h2o2inputrep2<- subset(h2o2inputrep1_h2o2inputrep2, h2o2inputrep2_q_mean!="NA")

	#Use join function to merge wt_f34rep1 and wt_f34rep2
	wt_f34rep1_wt_f34rep2<- join(wt_f34rep1_st, wt_f34rep2_st, by="position")
	#Remove the columns that contain NA
	wt_f34rep1_wt_f34rep2<- subset(wt_f34rep1_wt_f34rep2, wt_f34rep2_q_mean!="NA")

	#Use join function to merge h2o2_f34rep1 and h2o2_f34rep2
	h2o2_f34rep1_h2o2_f34rep2<- join(h2o2_f34rep1_st, h2o2_f34rep2_st, by="position")
	#Remove the columns that contain NA
	h2o2_f34rep1_h2o2_f34rep2<- subset(h2o2_f34rep1_h2o2_f34rep2, h2o2_f34rep2_q_mean!="NA")




	#Difference
	wtinputrep1_h2o2_f34rep1$diff_mis <-  abs(wtinputrep1_h2o2_f34rep1$wtinputrep1_mis-wtinputrep1_h2o2_f34rep1$h2o2_f34rep1_mis)
	wtinputrep1_wt_f34rep1$diff_mis <-  abs(wtinputrep1_wt_f34rep1$wtinputrep1_mis-wtinputrep1_wt_f34rep1$wt_f34rep1_mis)
	wtinputrep1_h2o2inputrep1$diff_mis <-  abs(wtinputrep1_h2o2inputrep1$wtinputrep1_mis-wtinputrep1_h2o2inputrep1$h2o2inputrep1_mis)
	h2o2inputrep1_h2o2_f34rep1$diff_mis <-  abs(h2o2inputrep1_h2o2_f34rep1$h2o2inputrep1_mis-h2o2inputrep1_h2o2_f34rep1$h2o2_f34rep1_mis)
	wt_f34rep1_h2o2_f34rep1$diff_mis <-  abs(wt_f34rep1_h2o2_f34rep1$wt_f34rep1_mis-wt_f34rep1_h2o2_f34rep1$h2o2_f34rep1_mis)
	h2o2inputrep1_wt_f34rep1$diff_mis <-  abs(h2o2inputrep1_wt_f34rep1$h2o2inputrep1_mis-h2o2inputrep1_wt_f34rep1$wt_f34rep1_mis)


	wtinputrep2_h2o2_f34rep2$diff_mis <-  abs(wtinputrep2_h2o2_f34rep2$wtinputrep2_mis-wtinputrep2_h2o2_f34rep2$h2o2_f34rep2_mis)
	wtinputrep2_wt_f34rep2$diff_mis <-  abs(wtinputrep2_wt_f34rep2$wtinputrep2_mis-wtinputrep2_wt_f34rep2$wt_f34rep2_mis)
	wtinputrep2_h2o2inputrep2$diff_mis <-  abs(wtinputrep2_h2o2inputrep2$wtinputrep2_mis-wtinputrep2_h2o2inputrep2$h2o2inputrep2_mis)
	h2o2inputrep2_h2o2_f34rep2$diff_mis <-  abs(h2o2inputrep2_h2o2_f34rep2$h2o2inputrep2_mis-h2o2inputrep2_h2o2_f34rep2$h2o2_f34rep2_mis)
	wt_f34rep2_h2o2_f34rep2$diff_mis <-  abs(wt_f34rep2_h2o2_f34rep2$wt_f34rep2_mis-wt_f34rep2_h2o2_f34rep2$h2o2_f34rep2_mis)
	h2o2inputrep2_wt_f34rep2$diff_mis <-  abs(h2o2inputrep2_wt_f34rep2$h2o2inputrep2_mis-h2o2inputrep2_wt_f34rep2$wt_f34rep2_mis)


	wtinputrep1_wtinputrep2$diff_mis <-  abs(wtinputrep1_wtinputrep2$wtinputrep1_mis-wtinputrep1_wtinputrep2$wtinputrep2_mis)
	h2o2inputrep1_h2o2inputrep2$diff_mis <-  abs(h2o2inputrep1_h2o2inputrep2$h2o2inputrep1_mis-h2o2inputrep1_h2o2inputrep2$h2o2inputrep2_mis)
	wt_f34rep1_wt_f34rep2$diff_mis <-  abs(wt_f34rep1_wt_f34rep2$wt_f34rep1_mis-wt_f34rep1_wt_f34rep2$wt_f34rep2_mis)
	h2o2_f34rep1_h2o2_f34rep2$diff_mis <-  abs(h2o2_f34rep1_h2o2_f34rep2$h2o2_f34rep1_mis-h2o2_f34rep1_h2o2_f34rep2$h2o2_f34rep2_mis)



		#Plotting the scatter plots using ggplot2
			#Rep1
				#For h2o2_f34rep1 and wtinputrep1
				#Export as PDF
				pdf(file="wtinputrep1_h2o2_f34rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep1_h2o2_f34rep1, aes(x=wtinputrep1_mis, y=h2o2_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep1_h2o2_f34rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep1_h2o2_f34rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("WT input vs H2O2 Fractions 3-4")+
			  			xlab("Mismatch Frequency (wtinputrep1)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep1)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()

	
				#For wtinputrep1 and wt_f34rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="wtinputrep1_wt_f34rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep1_wt_f34rep1, aes(x=wtinputrep1_mis, y=wt_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep1_wt_f34rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep1_wt_f34rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs wt Fractions 3-4")+
			  			xlab("Mismatch Frequency (wtinputrep1)")+
			  			ylab("Mismatch Frequency (wt_f34rep1)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep1 and wt_f34rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="h2o2inputrep1_wt_f34rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep1_wt_f34rep1, aes(x=h2o2inputrep1_mis, y=wt_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(h2o2inputrep1_wt_f34rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(h2o2inputrep1_wt_f34rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("H2o2 input vs WT Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2inputrep1)")+
			  			ylab("Mismatch Frequency (wt_f34rep1)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep1 and h2o2inputrep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="wtinputrep1_h2o2inputrep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep1_h2o2inputrep1, aes(x=wtinputrep1_mis, y=h2o2inputrep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep1_h2o2inputrep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep1_h2o2inputrep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs h2O2 Input")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wtinputrep1)")+
			  			ylab("Mismatch Frequency (h2o2inputrep1)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For h2o2inputrep1 and h2o2_f34rep1
				#Export as PDF
				pdf(file="h2o2inputrep1_h2o2_f34rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep1_h2o2_f34rep1, aes(x=h2o2inputrep1_mis, y=h2o2_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(h2o2inputrep1_h2o2_f34rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(h2o2inputrep1_h2o2_f34rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("H2O2 input vs H2O2 Fractions 3-4")+
			  			xlab("Mismatch Frequency (h2o2inputrep1)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep1)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wt_f34rep1 and h2o2_f34rep1
				#Export as PDF
				pdf(file="wt_f34rep1_h2o2_f34rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_f34rep1_h2o2_f34rep1, aes(x=wt_f34rep1_mis, y=h2o2_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wt_f34rep1_h2o2_f34rep1, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_f34rep1_h2o2_f34rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("H2o2 fractions 3-4 vs WT Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wt_f34rep1)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep1)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()







			### REP2 
				#For h2o2_f34rep2 and wtinputrep2
				#Export as PDF
				pdf(file="wtinputrep2_h2o2_f34rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep2_h2o2_f34rep2, aes(x=wtinputrep2_mis, y=h2o2_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep2_h2o2_f34rep2, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep2_h2o2_f34rep2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs h2O2 Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wtinputrep2)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()

	
				#For wtinputrep2 and wt_f34rep2
				#Export as PDF
				pdf(file="wtinputrep2_wt_f34rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep2_wt_f34rep2, aes(x=wtinputrep2_mis, y=wt_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep2_wt_f34rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep2_wt_f34rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("WT input vs WT Fractions 3-4")+
			  			xlab("Mismatch Frequency (wtinputrep2)")+
			  			ylab("Mismatch Frequency (wt_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep2 and wt_f34rep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="h2o2inputrep2_wt_f34rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep2_wt_f34rep2, aes(x=h2o2inputrep2_mis, y=wt_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(h2o2inputrep2_wt_f34rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(h2o2inputrep2_wt_f34rep2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("H2o2 input vs WT Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2inputrep2)")+
			  			ylab("Mismatch Frequency (wt_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep2 and h2o2inputrep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="wtinputrep2_h2o2inputrep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep2_h2o2inputrep2, aes(x=wtinputrep2_mis, y=h2o2inputrep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep2_h2o2inputrep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep2_h2o2inputrep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs h2O2 Input")+
			  			xlab("Mismatch Frequency (wtinputrep2)")+
			  			ylab("Mismatch Frequency (h2o2inputrep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For h2o2inputrep2 and h2o2_f34rep2
				#Export as PDF
				pdf(file="h2o2inputrep2_h2o2_f34rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep2_h2o2_f34rep2, aes(x=h2o2inputrep2_mis, y=h2o2_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(h2o2inputrep2_h2o2_f34rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(h2o2inputrep2_h2o2_f34rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("H2o2 input vs H2O2 Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2inputrep2)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wt_f34rep2 and h2o2_f34rep2
				#Export as PDF
				pdf(file="wt_f34rep2_h2o2_f34rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_f34rep2_h2o2_f34rep2, aes(x=wt_f34rep2_mis, y=h2o2_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wt_f34rep2_h2o2_f34rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_f34rep2_h2o2_f34rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("H2o2 fractions 3-4 vs WT Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wt_f34rep2)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()




			#REP1 VS REP2 

				#For h2o2_f34rep1 and h2o2_f34rep2
				#Export as PDF
				pdf(file="h2o2_f34rep1_h2o2_f34rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2_f34rep1_h2o2_f34rep2, aes(x=h2o2_f34rep1_mis, y=h2o2_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black",se=FALSE)+
			 			geom_point(data=subset(h2o2_f34rep1_h2o2_f34rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(h2o2_f34rep1_h2o2_f34rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("H2o2 fractioN 3-4  vs H2O2 Fraction 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2_f34rep1)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep1 and wtinputrep2
				#Export as PDF
				pdf(file="wtinputrep1_wtinputrep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep1_wtinputrep2, aes(x=wtinputrep1_mis, y=wtinputrep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep1_wtinputrep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep1_wtinputrep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs WT Input")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wtinputrep1)")+
			  			ylab("Mismatch Frequency (wtinputrep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For h2o2inputrep1 and h2o2inputrep2
				#Export as PDF
				pdf(file="h2o2inputrep1_h2o2inputrep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep1_h2o2inputrep2, aes(x=h2o2inputrep1_mis, y=h2o2inputrep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(h2o2inputrep1_h2o2inputrep2, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(h2o2inputrep1_h2o2inputrep2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("H2o2 input vS H2O2 Input")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2inputrep1)")+
			  			ylab("Mismatch Frequency (h2o2inputrep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()



				#For wt_f34rep1 and wt_f34rep2
				#Export as PDF
				pdf(file="wt_f34rep1_wt_f34rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_f34rep1_wt_f34rep2, aes(x=wt_f34rep1_mis, y=wt_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wt_f34rep1_wt_f34rep2, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(wt_f34rep1_wt_f34rep2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt fraction 3-4  vs WT Fraction 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wt_f34rep1)")+
			  			ylab("Mismatch Frequency (wt_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()






		# VENN DIAGRAMS BETWEEN REPLICATES


					#Rep1 Only

					wtinputrep1_h2o2_f34rep1.sig<- subset(wtinputrep1_h2o2_f34rep1, diff_mis>0.1)
					h2o2inputrep1_h2o2_f34rep1.sig<- subset(h2o2inputrep1_h2o2_f34rep1, diff_mis>0.1)
					wt_f34rep1_h2o2_f34rep1.sig<- subset(wt_f34rep1_h2o2_f34rep1, diff_mis>0.1)
					h2o2inputrep1_wt_f34rep1.sig<- subset(h2o2inputrep1_wt_f34rep1, diff_mis>0.1)
					wtinputrep1_wt_f34rep1.sig<- subset(wtinputrep1_wt_f34rep1, diff_mis>0.1)
					wtinputrep1_h2o2inputrep1.sig<- subset(wtinputrep1_h2o2inputrep1, diff_mis>0.1)	
					wtinputrep2_h2o2_f34rep2.sig<- subset(wtinputrep2_h2o2_f34rep2, diff_mis>0.1)
					h2o2inputrep2_h2o2_f34rep2.sig<- subset(h2o2inputrep2_h2o2_f34rep2, diff_mis>0.1)
					wt_f34rep2_h2o2_f34rep2.sig<- subset(wt_f34rep2_h2o2_f34rep2, diff_mis>0.1)
					h2o2inputrep2_wt_f34rep2.sig<- subset(h2o2inputrep2_wt_f34rep2, diff_mis>0.1)
					wtinputrep2_wt_f34rep2.sig<- subset(wtinputrep2_wt_f34rep2, diff_mis>0.1)
					wtinputrep2_h2o2inputrep2.sig<- subset(wtinputrep2_h2o2inputrep2, diff_mis>0.1)



					wtinput_h2o2_f34.common<- intersect(wtinputrep1_h2o2_f34rep1.sig$position,wtinputrep2_h2o2_f34rep2.sig$position)
					h2o2input_h2o2_f34.common<- intersect(h2o2inputrep1_h2o2_f34rep1.sig$position,h2o2inputrep2_h2o2_f34rep2.sig$position)
					wt_f34_h2o2_f34.common<- intersect(wt_f34rep1_h2o2_f34rep1.sig$position,wt_f34rep2_h2o2_f34rep2.sig$position)
					h2o2input_wt_f34.common<- intersect(h2o2inputrep1_wt_f34rep1.sig$position,h2o2inputrep2_wt_f34rep2.sig$position)
					wtinput_wt_f34.common<- intersect(wtinputrep1_wt_f34rep1.sig$position,wtinputrep2_wt_f34rep2.sig$position)
					wtinput_h2o2input.common<- intersect(wtinputrep1_h2o2inputrep1.sig$position,wtinputrep2_h2o2inputrep2.sig$position)


					write.table(wtinput_h2o2_f34.common, file="wtinput_h2o2_f34.common.txt")
					write.table(h2o2input_h2o2_f34.common, file="h2o2input_h2o2_f34.common.txt")
					write.table(wt_f34_h2o2_f34.common, file="wt_f34_h2o2_f34.common.txt")
					write.table(h2o2input_wt_f34.common, file="h2o2input_wt_f34.common.txt")
					write.table(wtinput_wt_f34.common, file="wtinput_wt_f34.common.txt")
					write.table(wtinput_h2o2input.common, file="wtinput_h2o2input.common.txt")



					#Venn diagram for common specific genes combination
					library(VennDiagram)
					pdf("Venn_Diagram_WtiNpUT_H2O2_F34.pdf",height=8,width=8)
					draw.pairwise.venn(length(wtinputrep1_h2o2_f34rep1.sig$position), length(wtinputrep2_h2o2_f34rep2.sig$position), length(wtinput_h2o2_f34.common), category = c("WTinPut vs H2O2 F3-4", "F1vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F2_F4.pdf",height=8,width=8)
					draw.pairwise.venn(length(h2o2inputrep1_h2o2_f34rep1.sig$position), length(h2o2inputrep2_h2o2_f34rep2.sig$position), length(h2o2input_h2o2_f34.common), category = c("H2o2InpUt vs H2O2 F3-4", "F2vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_Wt_F34_H2O2_F34.pdf",height=8,width=8)
					draw.pairwise.venn(length(wt_f34rep1_h2o2_f34rep1.sig$position), length(wt_f34rep2_h2o2_f34rep2.sig$position), length(wt_f34_h2o2_f34.common), category = c("WT f3-4 vs H2O2 F3-4", "F3vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_H2o2iNPUT_WT_F34.pdf",height=8,width=8)
					draw.pairwise.venn(length(h2o2inputrep1_wt_f34rep1.sig$position), length(h2o2inputrep2_wt_f34rep2.sig$position), length(h2o2input_wt_f34.common), category = c("H2o2InpUt vs WT F3-4", "F2vsF3 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_WtiNpUT_WT_F34.pdf",height=8,width=8)
					draw.pairwise.venn(length(wtinputrep1_wt_f34rep1.sig$position), length(wtinputrep2_wt_f34rep2.sig$position), length(wtinput_wt_f34.common), category = c("WTinPut vs WT F3-4", "F1vsF3 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()


					pdf("Venn_Diagram_WtiNpUT_H2O2INPUT.pdf",height=8,width=8)
					draw.pairwise.venn(length(wtinputrep1_h2o2inputrep1.sig$position), length(wtinputrep2_h2o2inputrep2.sig$position), length(wtinput_h2o2input.common), category = c("WT iNpuT vs H2O2 Input", "F1vsF2 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()












					#ONLY Y positions

					wtinputrep1_h2o2_f34rep1.sig<- subset(wtinputrep1_h2o2_f34rep1, diff_mis>0.1 & ModStatus == "Y")
					h2o2inputrep1_h2o2_f34rep1.sig<- subset(h2o2inputrep1_h2o2_f34rep1, diff_mis>0.1 & ModStatus == "Y")
					wt_f34rep1_h2o2_f34rep1.sig<- subset(wt_f34rep1_h2o2_f34rep1, diff_mis>0.1 & ModStatus == "Y")
					h2o2inputrep1_wt_f34rep1.sig<- subset(h2o2inputrep1_wt_f34rep1, diff_mis>0.1 & ModStatus == "Y")
					wtinputrep1_wt_f34rep1.sig<- subset(wtinputrep1_wt_f34rep1, diff_mis>0.1 & ModStatus == "Y")
					wtinputrep1_h2o2inputrep1.sig<- subset(wtinputrep1_h2o2inputrep1, diff_mis>0.1 & ModStatus == "Y")	
					wtinputrep2_h2o2_f34rep2.sig<- subset(wtinputrep2_h2o2_f34rep2, diff_mis>0.1 & ModStatus == "Y")
					h2o2inputrep2_h2o2_f34rep2.sig<- subset(h2o2inputrep2_h2o2_f34rep2, diff_mis>0.1 & ModStatus == "Y")
					wt_f34rep2_h2o2_f34rep2.sig<- subset(wt_f34rep2_h2o2_f34rep2, diff_mis>0.1 & ModStatus == "Y")
					h2o2inputrep2_wt_f34rep2.sig<- subset(h2o2inputrep2_wt_f34rep2, diff_mis>0.1 & ModStatus == "Y")
					wtinputrep2_wt_f34rep2.sig<- subset(wtinputrep2_wt_f34rep2, diff_mis>0.1 & ModStatus == "Y")
					wtinputrep2_h2o2inputrep2.sig<- subset(wtinputrep2_h2o2inputrep2, diff_mis>0.1 & ModStatus == "Y")



					wtinput_h2o2_f34.common<- intersect(wtinputrep1_h2o2_f34rep1.sig$position,wtinputrep2_h2o2_f34rep2.sig$position)
					h2o2input_h2o2_f34.common<- intersect(h2o2inputrep1_h2o2_f34rep1.sig$position,h2o2inputrep2_h2o2_f34rep2.sig$position)
					wt_f34_h2o2_f34.common<- intersect(wt_f34rep1_h2o2_f34rep1.sig$position,wt_f34rep2_h2o2_f34rep2.sig$position)
					h2o2input_wt_f34.common<- intersect(h2o2inputrep1_wt_f34rep1.sig$position,h2o2inputrep2_wt_f34rep2.sig$position)
					wtinput_wt_f34.common<- intersect(wtinputrep1_wt_f34rep1.sig$position,wtinputrep2_wt_f34rep2.sig$position)
					wtinput_h2o2input.common<- intersect(wtinputrep1_h2o2inputrep1.sig$position,wtinputrep2_h2o2inputrep2.sig$position)


					write.table(wtinput_h2o2_f34.common, file="wtinput_h2o2_f34.common_Y.txt")
					write.table(h2o2input_h2o2_f34.common, file="h2o2input_h2o2_f34.common_Y.txt")
					write.table(wt_f34_h2o2_f34.common, file="wt_f34_h2o2_f34.common_Y.txt")
					write.table(h2o2input_wt_f34.common, file="h2o2input_wt_f34.common_Y.txt")
					write.table(wtinput_wt_f34.common, file="wtinput_wt_f34.common_Y.txt")
					write.table(wtinput_h2o2input.common, file="wtinput_h2o2input.common_Y.txt")



					#Venn diagram for common specific genes combination
					library(VennDiagram)
					pdf("Venn_Diagram_WtiNpUT_H2O2_F34_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(wtinputrep1_h2o2_f34rep1.sig$position), length(wtinputrep2_h2o2_f34rep2.sig$position), length(wtinput_h2o2_f34.common), category = c("WTinPut vs H2O2 F3-4", "F1vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F2_F4_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(h2o2inputrep1_h2o2_f34rep1.sig$position), length(h2o2inputrep2_h2o2_f34rep2.sig$position), length(h2o2input_h2o2_f34.common), category = c("H2o2InpUt vs H2O2 F3-4", "F2vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_Wt_F34_H2O2_F34_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(wt_f34rep1_h2o2_f34rep1.sig$position), length(wt_f34rep2_h2o2_f34rep2.sig$position), length(wt_f34_h2o2_f34.common), category = c("WT f3-4 vs H2O2 F3-4", "F3vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_H2o2iNPUT_WT_F34_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(h2o2inputrep1_wt_f34rep1.sig$position), length(h2o2inputrep2_wt_f34rep2.sig$position), length(h2o2input_wt_f34.common), category = c("H2o2InpUt vs WT F3-4", "F2vsF3 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_WtiNpUT_WT_F34_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(wtinputrep1_wt_f34rep1.sig$position), length(wtinputrep2_wt_f34rep2.sig$position), length(wtinput_wt_f34.common), category = c("WTinPut vs WT F3-4", "F1vsF3 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()


					pdf("Venn_Diagram_WtiNpUT_H2O2INPUT_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(wtinputrep1_h2o2inputrep1.sig$position), length(wtinputrep2_h2o2inputrep2.sig$position), length(wtinput_h2o2input.common), category = c("WT iNpuT vs H2O2 Input", "F1vsF2 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()















		### PLOT THE SCATTER WITH ALL SIGNIFICANT POSITIONS (NOT ONLY Y)

		#Plotting the scatter plots using ggplot2
			#Rep1
				#For h2o2_f34rep1 and wtinputrep1
				#Export as PDF
				pdf(file="wtinputrep1_h2o2_f34rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep1_h2o2_f34rep1, aes(x=wtinputrep1_mis, y=h2o2_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep1_h2o2_f34rep1, position == "25s 776"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep1_h2o2_f34rep1, position == "25s 776"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs h2O2 Fractions 3-4")+
			  			xlab("Mismatch Frequency (wtinputrep1)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep1)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()

	
				#For wtinputrep1 and wt_f34rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="wtinputrep1_wt_f34rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep1_wt_f34rep1, aes(x=wtinputrep1_mis, y=wt_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep1_wt_f34rep1, position == "25s 776" | position == "25s 1051" | position == "25s 1335"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep1_wt_f34rep1, position == "25s 776" | position == "25s 1051" | position == "25s 1335"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs wt Fractions 3-4")+
			  			xlab("Mismatch Frequency (wtinputrep1)")+
			  			ylab("Mismatch Frequency (wt_f34rep1)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep1 and wt_f34rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="h2o2inputrep1_wt_f34rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep1_wt_f34rep1, aes(x=h2o2inputrep1_mis, y=wt_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("H2o2 input vs WT Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2inputrep1)")+
			  			ylab("Mismatch Frequency (wt_f34rep1)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep1 and h2o2inputrep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="wtinputrep1_h2o2inputrep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep1_h2o2inputrep1, aes(x=wtinputrep1_mis, y=h2o2inputrep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep1_h2o2inputrep1, position == "25s 776" | position == "25s 1051" | position == "25s 1449"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep1_h2o2inputrep1, position == "25s 776" | position == "25s 1051" | position == "25s 1449"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs h2O2 Input")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wtinputrep1)")+
			  			ylab("Mismatch Frequency (h2o2inputrep1)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For h2o2inputrep1 and h2o2_f34rep1
				#Export as PDF
				pdf(file="h2o2inputrep1_h2o2_f34rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep1_h2o2_f34rep1, aes(x=h2o2inputrep1_mis, y=h2o2_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("H2o2 input vs H2O2 Fractions 3-4")+
			  			xlab("Mismatch Frequency (h2o2inputrep1)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep1)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wt_f34rep1 and h2o2_f34rep1
				#Export as PDF
				pdf(file="wt_f34rep1_h2o2_f34rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_f34rep1_h2o2_f34rep1, aes(x=wt_f34rep1_mis, y=h2o2_f34rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("H2o2 fractions 3-4 vs WT Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wt_f34rep1)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep1)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()







			### REP2 
				#For h2o2_f34rep2 and wtinputrep2
				#Export as PDF
				pdf(file="wtinputrep2_h2o2_f34rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep2_h2o2_f34rep2, aes(x=wtinputrep2_mis, y=h2o2_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep2_h2o2_f34rep2, position == "25s 776"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep2_h2o2_f34rep2, position == "25s 776"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs h2O2 Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wtinputrep2)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep2)") +
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()

	
				#For wtinputrep2 and wt_f34rep2
				#Export as PDF
				pdf(file="wtinputrep2_wt_f34rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep2_wt_f34rep2, aes(x=wtinputrep2_mis, y=wt_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep2_wt_f34rep2, position == "25s 776" | position == "25s 1051" | position == "25s 1335"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep2_wt_f34rep2, position == "25s 776" | position == "25s 1051" | position == "25s 1335"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs wt Fractions 3-4")+
			  			xlab("Mismatch Frequency (wtinputrep2)")+
			  			ylab("Mismatch Frequency (wt_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep2 and wt_f34rep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="h2o2inputrep2_wt_f34rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep2_wt_f34rep2, aes(x=h2o2inputrep2_mis, y=wt_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("H2o2 input vs WT Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2inputrep2)")+
			  			ylab("Mismatch Frequency (wt_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep2 and h2o2inputrep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="wtinputrep2_h2o2inputrep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep2_h2o2inputrep2, aes(x=wtinputrep2_mis, y=h2o2inputrep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(wtinputrep2_h2o2inputrep2, position == "25s 776" | position == "25s 1051" | position == "25s 1449"), size=2, color="red")+
			  			geom_text_repel(data=subset(wtinputrep2_h2o2inputrep2, position == "25s 776" | position == "25s 1051" | position == "25s 1449" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Wt input vs h2O2 Input")+
			  			xlab("Mismatch Frequency (wtinputrep2)")+
			  			ylab("Mismatch Frequency (h2o2inputrep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For h2o2inputrep2 and h2o2_f34rep2
				#Export as PDF
				pdf(file="h2o2inputrep2_h2o2_f34rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep2_h2o2_f34rep2, aes(x=h2o2inputrep2_mis, y=h2o2_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("H2o2 input vs H2O2 Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2inputrep2)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wt_f34rep2 and h2o2_f34rep2
				#Export as PDF
				pdf(file="wt_f34rep2_h2o2_f34rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_f34rep2_h2o2_f34rep2, aes(x=wt_f34rep2_mis, y=h2o2_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("H2o2 fractions 3-4 vs WT Fractions 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wt_f34rep2)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()




			#REP1 VS REP2 

				#For h2o2_f34rep1 and h2o2_f34rep2
				#Export as PDF
				pdf(file="h2o2_f34rep1_h2o2_f34rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2_f34rep1_h2o2_f34rep2, aes(x=h2o2_f34rep1_mis, y=h2o2_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black",se=FALSE)+
			  			ggtitle("H2o2 fractioN 3-4  vs H2O2 Fraction 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2_f34rep1)")+
			  			ylab("Mismatch Frequency (h2o2_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For wtinputrep1 and wtinputrep2
				#Export as PDF
				pdf(file="wtinputrep1_wtinputrep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wtinputrep1_wtinputrep2, aes(x=wtinputrep1_mis, y=wtinputrep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Wt input vs WT Input")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wtinputrep1)")+
			  			ylab("Mismatch Frequency (wtinputrep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


				#For h2o2inputrep1 and h2o2inputrep2
				#Export as PDF
				pdf(file="h2o2inputrep1_h2o2inputrep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(h2o2inputrep1_h2o2inputrep2, aes(x=h2o2inputrep1_mis, y=h2o2inputrep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			  			ggtitle("H2o2 input vS H2O2 Input")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (h2o2inputrep1)")+
			  			ylab("Mismatch Frequency (h2o2inputrep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()



				#For wt_f34rep1 and wt_f34rep2
				#Export as PDF
				pdf(file="wt_f34rep1_wt_f34rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(wt_f34rep1_wt_f34rep2, aes(x=wt_f34rep1_mis, y=wt_f34rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Wt fraction 3-4  vs WT Fraction 3-4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (wt_f34rep1)")+
			  			ylab("Mismatch Frequency (wt_f34rep2)") +
			  			theme_bw()+
			  			theme(axis.text.x = element_text(face="bold", color="black",size=11),
	         				 axis.text.y = element_text(face="bold", color="black", size=11),
							plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
							axis.title.x = element_text(color="black", size=15, face="bold"),
							axis.title.y = element_text(color="black", size=15, face="bold"),
							panel.background = element_blank(),
							axis.line = element_line(colour = "black", size=0.5),
							legend.title = element_text(color = "black", size = 20,face="bold"),
	           				legend.text = element_text(color = "black", size=20))
					p
				dev.off()


```