# Yeast Ribosomal RNA 

## Epinano Analysis of  Run (Fractions in Normal Condition only)

Here we only explain the analysis of Epinano output of the Graphmap Default mapped reads

For this, we created a table that contains positions and modifications, called rrna_mod_status.tsv



```R

###FINAL 
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
	#For f1rep1 (Rep 1)
		#Read the Epinano per site table
		f1rep1 <- read.delim("f1_rep1.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		f1rep1<-subset(f1rep1, cov>30)
		#Add a column replicating f1rep1
		f1rep1$sample <- rep("f1rep1",nrow(f1rep1)) 
		#Create a column for unique positions (18s 145)
		f1rep1$position<- paste(f1rep1$X.Ref,f1rep1$pos)
		#Use join function to add the modification information to the Epinano outputs
		f1rep1_st<- join(f1rep1, status, by="position")
		#Remove redundant columns
		f1rep1_st<- f1rep1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(f1rep1_st)<- c("base","f1rep1_q_mean", "f1rep1_q_median", "f1rep1_q_std", "f1rep1_mis", "f1rep1_del", "f1rep1_ins", "position", "ModStatus","Status") #Rename the columns

	#For f2rep1 (Rep 2)
		#Read the Epinano per site table
		f2rep1 <- read.delim("f2_rep1.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		f2rep1<-subset(f2rep1, cov>30)
		#Add a column replicating f2rep1
		f2rep1$sample <- rep("f2rep1",nrow(f2rep1)) 
		#Create a column for unique positions (18s 145)
		f2rep1$position<- paste(f2rep1$X.Ref,f2rep1$pos)
		#Use join function to add the modification information to the Epinano outputs
		f2rep1_st<- join(f2rep1, status, by="position")
		#Remove redundant columns
		f2rep1_st<- f2rep1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(f2rep1_st)<- c("base","f2rep1_q_mean", "f2rep1_q_median", "f2rep1_q_std", "f2rep1_mis", "f2rep1_del", "f2rep1_ins", "position", "ModStatus","Status") #Rename the columns

	#For f3rep1 (Rep 3)
		#Read the Epinano per site table
		f3rep1 <- read.delim("f3_rep1.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		f3rep1<-subset(f3rep1, cov>30)
		#Add a column replicating f3rep1
		f3rep1$sample <- rep("f3rep1",nrow(f3rep1)) 
		#Create a column for unique positions (18s 145)
		f3rep1$position<- paste(f3rep1$X.Ref,f3rep1$pos)
		#Use join function to add the modification information to the Epinano outputs
		f3rep1_st<- join(f3rep1, status, by="position")
		#Remove redundant columns
		f3rep1_st<- f3rep1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(f3rep1_st)<- c("base","f3rep1_q_mean", "f3rep1_q_median", "f3rep1_q_std", "f3rep1_mis", "f3rep1_del", "f3rep1_ins", "position", "ModStatus","Status") #Rename the columns

	#For f4rep1 (Rep 4)
		#Read the Epinano per site table
		f4rep1 <- read.delim("f4_rep1.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		f4rep1<-subset(f4rep1, cov>30)
		#Add a column replicating f4rep1
		f4rep1$sample <- rep("f4rep1",nrow(f4rep1)) 
		#Create a column for unique positions (18s 145)
		f4rep1$position<- paste(f4rep1$X.Ref,f4rep1$pos)
		#Use join function to add the modification information to the Epinano outputs
		f4rep1_st<- join(f4rep1, status, by="position")
		#Remove redundant columns
		f4rep1_st<- f4rep1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(f4rep1_st)<- c("base","f4rep1_q_mean", "f4rep1_q_median", "f4rep1_q_std", "f4rep1_mis", "f4rep1_del", "f4rep1_ins", "position", "ModStatus","Status") #Rename the columns




#Importing and manipulating the Epinano outputs
	#For f1rep2 (Rep 2)
		#Read the Epinano per site table
		f1rep2 <- read.delim("f1_rep2.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		f1rep2<-subset(f1rep2, cov>30)
		#Add a column replicating f1rep2
		f1rep2$sample <- rep("f1rep2",nrow(f1rep2)) 
		#Create a column for unique positions (18s 145)
		f1rep2$position<- paste(f1rep2$X.Ref,f1rep2$pos)
		#Use join function to add the modification information to the Epinano outputs
		f1rep2_st<- join(f1rep2, status, by="position")
		#Remove redundant columns
		f1rep2_st<- f1rep2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(f1rep2_st)<- c("base","f1rep2_q_mean", "f1rep2_q_median", "f1rep2_q_std", "f1rep2_mis", "f1rep2_del", "f1rep2_ins", "position", "ModStatus","Status") #Rename the columns

	#For f2rep2 (Rep 2)
		#Read the Epinano per site table
		f2rep2 <- read.delim("f2_rep2.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		f2rep2<-subset(f2rep2, cov>30)
		#Add a column replicating f2rep2
		f2rep2$sample <- rep("f2rep2",nrow(f2rep2)) 
		#Create a column for unique positions (18s 145)
		f2rep2$position<- paste(f2rep2$X.Ref,f2rep2$pos)
		#Use join function to add the modification information to the Epinano outputs
		f2rep2_st<- join(f2rep2, status, by="position")
		#Remove redundant columns
		f2rep2_st<- f2rep2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(f2rep2_st)<- c("base","f2rep2_q_mean", "f2rep2_q_median", "f2rep2_q_std", "f2rep2_mis", "f2rep2_del", "f2rep2_ins", "position", "ModStatus","Status") #Rename the columns

	#For f3rep2 (Rep 3)
		#Read the Epinano per site table
		f3rep2 <- read.delim("f3_rep2.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		f3rep2<-subset(f3rep2, cov>30)
		#Add a column replicating f3rep2
		f3rep2$sample <- rep("f3rep2",nrow(f3rep2)) 
		#Create a column for unique positions (18s 145)
		f3rep2$position<- paste(f3rep2$X.Ref,f3rep2$pos)
		#Use join function to add the modification information to the Epinano outputs
		f3rep2_st<- join(f3rep2, status, by="position")
		#Remove redundant columns
		f3rep2_st<- f3rep2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(f3rep2_st)<- c("base","f3rep2_q_mean", "f3rep2_q_median", "f3rep2_q_std", "f3rep2_mis", "f3rep2_del", "f3rep2_ins", "position", "ModStatus","Status") #Rename the columns

	#For f4rep2 (Rep 4)
		#Read the Epinano per site table
		f4rep2 <- read.delim("f4_rep2.bam.tsv.per.site.var.csv",sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		f4rep2<-subset(f4rep2, cov>30)
		#Add a column replicating f4rep2
		f4rep2$sample <- rep("f4rep2",nrow(f4rep2)) 
		#Create a column for unique positions (18s 145)
		f4rep2$position<- paste(f4rep2$X.Ref,f4rep2$pos)
		#Use join function to add the modification information to the Epinano outputs
		f4rep2_st<- join(f4rep2, status, by="position")
		#Remove redundant columns
		f4rep2_st<- f4rep2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(f4rep2_st)<- c("base","f4rep2_q_mean", "f4rep2_q_median", "f4rep2_q_std", "f4rep2_mis", "f4rep2_del", "f4rep2_ins", "position", "ModStatus","Status")  

#Creating comparison files
	#Rep1 Only
	#Use join function to merge f1rep1 and f4rep1
	f1rep1_f4rep1<- join(f1rep1_st, f4rep1_st, by="position")
	#Remove the columns that contain NA
	f1rep1_f4rep1<- subset(f1rep1_f4rep1, f4rep1_q_mean!="NA")

	#Use join function to merge f2rep1 and f4rep1
	f2rep1_f4rep1<- join(f2rep1_st, f4rep1_st, by="position")
	#Remove the columns that contain NA
	f2rep1_f4rep1<- subset(f2rep1_f4rep1, f4rep1_q_mean!="NA")

	#Use join function to merge f3rep1 and f4rep1
	f3rep1_f4rep1<- join(f3rep1_st, f4rep1_st, by="position")
	#Remove the columns that contain NA
	f3rep1_f4rep1<- subset(f3rep1_f4rep1, f4rep1_q_mean!="NA")

	#Use join function to merge f2rep1 and f4rep1
	f2rep1_f3rep1<- join(f2rep1_st, f3rep1_st, by="position")
	#Remove the columns that contain NA
	f2rep1_f3rep1<- subset(f2rep1_f3rep1, f3rep1_q_mean!="NA")
	
	#Use join function to merge f1rep1 and f4rep1
	f1rep1_f3rep1<- join(f1rep1_st, f3rep1_st, by="position")
	#Remove the columns that contain NA
	f1rep1_f3rep1<- subset(f1rep1_f3rep1, f3rep1_q_mean!="NA")
	
	#Use join function to merge f1rep1 and f4rep1
	f1rep1_f2rep1<- join(f1rep1_st, f2rep1_st, by="position")
	#Remove the columns that contain NA
	f1rep1_f2rep1<- subset(f1rep1_f2rep1, f2rep1_q_mean!="NA")

	#Rep2 Only
	#Use join function to merge f1rep2 and f4rep2
	f1rep2_f4rep2<- join(f1rep2_st, f4rep2_st, by="position")
	#Remove the columns that contain NA
	f1rep2_f4rep2<- subset(f1rep2_f4rep2, f4rep2_q_mean!="NA")

	#Use join function to merge f2rep2 and f4rep2
	f2rep2_f4rep2<- join(f2rep2_st, f4rep2_st, by="position")
	#Remove the columns that contain NA
	f2rep2_f4rep2<- subset(f2rep2_f4rep2, f4rep2_q_mean!="NA")

	#Use join function to merge f3rep2 and f4rep2
	f3rep2_f4rep2<- join(f3rep2_st, f4rep2_st, by="position")
	#Remove the columns that contain NA
	f3rep2_f4rep2<- subset(f3rep2_f4rep2, f4rep2_q_mean!="NA")

	#Use join function to merge f2rep2 and f4rep2
	f2rep2_f3rep2<- join(f2rep2_st, f3rep2_st, by="position")
	#Remove the columns that contain NA
	f2rep2_f3rep2<- subset(f2rep2_f3rep2, f3rep2_q_mean!="NA")
	
	#Use join function to merge f1rep2 and f4rep2
	f1rep2_f3rep2<- join(f1rep2_st, f3rep2_st, by="position")
	#Remove the columns that contain NA
	f1rep2_f3rep2<- subset(f1rep2_f3rep2, f3rep2_q_mean!="NA")
	
	#Use join function to merge f1rep2 and f4rep2
	f1rep2_f2rep2<- join(f1rep2_st, f2rep2_st, by="position")
	#Remove the columns that contain NA
	f1rep2_f2rep2<- subset(f1rep2_f2rep2, f2rep2_q_mean!="NA")


	#Rep1 VS Rep2
	#Use join function to merge f1rep1 and f1rep2
	f1rep1_f1rep2<- join(f1rep1_st, f1rep2_st, by="position")
	#Remove the columns that contain NA
	f1rep1_f1rep2<- subset(f1rep1_f1rep2, f1rep2_q_mean!="NA")

	#Use join function to merge f2rep1 and f2rep2
	f2rep1_f2rep2<- join(f2rep1_st, f2rep2_st, by="position")
	#Remove the columns that contain NA
	f2rep1_f2rep2<- subset(f2rep1_f2rep2, f2rep2_q_mean!="NA")

	#Use join function to merge f3rep1 and f3rep2
	f3rep1_f3rep2<- join(f3rep1_st, f3rep2_st, by="position")
	#Remove the columns that contain NA
	f3rep1_f3rep2<- subset(f3rep1_f3rep2, f3rep2_q_mean!="NA")

	#Use join function to merge f4rep1 and f4rep2
	f4rep1_f4rep2<- join(f4rep1_st, f4rep2_st, by="position")
	#Remove the columns that contain NA
	f4rep1_f4rep2<- subset(f4rep1_f4rep2, f4rep2_q_mean!="NA")




	#Difference
	f1rep1_f4rep1$diff_mis <-  abs(f1rep1_f4rep1$f1rep1_mis-f1rep1_f4rep1$f4rep1_mis)
	f1rep1_f3rep1$diff_mis <-  abs(f1rep1_f3rep1$f1rep1_mis-f1rep1_f3rep1$f3rep1_mis)
	f1rep1_f2rep1$diff_mis <-  abs(f1rep1_f2rep1$f1rep1_mis-f1rep1_f2rep1$f2rep1_mis)
	f2rep1_f4rep1$diff_mis <-  abs(f2rep1_f4rep1$f2rep1_mis-f2rep1_f4rep1$f4rep1_mis)
	f3rep1_f4rep1$diff_mis <-  abs(f3rep1_f4rep1$f3rep1_mis-f3rep1_f4rep1$f4rep1_mis)
	f2rep1_f3rep1$diff_mis <-  abs(f2rep1_f3rep1$f2rep1_mis-f2rep1_f3rep1$f3rep1_mis)


	f1rep2_f4rep2$diff_mis <-  abs(f1rep2_f4rep2$f1rep2_mis-f1rep2_f4rep2$f4rep2_mis)
	f1rep2_f3rep2$diff_mis <-  abs(f1rep2_f3rep2$f1rep2_mis-f1rep2_f3rep2$f3rep2_mis)
	f1rep2_f2rep2$diff_mis <-  abs(f1rep2_f2rep2$f1rep2_mis-f1rep2_f2rep2$f2rep2_mis)
	f2rep2_f4rep2$diff_mis <-  abs(f2rep2_f4rep2$f2rep2_mis-f2rep2_f4rep2$f4rep2_mis)
	f3rep2_f4rep2$diff_mis <-  abs(f3rep2_f4rep2$f3rep2_mis-f3rep2_f4rep2$f4rep2_mis)
	f2rep2_f3rep2$diff_mis <-  abs(f2rep2_f3rep2$f2rep2_mis-f2rep2_f3rep2$f3rep2_mis)


	f1rep1_f1rep2$diff_mis <-  abs(f1rep1_f1rep2$f1rep1_mis-f1rep1_f1rep2$f1rep2_mis)
	f2rep1_f2rep2$diff_mis <-  abs(f2rep1_f2rep2$f2rep1_mis-f2rep1_f2rep2$f2rep2_mis)
	f3rep1_f3rep2$diff_mis <-  abs(f3rep1_f3rep2$f3rep1_mis-f3rep1_f3rep2$f3rep2_mis)
	f4rep1_f4rep2$diff_mis <-  abs(f4rep1_f4rep2$f4rep1_mis-f4rep1_f4rep2$f4rep2_mis)



		#Plotting the scatter plots using ggplot2
			#Rep1
				#For f4rep1 and f1rep1
				#Export as PDF
				pdf(file="f1rep1_f4rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f4rep1, aes(x=f1rep1_mis, y=f4rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f4rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f4rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 4")+
			  			xlab("Mismatch Frequency (f1rep1)")+
			  			ylab("Mismatch Frequency (f4rep1)") +
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

	
				#For f1rep1 and f3rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep1_f3rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f3rep1, aes(x=f1rep1_mis, y=f3rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f3rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f3rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 3")+
			  			xlab("Mismatch Frequency (f1rep1)")+
			  			ylab("Mismatch Frequency (f3rep1)") +
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


				#For f1rep1 and f3rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f2rep1_f3rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f3rep1, aes(x=f2rep1_mis, y=f3rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep1_f3rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep1_f3rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 vs 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f2rep1)")+
			  			ylab("Mismatch Frequency (f3rep1)") +
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


				#For f1rep1 and f2rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep1_f2rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f2rep1, aes(x=f1rep1_mis, y=f2rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f2rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f2rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 2")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f1rep1)")+
			  			ylab("Mismatch Frequency (f2rep1)") +
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


				#For f2rep1 and f4rep1
				#Export as PDF
				pdf(file="f2rep1_f4rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f4rep1, aes(x=f2rep1_mis, y=f4rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep1_f4rep1, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep1_f4rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 vs 4")+
			  			xlab("Mismatch Frequency (f2rep1)")+
			  			ylab("Mismatch Frequency (f4rep1)") +
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


				#For f3rep1 and f4rep1
				#Export as PDF
				pdf(file="f3rep1_f4rep1_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep1_f4rep1, aes(x=f3rep1_mis, y=f4rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f3rep1_f4rep1, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f3rep1_f4rep1, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 3 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f3rep1)")+
			  			ylab("Mismatch Frequency (f4rep1)") +
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







			### REP2 
				#For f4rep2 and f1rep2
				#Export as PDF
				pdf(file="f1rep2_f4rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f4rep2, aes(x=f1rep2_mis, y=f4rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f4rep2, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f4rep2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f1rep2)")+
			  			ylab("Mismatch Frequency (f4rep2)") +
			  			
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

	
				#For f1rep2 and f3rep2
				#Export as PDF
				pdf(file="f1rep2_f3rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f3rep2, aes(x=f1rep2_mis, y=f3rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f3rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f3rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 3")+
			  			xlab("Mismatch Frequency (f1rep2)")+
			  			ylab("Mismatch Frequency (f3rep2)")+
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


				#For f1rep2 and f3rep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f2rep2_f3rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep2_f3rep2, aes(x=f2rep2_mis, y=f3rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep2_f3rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep2_f3rep2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 vs 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f2rep2)")+
			  			ylab("Mismatch Frequency (f3rep2)") +
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


				#For f1rep2 and f2rep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep2_f2rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f2rep2, aes(x=f1rep2_mis, y=f2rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f2rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f2rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 2")+
			  			xlab("Mismatch Frequency (f1rep2)")+
			  			ylab("Mismatch Frequency (f2rep2)") +
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


				#For f2rep2 and f4rep2
				#Export as PDF
				pdf(file="f2rep2_f4rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep2_f4rep2, aes(x=f2rep2_mis, y=f4rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep2_f4rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep2_f4rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f2rep2)")+
			  			ylab("Mismatch Frequency (f4rep2)") +
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


				#For f3rep2 and f4rep2
				#Export as PDF
				pdf(file="f3rep2_f4rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep2_f4rep2, aes(x=f3rep2_mis, y=f4rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f3rep2_f4rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f3rep2_f4rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 3 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f3rep2)")+
			  			ylab("Mismatch Frequency (f4rep2)") +
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




			#REP1 VS REP2 

				#For f4rep1 and f4rep2
				#Export as PDF
				pdf(file="f4rep1_f4rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f4rep1_f4rep2, aes(x=f4rep1_mis, y=f4rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black",se=FALSE)+
			 			geom_point(data=subset(f4rep1_f4rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f4rep1_f4rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 4 VS 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f4rep1)")+
			  			ylab("Mismatch Frequency (f4rep2)") +
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


				#For f1rep1 and f1rep2
				#Export as PDF
				pdf(file="f1rep1_f1rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f1rep2, aes(x=f1rep1_mis, y=f1rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f1rep2, diff_mis>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f1rep2, diff_mis>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 VS 1")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f1rep1)")+
			  			ylab("Mismatch Frequency (f1rep2)") +
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


				#For f2rep1 and f2rep2
				#Export as PDF
				pdf(file="f2rep1_f2rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f2rep2, aes(x=f2rep1_mis, y=f2rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep1_f2rep2, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep1_f2rep2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 VS 2")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f2rep1)")+
			  			ylab("Mismatch Frequency (f2rep2)") +
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



				#For f3rep1 and f3rep2
				#Export as PDF
				pdf(file="f3rep1_f3rep2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep1_f3rep2, aes(x=f3rep1_mis, y=f3rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f3rep1_f3rep2, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f3rep1_f3rep2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 3 VS 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f3rep1)")+
			  			ylab("Mismatch Frequency (f3rep2)") +
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






		# VENN DIAGRAMS BETWEEN REPLICATES


					#Rep1 Only

					f1rep1_f4rep1.sig<- subset(f1rep1_f4rep1, diff_mis>0.1)
					f2rep1_f4rep1.sig<- subset(f2rep1_f4rep1, diff_mis>0.1)
					f3rep1_f4rep1.sig<- subset(f3rep1_f4rep1, diff_mis>0.1)
					f2rep1_f3rep1.sig<- subset(f2rep1_f3rep1, diff_mis>0.1)
					f1rep1_f3rep1.sig<- subset(f1rep1_f3rep1, diff_mis>0.1)
					f1rep1_f2rep1.sig<- subset(f1rep1_f2rep1, diff_mis>0.1)	
					f1rep2_f4rep2.sig<- subset(f1rep2_f4rep2, diff_mis>0.1)
					f2rep2_f4rep2.sig<- subset(f2rep2_f4rep2, diff_mis>0.1)
					f3rep2_f4rep2.sig<- subset(f3rep2_f4rep2, diff_mis>0.1)
					f2rep2_f3rep2.sig<- subset(f2rep2_f3rep2, diff_mis>0.1)
					f1rep2_f3rep2.sig<- subset(f1rep2_f3rep2, diff_mis>0.1)
					f1rep2_f2rep2.sig<- subset(f1rep2_f2rep2, diff_mis>0.1)



					f1_f4.common<- intersect(f1rep1_f4rep1.sig$position,f1rep2_f4rep2.sig$position)
					f2_f4.common<- intersect(f2rep1_f4rep1.sig$position,f2rep2_f4rep2.sig$position)
					f3_f4.common<- intersect(f3rep1_f4rep1.sig$position,f3rep2_f4rep2.sig$position)
					f2_f3.common<- intersect(f2rep1_f3rep1.sig$position,f2rep2_f3rep2.sig$position)
					f1_f3.common<- intersect(f1rep1_f3rep1.sig$position,f1rep2_f3rep2.sig$position)
					f1_f2.common<- intersect(f1rep1_f2rep1.sig$position,f1rep2_f2rep2.sig$position)


					write.table(f1_f4.common, file="f1_f4.common.txt")
					write.table(f2_f4.common, file="f2_f4.common.txt")
					write.table(f3_f4.common, file="f3_f4.common.txt")
					write.table(f2_f3.common, file="f2_f3.common.txt")
					write.table(f1_f3.common, file="f1_f3.common.txt")
					write.table(f1_f2.common, file="f1_f2.common.txt")



					#Venn diagram for common specific genes combination
					library(VennDiagram)
					pdf("Venn_Diagram_F1_F4_humanvsmouse.pdf",height=8,width=8)
					draw.pairwise.venn(length(f1rep1_f4rep1.sig$position), length(f1rep2_f4rep2.sig$position), length(f1_f4.common), category = c("F1vsF4 Rep1", "F1vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F2_F4_humanvsmouse.pdf",height=8,width=8)
					draw.pairwise.venn(length(f2rep1_f4rep1.sig$position), length(f2rep2_f4rep2.sig$position), length(f2_f4.common), category = c("F2vsF4 Rep1", "F2vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F3_F4_humanvsmouse.pdf",height=8,width=8)
					draw.pairwise.venn(length(f3rep1_f4rep1.sig$position), length(f3rep2_f4rep2.sig$position), length(f3_f4.common), category = c("F3vsF4 Rep1", "F3vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F2_F3_humanvsmouse.pdf",height=8,width=8)
					draw.pairwise.venn(length(f2rep1_f3rep1.sig$position), length(f2rep2_f3rep2.sig$position), length(f2_f3.common), category = c("F2vsF3 Rep1", "F2vsF3 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F1_F3_humanvsmouse.pdf",height=8,width=8)
					draw.pairwise.venn(length(f1rep1_f3rep1.sig$position), length(f1rep2_f3rep2.sig$position), length(f1_f3.common), category = c("F1vsF3 Rep1", "F1vsF3 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()


					pdf("Venn_Diagram_F1_F2_humanvsmouse.pdf",height=8,width=8)
					draw.pairwise.venn(length(f1rep1_f2rep1.sig$position), length(f1rep2_f2rep2.sig$position), length(f1_f2.common), category = c("F1vsF2 Rep1", "F1vsF2 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()












					#ONLY Y positions

					f1rep1_f4rep1.sig<- subset(f1rep1_f4rep1, diff_mis>0.1 & ModStatus == "Y")
					f2rep1_f4rep1.sig<- subset(f2rep1_f4rep1, diff_mis>0.1 & ModStatus == "Y")
					f3rep1_f4rep1.sig<- subset(f3rep1_f4rep1, diff_mis>0.1 & ModStatus == "Y")
					f2rep1_f3rep1.sig<- subset(f2rep1_f3rep1, diff_mis>0.1 & ModStatus == "Y")
					f1rep1_f3rep1.sig<- subset(f1rep1_f3rep1, diff_mis>0.1 & ModStatus == "Y")
					f1rep1_f2rep1.sig<- subset(f1rep1_f2rep1, diff_mis>0.1 & ModStatus == "Y")	
					f1rep2_f4rep2.sig<- subset(f1rep2_f4rep2, diff_mis>0.1 & ModStatus == "Y")
					f2rep2_f4rep2.sig<- subset(f2rep2_f4rep2, diff_mis>0.1 & ModStatus == "Y")
					f3rep2_f4rep2.sig<- subset(f3rep2_f4rep2, diff_mis>0.1 & ModStatus == "Y")
					f2rep2_f3rep2.sig<- subset(f2rep2_f3rep2, diff_mis>0.1 & ModStatus == "Y")
					f1rep2_f3rep2.sig<- subset(f1rep2_f3rep2, diff_mis>0.1 & ModStatus == "Y")
					f1rep2_f2rep2.sig<- subset(f1rep2_f2rep2, diff_mis>0.1 & ModStatus == "Y")



					f1_f4.common<- intersect(f1rep1_f4rep1.sig$position,f1rep2_f4rep2.sig$position)
					f2_f4.common<- intersect(f2rep1_f4rep1.sig$position,f2rep2_f4rep2.sig$position)
					f3_f4.common<- intersect(f3rep1_f4rep1.sig$position,f3rep2_f4rep2.sig$position)
					f2_f3.common<- intersect(f2rep1_f3rep1.sig$position,f2rep2_f3rep2.sig$position)
					f1_f3.common<- intersect(f1rep1_f3rep1.sig$position,f1rep2_f3rep2.sig$position)
					f1_f2.common<- intersect(f1rep1_f2rep1.sig$position,f1rep2_f2rep2.sig$position)


					write.table(f1_f4.common, file="f1_f4.common_Y.txt")
					write.table(f2_f4.common, file="f2_f4.common_Y.txt")
					write.table(f3_f4.common, file="f3_f4.common_Y.txt")
					write.table(f2_f3.common, file="f2_f3.common_Y.txt")
					write.table(f1_f3.common, file="f1_f3.common_Y.txt")
					write.table(f1_f2.common, file="f1_f2.common_Y.txt")



					#Venn diagram for common specific genes combination
					library(VennDiagram)
					pdf("Venn_Diagram_F1_F4_humanvsmouse_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(f1rep1_f4rep1.sig$position), length(f1rep2_f4rep2.sig$position), length(f1_f4.common), category = c("F1vsF4 Rep1", "F1vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F2_F4_humanvsmouse_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(f2rep1_f4rep1.sig$position), length(f2rep2_f4rep2.sig$position), length(f2_f4.common), category = c("F2vsF4 Rep1", "F2vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F3_F4_humanvsmouse_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(f3rep1_f4rep1.sig$position), length(f3rep2_f4rep2.sig$position), length(f3_f4.common), category = c("F3vsF4 Rep1", "F3vsF4 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F2_F3_humanvsmouse_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(f2rep1_f3rep1.sig$position), length(f2rep2_f3rep2.sig$position), length(f2_f3.common), category = c("F2vsF3 Rep1", "F2vsF3 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()

					pdf("Venn_Diagram_F1_F3_humanvsmouse_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(f1rep1_f3rep1.sig$position), length(f1rep2_f3rep2.sig$position), length(f1_f3.common), category = c("F1vsF3 Rep1", "F1vsF3 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()


					pdf("Venn_Diagram_F1_F2_humanvsmouse_Y.pdf",height=8,width=8)
					draw.pairwise.venn(length(f1rep1_f2rep1.sig$position), length(f1rep2_f2rep2.sig$position), length(f1_f2.common), category = c("F1vsF2 Rep1", "F1vsF2 Rep2"), lty = rep("blank", 
					    2), fill = c("#004a2f", "#ffa323"), alpha = rep(0.5, 2), cat.pos = c(0, 
					    0), cat.dist = rep(0.025, 2))
					dev.off()















		### PLOT THE SCATTER WITH ALL SIGNIFICANT POSITIONS (NOT ONLY Y)

		#Plotting the scatter plots using ggplot2
			#Rep1
				#For f4rep1 and f1rep1
				#Export as PDF
				pdf(file="f1rep1_f4rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f4rep1, aes(x=f1rep1_mis, y=f4rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f4rep1, position == "25s 776"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f4rep1, position == "25s 776"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 4")+
			  			xlab("Mismatch Frequency (f1rep1)")+
			  			ylab("Mismatch Frequency (f4rep1)") +
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

	
				#For f1rep1 and f3rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep1_f3rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f3rep1, aes(x=f1rep1_mis, y=f3rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f3rep1, position == "25s 776" | position == "25s 1051" | position == "25s 1335"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f3rep1, position == "25s 776" | position == "25s 1051" | position == "25s 1335"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 3")+
			  			xlab("Mismatch Frequency (f1rep1)")+
			  			ylab("Mismatch Frequency (f3rep1)") +
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


				#For f1rep1 and f3rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f2rep1_f3rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f3rep1, aes(x=f2rep1_mis, y=f3rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 2 vs 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f2rep1)")+
			  			ylab("Mismatch Frequency (f3rep1)") +
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


				#For f1rep1 and f2rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep1_f2rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f2rep1, aes(x=f1rep1_mis, y=f2rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f2rep1, position == "25s 776" | position == "25s 1051" | position == "25s 1449"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f2rep1, position == "25s 776" | position == "25s 1051" | position == "25s 1449"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 2")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f1rep1)")+
			  			ylab("Mismatch Frequency (f2rep1)") +
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


				#For f2rep1 and f4rep1
				#Export as PDF
				pdf(file="f2rep1_f4rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f4rep1, aes(x=f2rep1_mis, y=f4rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 2 vs 4")+
			  			xlab("Mismatch Frequency (f2rep1)")+
			  			ylab("Mismatch Frequency (f4rep1)")+
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


				#For f3rep1 and f4rep1
				#Export as PDF
				pdf(file="f3rep1_f4rep1_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep1_f4rep1, aes(x=f3rep1_mis, y=f4rep1_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 3 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f3rep1)")+
			  			ylab("Mismatch Frequency (f4rep1)") +
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







			### REP2 
				#For f4rep2 and f1rep2
				#Export as PDF
				pdf(file="f1rep2_f4rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f4rep2, aes(x=f1rep2_mis, y=f4rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f4rep2, position == "25s 776"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f4rep2, position == "25s 776"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f1rep2)")+
			  			ylab("Mismatch Frequency (f4rep2)") +
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

	
				#For f1rep2 and f3rep2
				#Export as PDF
				pdf(file="f1rep2_f3rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f3rep2, aes(x=f1rep2_mis, y=f3rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f3rep2, position == "25s 776" | position == "25s 1051" | position == "25s 1335"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f3rep2, position == "25s 776" | position == "25s 1051" | position == "25s 1335"), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 3")+
			  			xlab("Mismatch Frequency (f1rep2)")+
			  			ylab("Mismatch Frequency (f3rep2)") +
			  			theme_bw()+
			  			xlim(0,1)+
			  			ylim(0,1)+
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


				#For f1rep2 and f3rep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f2rep2_f3rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep2_f3rep2, aes(x=f2rep2_mis, y=f3rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 2 vs 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f2rep2)")+
			  			ylab("Mismatch Frequency (f3rep2)") +
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


				#For f1rep2 and f2rep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep2_f2rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f2rep2, aes(x=f1rep2_mis, y=f2rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f2rep2, position == "25s 776" | position == "25s 1051" | position == "25s 1449"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f2rep2, position == "25s 776" | position == "25s 1051" | position == "25s 1449" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 2")+
			  			xlab("Mismatch Frequency (f1rep2)")+
			  			ylab("Mismatch Frequency (f2rep2)") +
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


				#For f2rep2 and f4rep2
				#Export as PDF
				pdf(file="f2rep2_f4rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep2_f4rep2, aes(x=f2rep2_mis, y=f4rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 2 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f2rep2)")+
			  			ylab("Mismatch Frequency (f4rep2)") +
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


				#For f3rep2 and f4rep2
				#Export as PDF
				pdf(file="f3rep2_f4rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep2_f4rep2, aes(x=f3rep2_mis, y=f4rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 3 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f3rep2)")+
			  			ylab("Mismatch Frequency (f4rep2)") +
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




			#REP1 VS REP2 

				#For f4rep1 and f4rep2
				#Export as PDF
				pdf(file="f4rep1_f4rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f4rep1_f4rep2, aes(x=f4rep1_mis, y=f4rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black",se=FALSE)+
			  			ggtitle("Fractions 4 VS 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f4rep1)")+
			  			ylab("Mismatch Frequency (f4rep2)") +
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


				#For f1rep1 and f1rep2
				#Export as PDF
				pdf(file="f1rep1_f1rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f1rep2, aes(x=f1rep1_mis, y=f1rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 1 VS 1")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f1rep1)")+
			  			ylab("Mismatch Frequency (f1rep2)") +
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


				#For f2rep1 and f2rep2
				#Export as PDF
				pdf(file="f2rep1_f2rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f2rep2, aes(x=f2rep1_mis, y=f2rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 2 VS 2")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f2rep1)")+
			  			ylab("Mismatch Frequency (f2rep2)") +
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



				#For f3rep1 and f3rep2
				#Export as PDF
				pdf(file="f3rep1_f3rep2_scatter_mismatch_allsignificant.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep1_f3rep2, aes(x=f3rep1_mis, y=f3rep2_mis)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			  			ggtitle("Fractions 3 VS 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Mismatch Frequency (f3rep1)")+
			  			ylab("Mismatch Frequency (f3rep2)") +
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
















# SCATTER PLOTS FOR INSERTION


	#Difference
	f1rep1_f4rep1$diff_ins <-  abs(f1rep1_f4rep1$f1rep1_ins-f1rep1_f4rep1$f4rep1_ins)
	f1rep1_f3rep1$diff_ins <-  abs(f1rep1_f3rep1$f1rep1_ins-f1rep1_f3rep1$f3rep1_ins)
	f1rep1_f2rep1$diff_ins <-  abs(f1rep1_f2rep1$f1rep1_ins-f1rep1_f2rep1$f2rep1_ins)
	f2rep1_f4rep1$diff_ins <-  abs(f2rep1_f4rep1$f2rep1_ins-f2rep1_f4rep1$f4rep1_ins)
	f3rep1_f4rep1$diff_ins <-  abs(f3rep1_f4rep1$f3rep1_ins-f3rep1_f4rep1$f4rep1_ins)
	f2rep1_f3rep1$diff_ins <-  abs(f2rep1_f3rep1$f2rep1_ins-f2rep1_f3rep1$f3rep1_ins)


	f1rep2_f4rep2$diff_ins <-  abs(f1rep2_f4rep2$f1rep2_ins-f1rep2_f4rep2$f4rep2_ins)
	f1rep2_f3rep2$diff_ins <-  abs(f1rep2_f3rep2$f1rep2_ins-f1rep2_f3rep2$f3rep2_ins)
	f1rep2_f2rep2$diff_ins <-  abs(f1rep2_f2rep2$f1rep2_ins-f1rep2_f2rep2$f2rep2_ins)
	f2rep2_f4rep2$diff_ins <-  abs(f2rep2_f4rep2$f2rep2_ins-f2rep2_f4rep2$f4rep2_ins)
	f3rep2_f4rep2$diff_ins <-  abs(f3rep2_f4rep2$f3rep2_ins-f3rep2_f4rep2$f4rep2_ins)
	f2rep2_f3rep2$diff_ins <-  abs(f2rep2_f3rep2$f2rep2_ins-f2rep2_f3rep2$f3rep2_ins)


	f1rep1_f1rep2$diff_ins <-  abs(f1rep1_f1rep2$f1rep1_ins-f1rep1_f1rep2$f1rep2_ins)
	f2rep1_f2rep2$diff_ins <-  abs(f2rep1_f2rep2$f2rep1_ins-f2rep1_f2rep2$f2rep2_ins)
	f3rep1_f3rep2$diff_ins <-  abs(f3rep1_f3rep2$f3rep1_ins-f3rep1_f3rep2$f3rep2_ins)
	f4rep1_f4rep2$diff_ins <-  abs(f4rep1_f4rep2$f4rep1_ins-f4rep1_f4rep2$f4rep2_ins)



		#Plotting the scatter plots using ggplot2
			#Rep1
				#For f4rep1 and f1rep1
				#Export as PDF
				pdf(file="f1rep1_f4rep1_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f4rep1, aes(x=f1rep1_ins, y=f4rep1_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f4rep1, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f4rep1, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 4")+
			  			xlab("Insertion Frequency (f1rep1)")+
			  			ylab("Insertion Frequency (f4rep1)") +
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

	
				#For f1rep1 and f3rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep1_f3rep1_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f3rep1, aes(x=f1rep1_ins, y=f3rep1_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f3rep1, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f3rep1, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 3")+
			  			xlab("Insertion Frequency (f1rep1)")+
			  			ylab("Insertion Frequency (f3rep1)") +
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


				#For f1rep1 and f3rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f2rep1_f3rep1_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f3rep1, aes(x=f2rep1_ins, y=f3rep1_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep1_f3rep1, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep1_f3rep1, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 vs 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f2rep1)")+
			  			ylab("Insertion Frequency (f3rep1)") +
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


				#For f1rep1 and f2rep1
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep1_f2rep1_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f2rep1, aes(x=f1rep1_ins, y=f2rep1_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f2rep1, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f2rep1, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 2")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f1rep1)")+
			  			ylab("Insertion Frequency (f2rep1)") +
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


				#For f2rep1 and f4rep1
				#Export as PDF
				pdf(file="f2rep1_f4rep1_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f4rep1, aes(x=f2rep1_ins, y=f4rep1_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep1_f4rep1, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep1_f4rep1, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 vs 4")+
			  			xlab("Insertion Frequency (f2rep1)")+
			  			ylab("Insertion Frequency (f4rep1)") +
			  			theme_bw()+
			  			xlim(0,1)+
			  			ylim(0,1)+
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


				#For f3rep1 and f4rep1
				#Export as PDF
				pdf(file="f3rep1_f4rep1_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep1_f4rep1, aes(x=f3rep1_ins, y=f4rep1_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f3rep1_f4rep1, diff_ins>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f3rep1_f4rep1, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 3 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f3rep1)")+
			  			ylab("Insertion Frequency (f4rep1)") +
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







			### REP2 
				#For f4rep2 and f1rep2
				#Export as PDF
				pdf(file="f1rep2_f4rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f4rep2, aes(x=f1rep2_ins, y=f4rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f4rep2, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f4rep2, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f1rep2)")+
			  			ylab("Insertion Frequency (f4rep2)") +
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

	
				#For f1rep2 and f3rep2
				#Export as PDF
				pdf(file="f1rep2_f3rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f3rep2, aes(x=f1rep2_ins, y=f3rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f3rep2, diff_ins>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f3rep2, diff_ins>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 3")+
			  			xlab("Insertion Frequency (f1rep2)")+
			  			ylab("Insertion Frequency (f3rep2)") +
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


				#For f1rep2 and f3rep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f2rep2_f3rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep2_f3rep2, aes(x=f2rep2_ins, y=f3rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep2_f3rep2, diff_ins>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep2_f3rep2, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 vs 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f2rep2)")+
			  			ylab("Insertion Frequency (f3rep2)") +
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


				#For f1rep2 and f2rep2
				#Create a subset of the significantly different positions
				#Export as PDF
				pdf(file="f1rep2_f2rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep2_f2rep2, aes(x=f1rep2_ins, y=f2rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep2_f2rep2, diff_ins>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep2_f2rep2, diff_ins>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 vs 2")+
			  			xlab("Insertion Frequency (f1rep2)")+
			  			ylab("Insertion Frequency (f2rep2)") +
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


				#For f2rep2 and f4rep2
				#Export as PDF
				pdf(file="f2rep2_f4rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep2_f4rep2, aes(x=f2rep2_ins, y=f4rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep2_f4rep2, diff_ins>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep2_f4rep2, diff_ins>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f2rep2)")+
			  			ylab("Insertion Frequency (f4rep2)") +
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


				#For f3rep2 and f4rep2
				#Export as PDF
				pdf(file="f3rep2_f4rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep2_f4rep2, aes(x=f3rep2_ins, y=f4rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f3rep2_f4rep2, diff_ins>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f3rep2_f4rep2, diff_ins>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 3 vs 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f3rep2)")+
			  			ylab("Insertion Frequency (f4rep2)") +
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




			#REP1 VS REP2 

				#For f4rep1 and f4rep2
				#Export as PDF
				pdf(file="f4rep1_f4rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f4rep1_f4rep2, aes(x=f4rep1_ins, y=f4rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			 			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f4rep1_f4rep2, diff_ins>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f4rep1_f4rep2, diff_ins>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 4 VS 4")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f4rep1)")+
			  			ylab("Insertion Frequency (f4rep2)") +
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


				#For f1rep1 and f1rep2
				#Export as PDF
				pdf(file="f1rep1_f1rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f1rep1_f1rep2, aes(x=f1rep1_ins, y=f1rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f1rep1_f1rep2, diff_ins>0.1& ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f1rep1_f1rep2, diff_ins>0.1& ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 1 VS 1")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f1rep1)")+
			  			ylab("Insertion Frequency (f1rep2)") +
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


				#For f2rep1 and f2rep2
				#Export as PDF
				pdf(file="f2rep1_f2rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f2rep1_f2rep2, aes(x=f2rep1_ins, y=f2rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f2rep1_f2rep2, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f2rep1_f2rep2, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 2 VS 2")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f2rep1)")+
			  			ylab("Insertion Frequency (f2rep2)") +
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



				#For f3rep1 and f3rep2
				#Export as PDF
				pdf(file="f3rep1_f3rep2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
					p<-ggplot(f3rep1_f3rep2, aes(x=f3rep1_ins, y=f3rep2_ins)) +
			 			geom_point(size=1, color="grey")+
			  			geom_abline(slope=1, intercept=0)+
			 			geom_point(data=subset(f3rep1_f3rep2, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
			  			geom_text_repel(data=subset(f3rep1_f3rep2, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			  			ggtitle("Fractions 3 VS 3")+
			  			xlim(0,1)+
			  			ylim(0,1)+
			  			xlab("Insertion Frequency (f3rep1)")+
			  			ylab("Insertion Frequency (f3rep2)") +
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



	#SCATTER FOR DELETION


```