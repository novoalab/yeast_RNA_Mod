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
		ko <- read.delim(input2,sep=",")
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