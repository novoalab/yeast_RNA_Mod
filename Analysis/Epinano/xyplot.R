#Executable for scatter plot
#Rscript xyplot.R wt.bam.tsv.per.site.var.csv sn3.bam.tsv.per.site.var.csv
#Load the libraries needed for this script
library(plyr)
library(ggplot2)
library(ggrepel) 

#Input
args <- commandArgs(trailingOnly = TRUE)
cond1 <- args[1] #1st variable
cond2 <- args[2] #2nd variable

#Importing and manipulating the rna mod status file
	#Read the table for modification positions		
	status<- read.delim("rrna_mod_status.tsv")
	#Create a column for unique positions (18s 145)
	status$position<- paste(status$Chr,status$Position)

#Importing and manipulating the Epinano outputs
	#For condition1 s
		#Read the Epinano per site table
		condition1 <- read.delim(cond1,sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		condition1<-subset(condition1, cov>50)
		#Add a column replicating condition1
		condition1$sample <- rep("condition1",nrow(condition1)) 
		#Create a column for unique positions (18s 145)
		condition1$position<- paste(condition1$X.Ref,condition1$pos)
		#Use join function to add the modification information to the Epinano outputs
		condition1_st<- join(condition1, status, by="position")
		#Remove redundant columns
		condition1_st<- condition1_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(condition1_st)<- c("base","condition1_q_mean", "condition1_q_median", "condition1_q_std", "condition1_mis", "condition1_del", "condition1_ins", "position", "ModStatus","Status") #Rename the columns


	#For condition2
		#Read the Epinano per site table
		condition2 <- read.delim(cond2,sep=",")
		#Coverage threshold (remove the positions with less than X coverage)
		condition2<-subset(condition2, cov>50)
		#Add a column replicating condition2
		condition2$sample <- rep("condition2",nrow(condition2)) 
		#Create a column for unique positions (18s 145)
		condition2$position<- paste(condition2$X.Ref,condition2$pos)
		#Use join function to add the modification information to the Epinano outputs
		condition2_st<- join(condition2, status, by="position")
		#Remove redundant columns
		condition2_st<- condition2_st[,c("base","q_mean", "q_median", "q_std", "mis", "del","ins","position", "ModStatus","Status")]
		#Rename columns
		colnames(condition2_st)<- c("base","condition2_q_mean", "condition2_q_median", "condition2_q_std", "condition2_mis", "condition2_del", "condition2_ins", "position", "ModStatus","Status") #Rename the columns


#Creating comparison files
	#Use join function to merge Condition 1 and SN3Condition 2
	condition1_condition2<- join(condition1_st, condition2_st, by="position")
	#Remove the columns that contain NA
	condition1_condition2<- subset(condition1_condition2, condition2_q_mean!="NA")

# Scatter plots of Mismatch frequencies
	#Calculate the thresholds for the significant positions (to be labeled)
	condition1_condition2$diff_mis<- abs(condition1_condition2$condition1_mis - condition1_condition2$condition2_mis)
	condition1_condition2$diff_ins<- abs(condition1_condition2$condition1_ins - condition1_condition2$condition2_ins)
	condition1_condition2$diff_del<- abs(condition1_condition2$condition1_del - condition1_condition2$condition2_del)


#Plotting the scatter plots using ggplot2
	#For Condition 1 and SN3Condition 2
		#Export as PDF
		pdf(file="condition1_condition2_scatter_mismatch.pdf",height=5,width=5,onefile=FALSE)
			p<-ggplot(condition1_condition2, aes(x=condition1_mis, y=condition2_mis)) +
	 			geom_point(size=1, color="grey")+
	 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
	 			geom_point(data=subset(condition1_condition2, diff_mis>0.1 & ModStatus == "Y"), size=2, color="red")+
	  			geom_text_repel(data=subset(condition1_condition2, diff_mis>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
	  			ggtitle("Conditions 1 vs 2")+
	  			xlim(0,1)+
	  			ylim(0,1)+
	  			xlab("Mismatch Frequency (Condition 1)")+
	  			ylab("Mismatch Frequency (Condition 2)") +
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
		#For Condition 1 and Condition 2
			#Export as PDF
			pdf(file="condition1_condition2_scatter_insertion.pdf",height=5,width=5,onefile=FALSE)
				p<-ggplot(condition1_condition2, aes(x=condition1_ins, y=condition2_ins)) +
		 			geom_point(size=1, color="grey")+
		 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
		 			geom_point(data=subset(condition1_condition2, diff_ins>0.1 & ModStatus == "Y"), size=2, color="red")+
		  			geom_text_repel(data=subset(condition1_condition2, diff_ins>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
		  			ggtitle("Conditions 1 vs 2")+
		  			xlab("Insertion Frequency (Condition 1)")+
		  			ylab("Insertion Frequency (Condition 2)") +
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
		#For Condition 1 and Condition 2
			#Export as PDF
			pdf(file="condition1_condition2_scatter_deletion.pdf",height=5,width=5,onefile=FALSE)
				p<-ggplot(condition1_condition2, aes(x=condition1_del, y=condition2_del)) +
		 			geom_point(size=1, color="grey")+
		 			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
		 			geom_point(data=subset(condition1_condition2, diff_del>0.1 & ModStatus == "Y"), size=2, color="red")+
		  			geom_text_repel(data=subset(condition1_condition2, diff_del>0.1 & ModStatus == "Y" ), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
		  			ggtitle("Conditions 1 vs 2")+
		  			xlab("Deletion Frequency (Condition 1)")+
		  			ylab("Deletion Frequency (Condition 2)") +
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
