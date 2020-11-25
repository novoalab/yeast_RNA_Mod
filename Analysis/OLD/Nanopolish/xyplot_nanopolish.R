  #Scatter Plots
    library(plyr)
    library(ggplot2)
    library(ggrepel)
 
   args <- commandArgs(trailingOnly = TRUE)
   cond1 <- args[1] #1st variable
   cond2 <- args[2] #2nd variable
   threshold<-args[3]

    #Import the data
      #Import the Per position table for Sn3
      condition1 <- read.delim(cond1, sep=" ")
      #Add 3 nt to each position since Nanopolish output is 0 based and falls behind 2 nts
      condition1$position<- condition1$position+3 
      #Create a column for unique positions
      condition1$pos<- paste(condition1$contig, condition1$position)
      #Add a column replicating condition1
      condition1$sample <- rep("condition1",nrow(condition1)) 
      #Rename the columns
      colnames(condition1)<- c("contig","position", "kmer","readidx", "condition1_event_level_mean_mean", "pos", "sample")

      #Import the Per position table for WT
      condition2 <- read.delim(cond2, sep=" ")
      #Add 3 nt to each position since Nanopolish output is 0 based and falls behind 2 nts
      condition2$position<- condition2$position+3 
      #Create a column for unique positions
      condition2$pos<- paste(condition2$contig, condition2$position)
      #Add a column replicating condition2
      condition2$sample <- rep("condition2",nrow(condition2)) 
      #Rename the columns
      colnames(condition2)<- c("contig","position", "kmer","readidx", "condition2_event_level_mean_mean", "pos", "sample")



      # Use join function to merge the strains
      condition2_sn3<- join(condition2, condition1, by="pos") 
      #Remove the columns that do not match
      condition2_sn3<- subset(condition2_sn3, condition1_event_level_mean_mean!="NA")
      #Remove the 20 nt from 5'UTR
      condition2_sn3<- subset(condition2_sn3, position>20)
      #Calculate the difference 
      condition2_sn3$diff<- abs(condition2_sn3$condition2_event_level_mean_mean - condition2_sn3$condition1_event_level_mean_mean)


      #Plot for sn3KO
      for (rna in unique(condition2_sn3$contig)) {
      subs<-subset(condition2_sn3, contig==rna)  
      pdf(file=paste(rna, "condition1_condition2_event_mean_xyplot.pdf", sep="_"),height=5,width=5,onefile=FALSE)
      print(ggplot(subs, aes(x=condition1_event_level_mean_mean, y=condition2_event_level_mean_mean)) +
        geom_point(size=1, color="grey")+
        geom_smooth(method='lm', formula= y~x, linetype="dashed", size=0.2, color= "black")+
        geom_point(data=subset(subs, diff>threshold) , size=2, color="red")+
        geom_text_repel(data=subset(subs, diff>threshold), aes(label=position), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
        ggtitle(paste(rna,"snR3-KO"))+
        xlab("Mean Current Intensity (Cond1)")+
        ylab("Mean Current Intensity (Cond2)") +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold", color="black",size=11),
          axis.text.y = element_text(face="bold", color="black", size=11),
          plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
          axis.title.x = element_text(color="black", size=15, face="bold"),
          axis.title.y = element_text(color="black", size=15, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title = element_text(color = "black", size = 20,face="bold"),
          legend.text = element_text(color = "black", size=20)))
      dev.off()
      }
