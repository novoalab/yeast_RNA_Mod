#We just need to ignore the 15mer sequence information because it is shifted 2 nucleotides
#Its because Nanopolish event means are based on 5 mers (base in the middle) but the position is reported for the starting base (position1)
library(stringr)
library(ggplot2)
library(reshape2)
library(data.table)



#Input
args <- commandArgs(trailingOnly = TRUE)
cond1 <- args[1] #1st variable
cond2 <- args[2] #2nd variable
cond3 <- args[3] #2nd variable
cond4 <- args[4] #2nd variable


condition1<- read.delim(cond1)
condition1_2<- condition1[!duplicated(condition1[,1:2]),]
ref <- str_split_fixed(condition1_2$windown, ":" , 15)
ref_2<- as.numeric(ref[,11])
condition1_2$windown<- paste(ref_2)
condition1_2$base<- substring(condition1_2$kmer, 10,10)
event_level_mean<- str_split_fixed(condition1_2$mean_current, ":" ,15)
condition1_3<- cbind(condition1_2, event_level_mean)
condition1_3$mean_current<-NULL
condition1_3$read<- NULL
condition1_3[,-c(1,2,3,4)] <- data.frame(sapply(condition1_3[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
colnames(condition1_3)<- c("ref", "windown", "kmer", "base", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1","2","3","4","5","6" , "7" )

condition1_final<- melt(condition1_3)
condition1_final$reference<-paste(condition1_final$ref, condition1_final$windown)
condition1_final$Strain<- rep("condition1", nrow(condition1_final))


condition2<- read.delim(cond2)
condition2_2<- condition2[!duplicated(condition2[,1:2]),]
ref <- str_split_fixed(condition2_2$windown, ":" , 15)
ref_2<- as.numeric(ref[,11])
condition2_2$windown<- paste(ref_2)
condition2_2$base<- substring(condition2_2$kmer, 10,10)
event_level_mean<- str_split_fixed(condition2_2$mean_current, ":" ,15)
condition2_3<- cbind(condition2_2, event_level_mean)
condition2_3$mean_current<-NULL
condition2_3$read<- NULL
condition2_3[,-c(1,2,3,4)] <- data.frame(sapply(condition2_3[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
colnames(condition2_3)<- c("ref", "windown", "kmer", "base", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1","2","3","4","5","6" , "7" )
condition2_final<- melt(condition2_3)
condition2_final$reference<-paste(condition2_final$ref, condition2_final$windown)
condition2_final$Strain<- rep("condition2", nrow(condition2_final))



condition3<- read.delim(cond3)
condition3_2<- condition3[!duplicated(condition3[,1:2]),]
ref <- str_split_fixed(condition3_2$windown, ":" , 15)
ref_2<- as.numeric(ref[,11])
condition3_2$windown<- paste(ref_2)
condition3_2$base<- substring(condition3_2$kmer, 10,10)
event_level_mean<- str_split_fixed(condition3_2$mean_current, ":" ,15)
condition3_3<- cbind(condition3_2, event_level_mean)
condition3_3$mean_current<-NULL
condition3_3$read<- NULL
condition3_3[,-c(1,2,3,4)] <- data.frame(sapply(condition3_3[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
colnames(condition3_3)<- c("ref", "windown", "kmer", "base", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1","2","3","4","5","6" , "7" )
condition3_final<- melt(condition3_3)
condition3_final$reference<-paste(condition3_final$ref, condition3_final$windown)
condition3_final$Strain<- rep("condition3", nrow(condition3_final))



condition4<- read.delim(cond4)
condition4_2<- condition4[!duplicated(condition4[,1:2]),]
ref <- str_split_fixed(condition4_2$windown, ":" , 15)
ref_2<- as.numeric(ref[,11])
condition4_2$windown<- paste(ref_2)
condition4_2$base<- substring(condition4_2$kmer, 10,10)
event_level_mean<- str_split_fixed(condition4_2$mean_current, ":" ,15)
condition4_3<- cbind(condition4_2, event_level_mean)
condition4_3$mean_current<-NULL
condition4_3$read<- NULL
condition4_3[,-c(1,2,3,4)] <- data.frame(sapply(condition4_3[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
colnames(condition4_3)<- c("ref", "windown", "kmer", "base", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1","2","3","4","5","6" , "7" )
condition4_final<- melt(condition4_3)
condition4_final$reference<-paste(condition4_final$ref, condition4_final$windown)
condition4_final$Strain<- rep("condition4", nrow(condition4_final))


all_final<- rbind(condition1_final, condition2_final, condition3_final,condition4_final)


for (i in seq_along(unique(all_final$reference))) { 
  subs<- subset(all_final, all_final$reference == unique(all_final$reference)[i])
  pdf(file=paste(unique(all_final$reference)[i],"15nt_window.pdf",sep="."),height=10,width=25,onefile=FALSE)
  print(ggplot(subs, aes(x=variable, y=value, group=Strain)) +
    geom_line(aes(color=Strain), size=2)+
    geom_point(aes(color=Strain))+
    ggtitle(paste(unique(all_final$reference)[i]))+
    xlab("Relative position") +
    ylab("Mean Event Level")+
    theme(axis.text.x = element_text(face="bold", color="black",size=40),
      axis.text.y = element_text(face="bold", color="black", size=40),
      plot.title = element_text(color="black", size=40, face="bold.italic",hjust = 0.5),
      axis.title.x = element_text(color="black", size=40, face="bold"),
      axis.title.y = element_text(color="black", size=40, face="bold"),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size=1),
      panel.grid.major = element_line(colour = "black", size=0.01, linetype="dashed"),
      legend.key.size = unit(1.5, "cm"),
      legend.title = element_text(color = "black", size = 30,face="bold"),
      legend.text = element_text(color = "black", size=25)))
    dev.off()
}



for (i in seq_along(unique(all_final$reference))) { 
  subs<- subset(all_final, all_final$reference == unique(all_final$reference)[i])
  png(file=paste(unique(all_final$reference)[i],"15nt_window.png",sep="."),height=800,width=2400)
  print(ggplot(subs, aes(x=variable, y=value, group=Strain)) +
    geom_line(aes(color=Strain), size=2)+
    geom_point(aes(color=Strain))+
    ggtitle(paste(unique(all_final$reference)[i]))+
    xlab("Relative position") +
    ylab("Mean Event Level")+
    theme(axis.text.x = element_text(face="bold", color="black",size=40),
          axis.text.y = element_text(face="bold", color="black", size=40),
      plot.title = element_text(color="black", size=40, face="bold.italic",hjust = 0.5),
      axis.title.x = element_text(color="black", size=40, face="bold"),
      axis.title.y = element_text(color="black", size=40, face="bold"),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size=1),
      panel.grid.major = element_line(colour = "black", size=0.1),
      legend.key.size = unit(1.5, "cm"),
      legend.title = element_text(color = "black", size = 30,face="bold"),
            legend.text = element_text(color = "black", size=25)))
    dev.off()
}