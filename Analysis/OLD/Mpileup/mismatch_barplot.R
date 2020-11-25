## Loading libraries
library(ggplot2)
library(reshape2)
library(plyr)

#Input
args <- commandArgs(trailingOnly = TRUE)
mismatches_input <- args[1] #1st variable
status_input <- args[2]



mod_status <- read.delim(status_input)
colnames(mod_status)<- c("chr", "pos" , "ModStatus", "Status")
yeast<- read.delim(mismatches_input, sep="\t")
yeast2<-merge(mod_status,yeast, by=c("chr","pos"))
yeast3<- na.omit(yeast2)
yeast4 <- subset(yeast3, ModStatus != "Unm")

clean_mis_input_agg<-function(data) {
data<- data[,c(3,5,6,7,8,9)]
data$sum<- data$A+data$T+data$C+data$G
data$A<- data$A/data$sum
data$T<- data$T/data$sum
data$C<- data$C/data$sum
data$G<- data$G/data$sum
data2<- aggregate(data[3:7],by=list(data$ModStatus,data$ref_nuc),FUN=mean, na.rm=TRUE)
data3<- data2[,c("Group.1", "Group.2","A", "T", "C", "G")]
data4<- melt(data3)
return(data4)
}

data<-clean_mis_input_agg(yeast4)


  ref_base_data<- vector()
  for (mod in unique(data$Group.1)){
      subs<- subset(data, Group.1==mod)
      subs$variable<- gsub(as.character(unique(subs$Group.2)),"Ref", subs$variable)
      ref_base_data<- rbind(ref_base_data, subs)
    }


mismatch_profile_barplot<- function(data,label) {
  data$variable <- factor(data$variable, levels = c("A", "T", "C", "G", "Ref"))
pdf(file= paste(label,"_mismatch_profile_barplot.pdf",sep=""),height=2,width=6, onefile=FALSE)
  print(ggplot(data, aes(x=Group.1, y=value, fill=variable), width=1) +
      scale_fill_manual(values=c("#1fab89","#eb4d55","#1e56a0", "#f0cf85","#888888"))+
        geom_bar(stat='identity', colour="black")+
        theme(axis.text.x = element_text(angle = 45,hjust = 1)))
        #facet_wrap(~unique,nrow=length(unique(data$unique)))+
        #theme_bw()
  dev.off()
}


mismatch_profile_barplot(ref_base_data, "yeastmods")



mismatch_profile_barplot_facet<- function(data,label) {
  data$variable <- factor(data$variable, levels = c("A", "T", "C", "G", "Ref"))
  pdf(file= paste(label,"_mismatch_profile_barplot_facet.pdf",sep=""),height=8,width=6, onefile=FALSE)
  print(ggplot(data, aes(x=Group.1, y=value, fill=variable), width=1) +
      scale_fill_manual(values=c("#1fab89","#eb4d55","#1e56a0", "#f0cf85","#888888"))+
        geom_bar(stat='identity', colour="black")+
        theme(axis.text.x = element_text(angle = 45,hjust = 1))+
        facet_wrap(~Group.1, scales="free"))
        #theme_bw()
  dev.off()
}


mismatch_profile_barplot_facet(ref_base_data, "yeastmods")

