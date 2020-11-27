### Dot-Plots of Modified vs Unmodified bases
  # We are using a status table which contains ModStatus and Status values for each position. 
  # ModStatus only indicates whether a position is modified or not
  # Status indicates whether a position is AFFECTED by a neighboruing position (i.e. they are within the 5-mer of another modification)
  # So we will plot ALL positions and only Unaffected positions
  #Rscript dotplot.R wt.bam.tsv.per.site.var.per_site_var.5mer.csv all_rrna_mod_status.tsv
  #Libraries needed
	library(plyr)
	library(stringr)
	library(reshape2)
	library(dplyr)
	library(ggplot2)
	library(ggbeeswarm)
	library(ggpubr)

#Input
args <- commandArgs(trailingOnly = TRUE)
input <- args[1] #1st variable
status_input <- args[2]

#Importing the files
	#Status File
	status<-read.delim(status_input)
		#Add a column for unique positions
		status$chr_pos<- paste(status$Chr, status$Position)
	#Input Epinano file
	yeast<-read.delim(input, sep=",")
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

#Merge status file with processed input
yeast_mod<- join(yeast,status, by="chr_pos")
#Order by position
yeast_mod<-yeast_mod[order(yeast_mod$Position),]
#Order by Chromosome
yeast_mod<-yeast_mod[order(yeast_mod$Ref),]
#Extract the base information
yeast_mod$Nuc<- substring(yeast_mod$X.Kmer, 3,3)


#Process the data for ALL positions analysis
	#Extract mod positions which are ALL by other mods
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

	#Extract unmod positions which are ALL by other mods
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




#Process the data for UNAFFECTED positions analysis
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