### Dot-Plots of modifications on 5 mer context 
  #Rscript dotplot_5mer wt.bam.tsv.per.site.var.per_site_var.5mer.csv rrna_mod_5mer.tsv


#Load the libraries needed for this script
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


#Importing and manipulating the rna mod 5mer file
	#Read the table for modification positions		
	mod_positions<- read.delim(status_input)
	#Rename the columns
	colnames(mod_positions)<- c("rRNA","Ref_Pos","Mod", "X.Kmer")
	#Create a column for unique positions (18s 145)
	mod_positions$chr_pos<- paste(mod_positions$rRNA, mod_positions$Ref_Pos)


	#For WT
		#Importing the epinano output
		wt<- read.delim(input, sep=",")
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


	#Dot Plots
	for (mod in unique(wt_mod2$Mod)) {
		subset_wt_mod2<- subset(wt_mod2, Mod==mod)
		pdf(file=paste(mod, "dotplot.5mer.pdf", sep="_"),height=3,width=16,onefile=FALSE)
			print(ggplot(subset_wt_mod2, aes(x=variable, y=value)) + 
				geom_quasirandom(varwidth = TRUE, aes(color=feature))+
				geom_boxplot()+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
  				facet_wrap(~feature,scales="free",nrow=1)+
				theme_bw()+
				xlab("5-mer positions")+
              	ylab("Features") +
				theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            		axis.title=element_text(size=17,face="bold"),
            		legend.title = element_text(size = 20),
            		legend.text = element_text(color = "black", size=15)))
		dev.off()
	}
