#####################################################################
######## SCRIPTS FOR KNOWN ncRNA SITES  ##########
#####################################################################
###################### BY OGUZHAN BEGIK #############################
#####################################################################
#Rscript ncRNA_stress_dotplot.R normal1 stress1 normal2 stress2
## Loading libraries
library(stringr)
library(ggtern)
library(ggplot2)
library(reshape2)
library(dplyr)
library(plyr)





# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
input1 <- read.delim(args[1],sep=",")  #1st variable
input2 <- read.delim(args[2],sep=",")  #2nd variable
input3 <- read.delim(args[3],sep=",")  #3rd variable
input4 <- read.delim(args[4],sep=",")  #4th variable





############################################
### PART1####
############################################
mod_ncRNA <- read.delim("yeast_ncRNA_Y_Positions.tsv")
### Importing the data

normal_rep1_input <- input1
stress_rep1_input <- input2

normal_rep2_input <- input3
stress_rep2_input <- input4


############################################
### PART1####
############################################
#PROCESS ncRNA DATA


data_manipulation_formismatch_ncrna <- function(data,label) {
	#Only Curlcake 3
	#Coverage filter
	#Position filter
	chrom <- str_split_fixed(data$X.Ref, n=2, pattern="_")
	data$RNA <- chrom[,2]
	data <- subset(data, cov > 30)
	data <- subset(data, pos > 14)
	data_s <- subset(data , RNA =="snRNA" | RNA == "snoRNA") 
	data_s$chr_pos <- paste(data_s$X.Ref, data_s$pos, sep="_")
	data_mod <- merge(data_s, mod_ncRNA, by.x=c("X.Ref", "pos"), by.y=c("chr", "gene_pos"))
	data_mod2 <- data_mod[,c("X.Ref", "pos", "chr_pos.x", "RNA", "base","Reference", "Heat", "Enzyme", "cov","mis")]
	colnames(data_mod2) <- c("X.Ref", "pos", "chr_pos", "RNA", "base", "Reference", "Heat", "Enzyme" , paste(label, "cov", sep="_"), paste(label , "mis", sep="_"))
	return(data_mod2)
}

normal_rep1 <- data_manipulation_formismatch_ncrna(normal_rep1_input, "Normal_Rep1")
stress_rep1 <- data_manipulation_formismatch_ncrna(stress_rep1_input, "Stress_Rep1")

normal_rep2 <- data_manipulation_formismatch_ncrna(normal_rep2_input, "Normal_Rep2")
stress_rep2 <- data_manipulation_formismatch_ncrna(stress_rep2_input, "Stress_Rep2")



### MERGE WITH NORMAL
Temp_REP1_ncRNA1 <-  merge(normal_rep1,stress_rep1, by.x=c("X.Ref", "pos","chr_pos","RNA", "base","Reference","Heat","Enzyme"), by.y=c("X.Ref", "pos","chr_pos", "RNA","base","Reference","Heat","Enzyme"))
Temp_REP1_ncRNA1$Stress_Score_REP1 <- Temp_REP1_ncRNA1$Stress_Rep1_mis - Temp_REP1_ncRNA1$Normal_Rep1_mis


Temp_REP2_ncRNA1 <-  merge(normal_rep2,stress_rep2, by.x=c("X.Ref", "pos","chr_pos","RNA", "base","Reference","Heat","Enzyme"), by.y=c("X.Ref", "pos","chr_pos", "RNA","base","Reference","Heat","Enzyme"))
Temp_REP2_ncRNA1$Stress_Score_REP2 <- Temp_REP2_ncRNA1$Stress_Rep2_mis - Temp_REP2_ncRNA1$Normal_Rep2_mis


both_reps <- merge(Temp_REP1_ncRNA1, Temp_REP2_ncRNA1, by.x=c("X.Ref", "pos","chr_pos","RNA", "base","Reference","Heat","Enzyme"), by.y=c("X.Ref", "pos","chr_pos","RNA", "base","Reference","Heat","Enzyme"))

######################################
#Plot for Heat Responsive Sites
data_to_plot <- both_reps[,c("X.Ref","pos","chr_pos", "RNA","Heat","Stress_Score_REP1", "Stress_Score_REP2")]
data_to_plot$Stress_Score_Mean <- rowMeans(data_to_plot[,c("Stress_Score_REP1","Stress_Score_REP2")])
data_to_plot <- data_to_plot[!is.na(data_to_plot$Stress_Score_Mean),]

data_to_plot2 <- melt(data_to_plot, id.vars=c("X.Ref","pos","chr_pos", "RNA","Heat","Stress_Score_REP1", "Stress_Score_REP2"))


data_to_plot2_HEAT <- subset(data_to_plot2, Heat=="Yes")
data_to_plot2_HEAT$label <- paste(data_to_plot2_HEAT$X.Ref,data_to_plot2_HEAT$pos)
data_to_plot2_HEAT$label <- gsub("_snoRNA", "", data_to_plot2_HEAT$label  )
data_to_plot2_HEAT$label <- gsub("_snRNA", "", data_to_plot2_HEAT$label  )
data_to_plot2_HEAT$label <- factor(data_to_plot2_HEAT$label, levels = unique(data_to_plot2_HEAT$label[order(-data_to_plot2_HEAT$value )]))


pdf(file= "Yeast_nRNA_Stress_Response_Dotplot_ReportedHeatSites.pdf",height=4,width=8,onefile=FALSE)
print(ggplot(data_to_plot2_HEAT, aes(x=label, y=value, fill=variable)) +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.1), dotsize= 2)+
  theme(axis.text.x = element_text(face="bold", color="black",size=11, angle = 45, hjust=1),
			axis.text.y = element_text(face="bold", color="black", size=11),
			plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
			axis.title.x = element_text(color="black", size=15, face="bold"),
			axis.title.y = element_text(color="black", size=15, face="bold"),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black", size=0.5),
			legend.title = element_text(color = "black", size = 15,face="bold"),
			legend.text = element_text(color = "black", size=15),
			panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()


data_to_plot2_NO <- subset(data_to_plot2, Heat=="No")
data_to_plot2_NO$label <- paste(data_to_plot2_NO$X.Ref,data_to_plot2_NO$pos)
data_to_plot2_NO$label <- gsub("_snoRNA", "", data_to_plot2_NO$label  )
data_to_plot2_NO$label <- gsub("_snRNA", "", data_to_plot2_NO$label  )

data_to_plot2_NO$Stress_Score_Mean2 <- rowMeans(data_to_plot2_NO[,c("Stress_Score_REP1", "Stress_Score_REP2")])
data_to_plot2_NO$label <- factor(data_to_plot2_NO$label, levels = unique(data_to_plot2_NO$label[order(-data_to_plot2_NO$Stress_Score_Mean2 )]))


pdf(file= "Yeast_nRNA_Stress_Response_Dotplot_NonHeatResponsiveSites.pdf",height=4,width=15,onefile=FALSE)
print(ggplot(data_to_plot2_NO, aes(x=label, y=value, fill=variable)) +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.1), dotsize= 2)+
  theme(axis.text.x = element_text(face="bold", color="black",size=11, angle = 45, hjust=1),
			axis.text.y = element_text(face="bold", color="black", size=11),
			plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
			axis.title.x = element_text(color="black", size=15, face="bold"),
			axis.title.y = element_text(color="black", size=15, face="bold"),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black", size=0.5),
			legend.title = element_text(color = "black", size = 15,face="bold"),
			legend.text = element_text(color = "black", size=15),
			panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()
