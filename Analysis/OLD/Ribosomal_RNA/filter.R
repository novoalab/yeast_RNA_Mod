#####################################################################
######## SCRIPTS FOR RANDOM FOREST ON RIBOSOMAL RNAS.      ##########
#####################################################################
###################### BY OGUZHAN BEGIK #############################
#####################################################################
## Loading libraries
library(stringr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(plyr)
library(data.table)
library(cutpointr)


############################################
### PART1####
############################################
### Mod File
mod_rRNA <- read.delim("yeast_all_rrna_mod_status.tsv")
mod_ncRNA <- read.delim("supp/yeast_ncRNA_Y_Positions.tsv")
### Importing the data

rRNA_REP1 <- read.delim("RNA814001_REP1_Yeast_WT_Normal_all_rRNA_gdef.plus_strand.per.site.var.csv" ,sep=",")
rRNA_REP2 <- read.delim("RNA814001_REP1_Yeast_Sn3KO_Normal_all_rRNA_gdef.plus_strand.per.site.var.csv" ,sep=",")
rRNA_REP3 <- read.delim("RNA814001_REP1_Yeast_Sn36KO_Normal_all_rRNA_gdef.plus_strand.per.site.var.csv" ,sep=",")




############################################
### PART2####
############################################
#PROCESS rRNA DATA
data_manipulation_forbases <- function(data,label) {
	#Only Curlcake 3
	#Coverage filter
	data<- subset(data, cov > 30)
	data$chr_pos <- paste(data$X.Ref, data$pos, sep="_")
	data <- subset(data, chr_pos != "18s_1187" )
	data <- subset(data, chr_pos != "25s_2129")
	data <- subset(data, chr_pos != "25s_2133" )
	data <- subset(data, chr_pos != "25s_2264")
	#data <- subset(data, chr_pos != "18s_120")
	#data <- subset(data, chr_pos != "18s_302")
	data <- subset(data, chr_pos != "25s_2880")
	bases <- str_split_fixed(data$ACGT_freq, n=4, pattern=":")
	colnames(bases) <- c("A", "C", "G", "T")
	data2 <- cbind(data, bases)
	data3 <- data2[,c("X.Ref", "pos", "chr_pos", "base", "cov", "q_mean", "q_median", "q_std", "ins", "del", "mis", "A", "T", "C", "G")]
	data3$A <- as.numeric(as.character(data3$A))
	data3$T <- as.numeric(as.character(data3$T))
	data3$G <- as.numeric(as.character(data3$G))
	data3$C <- as.numeric(as.character(data3$C))
	#Calculate mismatches
	mismatches <- c()  
	for (i in c(1:nrow(data3))){   
		base <- data3[i,c("base")]
		a <- sum(data3[i,c("A","T", "C", "G")])-data3[i,toString(base)]
		mismatches <- c(mismatches, a)
	}
	data3 <- cbind(data3, mismatches)
 	data3$count <- data3$A+ data3$C + data3$G + data3$T
	data3$mis_freq <- data3$mismatches/data3$count
	#Position filter
	final<- vector()
	for (chro in unique(data3$X.Ref)) {
		subs<- subset(data3, X.Ref==chro)
		subs2<- subset(subs, pos >30)
		subs3<- subset(subs2, pos <(max(subs$pos)-30))
		final<- rbind(final, subs3)
	}
	final_U <- subset(final, base=="T")
	final_U$A_freq <- final_U$A / final_U$mismatches
	final_U$C_freq <- final_U$C / final_U$mismatches
	final_U$G_freq <- final_U$G / final_U$mismatches
	final_U <- final_U %>% mutate(A_freq = coalesce(A_freq, 0))
	final_U <- final_U %>% mutate(C_freq = coalesce(C_freq, 0))
	final_U <- final_U %>% mutate(G_freq = coalesce(G_freq, 0))
	final_U$sample <- label
	final_U2 <- merge(final_U, mod_rRNA, by.x=c("X.Ref", "pos"), by.y=c("Chr", "Position"))
	final_U3 <- final_U2[,c("X.Ref", "pos", "chr_pos", "base","sample", "ModStatus", "Status", "cov", "q_mean", "q_median", "q_std", "ins", "del", "mis", "A", "T", "C", "G","mis_freq","A_freq","C_freq","G_freq")]
	return(final_U3)
}



rRNA_Rep1_bases <- data_manipulation_forbases(rRNA_REP1, "rRNA_Rep1")

rRNA_Rep2_bases <- data_manipulation_forbases(rRNA_REP2, "rRNA_Rep2")

rRNA_Rep3_bases <- data_manipulation_forbases(rRNA_REP3, "rRNA_Rep3")


#Merge all the reps to gret higher statistical power
rRNA_All <- rbind(rRNA_Rep1_bases,rRNA_Rep2_bases,rRNA_Rep3_bases)



############################################
### PART 3#### CutpointR for filtering
 ############################################
### Using cutpoint in order to find best cutpoints 

#Extract ONLY Unmodified/Unaffected and Y positions
rRNA_All_UNM_Y<- subset(rRNA_All, Status=="Y" | Status=="Unm")

#Remove other levels
rRNA_All_UNM_Y$Status <- factor(rRNA_All_UNM_Y$Status)

#Use ONLY Cytosolic RNAs for training
rRNA_All_UNM_Y_TRAINING <- subset(rRNA_All_UNM_Y, X.Ref !="15s" & X.Ref !="21s")


opt_cut_mis <-cutpointr(rRNA_All_UNM_Y_TRAINING, mis, Status, method = oc_youden_kernel)
opt_cut_c <-cutpointr(rRNA_All_UNM_Y_TRAINING, C_freq, Status, method = oc_youden_kernel)
opt_cut_g <-cutpointr(rRNA_All_UNM_Y_TRAINING, G_freq,  pos_class = "Y", Status, method = oc_youden_kernel)
opt_cut_a <-cutpointr(rRNA_All_UNM_Y_TRAINING, A_freq,  pos_class = "Y", Status, method = oc_youden_kernel)


pdf(file= "opt_cut_mismatch.pdf",height=5,width=10,onefile=FALSE)
plot(opt_cut_mis)
#Accuracy : 0.967596
#AUC 0.988284
dev.off()
pdf(file= "opt_cut_c.pdf",height=5,width=10,onefile=FALSE)
plot(opt_cut_c)
#Accuracy : 0.875
#AUC : 0.828660
dev.off()
pdf(file= "opt_cut_g.pdf",height=5,width=10,onefile=FALSE)
plot(opt_cut_g)
dev.off()
pdf(file= "opt_cut_a.pdf",height=5,width=10,onefile=FALSE)
plot(opt_cut_a)
dev.off()




rRNA_All_UNM_Y_TESTING <- subset(rRNA_All_UNM_Y, X.Ref =="15s")
#rRNA_All_UNM_Y_TESTING_mis_filter <- subset(rRNA_All_UNM_Y_TESTING, mis >opt_cut_mis$optimal_cutpoint)
rRNA_All_UNM_Y_TESTING_all_filter <- subset(rRNA_All_UNM_Y_TESTING, mis > opt_cut_mis$optimal_cutpoint & C_freq > opt_cut_c$optimal_cutpoint & G_freq < opt_cut_g$optimal_cutpoint )
table(rRNA_All_UNM_Y_TESTING_all_filter$chr_pos)

mis_thr <- 0.1374356
c_thr <- 0.5782206
g_thr <- 0.03320507
