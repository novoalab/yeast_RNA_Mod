#####################################################################
######## SCRIPTS FOR THE TERNARY PLOTS OF RRNA MODIFICATION #########
#####################################################################
###################### BY OGUZHAN BEGIK #############################
#####################################################################
# Rscript mods_ternary.R wt_epinano.csv all_rrna_mod_status.tsv
## Loading libraries
library(stringr)
library(ggtern)
library(ggplot2)
library(reshape2)
library(ggbeeswarm)
library(ggpubr)
library(dplyr)


#Input
args <- commandArgs(trailingOnly = TRUE)
input <- args[1] #1st variable
status_input <- read.delim(args[2])


data <- read.delim(input,sep=",")

############################################
### PART1####
############################################
### Imporing the data



data_manipulation <- function(data,label) {
	#Coverage filter
	data<- subset(data, cov > 30)
	data$chr_pos <- paste(data$X.Ref, data$pos, sep="_")
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
		subs2<- subset(subs, pos >25)
		subs3<- subset(subs2, pos <(max(subs$pos)-25))
		final<- rbind(final, subs3)
	}
	status_input$chr_pos <- paste(status_input$Chr, status_input$Position, sep="_")
	final_mod <- merge(final, status_input , by.x=c("chr_pos"), by.y=c("chr_pos"))
	return(final_mod)
}


data_processed <- data_manipulation(data,"wt")



##TERNARY PLOT for Base A
	a_bases<- subset(data_processed, base=="A")
 for (mod in unique(a_bases$ModStatus)) {
 subs <- subset(a_bases, ModStatus==mod)
	 pdf(file= paste(mod,"A_ternary_plot.pdf",sep="_"),height=3,width=5,onefile=FALSE)
	    print(ggtern(subs, aes(G,C, T)) +
	    geom_mask()+
	    geom_point(size = 2, colour="#f67e7d") +
	    labs(color = "Mismatch Frequency") +
	    ggtitle( paste("Ternary Diagram for",mod,"A positions in  Ribosomal RNA")) + 
	    theme_bw() +
	    theme(plot.title = element_text(size =10)))
	 dev.off()
}

##TERNARY PLOT for Base U
	u_bases<- subset(data_processed, base=="T")
 for (mod in unique(u_bases$ModStatus)) {
 subs <- subset(u_bases, ModStatus==mod)
	 pdf(file= paste(mod,"U_ternary_plot.pdf",sep="_"),height=3,width=5,onefile=FALSE)
	    print(ggtern(subs, aes(A,G,C)) +
	    geom_mask()+
	    geom_point(size = 2, colour="#f67e7d") +
	    labs(color = "Mismatch Frequency") +
	    ggtitle( paste("Ternary Diagram for",mod,"U positions in  Ribosomal RNA")) + 
	    theme_bw() +
	    theme(plot.title = element_text(size =10)))
	 dev.off()
}

##TERNARY PLOT for Base G
	g_bases<- subset(data_processed, base=="G")
 for (mod in unique(g_bases$ModStatus)) {
 subs <- subset(g_bases, ModStatus==mod)
	 pdf(file= paste(mod,"G_ternary_plot.pdf",sep="_"),height=3,width=5,onefile=FALSE)
	    print(ggtern(subs, aes(A,C, T)) +
	    geom_mask()+
	    geom_point(size = 2, colour="#f67e7d") +
	    labs(color = "Mismatch Frequency") +
	    ggtitle( paste("Ternary Diagram for",mod,"G positions in  Ribosomal RNA")) + 
	    theme_bw() +
	    theme(plot.title = element_text(size =10)))
	 dev.off()
}


##TERNARY PLOT for Base C
	c_bases<- subset(data_processed, base=="C")
 for (mod in unique(c_bases$ModStatus)) {
 subs <- subset(c_bases, ModStatus==mod)
	 pdf(file= paste(mod,"C_ternary_plot.pdf",sep="_"),height=3,width=5,onefile=FALSE)
	    print(ggtern(subs, aes(A,G, T)) +
	    geom_mask()+
	    geom_point(size = 2, colour="#f67e7d") +
	    labs(color = "Mismatch Frequency") +
	    ggtitle( paste("Ternary Diagram for",mod,"C positions in  Ribosomal RNA")) + 
	    theme_bw() +
	    theme(plot.title = element_text(size =10)))
	 dev.off()
}

