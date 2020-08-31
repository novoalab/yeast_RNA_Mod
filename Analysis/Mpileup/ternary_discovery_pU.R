#Ternary plot using mpileuptostats script to discover novel sites showing pU signature
#Rscript ternary_discovery_pU.R wt_allrRNA.mismatches wt_allrRNA.STATS all_rrna_mod_status.tsv

## Loading libraries
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggtern)
library(lattice)
library(latticeExtra)
library(plyr)
library(ggrepel)

#Input
args <- commandArgs(trailingOnly = TRUE)
stats_input <- args[1] #1st variable
mismatches_input <- args[2] #1st variable
status_input <- args[3]



mod_status <- read.delim(status_input)
colnames(mod_status)<- c("chr", "pos" , "ModStatus", "Status")
mod_status$chr_pos<- paste(mod_status$chr, mod_status$pos)
yeast_S <- read.delim(stats_input, sep="\t")
yeast_S$chr_pos<- paste(yeast_S$chr, yeast_S$pos)
yeast_M <- read.delim(mismatches_input, sep="\t")
yeast_M$chr_pos<- paste(yeast_M$chr, yeast_M$pos)



yeast<-join(yeast_M,mod_status, by="chr_pos")
yeast2<-join(yeast,yeast_S, by="chr_pos")
yeast3<-subset(yeast2, coverage>30)

final<- vector()
for (chro in unique(yeast3$chr)) {
	subs<- subset(yeast3, chr==chro)
	subs2<- subset(subs, pos >30)
	subs3<- subset(subs2, pos <(max(subs$pos)-30))
	final<- rbind(final, subs3)
}


unm<- subset(final, Status=="Unm")
psu<- subset(final, Status=="Y")

# We create a function for computing the total number of mismatches
l <- list(unm, psu)
mismatches_function <- function(k){ 
  mismatches <- c()  
  for (i in c(1:nrow(k))){   
    base <- k[i,3]
    a <- sum(k[i,4:7])-k[i,toString(base)]
    mismatches <- c(mismatches, a)}
  k <- cbind(k, mismatches)}
l <- lapply(l, mismatches_function)

UNM <- l[[1]]
PSU <- l[[2]]

# We compute the mismatch frequency
UNM$mismatches_freq <- UNM$mismatches/UNM$coverage
PSU$mismatches_freq <- PSU$mismatches/PSU$coverage

# We add a variable to identify the modification
UNM$modification <- "UNM"
PSU$modification <- "Y"


U_UNM <- subset(UNM, ref_nuc=="T")
U_PSU <- subset(PSU, ref_nuc=="T")


## Plot for different chromosomes (rRNAs)
  for (ref in unique(U_UNM$chr)){
  	subset_mt_unm_u <- subset(U_UNM, chr==ref)
    pdf(file= paste(ref,"u.ternary_colored_by_mismatchfreq.pdf", sep="_"),height=5,width=5,onefile=FALSE)
      print(ggtern(subset_mt_unm_u, aes(A,G,C)) +
      geom_mask()+
      geom_point(size = 2,aes(colour=mismatches_freq)) +
      scale_colour_gradientn(limits = c(0,1), colours=c("blue2", "red"))+
      labs(color = "Mismatch Frequency") +
      ggtitle(paste("Ternary Diagram for T positions on", ref, sep=" ")) + 
      theme_bw() +
      theme(plot.title = element_text(size =10)))
   dev.off()
 }
