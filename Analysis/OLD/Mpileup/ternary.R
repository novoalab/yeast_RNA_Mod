#Ternary plot using mpileuptostats script
#Rscript ternary.R wt.mismatches wt.STATS all_rrna_mod_status.tsv

## Loading libraries
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggtern)
library(lattice)
library(latticeExtra)
library(plyr)

#Input
args <- commandArgs(trailingOnly = TRUE)
stats_input <- args[1] #1st variable
mismatches_input <- args[2] #1st variable
status_input <- args[3]



#Importing and manipulating the rna mod 5mer file
  #Read the table for modification positions  
  mod_status <- read.delim(status_input)
  colnames(mod_status)<- c("chr", "pos" , "ModStatus", "Status")
  yeast_S <- read.delim(stats_input, sep="\t")
  yeast_M <- read.delim(mismatches_input, sep="\t")

  yeast<-merge(mod_status,yeast_M, by=c("chr","pos"))
  yeast2<-merge(yeast,yeast_S, by=c("chr","pos"))
  yeast3<-subset(yeast2, coverage>50) 
  yeast4<- na.omit(yeast3)


  unm<- subset(yeast4, ModStatus=="Unm")
  psu<- subset(yeast4, ModStatus=="Y")
  Am<- subset(yeast4, ModStatus=="Am")
  Gm<- subset(yeast4, ModStatus=="Gm")
  Um<- subset(yeast4, ModStatus=="Um")
  Cm<- subset(yeast4, ModStatus=="Cm")

# We create a function for computing the total number of mismatches
l <- list(unm, psu, Am, Gm, Um, Cm)
mismatches_function <- function(k){ 
  mismatches <- c()  
  for (i in c(1:nrow(k))){   
    base <- k[i,5]
    a <- sum(k[i,6:9])-k[i,toString(base)]
    mismatches <- c(mismatches, a)}
  k <- cbind(k, mismatches)}
l <- lapply(l, mismatches_function)

UNM <- l[[1]]
PSU <- l[[2]]
AM <- l[[3]]
GM <- l[[4]]
UM <- l[[5]]
CM <- l[[6]]

# We compute the mismatch frequency
UNM$mismatches_freq <- UNM$mismatches/UNM$coverage
PSU$mismatches_freq <- PSU$mismatches/PSU$coverage
AM$mismatches_freq <- AM$mismatches/AM$coverage
GM$mismatches_freq <- GM$mismatches/GM$coverage
UM$mismatches_freq <- UM$mismatches/UM$coverage
CM$mismatches_freq <- CM$mismatches/CM$coverage

# We add a variable to identify the modification
UNM$modification <- "UNM"
PSU$modification <- "Y"
AM$modification <- "Am"
GM$modification <- "Gm"
UM$modification <- "Um"
CM$modification <- "Cm"

U_PSU <- subset(PSU, ref_nuc.x=="T")
A_AM <- subset(AM, ref_nuc.x=="A")
G_GM <- subset(GM, ref_nuc.x=="G")
U_UM <- subset(UM, ref_nuc.x=="T")
C_CM <- subset(CM, ref_nuc.x=="C")


  pdf(file= "psu.ternary.pdf",height=2,width=2,onefile=FALSE)
    print(ggtern(U_PSU, aes(A,G,C)) +
    geom_point(size = 1,aes(colour="red")) +
    geom_mask()+
    ggtitle("Ternary Diagram for Pseudouridylation") + theme_bw(base_size = 6) +
    theme(plot.title = element_text(size =5),legend.position = "none"))
  dev.off()


  pdf(file= "am.ternary.pdf",height=2,width=2,onefile=FALSE)
    print(ggtern(A_AM, aes(T,G,C)) +
    geom_point(size = 1,aes(colour="red")) +
    geom_mask()+
    ggtitle("Ternary Diagram for Am") + theme_bw(base_size = 6) +
    theme(plot.title = element_text(size = 5),legend.position = "none"))
  dev.off()



  pdf(file= "gm.ternary.pdf",height=2,width=2,onefile=FALSE)
    print(ggtern(G_GM, aes(A,T,C)) +
    geom_point(size = 1,aes(colour="red")) +
    geom_mask()+
    ggtitle("Ternary Diagram for Gm") + theme_bw(base_size = 6) +
    theme(plot.title = element_text(size = 5),legend.position = "none"))
  dev.off()


  pdf(file= "um.ternary.pdf",height=2,width=2,onefile=FALSE)
    print(ggtern(U_UM, aes(A,G,C)) +
    geom_point(size = 1,aes(colour="red")) +
    geom_mask()+
    ggtitle("Ternary Diagram for Um") + theme_bw(base_size = 6) +
    theme(plot.title = element_text(size = 5),legend.position = "none"))
  dev.off()


  pdf(file= "cm.ternary.pdf",height=2,width=2,onefile=FALSE)
    print(ggtern(C_CM, aes(A,G,T)) +
    geom_point(size = 1,aes(colour="red")) +
    geom_mask()+
    ggtitle("Ternary Diagram for Cm") + theme_bw(base_size = 6) +
    theme(plot.title = element_text(size = 5),legend.position = "none"))
  dev.off()
