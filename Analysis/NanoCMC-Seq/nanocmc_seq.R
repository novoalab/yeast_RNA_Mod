#########################################################################
## Processing the STATS output for NanoCMC-Seq
## April 2020, Oguzhan Begik
##########################################################################
#Rscript nanocmc_seq.R cmc_positions.tsv nocmc.STATS cmc.STATS

# 1. Read data and clean up (orignal vs custom)
####
####################################

library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(grid)
library(gridExtra)
library(ggExtra)
library(ggsci)
library(ggrepel)

#Input
args <- commandArgs(trailingOnly = TRUE)
cmc_positions <- args[1]#1st variable
nocmc <- args[2]  #2nd variable
cmc <- args[3] #3rd variable



mod<- read.delim(cmc_positions)
mod$position <- paste(mod$Chr, mod$Position, sep="_") 


no_cmc_input<- read.delim(nocmc)
cmc_input<- read.delim(cmc)




clean_input<-function(data,label) {
data$pos <- data$pos-13	
data<- subset(data, pos > 40)
data$Enzyme<- rep(label, nrow(data))
data$position<- paste(data$chr, data$pos, sep="_")
	merged<-vector() #To normalize coverage by chromosome
	for (ref in (unique(data$chr))){
		subs<- subset(data, chr==ref)
		subs$norm_cov<- subs$coverage/max(subs$coverage)
		merged<- rbind(merged, subs)
	}
merged$norm_rt<- merged$rtstop/merged$coverage #to normalize RTdrop by the coverage at that position
joined<- join(merged,mod, by="position")
return(joined)
}


no_cmc<-clean_input(no_cmc_input,"No CMC")
cmc<-clean_input(cmc_input,"CMC")

all<- rbind(no_cmc,cmc)

palette_all <- c("#e84a5f","#3f72af")
palette_all_delta <- c("#3f72af")


# 2. Plotting  RT Drop
########################################
plot_area_rtdrop_facet<- function(merged, name,palette) {
	for (chro in unique(merged$chr)){
	subs<- subset(merged, chr==chro)
	pdf(file= paste(chro,name,"_rtdrop_facet.pdf",sep=""),height=15,width=15,onefile=FALSE)
	print(ggplot(subs, aes(x=pos, y=norm_rt,colour=Enzyme)) +
    geom_line()
	+scale_fill_manual(values=palette)
	+scale_colour_manual(values=palette)
  	+geom_area(aes(fill = Enzyme, group = Enzyme),alpha = 1/10, position = 'identity')
  	+geom_text_repel(data=subset(subs, norm_rt>0.2),aes(label=pos),colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)
	+theme_classic()
	+geom_vline(data=subset(subs, ModStatus=="Y"), aes(xintercept=pos),linetype="dashed")
	+facet_wrap(~Enzyme, nrow=length(unique(subs$Enzyme))))
	dev.off()
	}
}


plot_area_rtdrop_facet(all,"All",palette_all)






## Delta RT-Drop

delta_rt <- function(data_nocmc, data_cmc,label) {
	data_nocmc2<- data_nocmc[,c("chr", "pos", "position", "norm_rt", "ModStatus")]
	data_cmc2<- data_cmc[,c("position", "norm_rt")]
	joined<- join(data_nocmc2,data_cmc2, by="position")
	joined$delta<- abs(joined[,4]- joined[,6])
	joined$Sample<- label
	return(joined)
}
 

delta <- delta_rt(no_cmc,cmc,"Delta_CMC")



plot_area_rtdrop_delta<- function(merged, name,palette) {
	for (chro in unique(merged$chr)){
	subs<- subset(merged, chr==chro)
	pdf(file= paste(chro,name,"_rtdrop.pdf",sep=""),height=2,width=15,onefile=FALSE)
	print(ggplot(subs, aes(x=pos, y=delta,colour=Sample)) +
    geom_line()
	+scale_fill_manual(values=palette)
	+scale_colour_manual(values=palette)
	+geom_text_repel(data=subset(subs, delta>0.3),aes(label=pos),colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)
  	+geom_area(aes(fill = Sample, group = Sample),alpha = 1/10, position = 'identity')
	+theme_classic()
	+geom_vline(data=subset(subs, ModStatus=="Y"), aes(xintercept=pos),linetype="dashed")
	)
	dev.off()
	}
}


plot_area_rtdrop_delta(delta,"Delta",palette_all_delta)
