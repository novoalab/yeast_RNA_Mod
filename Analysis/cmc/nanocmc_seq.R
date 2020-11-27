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


#Input
no_cmc_input<- read.delim(nocmc, sep="")
cmc30_input<- read.delim(cmc, sep="")




clean_input<-function(data,label) {
data$pos <- data$pos-13	#position fixing
data<- subset(data, pos > 10) #Cut first 40 nucleotide
data$position<- paste(data$chr, data$pos, sep="_")
data$Enzyme<- rep(label, nrow(data))
	merged<-vector() #To remove the last 40 nt
	for (ref in (unique(data$chr))){
		subs<- subset(data, chr==ref)
		subs2<- subset(subs, pos< max(subs$pos)-10)
		merged<- rbind(merged, subs2)
	}
# Take the coverage from +1 position
	final <- vector()
	for (upos in unique(merged$position)) {
		subs <- subset(merged, position==upos)
		numb <- subs[,"pos"] + 1 
		plus_one <- paste(subs$chr, numb, sep="_")
		subs_2 <- subset(data, position==plus_one)
		subs$cov_1 <- subs_2[,"coverage"] 
		final<- rbind(final, subs)
	}

final$norm_rt<- final$rtstop/final$cov_1 #to normalize RTdrop by the coverage at that position
joined<- join(final,mod, by="position")
return(joined)
}



no_cmc<-clean_input(no_cmc_input,"No_Treatment")
cmc30<-clean_input(cmc30_input,"CMC")


final_30 <- rbind(no_cmc,cmc30)



## Delta RT-Drop
delta_rt <- function(data_nocmc, data_cmc,label) {
	data_nocmc2<- data_nocmc[,c("chr", "pos", "position", "norm_rt", "CMC", "Mod", "ModStatus")]
	data_cmc2<- data_cmc[,c("position", "norm_rt")]
	joined<- join(data_nocmc2,data_cmc2, by="position")
	joined$delta<- abs(joined[,4]- joined[,8])
	joined$Sample<- label
	colnames(joined) <- c("chr","pos","position","untreated_norm_rt", "CMC", "Mod","ModStatus","treated_norm_rt", "delta", "Sample")
	final <- vector()
	for (ref in unique(joined$chr)){
		subs <- subset(joined, chr==ref)
		subs$CMC_Score <- subs$delta / median(subs$delta)
		final <- rbind(final, subs)
	}
	return(final)
}
 

delta_30<- delta_rt(no_cmc,cmc30,"30CMC_Delta")
write.table(delta_30, file="CMC_Scores_rRNA.tsv", sep="\t", quote=FALSE, row.names=FALSE)




palette_15s <-c("red","gray")
palette_21s <-c("gray")
palette_25s <-c("gray", "steelblue")
palette_18s <-c("gray", "steelblue")
plot_area_rtdrop_delta_facet<- function(merged, name, ref, ylim, palette) {
	subs <- subset(merged, chr==ref)
	pdf(file= paste(ref,name,"_rtdrop_facet.pdf",sep=""),height=3,width=10,onefile=FALSE)
	print(ggplot(subs, aes(x=pos, y=CMC_Score,color=Mod, fill=Mod)) +
    geom_bar(stat="identity")
    +ggtitle(paste(ref,name, sep="_"))
    +ylim(0, as.numeric(ylim))
    +scale_fill_manual(values=palette)
    +scale_colour_manual(values=palette)
	+geom_text_repel(data=subset(subs, CMC_Score > 25 & Mod== "Maybe"),aes(label=pos),colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)
	#+geom_text_repel(data=subset(merged, CMC_Score>25 & Mod== "No"),aes(label=pos),colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)
	+theme_classic()
	+geom_hline(yintercept = 25, linetype="dashed")	
	+facet_wrap(~chr, scales="free", nrow=length(unique(subs$chr))))
	dev.off()
}

plot_area_rtdrop_delta_facet(delta_30,"Delta_Merged_30","15s", "50", palette_15s )
plot_area_rtdrop_delta_facet(delta_30,"Delta_Merged_30","21s", "50", palette_21s )
plot_area_rtdrop_delta_facet(delta_30,"Delta_Merged_30","25s", "210" , palette_25s)
plot_area_rtdrop_delta_facet(delta_30,"Delta_Merged_30","18s", "210" , palette_18s)





























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
