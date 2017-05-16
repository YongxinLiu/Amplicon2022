if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("Biobase","edgeR","ggplot2","gplots","grid","RColorBrewer","reshape2","VennDiagram"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd("result_k1-c")
library("Biobase")
library("edgeR")
library("ggplot2")
library("gplots")
library("grid")
library("RColorBrewer")
library("reshape2")
library("VennDiagram")

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(size=.5, colour="black"),
                    axis.line.y=element_line(size=.5, colour="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black", size=7),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=7),
                    text=element_text(family="sans", size=7))

# Public file 1. "design.txt"  Design of experiment
design = read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/design.txt", header=T, row.names= 1, sep="\t") 

# Public file 2. "otu_table.txt"  raw reads count of each OTU in each sample
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

# Public file 3. "rep_seqs_tax.txt"  taxonomy for each OTU, tab seperated
taxonomy = read.delim("rep_seqs_tax.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","evalue")

# Public file 4. "otu_table_css.txt", css normalization of OTU table
cssnorm = read.delim("otu_table_css.txt", row.names= 1,  header=T, sep="\t")

# setting subset design
if (TRUE){
	sub_design = subset(design,genotype %in% c("WT","DM1","DM2","DO1","DO2") ) # select group1
}else{
	sub_design = design
}

if (TRUE){
	sub_design = subset(sub_design,batch %in% c("3") ) # select group2
}

if (FALSE){
	sub_design$group=paste(sub_design$genotype,sub_design$batch,sep = ".")
}else{
	sub_design$group=sub_design$genotype
}

# Set group order
if ("TRUE" == "TRUE") {
    sub_design$group  = factor(sub_design$group, levels=c("WT","DM1","DM2","DO1","DO2"))   # set group order
}

print(paste("Number of group: ",length(unique(sub_design$group)),sep="")) # show group numbers

# sub and reorder subdesign and otu_table
idx = rownames(sub_design) %in% colnames(otu_table)
sub_design = sub_design[idx,]
count = otu_table[, rownames(sub_design)]
norm = t(t(count)/colSums(count,na=T)) * 1000 # normalization to total 1000

# Pearson correlation among samples
sim=cor(norm,method="pearson")
sim=1-sim # dissimilarity
sim=round(sim,3) # Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures),  : invalid graphics state
write.table(sim, file=paste("sim.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
sim = read.table("sim.txt", header=T, row.names= 1, sep="\t") 
sim=as.matrix(sim)
pdf(file=paste("heat_cor_samples.pdf", sep=""), height = 16, width = 16)
heatmap.2(sim, Rowv=FALSE, Colv=FALSE, dendrogram='none', trace='none', margins=c(6,6), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),density.info="none") 
dev.off()
png(file=paste("heat_cor_samples.png", sep=""), height = 16, width = 16, units = "in", res = 300)
heatmap.2(sim, Rowv=FALSE, Colv=FALSE, dendrogram='none', trace='none', margins=c(6,6), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),density.info="none") 
dev.off()


# Pearson correlation among groups
sampFile = as.data.frame(sub_design$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = norm
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

# options(digits=3) # set digits length, but heatmap.2 cellnote still long digits
sim=cor(mat_mean_final,method="pearson")
sim=1-sim # dissimilarity
sim=round(sim,3) # Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures),  : invalid graphics state
write.table(sim, file=paste("sim.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
sim = read.table("sim.txt", header=T, row.names= 1, sep="\t") 
sim=as.matrix(sim)
pdf(file=paste("heat_cor_groups.pdf", sep=""), height = 8, width = 8)
heatmap.2(sim, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=sim, notecol="black", trace='none', margins=c(6,6), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),density.info="none")
dev.off()
png(file=paste("heat_cor_groups.png", sep=""), height = 8, width = 8, units = "in", res = 300)
heatmap.2(sim, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=sim, notecol="black", trace='none', margins=c(6,6), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),density.info="none")
dev.off()

# High abundance phylum
per= read.delim ("sum_taxa/otu_table_tax_L2.txt", sep = "\t", row.names=1, header=T)
per = per[(order(-rowSums(per))), ] # decrease sort
per$mean=rowMeans(per)
top_phylum=rownames(per[per$mean>0.005,]) 
print(paste("Number of phylum > 0.005: ",length(top_phylum), sep=""))
top_phylum=gsub("[\\w;_]+p__","",top_phylum,perl=TRUE) # regexp in perl mode
top_phylum

#############################################################
# Title: Differentially abundance OTU by edgeR - volcano plot, manhattan plot, heatmap
# Author: Yong-Xin Liu
# E-mail: yxliu@genetics.ac.cn
# Date: 3/14/2017
# Version: 1.1
# Enviroment: R 3.2.1 x64, OS Win10 x64
# Description: Script to draw scatterplot or heatmap show differentially abudance OTU
#############################################################
# create DGE list
groups = sub_design$group
d = DGEList(counts=count, group=groups)
d = calcNormFactors(d)

# fit the GLM
design.mat = model.matrix(~ 0 + d$samples$group)
colnames(design.mat)=levels(groups)
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)
write.table(paste("SampAvsB\tenriched\tdepleted\tnosig",sep="\t"), file=paste("otu_sum.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

# function DA limma
DAOTU_edgeR <- function(sampleV){
	sampleA <- as.vector(sampleV$sampA)
	sampleB <- as.vector(sampleV$sampB)
#	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	print(paste("Start DA OTU", sampleA, sampleB, ":",sep=" "))
	SampAvsB=paste(sampleA,"-", sampleB, sep="")
	print(SampAvsB)
	BvsA <- makeContrasts(contrasts = SampAvsB, levels=design.mat)
	lrt = glmLRT(fit,contrast=BvsA)
	de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

	x=lrt$table
	x$sig=de_lrt
	# No significant FDR OTU, change to pvalue
	if (dim(x[x$sig==1,])[1]==0 |dim(x[x$sig==-1,])[1]==0){
		x$level = ifelse(x$PValue<0.05 & x$logFC>0, "enriched",ifelse(x$PValue<0.05 & x$logFC<0, "depleted","nosig"))
	}else{
		x$level = ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig"))
	}
	x$otu = rownames(x)
	x$neglogp = -log(x$PValue)

	# Classify OTU by FDR < 0.05
	enriched = row.names(subset(x,level=="enriched"))
	nosig = row.names(subset(x,level== "nosig"))
	depleted = row.names(subset(x,level=="depleted"))

	# Order OTUs according to taxonomy phylum
	taxonomy = taxonomy[order(taxonomy[, 2]), ]
	idx = rownames(taxonomy) %in% x$otu
	tax = taxonomy[idx, ] # subset taxonomy from used OTU

	#idx = match(rownames(tax), x$otu)
	x = x[rownames(tax), ] # reorder according to tax
	x$tax = gsub("p__","",tax$phylum,perl=TRUE) 
	x$phylum = gsub("p__","",tax$phylum,perl=TRUE) 
	x$class = gsub("c__","",tax$class,perl=TRUE)
	x$order = gsub("o__","",tax$order,perl=TRUE)
	x$family = gsub("f__","",tax$family,perl=TRUE)
	x$genus = gsub("g__","",tax$genus,perl=TRUE)
	levels(x$tax)=c(levels(x$tax),"Low Abundance")
	x[!(x$tax %in% top_phylum),]$tax = "Low Abundance" # no level can get value


	# plot
	if (max(x$neglogp)>15){
	  x[x$neglogp>15,]$neglogp  = 15
	}

	# Draw scatter plot
	if (max(x$logFC)>4){x[x$logFC>4,]$logFC = 4} # norm x axis
	if (min(x$logFC)< -4){x[x$logFC< -4,]$logFC = -4} # norm x axis
	x$otu = factor(x$otu, levels=x$otu)   # set x order
	x$level = factor(x$level, levels=c("enriched","depleted","nosig"))
	x$tax = factor(x$tax, levels=c(top_phylum,"Low Abundance"))

	# Volcanol plot of fold change vs abundance plot
#	p = ggplot(x, aes(x=logFC, y=logCPM, color=level, size=logCPM, shape=tax)) + geom_point()  + # size by abundance
	p = ggplot(x, aes(x=logFC, y=logCPM, color=level)) + geom_point()  +
	  scale_colour_manual(values=c("red","green","grey"))+ xlim(-4, 4)+
	  #  ylim(5, 20)+
	  labs(x="log2(fold change)",y="log2(count per million)", title=paste(sampleA, "vs", sampleB, sep=" "))+main_theme
	p
	ggsave(file=paste("vol_otu_", sampleA, "vs", sampleB, ".pdf", sep=""), p, width = 5, height = 3)
	ggsave(file=paste("vol_otu_", sampleA, "vs", sampleB, ".png", sep=""), p, width = 5, height = 3)

	# Manhattan plot
	FDR = min(x$neglogp[x$level=="enriched"])
	p = ggplot(x, aes(x=otu, y=neglogp, color=tax, size=logCPM, shape=level)) +
	  geom_point(alpha=.7) + 
	  geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
	  #scale_color_manual(values=colors$colors) +
	  scale_shape_manual(values=c(17, 25, 20))+
	  scale_size(breaks=c(5, 10, 15)) +
#	  scale_size(breaks=c(1, 2, 3)) +
	  labs(x="OTU", y="-loge(P)", title=paste(sampleA, "vs", sampleB, sep=" ")) +main_theme +
	  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="right")
	p
	ggsave(file=paste("man_otu_", sampleA, "vs", sampleB, ".pdf", sep=""), p, width = 5*2, height = 3*1.5, useDingbats=F)
	ggsave(file=paste("man_otu_", sampleA, "vs", sampleB, ".png", sep=""), p, width = 5*2, height = 3*1.5)

	## Heatmap
	sub_group = subset(sub_design, group %in% c(sampleA, sampleB))
#    # All OTU in two group
#	all=c(enriched, nosig, depleted)
#	sub_norm = as.matrix(cssnorm[all, rownames(sub_group)])
#	pdf(file=paste("heat_otu_", sampleA, "vs", sampleB, "_all.pdf", sep=""), height = 8, width = 8)
#	heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), 
#				cexCol=1,keysize=1,density.info="none",main=NULL,trace="none")
#	dev.off() # if OTU > 1k, great slow in this step
#	png(file=paste("heat_otu_", sampleA, "vs", sampleB, "_all.png", sep=""), height = 8, width = 8, units = "in", res = 300)
#	heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), 
#				cexCol=1,keysize=1,density.info="none",main=NULL,trace="none")
#	dev.off() # if OTU > 1k, great slow in this step
    # Sig OTU in two group
	DE=c(enriched,depleted)
	sub_norm = as.matrix(cssnorm[DE, rownames(sub_group)])
    #colnames(sub_norm)=gsub("DM","KO",colnames(sub_norm),perl=TRUE) # rename samples ID
	pdf(file=paste("heat_otu_", sampleA, "vs", sampleB, "_sig.pdf", sep=""), height = 8, width = 8)
	# scale in row, dendrogram only in row, not cluster in column
	heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none")
	dev.off()
	png(file=paste("heat_otu_", sampleA, "vs", sampleB, "_sig.png", sep=""), height = 8, width = 8, units = "in", res = 300)
	# scale in row, dendrogram only in row, not cluster in column
	heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none")
	dev.off()
#    # Sig OTU in all group
#	sub_norm = as.matrix(cssnorm[DE,rownames(sub_design)])
#	pdf(file=paste("heat_otu_", sampleA, "vs", sampleB, "_data.pdf", sep=""), height = 8, width = 8)
#	# scale in row, dendrogram only in row, not cluster in column
#	heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none")
#	dev.off()
#	png(file=paste("heat_otu_", sampleA, "vs", sampleB, "_data.png", sep=""), height = 8, width = 8, units = "in", res = 300)
#	# scale in row, dendrogram only in row, not cluster in column
#	heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none")
#	dev.off()
 
    style = "none"
	x$percentage=(2^(x$logCPM))/10000
	x$fold=2^(x$logFC)
	if (style == "css"){
		x=cbind(x,cssnorm[rownames(x), rownames(sub_group)])
	}else if (style == "percentage"){
		x=cbind(x,norm[rownames(x), rownames(sub_group)])
	}
    # save each group DA taxonomy summary
    write.table(paste(SampAvsB,length(enriched),length(depleted),length(nosig),sep="\t"), file=paste("otu_sum.txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    # save each group DA OTU detail for plot_pie
	write.table(x[enriched,], file=paste("otu_", sampleA, "vs", sampleB, "_enriched.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
	write.table(x[depleted,], file=paste("otu_", sampleA, "vs", sampleB, "_depleted.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
    # save each group DA OTU list for venndiagram
    write.table(cbind(rownames(x[enriched,]),rep(paste(sampleA, "vs", sampleB, "_enriched", sep=""),length(rownames(x[enriched,]))),x[enriched,]$PValue), file=paste("otu", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    write.table(cbind(rownames(x[depleted,]),rep(paste(sampleA, "vs", sampleB, "_depleted", sep=""),length(rownames(x[depleted,]))),x[depleted,]$PValue), file=paste("otu", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
}

if ("TRUE" == "FALSE") {
	compare_data <- as.vector(unique(sub_design$group))
	len_compare_data <- length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare <- as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			DAOTU_edgeR(tmp_compare)
		}
	}
}else {
	compare_data <- read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) <- c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){DAOTU_edgeR(compare_data[i,])}
}
