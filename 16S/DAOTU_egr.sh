#!/bin/bash
set -e

### Default parameter
otu_table_css=otu_table_css.txt
compare_group=compare_group.txt
design=design.txt
execute='TRUE'
otu_table=otu_table.txt
ist='FALSE' # install package, default FALSE
merge_group='FALSE' # combine group1.group2 as new group name, not allow TRUE with group_order
output='result_k1-c' # default work directory, remove low abundance < .1% and p__Cyanobacteria,p__Chloroflexi
pair_group='TRUE' # when pair_group have value turn TRUE
short_tax='TRUE'
taxonomy=rep_seqs_tax.txt
#g1=genotype
#g1_list='"WT","DM1","DM2","DO1","DO2"'
#g2=batch
#g2_list=1
group_order='FALSE' # order group by input list, not allow TRUE with merge_group
select1='FALSE' # filter first group1 info
select2='FALSE' # filter first group2 info
tax_number=5 # set taxonomy number in stackplot
width=4
height=2.5
text_size=8

pvalue=0.05
thre_per=0.005 # threshold for show phylum/family
ymax=15 # manhattan plot ylab max size
fold_max=4 # fold change max for volcano plot
# style include css, percentage, default none
style=none

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    DAOTU_egr.sh
Revision:    1.3
Date:        2017/5/10
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: This script is used to perform calculate different abundance OTUs, and statistic by edgeR.
Notes:       Visualize in volcano plot, manhattan plot and heatmap
-------------------------------------------------------------------------------
Copyright:   2017 (c) Yong-Xin Liu
License:     GPL
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
If any changes are made to this script, please mail me a copy of the changes
-------------------------------------------------------------------------------
Version 1.0 2017/3/27
The first one , can different abundance OTUs, and statistic by edgeR.
Version 1.1 2017/4/8
Add detail usage, add order according with input list group, adjust font size to 8, add heatmap show DA OTU in all samples
Version 1.2 2017/4/11
Change hyper/hypo to enriched/depleted, add levels of taxonomy and significant levels
Version 1.3 2017/5/10
Add png format figure, for html report; if FDR no sig draw wrong color to vol and man plot, auto change to Pvalue; forbidden heatmap of all OTU and samples

# All input and output should be in -o directory, or give relative -o path, or absolute path
# Input files: design.txt, otu_table.txt, rep_seqs_tax.txt, compare_group.txt and otu_table_css.txt

# 1. design.txt, grouping samples, must have SampleID and group info, group1/2 can give to parameter g1 and g2, manually design
SampleID	BarcodeSequence	group	gene	batch	description
WT.1	TAGCTT	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down
WT.2	GGCTAC	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down
WT.3	CGCGCG	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down

# 2. otu_table.txt, reads count of each samples in each OTU
DM2.B1.4        DM1.B1.9        DM1.B1.6        WT.B1.8 DM2.B1.2      
OTU_28  2452.0  2139.0  1313.0  867.0   1883.0  1286.0  1442.0  1157.0
OTU_85  795.0   1034.0  834.0   130.0   942.0   1010.0  505.0   306.0 
OTU_88  206.0   232.0   196.0   36.0    337.0   207.0   108.0   72.0  

#3. rep_seqs_tax.txt, include OTU and kindom to species, and evalue
OTU_28  k__Bacteria     p__Bacteroidetes        c__Cytophagia   o__Cytophagales f__Cytophagaceae        g__     s__     1.000
OTU_85  k__Bacteria     p__Actinobacteria       c__Actinobacteria       o__Actinomycetales      f__Pseudonocardiaceae   g__Pseudonocardia       s__     1.000
OTU_88  k__Bacteria     p__Proteobacteria       c__Alphaproteobacteria  o__Rhodobacterales      f__Hyphomonadaceae      g__     s__     1.000

#4. compare_group.txt, must have SampleA and SampleB in one line and seperated by tab, manually design
OE	WT
KO	WT
OE	KO

# 5. otu_table_css.txt, normalize_table.py css normalized otu table, for heatmap double check edgeR DAOTU
DM2.B1.4        DM1.B1.9        DM1.B1.6        WT.B1.8 
OTU_28  8.6695  8.683   8.3577  9.1612  8.6105  7.9944  
OTU_85  7.0519  7.6381  7.7054  6.4379  7.6149  7.6474  
OTU_88  5.1343  5.5068  5.6386  4.6283  6.1451  5.3884  

# Output file
1. Volcano plot: vol_otu_SampleAvsSampleB.pdf
2. Manhattan plot: man_otu_SampleAvsSampleB.pdf
3. Heatmap: heat_otu_SampleAvsSampleB_all.pdf: all OTU; heat_otu_SampleAvsSampleB_sig.pdf: significant OTU; heat_Sotu_ampleAvsSampleB_data.pdf: sig OTU in all samples
4. DA OTU summary: otu_sum.txt
5. DA OTU detail for plot_pie: otu_SampleAvsSampleB_enriched/depleted.pdf
6. DA OTU list for venndiagram: otu.txt
7. Dissimilarity: 1- pearson correlation: heat_cor_groups/samples.pdf

OPTIONS:
	-a css OTU table, default in otu_table_css.txt
	-b pvalue threshold, default=0.05
	-c pairs needed to do comparasion, compare_group.txt
	-d design for each samples, default design.txt
	-e execuate Rscript, TRUE or FALSE
	-f OTU table, default in otu_table.txt
	-g group order by group1 list, default FALSE, not allow both TRUE with -m merge group, when B exist is TRUE
	-h figure height, default 2.5
	-i install package TRUE or FALSE, default FALSE
	-m merge group 1 and group 2 as new group name, default FALSE, not allow both TRUE with -g group order
	-n show top N taxonomy number, default 5, recommend phylum is 5-10
	-o output director, default result_k1-c/
	-p default TRUE for compare group, set FALSE to loop each pair group
	-s text size, default 8
#	-s default=TRUE, shorted taxonomy, only for view taxonomy, for pie chart need FALSE
	-t default rep_seqs_tax.txt, taxonomy of all OTU
	-w figure width, default 4
	-A group1 name
	-B group1 selected list
	-C group2 name
	-D group2 selected list
	-S sample abundance style, default none, alternative css and percentage
	-h/? show help of script

Example:
	DAOTU_egr.sh -o result_k1-c -c compare_group.txt -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D '"2"' # select group1/2
	DAOTU_egr.sh -o result_k1-c -c compare_group.txt -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D '"2"' -d design.txt -e TRUE -f otu_table.txt -g TRUE -i FALSE -m FALSE -p TRUE -t rep_seqs_tax.txt # all parameter

EOF
}


# Analysis parameter
while getopts "a:b:c:d:e:f:g:h:i:m:o:p:s:t:w:A:B:C:D:" OPTION
do
	case $OPTION in
		a)
			otu_table_css=$OPTARG
			;;
		b)
			pvalue=$OPTARG
			;;
		c)
			compare_group=$OPTARG
			;;
		d)
			design=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		f)
			norm_otu=$OPTARG
			;;
		g)
			group_order=$OPTARG
			;;
		h)
			height=$OPTARG
			;;
		i)
			ist=$OPTARG
			;;
		m)
			merge_group=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		p)
			pair_group=$OPTARG
			;;
		s)
			text_size=$OPTARG
			;;
		t)
			taxonomy=$OPTARG
			;;
		w)
			width=$OPTARG
			;;
		A)
			g1=$OPTARG
			;;
		B)
			g1_list=$OPTARG
			group_order=TRUE
			select1=TRUE
			;;
		C)
			g2=$OPTARG
			;;
		D)
			g2_list=$OPTARG
			select2=TRUE
			;;
		S)
			style=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done



cat <<END >DAOTU_egr.r
if ($ist){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("Biobase","edgeR","ggplot2","gplots","grid","RColorBrewer","reshape2","VennDiagram"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd("${output}")
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
                    axis.text=element_text(color="black", size=${text_size}),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=${text_size}),
                    text=element_text(family="sans", size=${text_size}))

# Public file 1. "design.txt"  Design of experiment
design = read.table("${design}", header=T, row.names= 1, sep="\t") 

# Public file 2. "otu_table.txt"  raw reads count of each OTU in each sample
otu_table = read.delim("${otu_table}", row.names= 1,  header=T, sep="\t")

# Public file 3. "rep_seqs_tax.txt"  taxonomy for each OTU, tab seperated
taxonomy = read.delim("${taxonomy}", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","evalue")

# Public file 4. "otu_table_css.txt", css normalization of OTU table
cssnorm = read.delim("${otu_table_css}", row.names= 1,  header=T, sep="\t")

# setting subset design
if ($select1){
	sub_design = subset(design,${g1} %in% c(${g1_list}) ) # select group1
}else{
	sub_design = design
}

if ($select2){
	sub_design = subset(sub_design,${g2} %in% c(${g2_list}) ) # select group2
}

if ($merge_group){
	sub_design\$group=paste(sub_design\$${g1},sub_design\$${g2},sep = ".")
}else{
	sub_design\$group=sub_design\$${g1}
}

# Set group order
if ("${group_order}" == "TRUE") {
    sub_design\$group  = factor(sub_design\$group, levels=c(${g1_list}))   # set group order
}

print(paste("Number of group: ",length(unique(sub_design\$group)),sep="")) # show group numbers

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
sampFile = as.data.frame(sub_design\$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = norm
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean\$group
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
per\$mean=rowMeans(per)
top_phylum=rownames(per[per\$mean>${thre_per},]) 
print(paste("Number of phylum > ${thre_per}: ",length(top_phylum), sep=""))
top_phylum=gsub("[\\\\w;_]+p__","",top_phylum,perl=TRUE) # regexp in perl mode
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
groups = sub_design\$group
d = DGEList(counts=count, group=groups)
d = calcNormFactors(d)

# fit the GLM
design.mat = model.matrix(~ 0 + d\$samples\$group)
colnames(design.mat)=levels(groups)
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)
write.table(paste("SampAvsB\tenriched\tdepleted\tnosig",sep="\t"), file=paste("otu_sum.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

# function DA limma
DAOTU_edgeR <- function(sampleV){
	sampleA <- as.vector(sampleV\$sampA)
	sampleB <- as.vector(sampleV\$sampB)
#	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	print(paste("Start DA OTU", sampleA, sampleB, ":",sep=" "))
	SampAvsB=paste(sampleA,"-", sampleB, sep="")
	print(SampAvsB)
	BvsA <- makeContrasts(contrasts = SampAvsB, levels=design.mat)
	lrt = glmLRT(fit,contrast=BvsA)
	de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=${pvalue})

	x=lrt\$table
	x\$sig=de_lrt
	# No significant FDR OTU, change to pvalue
	if (dim(x[x\$sig==1,])[1]==0 |dim(x[x\$sig==-1,])[1]==0){
		x\$level = ifelse(x\$PValue<${pvalue} & x\$logFC>0, "enriched",ifelse(x\$PValue<${pvalue} & x\$logFC<0, "depleted","nosig"))
	}else{
		x\$level = ifelse(x\$sig==1, "enriched",ifelse(x\$sig==-1, "depleted","nosig"))
	}
	x\$otu = rownames(x)
	x\$neglogp = -log(x\$PValue)

	# Classify OTU by FDR < 0.05
	enriched = row.names(subset(x,level=="enriched"))
	nosig = row.names(subset(x,level== "nosig"))
	depleted = row.names(subset(x,level=="depleted"))

	# Order OTUs according to taxonomy phylum
	taxonomy = taxonomy[order(taxonomy[, 2]), ]
	idx = rownames(taxonomy) %in% x\$otu
	tax = taxonomy[idx, ] # subset taxonomy from used OTU

	#idx = match(rownames(tax), x\$otu)
	x = x[rownames(tax), ] # reorder according to tax
	x\$tax = gsub("p__","",tax\$phylum,perl=TRUE) 
	x\$phylum = gsub("p__","",tax\$phylum,perl=TRUE) 
	x\$class = gsub("c__","",tax\$class,perl=TRUE)
	x\$order = gsub("o__","",tax\$order,perl=TRUE)
	x\$family = gsub("f__","",tax\$family,perl=TRUE)
	x\$genus = gsub("g__","",tax\$genus,perl=TRUE)
	levels(x\$tax)=c(levels(x\$tax),"Low Abundance")
	x[!(x\$tax %in% top_phylum),]\$tax = "Low Abundance" # no level can get value


	# plot
	if (max(x\$neglogp)>${ymax}){
	  x[x\$neglogp>${ymax},]\$neglogp  = ${ymax}
	}

	# Draw scatter plot
	if (max(x\$logFC)>${fold_max}){x[x\$logFC>${fold_max},]\$logFC = ${fold_max}} # norm x axis
	if (min(x\$logFC)< -${fold_max}){x[x\$logFC< -${fold_max},]\$logFC = -${fold_max}} # norm x axis
	x\$otu = factor(x\$otu, levels=x\$otu)   # set x order
	x\$level = factor(x\$level, levels=c("enriched","depleted","nosig"))
	x\$tax = factor(x\$tax, levels=c(top_phylum,"Low Abundance"))

	# Volcanol plot of fold change vs abundance plot
#	p = ggplot(x, aes(x=logFC, y=logCPM, color=level, size=logCPM, shape=tax)) + geom_point()  + # size by abundance
	p = ggplot(x, aes(x=logFC, y=logCPM, color=level)) + geom_point()  +
	  scale_colour_manual(values=c("red","green","grey"))+ xlim(-4, 4)+
	  #  ylim(5, 20)+
	  labs(x="log2(fold change)",y="log2(count per million)", title=paste(sampleA, "vs", sampleB, sep=" "))+main_theme
	p
	ggsave(file=paste("vol_otu_", sampleA, "vs", sampleB, ".pdf", sep=""), p, width = ${width}, height = ${height})
	ggsave(file=paste("vol_otu_", sampleA, "vs", sampleB, ".png", sep=""), p, width = ${width}, height = ${height})

	# Manhattan plot
	FDR = min(x\$neglogp[x\$level=="enriched"])
	p = ggplot(x, aes(x=otu, y=neglogp, color=tax, size=logCPM, shape=level)) +
	  geom_point(alpha=.7) + 
	  geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
	  #scale_color_manual(values=colors\$colors) +
	  scale_shape_manual(values=c(17, 25, 20))+
	  scale_size(breaks=c(5, 10, 15)) +
#	  scale_size(breaks=c(1, 2, 3)) +
	  labs(x="OTU", y="-loge(P)", title=paste(sampleA, "vs", sampleB, sep=" ")) +main_theme +
	  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="right")
	p
	ggsave(file=paste("man_otu_", sampleA, "vs", sampleB, ".pdf", sep=""), p, width = ${width}*2, height = ${height}*1.5, useDingbats=F)
	ggsave(file=paste("man_otu_", sampleA, "vs", sampleB, ".png", sep=""), p, width = ${width}*2, height = ${height}*1.5)

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
 
    style = "${style}"
	x\$percentage=(2^(x\$logCPM))/10000
	x\$fold=2^(x\$logFC)
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
    write.table(cbind(rownames(x[enriched,]),rep(paste(sampleA, "vs", sampleB, "_enriched", sep=""),length(rownames(x[enriched,]))),x[enriched,]\$PValue), file=paste("otu", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    write.table(cbind(rownames(x[depleted,]),rep(paste(sampleA, "vs", sampleB, "_depleted", sep=""),length(rownames(x[depleted,]))),x[depleted,]\$PValue), file=paste("otu", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
}

if ("${pair_group}" == "FALSE") {
	compare_data <- as.vector(unique(sub_design\$group))
	len_compare_data <- length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare <- as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			DAOTU_edgeR(tmp_compare)
		}
	}
}else {
	compare_data <- read.table("${compare_group}", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) <- c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){DAOTU_edgeR(compare_data[i,])}
}
END



if test "${execute}" == "TRUE";
then
	Rscript DAOTU_egr.r
fi
