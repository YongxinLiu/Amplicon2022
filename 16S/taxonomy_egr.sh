#!/bin/bash
set -e

### Default parameter
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
pvalue=0.05
width=4
height=2.5
text_size=8

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    taxonomy_egr.sh
Revision:    1.6
Date:        2017/5/10
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: This script is used to perform calculate taxonomy in phylum, order, class, family, genus in bar and stack plot, and statistic by edgeR.
Notes:       Statistic count by edgeR/DESeq2 negtive binomial is better than percentage/rpm by limma
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
The first one , can calculate taxonomy bar or stack plot, and statistics by GLM LRT of edgeR
Version 1.1 2017/4/7
Add detail usage, add order according with input list group, adjust font size to 8 and figure size to 4x2.5. 
Version 1.2 2017/4/8
Debug: Input group level introduce wrong with no order vector, modify design.mat and geno order; 
only output phylum and faimly, too much result can confuse; remove e_stat of filename
Version 1.3 2017/4/10
Heatmap add margin c(5,5) to c(5,15), size from 8x8 to 8x10, remove taxonomy k__xxx;p__ for shorten
Version 1.4 2017/4/11
Change hyper to enriched, hypo to depleted; adjust summary DA taxonomy nosig to end
Version 1.5 2017/5/3
design and compare file need full or relative(output) directory; add if to compare_file exclude error
Version 1.6 2017/5/10
Add png format figure, for html report

# All input and output should be in -o directory, or give relative -o path, or absolute path
# Input files: design.txt, otu_table.txt, rep_seqs_tax.txt and compare_group.txt

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

# Output file
1. bar plot: tax_bar_\${tax[$level]}_top5.pdf
2. stack plot: tax_stack_\${tax[$level]}_top9.pdf
3. TopN ordered taxonomy name: tax_\${tax[$level]}.topN, for pie plot same color with taxonomy
4. DA taxonomy summary: tax_sum.txt
5. statistics of each group DA list: \${tax[$level]}_sampleAvssampleB_enriched/depleted.txt
6. All DA list and group info, for venndiagram: \${tax[$level]}.txt

OPTIONS:
	-a default=TRUE, shorted taxonomy, only for view taxonomy, for pie chart need FALSE
	-b pvalue threshold, default=0.05
	-c pairs needed to do comparasion, compare_group.txt
	-d design for each samples, default design.txt
	-e execuate Rscript, TRUE or FALSE
	-f OTU table, default in result/otu_table_css.txt
	-g group order by group1 list, default FALSE, not allow both TRUE with -m merge group, when B exist is TRUE
	-h figure height, default 2.5
	-i install package TRUE or FALSE, default FALSE
	-m merge group 1 and group 2 as new group name, default FALSE, not allow both TRUE with -g group order
	-n show top N taxonomy number, default 5, recommend phylum is 5-10
	-o output director, default result_k1-c/
	-p default TRUE for compare group, set FALSE to loop each pair group
	-s text size, default 8
	-t default rep_seqs_tax.txt, taxonomy of all OTU
	-w figure width, default 4
	-A group1 name
	-B group1 selected list
	-C group2 name
	-D group2 selected list

Example:
	taxonomy_egr.sh -o result_k1-c -c compare_group.txt -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D '"2"' -n 5 -s TRUE # select group1/2
	taxonomy_egr.sh -o result_k1-c -c compare_group.txt -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D '"2"' -d design.txt -e TRUE -f otu_table.txt -g TRUE -i FALSE -m FALSE -n 5 -p TRUE -s TRUE -t rep_seqs_tax.txt # all parameter

EOF
}


# Analysis parameter
while getopts "a:b:c:d:e:f:g:h:i:m:n:o:p:s:t:w:A:B:C:D:" OPTION
do
	case $OPTION in
		a)
            short_tax=$OPTARG
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
		n)
			tax_number=$OPTARG
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
		?)
			usage
			exit 1
			;;
	esac
done



cat <<END >taxonomy_egr.r
# Install related packages
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
taxonomy\$full=taxonomy\$kingdom

# setting subset design
if ($select1){
	sub_design = subset(design,${g1} %in% c(${g1_list}) ) # select group1
}else{
	sub_design = design
}
if ($select2){
	sub_design = subset(sub_design,${g2} %in% c(${g2_list}) ) # select group2
}

# Set group style, single or combine
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

idx = rownames(sub_design) %in% colnames(otu_table) 
sub_design = sub_design[idx,]
sub_counts = otu_table[, rownames(sub_design)] # reorder according to design

#############################################################
# Title: Taxonomy barplot + error bar, stackplot scripts
# Author: Yong-Xin Liu
# E-mail: yxliu@genetics.ac.cn
# Date: 3/13/2016
# Description: Script to draw barplots/stackplot of top N taxonomy level, such as phylum and family
# Version 1.3
# Run enviroment: R3.3.2, ggplot2, reshape2
#############################################################
# Define useful functions, standard error (se), and plot se
se = function(x) sd(x)/sqrt(length(x)) # function for Standard Error
# Define function of plotting error bars
error.barsDB = function(x,y,z){g = (max(y)-min(y))/(3*length(y))
for (i in 1:length(y)){lines(c(x[i]+z[i],x[i]-0),c(y[i],y[i]))
  lines(c(x[i]+z[i],x[i]+z[i]),c(y[i]+g,y[i]-g))}}

write.table(paste("taxonomy\tSampAvsB\tenriched\tdepleted\tnosig",sep="\t"), file=paste("tax_sum.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)
END


## Plot each level taxonomy and stat by edgeR, too many result, mainly focus on phylum and family

tax=(domain kingdom phylum class order family genus species)
for level in 2 3 4 5 6
#for level in 2 5
	do

cat <<END >>taxonomy_egr.r
##############################
# ${tax[$level]} 
##############################

write.table(paste("taxonomy\tSampAvsB\tPvalue",sep="\t"), file=paste("${tax[$level]}.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

taxonomy\$full=paste(taxonomy\$full,taxonomy\$${tax[$level]},sep=";")
tax_count = merge(taxonomy, sub_counts, by="row.names")

tax_count_sum = aggregate(tax_count[,-(1:10)], by=tax_count[10], FUN=sum) # mean
rownames(tax_count_sum) = tax_count_sum\$full
tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100


sampFile = as.data.frame(sub_design\$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = per
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean\$group
colnames(mat_mean_final) = geno

mat_se = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=se) # se
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

mean_sort = mat_mean_final[(order(-rowSums(mat_mean_final))), ] # decrease sort
colSums(mat_mean_final)
mean_topRank = mean_sort[1:5, ] # get top 1-5 line
se_topRank = as.matrix(mat_se_final[rownames(mean_topRank), ]) # get same taxonomy line with mean
if (${short_tax}){
  rownames(mean_topRank) = gsub("[\\\\w;_]+__","",rownames(mean_topRank),perl=TRUE) # rowname unallowed same name
  rownames(se_topRank) = gsub("[\\\\w;_]+__","",rownames(se_topRank),perl=TRUE) # rowname unallowed same name
}

# Plotting # par(mfrow=c(2,1))
color = rainbow(length(geno))
pdf(file="tax_bar_${tax[$level]}_top5.pdf", width=5, height=10) # output to PDF or screen
# modify xlim scale, ${tax[$level]}  recommand 70, and family usually 50
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="${tax[$level]}", axis.lty=1, xlim = c(0,70), main="${tax[$level]} distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
png(file="tax_bar_${tax[$level]}_top5.png", width=5, height=10, units = "in", res = 300) # output to png
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="${tax[$level]}", axis.lty=1, xlim = c(0,70), main="${tax[$level]} distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
print("tax_bar_${tax[$level]}_top5.pdf finished.")

# Stackplot
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[${tax_number}:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(${tax_number}-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[${tax_number}] = c("Low Abundance")
# ordered taxonomy
write.table(rownames(mean_sort), file="tax_${tax[$level]}.topN", append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

if (${short_tax}){
    rownames(mean_sort) = gsub("[\\\\w;_]+__","",rownames(mean_sort),perl=TRUE) # rowname unallowed same name
}
#colSums(mean_sort)

mean_sort\$${tax[$level]} = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("${tax[$level]}")))
data_all\$${tax[$level]}  = factor(data_all\$${tax[$level]}, levels=rownames(mean_sort))   # set taxonomy order
if ("${group_order}" == "TRUE") {
    data_all\$variable  = factor(data_all\$variable, levels=c(${g1_list}))   # set group order
}
p = ggplot(data_all, aes(x=variable, y = value, fill = ${tax[$level]} )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme
p
ggsave("tax_stack_${tax[$level]}_top9.pdf", p, width = ${width}, height = ${height})
ggsave("tax_stack_${tax[$level]}_top9.png", p, width = ${width}, height = ${height})
print("tax_stack_${tax[$level]}_top9.pdf finished.")

## Statistics pair group by edgeR
# create DGE list
g = sub_design\$group
d = DGEList(counts=tax_count_sum, group=g)
d = calcNormFactors(d)

# fit the GLM
design.mat = model.matrix(~ 0 + group,data=d\$samples)
colnames(design.mat) = gsub("group","",colnames(design.mat))
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# function DA edgeR
da_edger = function(sampleV){
	sampleA = as.vector(sampleV\$sampA)
	sampleB = as.vector(sampleV\$sampB)
	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	# manual setting group1 vs group2, DM1/DO2 vs WT
	SampAvsB=paste(sampleA,"-", sampleB, sep="")
    print(paste("Start DA OTU", SampAvsB, ":",sep=" "))
	BvsA = makeContrasts(contrasts = SampAvsB, levels=design.mat)
	lrt = glmLRT(fit,contrast=BvsA)
    #topTags(lrt) # show top 10 significant pvalue taxonomy or OTU
	de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=${pvalue})

	x=lrt\$table
	x\$sig=de_lrt
	x\$level = ifelse(x\$sig==1, "enriched",ifelse(x\$sig==-1, "depleted","nosig")) # sig by FDR
	#x\$level = ifelse(x\$PValue<pvalue & x\$logFC>0, "enriched",ifelse(x\$PValue<pvalue & x\$logFC<0, "depleted","nosig")) # sig by pvalue

	# Classify taxonomy by FDR < ${pvalue}
	enriched = row.names(subset(x,sig==1))
	print(paste("enriched",length(enriched),sep=" "))
	nosig = row.names(subset(x,sig== 0))
	print(paste("nosig",length(nosig),sep=" "))
	depleted = row.names(subset(x,sig== -1))
	print(paste("depleted",length(depleted),sep=" "))

	## Heatmap
	DE=c(enriched,depleted)
    if (length(DE)>1){
        sub_norm = as.matrix(per[DE, rownames(design2)])
        rownames(sub_norm) = gsub("[\\\\w;_]+p__","",rownames(sub_norm),perl=TRUE) # rowname unallowed same name
        pdf(file=paste("heat_${tax[$level]}_", sampleA, "vs", sampleB, "_sig.pdf", sep=""), width = 8, height = 8)
        # scale in row, dendrogram only in row, not cluster in column
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
        png(file=paste("heat_${tax[$level]}_", sampleA, "vs", sampleB, "_sig.png", sep=""), width = 8, height = 8, units = "in", res = 300)
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
    }
	x\$percentage=(2^(x\$logCPM))/10000
	x\$fold=2^(x\$logFC)
	x=cbind(x,per[rownames(x), rownames(design2)])
    # save each group DA taxonomy summary
   	write.table(paste("${tax[$level]}",SampAvsB,length(enriched),length(depleted),length(nosig),sep="\t"), file=paste("tax_sum.txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    # save each group DA taxonomy detail for plot_pie
    write.table(x[enriched,], file=paste("${tax[$level]}_", sampleA, "vs", sampleB, "_enriched.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
	write.table(x[depleted,], file=paste("${tax[$level]}_", sampleA, "vs", sampleB, "_depleted.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
    # save each group DA taxonomy list for venndiagram
    write.table(cbind(rownames(x[enriched,]),rep(paste(sampleA, "vs", sampleB, "_enriched", sep=""),length(rownames(x[enriched,]))),x[enriched,]\$PValue), file=paste("${tax[$level]}", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    write.table(cbind(rownames(x[depleted,]),rep(paste(sampleA, "vs", sampleB, "_depleted", sep=""),length(rownames(x[depleted,]))),x[depleted,]\$PValue), file=paste("${tax[$level]}", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
	print(paste("Statistics significant ${tax[$level]}", sampleA, sampleB, "finished!",sep=" "))

}

if ("${compare_pair}" == "FALSE") {
	compare_data = as.vector(unique(sub_design\$group))
	len_compare_data = length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_edger(tmp_compare)
		}
	}
}else {
	compare_data = read.table("${compare_group}", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) = c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_edger(compare_data[i,])}
}	


END
done



if test "${execute}" == "TRUE";
then
	Rscript taxonomy_egr.r
fi
