# Install related packages
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
taxonomy$full=taxonomy$kingdom

# setting subset design
if (TRUE){
	sub_design = subset(design,genotype %in% c("WT","DM1","DM2","DO1","DO2") ) # select group1
}else{
	sub_design = design
}
if (TRUE){
	sub_design = subset(sub_design,batch %in% c("3") ) # select group2
}

# Set group style, single or combine
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
##############################
# phylum 
##############################

write.table(paste("taxonomy\tSampAvsB\tPvalue",sep="\t"), file=paste("phylum.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

taxonomy$full=paste(taxonomy$full,taxonomy$phylum,sep=";")
tax_count = merge(taxonomy, sub_counts, by="row.names")

tax_count_sum = aggregate(tax_count[,-(1:10)], by=tax_count[10], FUN=sum) # mean
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100


sampFile = as.data.frame(sub_design$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = per
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

mat_se = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=se) # se
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

mean_sort = mat_mean_final[(order(-rowSums(mat_mean_final))), ] # decrease sort
colSums(mat_mean_final)
mean_topRank = mean_sort[1:5, ] # get top 1-5 line
se_topRank = as.matrix(mat_se_final[rownames(mean_topRank), ]) # get same taxonomy line with mean
if (TRUE){
  rownames(mean_topRank) = gsub("[\\w;_]+__","",rownames(mean_topRank),perl=TRUE) # rowname unallowed same name
  rownames(se_topRank) = gsub("[\\w;_]+__","",rownames(se_topRank),perl=TRUE) # rowname unallowed same name
}

# Plotting # par(mfrow=c(2,1))
color = rainbow(length(geno))
pdf(file="tax_bar_phylum_top5.pdf", width=5, height=10) # output to PDF or screen
# modify xlim scale, phylum  recommand 70, and family usually 50
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="phylum", axis.lty=1, xlim = c(0,70), main="phylum distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
png(file="tax_bar_phylum_top5.png", width=5, height=10, units = "in", res = 300) # output to png
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="phylum", axis.lty=1, xlim = c(0,70), main="phylum distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
print("tax_bar_phylum_top5.pdf finished.")

# Stackplot
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[5:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(5-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[5] = c("Low Abundance")
# ordered taxonomy
write.table(rownames(mean_sort), file="tax_phylum.topN", append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

if (TRUE){
    rownames(mean_sort) = gsub("[\\w;_]+__","",rownames(mean_sort),perl=TRUE) # rowname unallowed same name
}
#colSums(mean_sort)

mean_sort$phylum = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("phylum")))
data_all$phylum  = factor(data_all$phylum, levels=rownames(mean_sort))   # set taxonomy order
if ("TRUE" == "TRUE") {
    data_all$variable  = factor(data_all$variable, levels=c("WT","DM1","DM2","DO1","DO2"))   # set group order
}
p = ggplot(data_all, aes(x=variable, y = value, fill = phylum )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme
p
ggsave("tax_stack_phylum_top9.pdf", p, width = 5, height = 3)
ggsave("tax_stack_phylum_top9.png", p, width = 5, height = 3)
print("tax_stack_phylum_top9.pdf finished.")

## Statistics pair group by edgeR
# create DGE list
g = sub_design$group
d = DGEList(counts=tax_count_sum, group=g)
d = calcNormFactors(d)

# fit the GLM
design.mat = model.matrix(~ 0 + group,data=d$samples)
colnames(design.mat) = gsub("group","",colnames(design.mat))
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# function DA edgeR
da_edger = function(sampleV){
	sampleA = as.vector(sampleV$sampA)
	sampleB = as.vector(sampleV$sampB)
	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	# manual setting group1 vs group2, DM1/DO2 vs WT
	SampAvsB=paste(sampleA,"-", sampleB, sep="")
    print(paste("Start DA OTU", SampAvsB, ":",sep=" "))
	BvsA = makeContrasts(contrasts = SampAvsB, levels=design.mat)
	lrt = glmLRT(fit,contrast=BvsA)
    #topTags(lrt) # show top 10 significant pvalue taxonomy or OTU
	de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

	x=lrt$table
	x$sig=de_lrt
	x$level = ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")) # sig by FDR
	#x$level = ifelse(x$PValue<pvalue & x$logFC>0, "enriched",ifelse(x$PValue<pvalue & x$logFC<0, "depleted","nosig")) # sig by pvalue

	# Classify taxonomy by FDR < 0.05
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
        rownames(sub_norm) = gsub("[\\w;_]+p__","",rownames(sub_norm),perl=TRUE) # rowname unallowed same name
        pdf(file=paste("heat_phylum_", sampleA, "vs", sampleB, "_sig.pdf", sep=""), width = 8, height = 8)
        # scale in row, dendrogram only in row, not cluster in column
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
        png(file=paste("heat_phylum_", sampleA, "vs", sampleB, "_sig.png", sep=""), width = 8, height = 8, units = "in", res = 300)
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
    }
	x$percentage=(2^(x$logCPM))/10000
	x$fold=2^(x$logFC)
	x=cbind(x,per[rownames(x), rownames(design2)])
    # save each group DA taxonomy summary
   	write.table(paste("phylum",SampAvsB,length(enriched),length(depleted),length(nosig),sep="\t"), file=paste("tax_sum.txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    # save each group DA taxonomy detail for plot_pie
    write.table(x[enriched,], file=paste("phylum_", sampleA, "vs", sampleB, "_enriched.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
	write.table(x[depleted,], file=paste("phylum_", sampleA, "vs", sampleB, "_depleted.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
    # save each group DA taxonomy list for venndiagram
    write.table(cbind(rownames(x[enriched,]),rep(paste(sampleA, "vs", sampleB, "_enriched", sep=""),length(rownames(x[enriched,]))),x[enriched,]$PValue), file=paste("phylum", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    write.table(cbind(rownames(x[depleted,]),rep(paste(sampleA, "vs", sampleB, "_depleted", sep=""),length(rownames(x[depleted,]))),x[depleted,]$PValue), file=paste("phylum", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
	print(paste("Statistics significant phylum", sampleA, sampleB, "finished!",sep=" "))

}

if ("" == "FALSE") {
	compare_data = as.vector(unique(sub_design$group))
	len_compare_data = length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_edger(tmp_compare)
		}
	}
}else {
	compare_data = read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) = c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_edger(compare_data[i,])}
}	


##############################
# class 
##############################

write.table(paste("taxonomy\tSampAvsB\tPvalue",sep="\t"), file=paste("class.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

taxonomy$full=paste(taxonomy$full,taxonomy$class,sep=";")
tax_count = merge(taxonomy, sub_counts, by="row.names")

tax_count_sum = aggregate(tax_count[,-(1:10)], by=tax_count[10], FUN=sum) # mean
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100


sampFile = as.data.frame(sub_design$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = per
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

mat_se = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=se) # se
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

mean_sort = mat_mean_final[(order(-rowSums(mat_mean_final))), ] # decrease sort
colSums(mat_mean_final)
mean_topRank = mean_sort[1:5, ] # get top 1-5 line
se_topRank = as.matrix(mat_se_final[rownames(mean_topRank), ]) # get same taxonomy line with mean
if (TRUE){
  rownames(mean_topRank) = gsub("[\\w;_]+__","",rownames(mean_topRank),perl=TRUE) # rowname unallowed same name
  rownames(se_topRank) = gsub("[\\w;_]+__","",rownames(se_topRank),perl=TRUE) # rowname unallowed same name
}

# Plotting # par(mfrow=c(2,1))
color = rainbow(length(geno))
pdf(file="tax_bar_class_top5.pdf", width=5, height=10) # output to PDF or screen
# modify xlim scale, class  recommand 70, and family usually 50
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="class", axis.lty=1, xlim = c(0,70), main="class distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
png(file="tax_bar_class_top5.png", width=5, height=10, units = "in", res = 300) # output to png
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="class", axis.lty=1, xlim = c(0,70), main="class distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
print("tax_bar_class_top5.pdf finished.")

# Stackplot
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[5:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(5-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[5] = c("Low Abundance")
# ordered taxonomy
write.table(rownames(mean_sort), file="tax_class.topN", append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

if (TRUE){
    rownames(mean_sort) = gsub("[\\w;_]+__","",rownames(mean_sort),perl=TRUE) # rowname unallowed same name
}
#colSums(mean_sort)

mean_sort$class = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("class")))
data_all$class  = factor(data_all$class, levels=rownames(mean_sort))   # set taxonomy order
if ("TRUE" == "TRUE") {
    data_all$variable  = factor(data_all$variable, levels=c("WT","DM1","DM2","DO1","DO2"))   # set group order
}
p = ggplot(data_all, aes(x=variable, y = value, fill = class )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme
p
ggsave("tax_stack_class_top9.pdf", p, width = 5, height = 3)
ggsave("tax_stack_class_top9.png", p, width = 5, height = 3)
print("tax_stack_class_top9.pdf finished.")

## Statistics pair group by edgeR
# create DGE list
g = sub_design$group
d = DGEList(counts=tax_count_sum, group=g)
d = calcNormFactors(d)

# fit the GLM
design.mat = model.matrix(~ 0 + group,data=d$samples)
colnames(design.mat) = gsub("group","",colnames(design.mat))
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# function DA edgeR
da_edger = function(sampleV){
	sampleA = as.vector(sampleV$sampA)
	sampleB = as.vector(sampleV$sampB)
	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	# manual setting group1 vs group2, DM1/DO2 vs WT
	SampAvsB=paste(sampleA,"-", sampleB, sep="")
    print(paste("Start DA OTU", SampAvsB, ":",sep=" "))
	BvsA = makeContrasts(contrasts = SampAvsB, levels=design.mat)
	lrt = glmLRT(fit,contrast=BvsA)
    #topTags(lrt) # show top 10 significant pvalue taxonomy or OTU
	de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

	x=lrt$table
	x$sig=de_lrt
	x$level = ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")) # sig by FDR
	#x$level = ifelse(x$PValue<pvalue & x$logFC>0, "enriched",ifelse(x$PValue<pvalue & x$logFC<0, "depleted","nosig")) # sig by pvalue

	# Classify taxonomy by FDR < 0.05
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
        rownames(sub_norm) = gsub("[\\w;_]+p__","",rownames(sub_norm),perl=TRUE) # rowname unallowed same name
        pdf(file=paste("heat_class_", sampleA, "vs", sampleB, "_sig.pdf", sep=""), width = 8, height = 8)
        # scale in row, dendrogram only in row, not cluster in column
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
        png(file=paste("heat_class_", sampleA, "vs", sampleB, "_sig.png", sep=""), width = 8, height = 8, units = "in", res = 300)
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
    }
	x$percentage=(2^(x$logCPM))/10000
	x$fold=2^(x$logFC)
	x=cbind(x,per[rownames(x), rownames(design2)])
    # save each group DA taxonomy summary
   	write.table(paste("class",SampAvsB,length(enriched),length(depleted),length(nosig),sep="\t"), file=paste("tax_sum.txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    # save each group DA taxonomy detail for plot_pie
    write.table(x[enriched,], file=paste("class_", sampleA, "vs", sampleB, "_enriched.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
	write.table(x[depleted,], file=paste("class_", sampleA, "vs", sampleB, "_depleted.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
    # save each group DA taxonomy list for venndiagram
    write.table(cbind(rownames(x[enriched,]),rep(paste(sampleA, "vs", sampleB, "_enriched", sep=""),length(rownames(x[enriched,]))),x[enriched,]$PValue), file=paste("class", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    write.table(cbind(rownames(x[depleted,]),rep(paste(sampleA, "vs", sampleB, "_depleted", sep=""),length(rownames(x[depleted,]))),x[depleted,]$PValue), file=paste("class", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
	print(paste("Statistics significant class", sampleA, sampleB, "finished!",sep=" "))

}

if ("" == "FALSE") {
	compare_data = as.vector(unique(sub_design$group))
	len_compare_data = length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_edger(tmp_compare)
		}
	}
}else {
	compare_data = read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) = c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_edger(compare_data[i,])}
}	


##############################
# order 
##############################

write.table(paste("taxonomy\tSampAvsB\tPvalue",sep="\t"), file=paste("order.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

taxonomy$full=paste(taxonomy$full,taxonomy$order,sep=";")
tax_count = merge(taxonomy, sub_counts, by="row.names")

tax_count_sum = aggregate(tax_count[,-(1:10)], by=tax_count[10], FUN=sum) # mean
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100


sampFile = as.data.frame(sub_design$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = per
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

mat_se = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=se) # se
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

mean_sort = mat_mean_final[(order(-rowSums(mat_mean_final))), ] # decrease sort
colSums(mat_mean_final)
mean_topRank = mean_sort[1:5, ] # get top 1-5 line
se_topRank = as.matrix(mat_se_final[rownames(mean_topRank), ]) # get same taxonomy line with mean
if (TRUE){
  rownames(mean_topRank) = gsub("[\\w;_]+__","",rownames(mean_topRank),perl=TRUE) # rowname unallowed same name
  rownames(se_topRank) = gsub("[\\w;_]+__","",rownames(se_topRank),perl=TRUE) # rowname unallowed same name
}

# Plotting # par(mfrow=c(2,1))
color = rainbow(length(geno))
pdf(file="tax_bar_order_top5.pdf", width=5, height=10) # output to PDF or screen
# modify xlim scale, order  recommand 70, and family usually 50
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="order", axis.lty=1, xlim = c(0,70), main="order distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
png(file="tax_bar_order_top5.png", width=5, height=10, units = "in", res = 300) # output to png
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="order", axis.lty=1, xlim = c(0,70), main="order distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
print("tax_bar_order_top5.pdf finished.")

# Stackplot
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[5:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(5-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[5] = c("Low Abundance")
# ordered taxonomy
write.table(rownames(mean_sort), file="tax_order.topN", append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

if (TRUE){
    rownames(mean_sort) = gsub("[\\w;_]+__","",rownames(mean_sort),perl=TRUE) # rowname unallowed same name
}
#colSums(mean_sort)

mean_sort$order = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("order")))
data_all$order  = factor(data_all$order, levels=rownames(mean_sort))   # set taxonomy order
if ("TRUE" == "TRUE") {
    data_all$variable  = factor(data_all$variable, levels=c("WT","DM1","DM2","DO1","DO2"))   # set group order
}
p = ggplot(data_all, aes(x=variable, y = value, fill = order )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme
p
ggsave("tax_stack_order_top9.pdf", p, width = 5, height = 3)
ggsave("tax_stack_order_top9.png", p, width = 5, height = 3)
print("tax_stack_order_top9.pdf finished.")

## Statistics pair group by edgeR
# create DGE list
g = sub_design$group
d = DGEList(counts=tax_count_sum, group=g)
d = calcNormFactors(d)

# fit the GLM
design.mat = model.matrix(~ 0 + group,data=d$samples)
colnames(design.mat) = gsub("group","",colnames(design.mat))
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# function DA edgeR
da_edger = function(sampleV){
	sampleA = as.vector(sampleV$sampA)
	sampleB = as.vector(sampleV$sampB)
	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	# manual setting group1 vs group2, DM1/DO2 vs WT
	SampAvsB=paste(sampleA,"-", sampleB, sep="")
    print(paste("Start DA OTU", SampAvsB, ":",sep=" "))
	BvsA = makeContrasts(contrasts = SampAvsB, levels=design.mat)
	lrt = glmLRT(fit,contrast=BvsA)
    #topTags(lrt) # show top 10 significant pvalue taxonomy or OTU
	de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

	x=lrt$table
	x$sig=de_lrt
	x$level = ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")) # sig by FDR
	#x$level = ifelse(x$PValue<pvalue & x$logFC>0, "enriched",ifelse(x$PValue<pvalue & x$logFC<0, "depleted","nosig")) # sig by pvalue

	# Classify taxonomy by FDR < 0.05
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
        rownames(sub_norm) = gsub("[\\w;_]+p__","",rownames(sub_norm),perl=TRUE) # rowname unallowed same name
        pdf(file=paste("heat_order_", sampleA, "vs", sampleB, "_sig.pdf", sep=""), width = 8, height = 8)
        # scale in row, dendrogram only in row, not cluster in column
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
        png(file=paste("heat_order_", sampleA, "vs", sampleB, "_sig.png", sep=""), width = 8, height = 8, units = "in", res = 300)
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
    }
	x$percentage=(2^(x$logCPM))/10000
	x$fold=2^(x$logFC)
	x=cbind(x,per[rownames(x), rownames(design2)])
    # save each group DA taxonomy summary
   	write.table(paste("order",SampAvsB,length(enriched),length(depleted),length(nosig),sep="\t"), file=paste("tax_sum.txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    # save each group DA taxonomy detail for plot_pie
    write.table(x[enriched,], file=paste("order_", sampleA, "vs", sampleB, "_enriched.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
	write.table(x[depleted,], file=paste("order_", sampleA, "vs", sampleB, "_depleted.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
    # save each group DA taxonomy list for venndiagram
    write.table(cbind(rownames(x[enriched,]),rep(paste(sampleA, "vs", sampleB, "_enriched", sep=""),length(rownames(x[enriched,]))),x[enriched,]$PValue), file=paste("order", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    write.table(cbind(rownames(x[depleted,]),rep(paste(sampleA, "vs", sampleB, "_depleted", sep=""),length(rownames(x[depleted,]))),x[depleted,]$PValue), file=paste("order", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
	print(paste("Statistics significant order", sampleA, sampleB, "finished!",sep=" "))

}

if ("" == "FALSE") {
	compare_data = as.vector(unique(sub_design$group))
	len_compare_data = length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_edger(tmp_compare)
		}
	}
}else {
	compare_data = read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) = c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_edger(compare_data[i,])}
}	


##############################
# family 
##############################

write.table(paste("taxonomy\tSampAvsB\tPvalue",sep="\t"), file=paste("family.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

taxonomy$full=paste(taxonomy$full,taxonomy$family,sep=";")
tax_count = merge(taxonomy, sub_counts, by="row.names")

tax_count_sum = aggregate(tax_count[,-(1:10)], by=tax_count[10], FUN=sum) # mean
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100


sampFile = as.data.frame(sub_design$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = per
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

mat_se = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=se) # se
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

mean_sort = mat_mean_final[(order(-rowSums(mat_mean_final))), ] # decrease sort
colSums(mat_mean_final)
mean_topRank = mean_sort[1:5, ] # get top 1-5 line
se_topRank = as.matrix(mat_se_final[rownames(mean_topRank), ]) # get same taxonomy line with mean
if (TRUE){
  rownames(mean_topRank) = gsub("[\\w;_]+__","",rownames(mean_topRank),perl=TRUE) # rowname unallowed same name
  rownames(se_topRank) = gsub("[\\w;_]+__","",rownames(se_topRank),perl=TRUE) # rowname unallowed same name
}

# Plotting # par(mfrow=c(2,1))
color = rainbow(length(geno))
pdf(file="tax_bar_family_top5.pdf", width=5, height=10) # output to PDF or screen
# modify xlim scale, family  recommand 70, and family usually 50
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="family", axis.lty=1, xlim = c(0,70), main="family distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
png(file="tax_bar_family_top5.png", width=5, height=10, units = "in", res = 300) # output to png
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="family", axis.lty=1, xlim = c(0,70), main="family distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
print("tax_bar_family_top5.pdf finished.")

# Stackplot
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[5:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(5-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[5] = c("Low Abundance")
# ordered taxonomy
write.table(rownames(mean_sort), file="tax_family.topN", append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

if (TRUE){
    rownames(mean_sort) = gsub("[\\w;_]+__","",rownames(mean_sort),perl=TRUE) # rowname unallowed same name
}
#colSums(mean_sort)

mean_sort$family = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("family")))
data_all$family  = factor(data_all$family, levels=rownames(mean_sort))   # set taxonomy order
if ("TRUE" == "TRUE") {
    data_all$variable  = factor(data_all$variable, levels=c("WT","DM1","DM2","DO1","DO2"))   # set group order
}
p = ggplot(data_all, aes(x=variable, y = value, fill = family )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme
p
ggsave("tax_stack_family_top9.pdf", p, width = 5, height = 3)
ggsave("tax_stack_family_top9.png", p, width = 5, height = 3)
print("tax_stack_family_top9.pdf finished.")

## Statistics pair group by edgeR
# create DGE list
g = sub_design$group
d = DGEList(counts=tax_count_sum, group=g)
d = calcNormFactors(d)

# fit the GLM
design.mat = model.matrix(~ 0 + group,data=d$samples)
colnames(design.mat) = gsub("group","",colnames(design.mat))
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# function DA edgeR
da_edger = function(sampleV){
	sampleA = as.vector(sampleV$sampA)
	sampleB = as.vector(sampleV$sampB)
	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	# manual setting group1 vs group2, DM1/DO2 vs WT
	SampAvsB=paste(sampleA,"-", sampleB, sep="")
    print(paste("Start DA OTU", SampAvsB, ":",sep=" "))
	BvsA = makeContrasts(contrasts = SampAvsB, levels=design.mat)
	lrt = glmLRT(fit,contrast=BvsA)
    #topTags(lrt) # show top 10 significant pvalue taxonomy or OTU
	de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

	x=lrt$table
	x$sig=de_lrt
	x$level = ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")) # sig by FDR
	#x$level = ifelse(x$PValue<pvalue & x$logFC>0, "enriched",ifelse(x$PValue<pvalue & x$logFC<0, "depleted","nosig")) # sig by pvalue

	# Classify taxonomy by FDR < 0.05
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
        rownames(sub_norm) = gsub("[\\w;_]+p__","",rownames(sub_norm),perl=TRUE) # rowname unallowed same name
        pdf(file=paste("heat_family_", sampleA, "vs", sampleB, "_sig.pdf", sep=""), width = 8, height = 8)
        # scale in row, dendrogram only in row, not cluster in column
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
        png(file=paste("heat_family_", sampleA, "vs", sampleB, "_sig.png", sep=""), width = 8, height = 8, units = "in", res = 300)
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
    }
	x$percentage=(2^(x$logCPM))/10000
	x$fold=2^(x$logFC)
	x=cbind(x,per[rownames(x), rownames(design2)])
    # save each group DA taxonomy summary
   	write.table(paste("family",SampAvsB,length(enriched),length(depleted),length(nosig),sep="\t"), file=paste("tax_sum.txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    # save each group DA taxonomy detail for plot_pie
    write.table(x[enriched,], file=paste("family_", sampleA, "vs", sampleB, "_enriched.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
	write.table(x[depleted,], file=paste("family_", sampleA, "vs", sampleB, "_depleted.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
    # save each group DA taxonomy list for venndiagram
    write.table(cbind(rownames(x[enriched,]),rep(paste(sampleA, "vs", sampleB, "_enriched", sep=""),length(rownames(x[enriched,]))),x[enriched,]$PValue), file=paste("family", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    write.table(cbind(rownames(x[depleted,]),rep(paste(sampleA, "vs", sampleB, "_depleted", sep=""),length(rownames(x[depleted,]))),x[depleted,]$PValue), file=paste("family", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
	print(paste("Statistics significant family", sampleA, sampleB, "finished!",sep=" "))

}

if ("" == "FALSE") {
	compare_data = as.vector(unique(sub_design$group))
	len_compare_data = length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_edger(tmp_compare)
		}
	}
}else {
	compare_data = read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) = c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_edger(compare_data[i,])}
}	


##############################
# genus 
##############################

write.table(paste("taxonomy\tSampAvsB\tPvalue",sep="\t"), file=paste("genus.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

taxonomy$full=paste(taxonomy$full,taxonomy$genus,sep=";")
tax_count = merge(taxonomy, sub_counts, by="row.names")

tax_count_sum = aggregate(tax_count[,-(1:10)], by=tax_count[10], FUN=sum) # mean
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100


sampFile = as.data.frame(sub_design$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
mat = per
mat_t = t(mat)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

mat_se = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=se) # se
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

mean_sort = mat_mean_final[(order(-rowSums(mat_mean_final))), ] # decrease sort
colSums(mat_mean_final)
mean_topRank = mean_sort[1:5, ] # get top 1-5 line
se_topRank = as.matrix(mat_se_final[rownames(mean_topRank), ]) # get same taxonomy line with mean
if (TRUE){
  rownames(mean_topRank) = gsub("[\\w;_]+__","",rownames(mean_topRank),perl=TRUE) # rowname unallowed same name
  rownames(se_topRank) = gsub("[\\w;_]+__","",rownames(se_topRank),perl=TRUE) # rowname unallowed same name
}

# Plotting # par(mfrow=c(2,1))
color = rainbow(length(geno))
pdf(file="tax_bar_genus_top5.pdf", width=5, height=10) # output to PDF or screen
# modify xlim scale, genus  recommand 70, and family usually 50
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="genus", axis.lty=1, xlim = c(0,70), main="genus distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
png(file="tax_bar_genus_top5.png", width=5, height=10, units = "in", res = 300) # output to png
bar_mean = barplot(t(mean_topRank), horiz=TRUE, beside=TRUE ,col=color, xlab = "Percentage (%)", ylab="genus", axis.lty=1, xlim = c(0,70), main="genus distribution")
error.barsDB(t(mean_topRank), bar_mean, t(se_topRank))
geno=as.vector(geno)
legend("topright", geno, cex=1, bty="n", fill=color)
dev.off()
print("tax_bar_genus_top5.pdf finished.")

# Stackplot
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[5:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(5-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[5] = c("Low Abundance")
# ordered taxonomy
write.table(rownames(mean_sort), file="tax_genus.topN", append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)

if (TRUE){
    rownames(mean_sort) = gsub("[\\w;_]+__","",rownames(mean_sort),perl=TRUE) # rowname unallowed same name
}
#colSums(mean_sort)

mean_sort$genus = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("genus")))
data_all$genus  = factor(data_all$genus, levels=rownames(mean_sort))   # set taxonomy order
if ("TRUE" == "TRUE") {
    data_all$variable  = factor(data_all$variable, levels=c("WT","DM1","DM2","DO1","DO2"))   # set group order
}
p = ggplot(data_all, aes(x=variable, y = value, fill = genus )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme
p
ggsave("tax_stack_genus_top9.pdf", p, width = 5, height = 3)
ggsave("tax_stack_genus_top9.png", p, width = 5, height = 3)
print("tax_stack_genus_top9.pdf finished.")

## Statistics pair group by edgeR
# create DGE list
g = sub_design$group
d = DGEList(counts=tax_count_sum, group=g)
d = calcNormFactors(d)

# fit the GLM
design.mat = model.matrix(~ 0 + group,data=d$samples)
colnames(design.mat) = gsub("group","",colnames(design.mat))
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

# function DA edgeR
da_edger = function(sampleV){
	sampleA = as.vector(sampleV$sampA)
	sampleB = as.vector(sampleV$sampB)
	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	# manual setting group1 vs group2, DM1/DO2 vs WT
	SampAvsB=paste(sampleA,"-", sampleB, sep="")
    print(paste("Start DA OTU", SampAvsB, ":",sep=" "))
	BvsA = makeContrasts(contrasts = SampAvsB, levels=design.mat)
	lrt = glmLRT(fit,contrast=BvsA)
    #topTags(lrt) # show top 10 significant pvalue taxonomy or OTU
	de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)

	x=lrt$table
	x$sig=de_lrt
	x$level = ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")) # sig by FDR
	#x$level = ifelse(x$PValue<pvalue & x$logFC>0, "enriched",ifelse(x$PValue<pvalue & x$logFC<0, "depleted","nosig")) # sig by pvalue

	# Classify taxonomy by FDR < 0.05
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
        rownames(sub_norm) = gsub("[\\w;_]+p__","",rownames(sub_norm),perl=TRUE) # rowname unallowed same name
        pdf(file=paste("heat_genus_", sampleA, "vs", sampleB, "_sig.pdf", sep=""), width = 8, height = 8)
        # scale in row, dendrogram only in row, not cluster in column
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
        png(file=paste("heat_genus_", sampleA, "vs", sampleB, "_sig.png", sep=""), width = 8, height = 8, units = "in", res = 300)
        heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1,keysize=1,density.info="none",main=NULL,trace="none", margins = c(5,15))
        dev.off()
    }
	x$percentage=(2^(x$logCPM))/10000
	x$fold=2^(x$logFC)
	x=cbind(x,per[rownames(x), rownames(design2)])
    # save each group DA taxonomy summary
   	write.table(paste("genus",SampAvsB,length(enriched),length(depleted),length(nosig),sep="\t"), file=paste("tax_sum.txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    # save each group DA taxonomy detail for plot_pie
    write.table(x[enriched,], file=paste("genus_", sampleA, "vs", sampleB, "_enriched.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
	write.table(x[depleted,], file=paste("genus_", sampleA, "vs", sampleB, "_depleted.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)
    # save each group DA taxonomy list for venndiagram
    write.table(cbind(rownames(x[enriched,]),rep(paste(sampleA, "vs", sampleB, "_enriched", sep=""),length(rownames(x[enriched,]))),x[enriched,]$PValue), file=paste("genus", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    write.table(cbind(rownames(x[depleted,]),rep(paste(sampleA, "vs", sampleB, "_depleted", sep=""),length(rownames(x[depleted,]))),x[depleted,]$PValue), file=paste("genus", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
	print(paste("Statistics significant genus", sampleA, sampleB, "finished!",sep=" "))

}

if ("" == "FALSE") {
	compare_data = as.vector(unique(sub_design$group))
	len_compare_data = length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_edger(tmp_compare)
		}
	}
}else {
	compare_data = read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) = c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_edger(compare_data[i,])}
}	


