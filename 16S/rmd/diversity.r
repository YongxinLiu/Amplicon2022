# Install related packages
if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","grid","scales","vegan"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("result_k1-c") # set work directory
library("ggplot2") # load related packages
library("grid")
library("scales")
library("vegan")

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

#############################################################
# Title: Alpha diversity - boxplot
# Author: Yong-Xin Liu
# E-mail: yxliu@genetics.ac.cn
# Date: 03/14/2017
# Description: Script to draw alpha diversity - boxplot
# Version 1.1
# Run enviroment: R3.3.1, ggplot2
# Input File required in the blow list: 
# 1."alpha.txt" output by alpha_diversity.py from qiime
#############################################################

# alpha diversity genotype and batch, read file and match with design
alpha = read.table("alpha.txt", header=T, row.names= 1, sep="\t")
idx = rownames(sub_design) %in% rownames(alpha) # match design with alpha
sub_design_alpha = sub_design [idx,] # sub design for alpha, due to rarefication can remove some samples
alpha = alpha [rownames(sub_design_alpha),] # reorder and subset alpha by design

# set colors by rainbow according with group
colors = data.frame(group=unique(sub_design_alpha$group), 
    color=rainbow(length(unique(sub_design_alpha$group)))) 

## shannon index
# add design to alpha
index = cbind(alpha$shannon, sub_design_alpha[match(rownames(alpha), rownames(sub_design_alpha)), ]) 
colnames(index)[1] = "shannon" # add shannon colname is value
p = ggplot(index, aes(x=group, y=shannon, color=group)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Groups", y="shannon index") + main_theme
p
ggsave(paste("alpha_shannon.pdf", sep=""), p, width = 5, height = 3)
ggsave(paste("alpha_shannon.png", sep=""), p, width = 5, height = 3)
print("alpha_shannon.pdf finished.")

# shannon Statistics
shannon_stats <- aov(shannon ~ group, data = index)
Tukey_HSD_shannon <- TukeyHSD(shannon_stats, ordered = FALSE, conf.level = 0.95)
Tukey_HSD_shannon_table <- as.data.frame(Tukey_HSD_shannon$group)
write.table(Tukey_HSD_shannon_table[order(Tukey_HSD_shannon_table$p, decreasing=FALSE), ], file="alpha_shannon_stats.txt",append = FALSE, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
print("alpha_shannon_stats.txt finished.")

## chao1 index
# add design to alpha
index = cbind(alpha$chao1, sub_design_alpha[match(rownames(alpha), rownames(sub_design_alpha)), ]) 
colnames(index)[1] = "chao1" # add chao1 colname is value
p = ggplot(index, aes(x=group, y=chao1, color=group)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Groups", y="chao1 index") + main_theme
p
ggsave(paste("alpha_chao1.pdf", sep=""), p, width = 5, height = 3)
ggsave(paste("alpha_chao1.png", sep=""), p, width = 5, height = 3)
print("alpha_chao1.pdf finished.")

# chao1 Statistics
chao1_stats <- aov(chao1 ~ group, data = index)
Tukey_HSD_chao1 <- TukeyHSD(chao1_stats, ordered = FALSE, conf.level = 0.95)
Tukey_HSD_chao1_table <- as.data.frame(Tukey_HSD_chao1$group)
write.table(Tukey_HSD_chao1_table[order(Tukey_HSD_chao1_table$p, decreasing=FALSE), ], file="alpha_chao1_stats.txt",append = FALSE, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
print("alpha_chao1_stats.txt finished.")

## observed_otus index
# add design to alpha
index = cbind(alpha$observed_otus, sub_design_alpha[match(rownames(alpha), rownames(sub_design_alpha)), ]) 
colnames(index)[1] = "observed_otus" # add observed_otus colname is value
p = ggplot(index, aes(x=group, y=observed_otus, color=group)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Groups", y="observed_otus index") + main_theme
p
ggsave(paste("alpha_observed_otus.pdf", sep=""), p, width = 5, height = 3)
ggsave(paste("alpha_observed_otus.png", sep=""), p, width = 5, height = 3)
print("alpha_observed_otus.pdf finished.")

# observed_otus Statistics
observed_otus_stats <- aov(observed_otus ~ group, data = index)
Tukey_HSD_observed_otus <- TukeyHSD(observed_otus_stats, ordered = FALSE, conf.level = 0.95)
Tukey_HSD_observed_otus_table <- as.data.frame(Tukey_HSD_observed_otus$group)
write.table(Tukey_HSD_observed_otus_table[order(Tukey_HSD_observed_otus_table$p, decreasing=FALSE), ], file="alpha_observed_otus_stats.txt",append = FALSE, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
print("alpha_observed_otus_stats.txt finished.")

## PD_whole_tree index
# add design to alpha
index = cbind(alpha$PD_whole_tree, sub_design_alpha[match(rownames(alpha), rownames(sub_design_alpha)), ]) 
colnames(index)[1] = "PD_whole_tree" # add PD_whole_tree colname is value
p = ggplot(index, aes(x=group, y=PD_whole_tree, color=group)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Groups", y="PD_whole_tree index") + main_theme
p
ggsave(paste("alpha_PD_whole_tree.pdf", sep=""), p, width = 5, height = 3)
ggsave(paste("alpha_PD_whole_tree.png", sep=""), p, width = 5, height = 3)
print("alpha_PD_whole_tree.pdf finished.")

# PD_whole_tree Statistics
PD_whole_tree_stats <- aov(PD_whole_tree ~ group, data = index)
Tukey_HSD_PD_whole_tree <- TukeyHSD(PD_whole_tree_stats, ordered = FALSE, conf.level = 0.95)
Tukey_HSD_PD_whole_tree_table <- as.data.frame(Tukey_HSD_PD_whole_tree$group)
write.table(Tukey_HSD_PD_whole_tree_table[order(Tukey_HSD_PD_whole_tree_table$p, decreasing=FALSE), ], file="alpha_PD_whole_tree_stats.txt",append = FALSE, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
print("alpha_PD_whole_tree_stats.txt finished.")

#############################################################
# Title: Beta diversity - PCoA
# Author: Yong-Xin Liu
# E-mail: yxliu@genetics.ac.cn
# Date: 03/14/2017
# Description: Script to draw beta diversity - scatterplot PCoA
# Version 1.1
# Run enviroment: R3.3.1, ggplot2
# Input File required in the blow list: 
# 1. "beta/bray_curtis_otu_table_css.txt" output by beta_diversity.py from qiime
# 2. "beta/unweighted_unifrac_otu_table_css.txt" 
# 3. "beta/weighted_unifrac_otu_table_css.txt"
#############################################################

# Set colors and shapes
colors = data.frame(group=unique(sub_design$group), 
                    color=rainbow(length(unique(sub_design$group)))) 
shapes = data.frame(group=unique(sub_design$group),
                     shape=c(1:length(unique(sub_design$group))))
write.table("Statistics of PCoA by Adonis of vegan", file="beta.stat", append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)
#  PCoA bray_curtis
bray_curtis = read.table("beta/bray_curtis_otu_table_css.txt", sep="\t", header=T, check.names=F)

# subset matrix and design
idx = rownames(sub_design) %in% colnames(bray_curtis) 
sub_design = sub_design[idx,]
bray_curtis = bray_curtis[rownames(sub_design), rownames(sub_design)] # subset and reorder distance matrix

# cmdscale {stats}, Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
pcoa = cmdscale(bray_curtis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
#points$group = factor(points$group, levels=colors$group)

# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group, shape=group)) +
  geom_point(alpha=.7, size=2) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="bray_curtis PCoA") + main_theme
p
ggsave("beta_pcoa_bray_curtis.pdf", p, width = 5, height = 3)
ggsave("beta_pcoa_bray_curtis.png", p, width = 5, height = 3)
print("beta_pcoa_bray_curtis.pdf finished.")

# Compare each group beta by vegan adonis in bray_curtis
da_adonis <- function(sampleV){
	  sampleA <- as.matrix(sampleV$sampA)
	  sampleB <- as.matrix(sampleV$sampB)
	  design2 = subset(sub_design, group %in% c(sampleA,sampleB))
      if (length(unique(design2$group))>1) {
	  sub_dis_table = dis_table[rownames(design2),rownames(design2)]
	  sub_dis_table <- as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
	  adonis_table = adonis(sub_dis_table~group, data=design2, permutations = 10000) 
	  adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
	  print(paste("In bray_curtis, pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
	  adonis_pvalue <- paste("bray_curtis", sampleA, sampleB, adonis_pvalue, sep="\t")
	  write.table(adonis_pvalue, file="beta.stat", append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
      }
}
dis_table <- as.matrix(bray_curtis)
if ("TRUE" == "FALSE") {
	compare_data <- as.vector(unique(sub_design$group))
	len_compare_data <- length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare <- as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_adonis(tmp_compare)
		}
	}
}else {
	compare_data <- read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) <- c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
}	

#  PCoA weighted_unifrac
weighted_unifrac = read.table("beta/weighted_unifrac_otu_table_css.txt", sep="\t", header=T, check.names=F)

# subset matrix and design
idx = rownames(sub_design) %in% colnames(weighted_unifrac) 
sub_design = sub_design[idx,]
weighted_unifrac = weighted_unifrac[rownames(sub_design), rownames(sub_design)] # subset and reorder distance matrix

# cmdscale {stats}, Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
pcoa = cmdscale(weighted_unifrac, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
#points$group = factor(points$group, levels=colors$group)

# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group, shape=group)) +
  geom_point(alpha=.7, size=2) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="weighted_unifrac PCoA") + main_theme
p
ggsave("beta_pcoa_weighted_unifrac.pdf", p, width = 5, height = 3)
ggsave("beta_pcoa_weighted_unifrac.png", p, width = 5, height = 3)
print("beta_pcoa_weighted_unifrac.pdf finished.")

# Compare each group beta by vegan adonis in weighted_unifrac
da_adonis <- function(sampleV){
	  sampleA <- as.matrix(sampleV$sampA)
	  sampleB <- as.matrix(sampleV$sampB)
	  design2 = subset(sub_design, group %in% c(sampleA,sampleB))
      if (length(unique(design2$group))>1) {
	  sub_dis_table = dis_table[rownames(design2),rownames(design2)]
	  sub_dis_table <- as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
	  adonis_table = adonis(sub_dis_table~group, data=design2, permutations = 10000) 
	  adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
	  print(paste("In weighted_unifrac, pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
	  adonis_pvalue <- paste("weighted_unifrac", sampleA, sampleB, adonis_pvalue, sep="\t")
	  write.table(adonis_pvalue, file="beta.stat", append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
      }
}
dis_table <- as.matrix(weighted_unifrac)
if ("TRUE" == "FALSE") {
	compare_data <- as.vector(unique(sub_design$group))
	len_compare_data <- length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare <- as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_adonis(tmp_compare)
		}
	}
}else {
	compare_data <- read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) <- c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
}	

#  PCoA unweighted_unifrac
unweighted_unifrac = read.table("beta/unweighted_unifrac_otu_table_css.txt", sep="\t", header=T, check.names=F)

# subset matrix and design
idx = rownames(sub_design) %in% colnames(unweighted_unifrac) 
sub_design = sub_design[idx,]
unweighted_unifrac = unweighted_unifrac[rownames(sub_design), rownames(sub_design)] # subset and reorder distance matrix

# cmdscale {stats}, Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
pcoa = cmdscale(unweighted_unifrac, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
#points$group = factor(points$group, levels=colors$group)

# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group, shape=group)) +
  geom_point(alpha=.7, size=2) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="unweighted_unifrac PCoA") + main_theme
p
ggsave("beta_pcoa_unweighted_unifrac.pdf", p, width = 5, height = 3)
ggsave("beta_pcoa_unweighted_unifrac.png", p, width = 5, height = 3)
print("beta_pcoa_unweighted_unifrac.pdf finished.")

# Compare each group beta by vegan adonis in unweighted_unifrac
da_adonis <- function(sampleV){
	  sampleA <- as.matrix(sampleV$sampA)
	  sampleB <- as.matrix(sampleV$sampB)
	  design2 = subset(sub_design, group %in% c(sampleA,sampleB))
      if (length(unique(design2$group))>1) {
	  sub_dis_table = dis_table[rownames(design2),rownames(design2)]
	  sub_dis_table <- as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
	  adonis_table = adonis(sub_dis_table~group, data=design2, permutations = 10000) 
	  adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
	  print(paste("In unweighted_unifrac, pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
	  adonis_pvalue <- paste("unweighted_unifrac", sampleA, sampleB, adonis_pvalue, sep="\t")
	  write.table(adonis_pvalue, file="beta.stat", append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
      }
}
dis_table <- as.matrix(unweighted_unifrac)
if ("TRUE" == "FALSE") {
	compare_data <- as.vector(unique(sub_design$group))
	len_compare_data <- length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare <- as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_adonis(tmp_compare)
		}
	}
}else {
	compare_data <- read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) <- c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
}	

#############################################################
# Title: Beta diversity - Constraind PCoA
# Author: Yong-Xin Liu
# E-mail: yxliu@genetics.ac.cn
# Date: 3/14/2016
# Description: Script to draw beta diversity - Constraind PCoA
# Version 1.1
# Run enviroment: R3.3.1, ggplot2
# Input File required in the blow list: 
# 1. "otu_table_css.txt" output by beta_diversity.py from qiime
#############################################################

# Function for analysis CPCoA/CCA result
variability_table = function(cca){
  chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table = cbind(chi, chi/chi[1])
  colnames(variability_table) = c("inertia", "proportion")
  rownames(variability_table) = c("total", "constrained", "unconstrained")
  return(variability_table)
}

# Load css OTU table, shape color same with beta diversity
otu_table = read.table("otu_table_css.txt", sep="\t", header=T, row.names= 1) # CSS norm otu table
idx = rownames(sub_design) %in% colnames(otu_table) 
sub_design = sub_design[idx,]
sub_otu_table = otu_table[, rownames(sub_design)] 

# Constrained analysis OTU table by genotype
capscale.gen = capscale(t(sub_otu_table) ~ group, data=sub_design, add=F, sqrt.dist=T, distance="bray") 

# ANOVA-like permutation analysis
perm_anova.gen = anova.cca(capscale.gen)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen = variability_table(capscale.gen)
eig = capscale.gen$CCA$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]

# extract the weighted average (sample) scores
points = capscale.gen$CCA$wa[, 1:2]
points = as.data.frame(points)
colnames(points) = c("x", "y")
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])

# plot CPCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group, shape=group)) +
  geom_point(alpha=.7, size=1.5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape)+
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",format(p.val, digits=2),sep="")) + main_theme
p
ggsave(paste( "CPCoA_genotype.pdf", sep=""), p, width = 5, height = 3)
ggsave(paste( "CPCoA_genotype.png", sep=""), p, width = 5, height = 3)
print("CPCoA_genotype.pdf finished.")
print("All done!!!")

