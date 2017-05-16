#!/bin/bash
set -e # stop when error

### Default parameter
alpha=alpha.txt
beta=beta
compare_group=compare_group.txt
design=design.txt
execute='TRUE'
norm_otu=otu_table_css.txt
ist='FALSE' # install package, default FALSE
merge_group='FALSE' # combine group1.group2 as new group name, not allow TRUE with group_order
output='result' # default work directory, remove low abundance < .1% and p__Cyanobacteria,p__Chloroflexi
pair_group='TRUE' # when pair_group have value turn TRUE
g1=genotype
g1_list='"WT","DM1","DM2","DO1","DO2"'
g2=batch
g2_list=1
group_order='FALSE' # order group by input list, not allow TRUE with merge_group
select1='FALSE' # filter first group1 info
select2='FALSE' # filter first group2 info
width=4 
height=2.5
text_size=6

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    diversity.sh
Revision:    1.3
Date:        2017/5/10
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: This script is used to perform calculate alpha and beta diversity by R.
Notes:       Statistic alpha by aov and TukeyHSD, beta by adonis of vegan. All input and output should be in -o directory.
-------------------------------------------------------------------------------
Copyright:   2017 (c) Yong-Xin Liu
License:     GPL
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
If any changes are made to this script, please mail me a copy of the changes
-------------------------------------------------------------------------------
Version 1.0 2017/3/26
The first one , can calculate alpha boxplot and beta PCoA, and statistics by aov and adonis
Version 1.1 2017/4/7
Add detail usage, add order according with input list group, adjust font size to 8 and figure size to 4x2.5
Version 1.2 2017/5/3
design and compare file need full or relative(output) directory; add if to compare_file exclude error
Version 1.3 2017/5/10
Add png format figure, for html report

# All input and output should be in -o directory, or give relative -o path, or absolute path
# Input files: design.txt, alpha.txt, beta/${method}_otu_table_css.txt, compare_group.txt, and otu_table_css.txt

# 1. design.txt, grouping samples, must have SampleID and group info, group1/2 can give to parameter g1 and g2, manually design
SampleID	BarcodeSequence	group	gene	batch	description
WT.1	TAGCTT	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down
WT.2	GGCTAC	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down
WT.3	CGCGCG	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down

# 2. alpha.txt, calculate by alpha_diversity.py, boxplot show alpha diversity
	shannon	chao1	observed_otus	PD_whole_tree
DM2.B1.4	6.66186782615	365.037037037	343.0	8.62306
DM1.B1.9	6.5488329159	352.64	323.0	8.46164
DM1.B1.6	6.71016796233	364.107142857	339.0	8.40676

# 3. beta/\${method}_otu_table_css.txt, distance calcuate by beta_diversity.py, scatterplot show beta diversity PCoA, 
method include: bray_curtis weighted_unifrac unweighted_unifrac
DM2.B1.4        DM1.B1.9        DM1.B1.6        WT.B1.8 DM2.B1.2        DM2.B1.6       
DM2.B1.4        0.0     0.0840606514582 0.0894695299217 0.121016152707  0.0897330554942
DM1.B1.9        0.0840606514582 0.0     0.0782910908055 0.123924239637  0.0671587931741
DM1.B1.6        0.0894695299217 0.0782910908055 0.0     0.143871025475  0.067721760309 

#4. compare_group.txt, must have SampleA and SampleB in one line and seperated by tab, manually design
OE	WT
KO	WT
OE	KO

# 5. otu_table_css.txt, normalize_table.py css normalized otu table, for constrained PCoA
DM2.B1.4        DM1.B1.9        DM1.B1.6        WT.B1.8 
OTU_28  8.6695  8.683   8.3577  9.1612  8.6105  7.9944  
OTU_85  7.0519  7.6381  7.7054  6.4379  7.6149  7.6474  
OTU_88  5.1343  5.5068  5.6386  4.6283  6.1451  5.3884  

# Output file
1. Alpha diversity boxplot: alpha_\${method}.pdf, method include shannon chao1 observed_otus PD_whole_tree
2. Alpha diversity statisitcs: alpha_\${method}_stats.txt, aov statistics, TukeyHSD test, conf.level = 0.95
3. Beta diversity scatterplot PCoA: beta_pcoa_\${method}.pdf
4. Beta diversity statisitcs: beta.txt, by adonis of vegan
5. Constrianed PCoA scatterplot: CPCoA_${g1}.pdf, constrianed by condition of group1

OPTIONS:
	-a alpha diversity, default alpha.txt
	-b beta diversity, default beta directory include ${method}_otu_table_css.txt
	-c pairs needed to do comparasion, compare_group.txt
	-d design for each samples, default design.txt
	-e execuate Rscript, TRUE or FALSE
	-f OTU table, default in result/otu_table_css.txt
	-g group order by group1 list, default FALSE, not allow both TRUE with -m merge group, when B exist is TRUE
	-h figure height, default 2.5
	-i install package TRUE or FALSE, default FALSE
	-m merge group 1 and group 2 as new group name, default FALSE, not allow both TRUE with -g group order
	-o output director, default result_k1-c/
	-p default TRUE for compare group, set FALSE to loop each pair group
	-s text size, default 8
	-w figure width, default 4
	-A group1 name
	-B group1 selected list
	-C group2 name
	-D group2 selected list
	-? show help of script

Example:
	diversity.sh -o result_k1-c -c compare_group.txt -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D '"2"' # select group1/2
	diversity.sh -o result_k1-c -c compare_group.txt -A genotype -B '"WT","DM1","DM2","DO1","DO2"' -C batch -D '"2"' -a alpha.txt -b beta -d design.txt -e TRUE -f otu_table_css.txt -g TRUE -i FALSE -m FALSE -p TRUE # all parameter
	diversity.sh -o result_k1-c -A genotype -C batch -D '"2"' -p FALSE # only select group2, compare all groups

EOF
}


# Analysis parameter
while getopts "a:b:c:d:e:f:g:h:i:m:o:p:s:w:A:B:C:D:" OPTION
do
	case $OPTION in
		a)
			alpha=$OPTARG
			;;
		b)
			beta=$OPTARG
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
		w)
            width=$OPTARG
			;;
		A)
			g1=$OPTARG
			;;
		B)
			g1_list=$OPTARG
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




cat <<END >diversity.r
# Install related packages
if ($ist){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","grid","scales","vegan"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("${output}") # set work directory
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
                    axis.text=element_text(color="black", size=${text_size}),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=${text_size}),
                    text=element_text(family="sans", size=${text_size}))

# Public file 1. "design.txt"  Design of experiment
design = read.table("${design}", header=T, row.names= 1, sep="\t") 

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
alpha = read.table("${alpha}", header=T, row.names= 1, sep="\t")
idx = rownames(sub_design) %in% rownames(alpha) # match design with alpha
sub_design_alpha = sub_design [idx,] # sub design for alpha, due to rarefication can remove some samples
alpha = alpha [rownames(sub_design_alpha),] # reorder and subset alpha by design

# set colors by rainbow according with group
colors = data.frame(group=unique(sub_design_alpha\$group), 
    color=rainbow(length(unique(sub_design_alpha\$group)))) 

END


## Plot each alpha diversity
for method in shannon chao1 observed_otus PD_whole_tree
	do

cat <<END >>diversity.r
## ${method} index
# add design to alpha
index = cbind(alpha\$${method}, sub_design_alpha[match(rownames(alpha), rownames(sub_design_alpha)), ]) 
colnames(index)[1] = "${method}" # add ${method} colname is value
p = ggplot(index, aes(x=group, y=${method}, color=group)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors\$color)) +
            labs(x="Groups", y="${method} index") + main_theme
p
ggsave(paste("alpha_${method}.pdf", sep=""), p, width = ${width}, height = ${height})
ggsave(paste("alpha_${method}.png", sep=""), p, width = ${width}, height = ${height})
print("alpha_${method}.pdf finished.")

# ${method} Statistics
${method}_stats <- aov(${method} ~ group, data = index)
Tukey_HSD_${method} <- TukeyHSD(${method}_stats, ordered = FALSE, conf.level = 0.95)
Tukey_HSD_${method}_table <- as.data.frame(Tukey_HSD_${method}\$group)
write.table(Tukey_HSD_${method}_table[order(Tukey_HSD_${method}_table\$p, decreasing=FALSE), ], file="alpha_${method}_stats.txt",append = FALSE, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
print("alpha_${method}_stats.txt finished.")

END
done


cat <<END >>diversity.r
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
colors = data.frame(group=unique(sub_design\$group), 
                    color=rainbow(length(unique(sub_design\$group)))) 
shapes = data.frame(group=unique(sub_design\$group),
                     shape=c(1:length(unique(sub_design\$group))))
write.table("Statistics of PCoA by Adonis of vegan", file="beta.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)
END

for method in bray_curtis weighted_unifrac unweighted_unifrac
	do
cat <<END >>diversity.r
#  PCoA ${method}
${method} = read.table("${beta}/${method}_otu_table_css.txt", sep="\t", header=T, check.names=F)

# subset matrix and design
idx = rownames(sub_design) %in% colnames(${method}) 
sub_design = sub_design[idx,]
${method} = ${method}[rownames(sub_design), rownames(sub_design)] # subset and reorder distance matrix

# cmdscale {stats}, Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
pcoa = cmdscale(${method}, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa\$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa\$eig
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
#points\$group = factor(points\$group, levels=colors\$group)

# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group, shape=group)) +
  geom_point(alpha=.7, size=2) +
  scale_colour_manual(values=as.character(colors\$color)) +
  scale_shape_manual(values=shapes\$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="${method} PCoA") + main_theme
p
ggsave("beta_pcoa_${method}.pdf", p, width = ${width}, height = ${height})
ggsave("beta_pcoa_${method}.png", p, width = ${width}, height = ${height})
print("beta_pcoa_${method}.pdf finished.")

# Compare each group beta by vegan adonis in ${method}
da_adonis <- function(sampleV){
	  sampleA <- as.matrix(sampleV\$sampA)
	  sampleB <- as.matrix(sampleV\$sampB)
	  design2 = subset(sub_design, group %in% c(sampleA,sampleB))
      if (length(unique(design2\$group))>1) {
	  sub_dis_table = dis_table[rownames(design2),rownames(design2)]
	  sub_dis_table <- as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
	  adonis_table = adonis(sub_dis_table~group, data=design2, permutations = 10000) 
	  adonis_pvalue = adonis_table\$aov.tab\$\`Pr(>F)\`[1]
	  print(paste("In ${method}, pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
	  adonis_pvalue <- paste("${method}", sampleA, sampleB, adonis_pvalue, sep="\t")
	  write.table(adonis_pvalue, file="beta.txt", append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
      }
}
dis_table <- as.matrix(${method})
if ("${pair_group}" == "FALSE") {
	compare_data <- as.vector(unique(sub_design\$group))
	len_compare_data <- length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare <- as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			da_adonis(tmp_compare)
		}
	}
}else {
	compare_data <- read.table("${compare_group}", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) <- c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
}	

END
done


cat <<END >>diversity.r
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
  chi = c(cca\$tot.chi, cca\$CCA\$tot.chi, cca\$CA\$tot.chi)
  variability_table = cbind(chi, chi/chi[1])
  colnames(variability_table) = c("inertia", "proportion")
  rownames(variability_table) = c("total", "constrained", "unconstrained")
  return(variability_table)
}

# Load css OTU table, shape color same with beta diversity
otu_table = read.table("${norm_otu}", sep="\t", header=T, row.names= 1) # CSS norm otu table
idx = rownames(sub_design) %in% colnames(otu_table) 
sub_design = sub_design[idx,]
sub_otu_table = otu_table[, rownames(sub_design)] 

# Constrained analysis OTU table by genotype
capscale.gen = capscale(t(sub_otu_table) ~ group, data=sub_design, add=F, sqrt.dist=T, distance="bray") 

# ANOVA-like permutation analysis
perm_anova.gen = anova.cca(capscale.gen)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen = variability_table(capscale.gen)
eig = capscale.gen\$CCA\$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]

# extract the weighted average (sample) scores
points = capscale.gen\$CCA\$wa[, 1:2]
points = as.data.frame(points)
colnames(points) = c("x", "y")
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])

# plot CPCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group, shape=group)) +
  geom_point(alpha=.7, size=1.5) +
  scale_colour_manual(values=as.character(colors\$color)) +
  scale_shape_manual(values=shapes\$shape)+
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",format(p.val, digits=2),sep="")) + main_theme
p
ggsave(paste( "CPCoA_${g1}.pdf", sep=""), p, width = ${width}, height = ${height})
ggsave(paste( "CPCoA_${g1}.png", sep=""), p, width = ${width}, height = ${height})
print("CPCoA_${g1}.pdf finished.")
print("All done!!!")

END




if test "${execute}" == "TRUE";
then
	Rscript diversity.r
fi
