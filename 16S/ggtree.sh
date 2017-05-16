#!/bin/bash
set -e

### Default parameter
execute='TRUE'
ist='FALSE' # install package, default FALSE
output='result' # default work directory, remove low abundance < .1% and p__Cyanobacteria,p__Chloroflexi
width=8
height=8
text_size=7

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    ggtree.sh
Revision:    1.0
Date:        2017/5/9
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: ggtree draw basic circle tree, color by taxonomy
Notes:       Default output order, family and genus colored tree
-------------------------------------------------------------------------------
Copyright:   2017 (c) Yong-Xin Liu
License:     GPL
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
If any changes are made to this script, please mail me a copy of the changes
-------------------------------------------------------------------------------
Version 1.0 2017/5/9
The first one ,ggtree draw basic circle tree, color by taxonomy, default output order, family and genus colored tree. First R in /mnt/bai/yongxin/bin/R/ggtree1.0.r
Version 1.1 2017/5/10
Debug: labels=unique can result wrong taxonomy lable, modify to labels=levels is correct

# All input and output should be in -o directory, or give relative -o path, or absolute path
# Input files: tax_rep_seqs.tree, tax_rep_seqs.tax

# 1. result/tax_rep_seqs.tree
((OTU_1353:0.03354,OTU_384:0.01985)0.849:0.01531,((OTU_2669:0.02448,((OTU_4528:0.02014,

# 2. result/tax_rep_seqs.tax
OTU_48	k__Bacteria	p__Proteobacteria	c__Betaproteobacteria	o__SC-I-84
OTU_384	k__Bacteria	p__Actinobacteria	c__Actinobacteria	o__Actinomycetales
OTU_399	k__Bacteria	p__Proteobacteria	c__Alphaproteobacteria	o__Rhizobiales

# Output file
1. bar plot: ggtree_order.pdf .png

OPTIONS:
	-e execuate Rscript, TRUE or FALSE
	-h figure height, default 2.5
	-i install package TRUE or FALSE, default FALSE
	-o output director, default result
	-s text size, default 7
	-w figure width, default 4

Example:
	ggtree.sh -o result

EOF
}


# Analysis parameter
while getopts "a:b:c:d:e:f:g:h:i:l:m:n:o:p:s:t:w:A:B:C:D:" OPTION
do
	case $OPTION in
		d)
			design=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		h)
			height=$OPTARG
			;;
		i)
			ist=$OPTARG
			;;
		l)
			library=$OPTARG
			;;
		m)
			merge_group=$OPTARG
			;;
		o)
			output=$OPTARG
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
		C)
			g2=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done



cat <<END >ggtree.r
# Install related packages
# Work well in R3.3.3
if ($ist){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggtree","colorspace"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("${output}")
library("ggtree")
library("colorspace")

tree <- read.tree("tax_rep_seqs.tree")
tax <- read.table("tax_rep_seqs.tax",row.names=1)
colnames(tax) = c("kingdom","phylum","class","order")

END


for tax in phylum class order
	do
cat <<END >>ggtree.r

groupInfo <- split(row.names(tax), tax\$${tax}) # OTU and ${tax} for group
tree <- groupOTU(tree, groupInfo)
pdf(file="ggtree_${tax}.pdf", width=${width}, height=${height})
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+
  scale_color_manual(values=c(rainbow_hcl(length(unique(tax\$${tax}))+1)), breaks=1:length(unique(tax\$${tax})), labels=levels(tax\$${tax}))+
  theme(legend.position = "right") +geom_tiplab2(size=3)
dev.off()
png(file="ggtree_${tax}.png", width=${width}, height=${height}, units = "in", res = 300)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+
  scale_color_manual(values=c(rainbow_hcl(length(unique(tax\$${tax}))+1)), breaks=1:length(unique(tax\$${tax})), labels=levels(tax\$${tax}))+
  theme(legend.position = "right") +geom_tiplab2(size=3)
dev.off()


END
done

if test "${execute}" == "TRUE";
then
	Rscript ggtree.r
fi
