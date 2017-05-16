#!/bin/bash
set -e

### Default parameter
design=../doc/design.txt
execute='TRUE'
ist='FALSE' # install package, default FALSE
merge_group='FALSE' # combine group1.group2 as new group name, not allow TRUE with group_order
output='result' # default work directory, remove low abundance < .1% and p__Cyanobacteria,p__Chloroflexi
group_order='FALSE' # order group by input list, not allow TRUE with merge_group
select1='FALSE' # filter first group1 info
select2='FALSE' # filter first group2 info
width=8
height=5
text_size=7
library=L1

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    stat_16s_lib_split.sh
Revision:    1.0
Date:        2017/5/9
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: Barplot of each samples counts after split_libraries_fastq.py
Notes:       Color by design group
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
The first one , Barplot of each samples counts after split_libraries_fastq.py, color by design group. First R in /mnt/bai/yongxin/bin/plot_16s_split.r

# All input and output should be in -o directory, or give relative -o path, or absolute path
# Input files: design.txt, L1_split.count

# 1. design.txt, grouping samples, must have SampleID and group info, group1/2 can give to parameter g1 and g2, manually design
SampleID	BarcodeSequence	group	gene	batch	description
WT.1	TAGCTT	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down
WT.2	GGCTAC	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down
WT.3	CGCGCG	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down

# 2. L1_split.count
WT.4    271108
DO1.5   256951
WT.5    245736


# Output file
1. bar plot: bar_split_library.pdf .png

OPTIONS:
	-d design for each samples, default design.txt
	-e execuate Rscript, TRUE or FALSE
	-h figure height, default 2.5
	-i install package TRUE or FALSE, default FALSE
	-l library name, must
	-m merge group 1 and group 2 as new group name, default FALSE, not allow both TRUE with -g group order
	-o output director, default result_k1-c/
	-s text size, default 7
	-w figure width, default 4
	-A group1 name
	-B group1 selected list
	-C group2 name
	-D group2 selected list

Example:
	stat_16s_lib_split.sh -o result -A genotype -C batch -d ../doc/design.txt -l L1

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



cat <<END >stat_16s_lib_split.r
# Install related packages
if ($ist){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","reshape2"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("${output}")
library("ggplot2")
library("reshape2")

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
sample_count = read.table("${library}_split.count", header=F, sep="\t") # dataframe
colnames(sample_count)=c("Sample","Count")

# Set group style, single or combine
if ($merge_group){
	design\$group=paste(design\$${g1},sub_design\$${g2},sep = ".")
}else{
	design\$group=design\$${g1}
}

sample_count_group = cbind(sample_count, design[match(sample_count\$Sample, rownames(design)), ]) 

# stat: count identify, position: stack dodge fill
p = ggplot(sample_count_group, aes(x=Sample, y = Count, fill=group))+ 
  geom_bar(stat = "identity",position="dodge", width=0.7)+ 
  xlab("Library")+ylab("Pair Reads count")+main_theme+ theme(axis.text.x = element_text(angle = 90))+labs(title="Library ${library}")
# + theme(axis.text.x = element_text(size = 15, family = "myFont", color = "green", face = "bold", vjust = 0.5, hjust = 0.5, angle = 45))
p
ggsave("stat_lib_split_${library}.pdf", p, width = ${width}, height = ${height})
ggsave("stat_lib_split_${library}.png", p, width = ${width}, height = ${height})

END




if test "${execute}" == "TRUE";
then
	Rscript stat_16s_lib_split.r
fi
