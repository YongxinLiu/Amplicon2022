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
#g1=genotype
#g1_list='"WT","DM1","DM2","DO1","DO2"'
#g2=batch
#g2_list=1
group_order='FALSE' # order group by input list, not allow TRUE with merge_group
select1='FALSE' # filter first group1 info
select2='FALSE' # filter first group2 info
width=8
height=2.5
text_size=8

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    stat_16s_lib.sh
Revision:    1.0
Date:        2017/4/25
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
Draw dodge barplot by ggplot2 of each library

# All input and output should be in -o directory, or give relative -o path, or absolute path
# Input files: qc.sum

# 1. qc.sum
Library Total   Paired  Oriented        Split(Q>20)     TrimP5  TrimP3(L>300)
L1      7082981 6628222 6567319 5515917 5499119 5400354
L2      8725459 8161467 8104636 6816820 6793780 6669955
L3      10888609        10164854        10084135        7833967 7821518 7482847
L4      7330556 6690543 6639366 5553840 5546780 5480698

# Output file
1. dodge barplot of each library in pdf and png

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
	diversity.sh -o result

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




cat <<END >plot_16s_lib.r
# Install related packages
if ($ist){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","reshape2"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("${output}") # set work directory
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

sum = read.table("qc.sum", header=T, sep="\t") # dataframe
data_all = as.data.frame(melt(sum, id.vars=c("Library")))
colnames(data_all)=c("Library","Process","value")

# stat: count identify, position: stack dodge fill
p = ggplot(data_all, aes(x=Library, y = value, fill = Process ))+ 
  geom_bar(stat = "identity",position="dodge", width=0.7)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("Library")+ylab("Pair Reads count")+main_theme
p
ggsave("stat_lib_qc_sum.pdf", p, width = ${width}, height = ${height})
ggsave("stat_lib_qc_sum.png", p, width = ${width}, height = ${height})


lib_len = read.table("length.txt", header=F, sep="\t") # dataframe
colnames(lib_len)=c("Library","Count","Length")

# stat: count identify, position: stack dodge fill
p = ggplot(lib_len, aes(x=Length, y = Count, color = Library ,shapes = Library))+  geom_line() +
  xlab("Length")+ylab("Pair Reads count")+main_theme
p
ggsave("stat_lib_length.pdf", p, width = ${width}, height = ${height})
ggsave("stat_lib_length.png", p, width = ${width}, height = ${height})




END




if test "${execute}" == "TRUE";
then
	Rscript plot_16s_lib.r
fi
