#!/bin/bash
set -e

### Default parameter
design=../doc/design.txt
execute='TRUE'
ist='FALSE' # install package, default FALSE
merge_group='FALSE' # combine group1.group2 as new group name, not allow TRUE with group_order
output='result/solit' # default work directory, remove low abundance < .1% and p__Cyanobacteria,p__Chloroflexi
group_order='FALSE' # order group by input list, not allow TRUE with merge_group
select1='FALSE' # filter first group1 info
select2='FALSE' # filter first group2 info
width=40
height=5
text_size=1
library=L1
g1=plate
g2=batch

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    stat_16s_lib_split2.sh
Revision:    2.0
Date:        2019/7/1
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
Version 1.1 2017/6/26
Add rscript library ID, debug rewrite script when batch run; remove rscript after run
Version 2.0 2019/7/1


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
	-o output director, default result/
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



cat <<END >script/stat_16s_lib_split_${library}.r
# Install related packages
if ($ist){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","reshape2","pheatmap"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("${output}")
library("ggplot2")
library("reshape2")
library("pheatmap")

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
sample_count = read.table("${library}.txt", header=F, sep="\t") # dataframe
colnames(sample_count)=c("Sample","Count")

# Set group style, single or combine
if ($merge_group){
	design\$group=paste(design\$${g1},design\$${g2},sep = "")
}else{
	design\$group=design\$${g1}
}

sample_count_group = cbind(sample_count, design[match(sample_count\$Sample, rownames(design)), ]) 

# stat: count identify, position: stack dodge fill
p = ggplot(sample_count_group, aes(x=Sample, y = Count, fill=group))+ 
  geom_bar(stat = "identity",position="dodge", width=0.7)+ 
  xlab("Library")+ylab("Pair Reads count")+main_theme+ theme(axis.text.x = element_text(angle = 90))+labs(title="Library ${library}")
# + theme(axis.text.x = element_text(size = 15, family = "myFont", color = "green", face = "bold", vjust = 0.5, hjust = 0.5, angle = 45))
#p
ggsave("stat_lib_split_${library}.pdf", p, width = ${width}, height = ${height})
ggsave("stat_lib_split_${library}.png", p, width = ${width}, height = ${height})

## 绘制板热图，只画前第一个示例
#plate_count = read.table("plate/${library}P01.plate", header=T, row.names= 1,sep="\t") # dataframe
#pheatmap(plate_count, cellwidth=30, cellheight=30, scale="none", cluster_rows=F, cluster_cols=F, labels_col=c(1:12), display_numbers=T, number_format="%.f", filename = "${library}P01.plate.pdf", width=9, height=5)  
#
## 绘制胶热图，只画前第一个示例
#gel_count = read.table("plate/${library}P01.gel", header=T, row.names= 1,sep="\t") # dataframe
#pheatmap(gel_count, cellwidth=20, cellheight=20, scale="none", cluster_rows=F, cluster_cols=F, labels_col=c(1:12), display_numbers=T, number_format="%.f", filename = "${library}P01.gel.pdf", width=9, height=5)  
#
### 根据变量"Count"作传统的直方图，窗宽为0.5 
#p=ggplot(sample_count, aes(x=Count)) + geom_histogram(binwidth=10) 
### 下面这行代码也能达到同样的效果： 
## qplot(sample_count$Count, binwidth=.5) 
#ggsave("stat_lib_density_all_${library}.pdf", p, width = 8, height = 4.5)
#
## 白色填充，黑色轮廓 
#ggplot(sample_count, aes(x=Count)) + 
#  geom_histogram(binwidth=10, colour="red", fill="white") 
#
## 密度曲线 
##ggplot(sample_count, aes(x=Count)) + geom_density() 
#
## 带核密度曲线的概率分布直方图 
#p=ggplot(sample_count, aes(x=Count)) + 
#  geom_histogram(aes(y=..density..), # y轴刻度对应的为概率密度而非计数 
#                 binwidth=10, 
#                 colour="red", fill="white") + xlim(20, 1200)+
#  geom_density(alpha=.3, fill="green") # 设置附加半透明的密度图 
#ggsave("stat_lib_density_20_1200_${library}.pdf", p, width = 8, height = 4.5)
#
#
#
## 总样本数
#length(sample_count\$Count)
#
## 大于20个reads的样本数
#norm=sample_count\$Count[sample_count\$Count>=20]
#length(norm)
## 大于40个reads的样本数
#norm=sample_count\$Count[sample_count\$Count>=40]
#length(norm)


END




if test "${execute}" == "TRUE";
then
	Rscript script/stat_16s_lib_split_${library}.r
#	rm script/stat_16s_lib_split_${library}.r
fi
