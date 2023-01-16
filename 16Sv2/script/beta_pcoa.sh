#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
compare=doc/compare.txt
input=result/beta/
method='"binary_jaccard","bray_curtis","unweighted_unifrac","weighted_unifrac"'
design=doc/design.txt
g1=groupID
g1_list=''
output=result/beta/
width=8 
height=5
text_size=7
execute=TRUE
ellipse=TRUE

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    beta_pcoa.sh
Version:     1.2
Date:        2018/4/6
Author:      Yong-Xin Liu
Email:       metagenome@126.com
Website:     https://blog.csdn.net/woodcorpse
Description: Based on distance matrix and design, draw beta pcoa and statistics 
Notes:       adnois: Permutational Multivariate Analysis of Variance Using Distance Matrices
-------------------------------------------------------------------------------
Copyright:   2016-2022 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8
-------------------------------------------------------------------------------
Version 1.0 2018/4/5
Based on QIIME beta_diversity.py and design, draw beta pcoa and statistics.
Version 1.1 2018/5/8
Output pcoa coordinate files
Version 1.2 2019/4/10
All comment is a new line, can filter by #


# All input and output should be in default directory, or give relative or absolute path by -i/-d

# Input files: design.txt, index.txt

# 1. 实验设计 doc/design.txt, SampleID and groupID column is needed
#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	groupID	genotype
GroupAr1	ACGCTCGACA	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT
GroupAr2	ATCAGACACG	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT

# 2. 距离矩阵 result/beta/index.txt, calculate by usearch10 -beta_div
        ACT1KDr1        ACT1KDr10       ACT1KDr11       ACT1KDr13       ACT1KDr14       ACT1KDr15       ACT1KDr16     
ACT1KDr1        0.0     0.3368211677    0.273395844824  0.456519556998  0.273334781229  0.290784172782  0.289115191987
ACT1KDr10       0.3368211677    0.0     0.230635393315  0.333813359232  0.281763151734  0.246929692467  0.30360096917 

# Output file
1. PCoA scatterplot: result/beta/${method}.pdf/png, default method include "binary_jaccard","bray_curtis","unweighted_unifrac","weighted_unifra"
2. Alpha diversity statisitcs: result/beta/${method}.stat, pvalue of adonis

OPTIONS:
	-c two groups compare list, default doc/compare.txt
	-d design for each samples, default doc/design.txt
	-e execuate Rscript, default TRUE
	-h figure height, default 5
	-i beta diversity index, default result/beta/index.txt
	-m index method, default "chao1","richness","shannon_e"; beta指数种类14种 "berger_parker","buzas_gibson","chao1","dominance","equitability","jost","jost1","reads","richness","robbins","simpson","shannon_e","shannon_2","shannon_10"
	-o output director, default result/beta/
	-s text size, default 7
	-w figure width, default 8
	-A group name 组名
	-B group selected list, empty will not select 组列表
	-E add ellipse for each group
	-? show help of script

Example:
beta_pcoa.sh -i ${input} -m '${method}' -d ${design} -A ${g1} -B '${g1_list}' -o ${output} -w ${width} -h ${height}

EOF
}


# 参数解析 Analysis parameter
while getopts "c:d:e:h:i:m:o:s:w:A:B:E:" OPTION
do
	case $OPTION in
		c)
			compare=$OPTARG
			;;
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
			input=$OPTARG
			;;
		m)
			method=$OPTARG
			;;
		o)
			output=$OPTARG
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
		E)
			ellipse=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

# 当选择列表为空时，关闭实验设计筛选
if [ ${g1_list} == ""]; then
select1=FALSE
fi

# 建立脚本目录
mkdir -p script

# 开始写R统计绘图脚本
cat <<END >script/beta_pcoa.R
#!/usr/bin/env Rscript
# 
# Copyright 2016-2022 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：Beta多样性主坐标轴分析及组间统计
# Functions: PCoA analysis of samples and groups comparing
# Main steps: 
# - Reads distance matrix input.txt
# - Calculate orrdinate by PCoA and show in scatter plot
# - Adonis calculate significant between groups distance and group inner distance

# 清空工作环境 Clean enviroment object
rm(list=ls()) 


# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","vegan")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list = c("digest","ggrepel")
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		source("https://bioconductor.org/biocLite.R")
		biocLite(p)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
	}
}

# 2.3 安装Github常用包
# 参数解析、数据变换、绘图和开发包安装
package_list = c("kassambara/ggpubr")
for(p in package_list){
	q=unlist(strsplit(p,split = "/"))[2]
	if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install_github(p)
		suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
	}
}

# 函数1：采用adonis对距离矩阵进行组间差异统计
# Compare each group distance matrix by vegan adonis in bray_curtis
da_adonis = function(sampleV){
	sampleA = as.matrix(sampleV\$sampA)
	sampleB = as.matrix(sampleV\$sampB)
	design2 = subset(sub_design, group %in% c(sampleA,sampleB))
	if (length(unique(design2\$group))>1) {
		sub_dis_table = dis_table[rownames(design2),rownames(design2)]
		sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
		adonis_table = adonis(sub_dis_table ~ group, data = design2, permutations = 1000) 
		adonis_pvalue = adonis_table\$aov.tab\$\`Pr(>F)\`[1]
		print(paste("In", m, "pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
		adonis_pvalue = paste(m, sampleA, sampleB, adonis_pvalue, sep="\t")
		return(adonis_pvalue)
	}
}


# 3. 读取输入文件

# 读取实验设计
design = read.table("${design}", header=T, row.names=1, sep="\t", quote = "")
# 统一改实验列为group
design\$group=design\$${g1}

# 按实验组筛选 Select by manual set group
if ($select1){
	sub_design = subset(design, group %in% c(${g1_list}))
# 调置组排序 Set group order
	sub_design\$group  = factor(sub_design\$group, levels=c(${g1_list}))
}else{
	sub_design = design
}



# 4. 循环对每种距离矩阵统计和绘图


method = c(${method})
for(m in method){
	# 读取usearch beta文件
	beta = read.table(paste("${input}",m,".txt",sep=""), header=T, row.names=1, sep="\t", comment.char="", quote = "") 
	write.table(paste("#Permutational Multivariate Analysis of Variance Using Distance Matrices\nDistanceMatrices\tGroupA\tGroupB\tP-value",  sep=""), 
				file=paste("$output", m, ".stat", sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=F)

	# 实验设计与输入文件交叉筛选
	idx = rownames(sub_design) %in% rownames(beta)
	sub_design=sub_design[idx,]
	sub_beta=beta[rownames(sub_design),rownames(sub_design)]

	# vegan:cmdscale计算矩阵矩阵中主坐标轴坐标，取前4维
    # # k is dimension, 3 is recommended; eig is eigenvalues
	pcoa = cmdscale(sub_beta, k=4, eig=T) 
    # get coordinate string, format to dataframe
	points = as.data.frame(pcoa\$points) 
	eig = pcoa\$eig
	points = cbind(points, sub_design\$group)
	colnames(points) = c("PC1", "PC2", "PC3", "PC4","group") 
	write.table("Samples\t", file=paste("${output}", m, "14.txt",sep = ""), append = F, sep="\t", quote=F,  eol = "",row.names=F, col.names=F)
	suppressWarnings(write.table(points[,c(1:4)], file=paste("${output}", m, "14.txt",sep = ""), append = T, sep="\t", quote=F, row.names=T, col.names=T))

	# plot PC 1 and 2
	p = ggplot(points, aes(x=PC1, y=PC2, color=group)) + geom_point(alpha=.7, size=2) +
		labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
		y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
		title=paste(m," PCoA",sep="")) + theme_classic()

	# 是否添加置信椭圆
	if (${ellipse}){
		p = p + stat_ellipse(level=0.68)
	}
	p
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("${output}", m, ".pdf", sep=""), p, width = $width, height = $height)
	ggsave(paste("${output}", m, ".png", sep=""), p, width = $width, height = $height)
	# 提示工作完成
	print(paste("Output in ${output}", m, ".pdf finished.", sep = ""))

	# 添加样品标签
	p = p + geom_text_repel(label = paste(rownames(points)), colour="black", size=3)
	p
	# 保存pdf和png格式方便查看和编辑
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("${output}", m, "_label.pdf", sep=""), p, width = $width*2, height = $height*2)
#	ggsave(paste("${output}", m, "_label.png", sep=""), p, width = $width, height = $height)
	# 提示工作完成
	print(paste("Output in ${output}", m, "_label.pdf finished.", sep = ""))

	# plot PC 3 and 4
	p = ggplot(points, aes(x=PC3, y=PC4, color=group)) + geom_point(alpha=.7, size=2) +
		labs(x=paste("PCo 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""),
		y=paste("PCo 4 (", format(100 * eig[4] / sum(eig), digits=4), "%)", sep=""),
		title=paste(m," PCoA",sep="")) + theme_classic()
	# 是否添加置信椭圆
	if (${ellipse}){p = p + stat_ellipse(level=0.68)}
	p
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("${output}", m, "_34.pdf", sep=""), p, width = $width, height = $height)
#	ggsave(paste("${output}", m, "_34.png", sep=""), p, width = $width, height = $height)
	# 提示工作完成
	print(paste("Output in ${output}", m, "_34.pdf finished.", sep = ""))

	# 循环统计各比较组或所有组 loop for each group pair
	dis_table = as.matrix(sub_beta)
	# 如果没有比较文件，则自动全循环
	if (!file.exists("${compare}")) {
		compare_data = as.vector(unique(sub_design\$group))
		len_compare_data = length(compare_data)
		for(i in 1:(len_compare_data-1)) {
			for(j in (i+1):len_compare_data) {
				tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
#				print(tmp_compare)
				adonis_pvalue = da_adonis(tmp_compare)
				write.table(adonis_pvalue, file=paste("$output", m, ".stat", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
			}
		}
	# 有比较文件，按设计比较
	}else {
		compare_data = read.table("${compare}", sep="\t", check.names=F, quote='', comment.char="")
		colnames(compare_data) = c("sampA", "sampB")
		for(i in 1:dim(compare_data)[1]){
			adonis_pvalue = da_adonis(compare_data[i,])
			write.table(adonis_pvalue, file=paste("$output", m, ".stat", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
		}
	}	 
	print(paste("Adnois statistics result in ", "$output", m, ".stat\n is finished.", sep = ""))
	print("")
}

END



# 执行脚本，脚本运行目录即工作目录(与脚本位置无关)
if test "${execute}" == "TRUE";
then
	Rscript script/beta_pcoa.R
fi
