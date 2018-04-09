#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
input=result/compare/
# 统计方法，默认edgeR基于负二项分布的检验，可选wilcoxon秩和检验，也叫‘Mann-Whitney’ test.
method="edgeR"
design=doc/design.txt
g1=groupID
g1_list=''
compare=doc/compare.txt
output=result/compare/
execute=TRUE
order=FALSE
pvaule=0.01
FDR=0.05
fold_change=1.3
abundance_threshold=0.0005
width=4
height=3

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    plot_volcano.sh
Version:     1.0
Date:        2018/4/9
Author:      Yong-Xin Liu
Email:       metagenome@126.com
Website:     https://blog.csdn.net/woodcorpse
Description: Group compare by edgeR or wilcon.test
Notes:       Input OTU table mustbe in raw reads counts
-------------------------------------------------------------------------------
Copyright:   2018 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
https://doi.org/10.1007/s11427-018-9284-4
-------------------------------------------------------------------------------
Version 1.0 2018/4/9
Group compare by edgeR or wilcon.test, input OTU table mustbe in raw reads counts
# All input and output should be in default directory, or give relative or absolute path by -i/-d

# Input files: design.txt, otutab.txt

# 1. 差异比较OTU，有logFC, logCPM, level三列即可
ACT2KO-Col      logFC   logCPM  PValue  FDR     level   MeanA   MeanB   ACT2KOr1        ACT2KOr
OTU_1   2.325   16.526  4.89025111048812e-21    2.60650384189017e-18    Enriched        13.208 
OTU_14  1.855   13.241  1.93349602816079e-15    5.1527669150485e-13     Enriched        1.31   

# Output file
1. volacno plot in pdf and png

OPTIONS:
	-c compare list file, default doc/compare.txt
	-d design for each samples, default doc/design.txt
	-e execuate Rscript, default TRUE
	-i OTU table in reads counts, default result/otutab.txt
	-m statistics method, default edgeR, alternative wilcon
	-o output director, default result/tax/
	-p pvaule, default 0.01
	-q FDR/qvalue, default 0.05
	-s text size, default 7
	-w figure width, default 8
	-A group name
	-B group selected list, empty will not select
	-F fold change, default 1.3
	-O order of legend, default FALSE alphabet, set TRUE abundance
	-? show help of script

Example:
compare.sh -i ${input} -m '${method}' -d ${design} -A ${g1} -B '${g1_list}' -o ${output} -O ${order} -w ${width} -h ${height}

EOF
}


# 参数解析 Analysis parameter
while getopts "c:d:e:h:i:m:n:o:p:q:s:t:w:A:B:F:O:" OPTION
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
		n)
			number=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		p)
			pvalue=$OPTARG
			;;
		q)
			FDR=$OPTARG
			;;
		s)
			text_size=$OPTARG
			;;
		t)
			abundance_threshold=$OPTARG
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
		F)
			foldchange=$OPTARG
			;;
		O)
			order=$OPTARG
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
cat <<END >script/plot_volcano.R
#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：高通量测序reads counts值的组间比较并筛选
# Functions: Calculate pvalue and FDR for each OTUs by edgeR or wilcon
# Main steps: 
# - Reads data matrix and design
# - Calculate pvalue and FDR
# - Save result table in *_all/sig.txt

# 清空工作环境 Clean enviroment object
rm(list=ls()) 


# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("limma","ggplot2","pheatmap","dplyr","devtools")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list = c("edgeR")
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


# 3. 读取输入文件

# 读取比较列表
input = read.table("${input}", header=T, row.names=1, sep="\t", comment.char="")
input\$level=factor(input\$level,levels = c("Enriched","Depleted","NotSig"))
NoE= dim(input[input\$level=="Enriched",])[1]
NoD= dim(input[input\$level=="Depleted",])[1]

# 绘制火山图
p = ggplot(input, aes(x=logFC, y=logCPM, color=level)) + 
	geom_point() + xlim(-4, 4) + theme_classic()+
	scale_colour_manual(values=c("red","green","grey")) + 
	labs(x="log2(fold change)", y="log2(count per million)")+ 
	annotate("text",x=-3,y=15,label=paste(NoD,sep=""))+ 
	annotate("text",x=3,y=15,label=paste(NoE,sep=""))
#p
suppressWarnings(ggsave(paste("${output}", "_volcano.pdf", sep=""), p, width = $width, height = $height))
suppressWarnings(ggsave(paste("${output}", "_volcano.png", sep=""), p, width = $width, height = $height))
# 提示工作完成
print(paste("Output in ${output}", "_volcano.pdf finished.", sep = ""))

END



# 执行脚本，脚本运行目录即工作目录(与脚本位置无关)
if test "${execute}" == "TRUE";
then
	mkdir -p ${output}
	Rscript script/plot_volcano.R
fi
