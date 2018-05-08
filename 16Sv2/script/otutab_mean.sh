#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
input=result/otutab.txt
design=doc/design.txt
g1=groupID
g1_list=''
output=temp/otutab.mean
execute=TRUE

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    otutab_mean.sh
Version:     1.0
Date:        2018/5/6
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
Version 1.0 2018/5/6
Calculate mean of OTU table 计算OTU表的均值

# All input and output should be in default directory, or give relative or absolute path by -i/-d

# Input files: design.txt, otutab.txt

# 1. 实验设计 doc/design.txt, SampleID and groupID column is needed
#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	groupID	genotype
GroupAr1	ACGCTCGACA	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT
GroupAr2	ATCAGACACG	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT

# 2. 标准化物种丰度表 result/tax/sum_*.txt, calculate by usearch10 -tax_div
#OTU ID ACT1KDr1        ACT1KDr10       ACT1KDr11       ACT1KDr13   
OTU_1   6898    4153    5775    1562    4774    4346    6469    4328
OTU_10  1085    524     948     349     1000    741     1214    739 

# Output file
1. OTUs with pvalue & FDR & fold change
2. Signifcantly abundance OTU.

OPTIONS:
	-c compare list file, default doc/compare.txt
	-d design for each samples, default doc/design.txt
	-e execuate Rscript, default TRUE
	-i OTU table in reads counts, default result/otutab.txt
	-m statistics method, default edgeR, alternative wilcon
	-o output director, default temp/otutab.mean
	-p pvaule, default 0.01
	-q FDR/qvalue, default 0.05
	-s text size, default 7
	-t taxonomy file, default 7
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
if [ ${g1_list} = ""]; then
	select1=FALSE
fi

# 建立脚本目录
mkdir -p script

# 开始写R统计绘图脚本
cat <<END >script/otutab_mean.R
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

# 3. 读取输入文件

# 读取实验设计
design = read.table("${design}", header=T, row.names=1, sep="\t", comment.char="")
# 统一改实验列为group
design\$group = design\$${g1}

# 按实验组筛选 Select by manual set group
if ($select1){
	design = subset(design, group %in% c(${g1_list}))
# 调置组排序 Set group order
	design\$group  = factor(design\$group, levels=c(${g1_list}))
}

# 读取OTU表
otutab = read.table(paste("${input}", sep=""), header=T, row.names=1, sep="\t", comment.char="") 

# 实验设计与输入文件交叉筛选
#idx = rownames(design) %in% colnames(otutab)
#design = design[idx,]
#otutab = otutab[,rownames(design)]

# 按丰度值按组中位数筛选OTU
# 标准化为比例，并转置
norm = t(otutab)/colSums(otutab,na=T)*100

mean = data.frame(OTUID = rownames(otutab), Mean = round(colMeans(norm),3))
# 筛选组信
#grp = design[, "group", drop=F]
## 按行名合并
#mat_t2 = merge(grp, norm, by="row.names")
#mat_t2 = mat_t2[,-1]
## 按组求中位数
#mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=median) # mean
#mat_mean_final = do.call(rbind, mat_mean)[-1,]
#geno = mat_mean\$group
#colnames(mat_mean_final) = geno
## 按丰度按组中位数筛选
#filtered = mat_mean_final[apply(mat_mean_final,1,max) > ${abundance_threshold}, ] # select OTU at least one sample > 0.1%
#otutab = otutab[rownames(filtered),]


# write.table(paste("OTUID", "\t",sep=""), file=paste("$output",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(mean, file=paste("temp/otutab.mean", sep="\t"), append = F, quote = F, row.names = F, col.names = T, sep="\t")


END



# 执行脚本，脚本运行目录即工作目录(与脚本位置无关)
if test "${execute}" == "TRUE";
then
	Rscript script/otutab_mean.R
fi
