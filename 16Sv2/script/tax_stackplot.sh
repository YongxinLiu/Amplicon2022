#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
input=result/tax/sum_
method='"p","c","o","f","g"'
design=doc/design.txt
g1=groupID
g1_list=''
output=result/tax/
number=10
width=8 
height=5
text_size=7
execute=TRUE
order=FALSE

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    tax_stackplot.sh
Version:     1.0
Date:        2018/4/7
Author:      Yong-Xin Liu
Email:       metagenome@126.com
Website:     https://blog.csdn.net/woodcorpse
Description: Stackplot of each taxonomy level
Notes:       
-------------------------------------------------------------------------------
Copyright:   2018 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
https://doi.org/10.1007/s11427-018-9284-4
-------------------------------------------------------------------------------
Version 1.0 2018/4/7
Based on usearch taxonomy summary, stackplot of each taxonomy level
# All input and output should be in default directory, or give relative or absolute path by -i/-d

# Input files: design.txt, index.txt

# 1. 实验设计 doc/design.txt, SampleID and groupID column is needed
#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	groupID	genotype
GroupAr1	ACGCTCGACA	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT
GroupAr2	ATCAGACACG	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT

# 2. 标准化物种丰度表 result/tax/sum_*.txt, calculate by usearch10 -tax_div
Phylum  ACT1KDr1        ACT1KDr10       ACT1KDr11       ACT1KDr13       ACT1KDr14   
Proteobacteria  66.3    60.1    69.8    75.1    76      69.4    75.6    64.5    34.1
Actinobacteria  10.5    21.7    12.6    11.2    8.51    11      7.2     16.5    17.5

# Output file
1. Taxonomy scatterplot: result/tax/sum_*.pdf/png, default level include "p","c","o","f","g"

OPTIONS:
	-d design for each samples, default doc/design.txt
	-e execuate Rscript, default TRUE
	-h figure height, default 5
	-i taxonomy summary, default result/tax/sum_*.txt
	-m taxonomy level, default "p","c","o","f","g"
	-o output director, default result/tax/
	-s text size, default 7
	-w figure width, default 8
	-A group name
	-B group selected list, empty will not select
	-O order of legend, default FALSE alphabet, set TRUE abundance
	-? show help of script

Example:
tax_stackplot.sh -i ${input} -m '${method}' -d ${design} -A ${g1} -B '${g1_list}' -o ${output} -O ${order} -w ${width} -h ${height}

EOF
}


# 参数解析 Analysis parameter
while getopts "c:d:e:h:i:m:n:o:s:w:A:B:O:" OPTION
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
cat <<END >script/tax_stackplot.R
#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：限制性主坐标轴分析及组间统计
# Functions: Stackplot of each taxonomy level
# Main steps: 
# - Reads taxonomy summary input.txt and design
# - Draw taxonomy stackplot by ggplot2

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


# 3. 读取输入文件

# 读取实验设计
design = read.table("${design}", header=T, row.names=1, sep="\t")
# 统一改实验列为group
design\$group = design\$${g1}

# 按实验组筛选 Select by manual set group
if ($select1){
	sub_design = subset(design, group %in% c(${g1_list}))
# 调置组排序 Set group order
	sub_design\$group  = factor(sub_design\$group, levels=c(${g1_list}))
}else{
	sub_design = design[,c("SampleID","group")]
}

# 4. 循环对每个分类级统计与绘图

method = c(${method})
for(m in method){
	# 读取usearch tax文件
	tax_sample = read.table(paste("${input}", m, ".txt", sep=""), header=T, row.names=1, sep="\t", comment.char="") 

	# 按丰度降序排序
	mean_sort = tax_sample[(order(-rowSums(tax_sample))), ]
	mean_sort = as.data.frame(mean_sort)
	# 筛选高丰度的 ${number} 类，其它归为低丰度Low abundance，可设置选择数量，即图例显示数量
	other = colSums(mean_sort[${number}:dim(mean_sort)[1], ])
	mean_sort = mean_sort[1:(${number} - 1), ]
	mean_sort = rbind(mean_sort,other)
	rownames(mean_sort)[${number}] = c("Low abundance")
	# 再次检验计算是否出错
	# colSums(mean_sort)

	# 实验设计与输入文件交叉筛选
	idx = rownames(sub_design) %in% colnames(mean_sort)
	sub_design=sub_design[idx,]
	mean_sort = mean_sort[,rownames(sub_design)]


	# 4.1 每个样品堆叠图 Stackplot for each samples

	# 保存变量备份，并输出至文件
	merge_tax=mean_sort
	write.table("Taxonomy\t", file=paste("${output}", m, "_sample.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
	suppressWarnings(write.table(merge_tax, file=paste("${output}", m, "_sample.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

	# 提取样品组信息,默认为genotype可指定
	sampFile = data.frame(sample=row.names(sub_design), group=sub_design\$group,row.names = row.names(sub_design))

	# 添加分类学列
	mean_sort\$tax = rownames(mean_sort)
	data_all = as.data.frame(melt(mean_sort, id.vars=c("tax")))
	# 按丰度顺序，默认为字母顺序 set taxonomy order by abundance, default by alphabet
	if (${order}){
		data_all\$tax  = factor(data_all\$tax, levels=rownames(mean_sort))
	}
	data_all = merge(data_all, sampFile, by.x="variable", by.y = "sample")

	p = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
		geom_bar(stat = "identity",position="fill", width=1)+ 
		scale_y_continuous(labels = scales::percent) + 
		# 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
		facet_grid( ~ group, scales = "free_x", switch = "x") +  theme(strip.background = element_blank())+
		# 关闭x轴刻度和标签
		theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
		xlab("Groups")+ylab("Percentage (%)")+ theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
	p
	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("${output}", m, "_sample.pdf", sep=""), p, width = $width, height = $height)
	ggsave(paste("${output}", m, "_sample.png", sep=""), p, width = $width, height = $height)
	print(paste("${output}", m, "_sample.pdf/png/txt finished.", sep = ""))


	# 4.2 按组均值绘制柱状图

	# 按组合并求均值
	# 转置样品名添加组名，并去除多余的两个样品列
	mat_t = t(merge_tax)
	mat_t2 = merge(sampFile, mat_t, by="row.names")
	mat_t2 = mat_t2[,c(-1,-2)]

	# 按组求均值，转置，再添加列名
	mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
	mat_mean_final = do.call(rbind, mat_mean)[-1,]
	geno = mat_mean\$group
	colnames(mat_mean_final) = geno

	# 保存变量备份，并输出至文件
	mean_sort=as.data.frame(mat_mean_final)
	write.table("Taxonomy\t", file=paste("${output}", m, "_group.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
	suppressWarnings(write.table(mean_sort, file=paste("${output}", m, "_group.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

	# 数据转换长表格并绘图
	mean_sort\$tax = rownames(mean_sort)
	data_all = as.data.frame(melt(mean_sort, id.vars=c("tax")))
	# 按丰度顺序，默认为字母顺序 set taxonomy order by abundance, default by alphabet
	if (${order}){
		data_all\$tax  = factor(data_all\$tax, levels=rownames(mean_sort))
	}

	p = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
	  geom_bar(stat = "identity",position="fill", width=0.7)+ 
	  scale_y_continuous(labels = scales::percent) + 
	  xlab("Groups")+ylab("Percentage (%)")+ theme_classic()
	if (length(unique(data_all\$variable))>3){
		p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
	}
	p

	# 保存pdf和png格式方便查看和编辑
	ggsave(paste("${output}", m, "_group.pdf", sep=""), p, width = $width, height = $height)
	ggsave(paste("${output}", m, "_group.png", sep=""), p, width = $width, height = $height)
	print(paste("${output}", m, "_group.pdf/txt finished.", sep = ""))

	print("")
}

END



# 执行脚本，脚本运行目录即工作目录(与脚本位置无关)
if test "${execute}" == "TRUE";
then
	Rscript script/tax_stackplot.R
fi
