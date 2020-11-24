#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
input=metadata.txt
g1=groupID
g1_list=''
output=metadata_filtered.txt
execute=TRUE

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    table_subset.sh
Version:     1.0
Date:        2020/3/5
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     https://blog.csdn.net/woodcorpse
Description: Select group by metadata
Notes:       Input table have unique row and column ID
-------------------------------------------------------------------------------
Copyright:   2020 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
    Jingying Zhang, Yong-Xin Liu, et al. 
    NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
    Nature Biotechnology. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4
-------------------------------------------------------------------------------
Version 1.0 2020/3/5


# All input and output should be in default directory, or give relative or absolute path by -i/-d

# Input files: metadata

# 1. 实验设计 metadata, SampleID and groupID column is needed
#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	groupID	genotype
GroupAr1	ACGCTCGACA	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT
GroupAr2	ATCAGACACG	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT

# Output file
1. OTUs with pvalue & FDR & fold change
2. Signifcantly abundance OTU.

OPTIONS:
	-i table, default metadata.txt
	-o filtered table, default metadata_filtered.txt
	-A group name
	-B group selected list, empty will not select
	-? show help of script

Example:
    table_subset.sh -i ${input} -A ${g1} -B ${g1_list} -o ${output}

EOF
}


# 参数解析 Analysis parameter
while getopts "e:i:o:A:B:" OPTION
do
	case $OPTION in
		e)
			execute=$OPTARG
			;;
		i)
			input=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		A)
			g1=$OPTARG
			;;
		B)
			g1_list=$OPTARG
			select1=TRUE
			;;
		?)
			usage
			exit 1
			;;
	esac
done

# 建立脚本目录
mkdir -p script

# 开始写R统计绘图脚本
cat <<END >script/table_subset.R
#!/usr/bin/env Rscript
#
# Copyright 2016-2020 Yong-Xin Liu <yxliu@genetics.ac.cn>

#----1. 准备#----

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1.1 简介和引文#----

# 程序功能：按列筛选表格
# Functions: Filter table by column
# 主要步骤 Main steps: 
# - 读取文件 Reads data matrix
# - 按列名和分组筛选 Filter by column name and value
# - 保存结果表 Save result table in *filtered.txt

# If used this script, please cited:
#    Jingying Zhang, Yong-Xin Liu, et al. 
#    NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#    Nature Biotechnology. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4

# 清空工作环境 Clean enviroment object
rm(list=ls()) 

#----1.2 加载包#----

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

#----2. 分析#----

#----2.1 读取输入文件#----

# 读取实验设计
input = read.table("${input}", header=T, row.names=1, sep="\t") # , comment.char=""
print(paste0("Input lines: ",dim(input)[1]))

#----2.2 按分组列和指定分组筛选#----

# 按实验组筛选 Select by manual set group
output = subset(input, ${g1} %in% c(${g1_list}))
print(paste0("Output lines: ",dim(output)[1]))

#----2.3 结果保存#----

write.table(paste("ID\t",sep=""), file=paste0("${output}"), append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(output, file=paste0("${output}"), append = T, quote = F, sep = '\t', row.names = T))

END

# 执行脚本，脚本运行目录即工作目录
if test "${execute}" == "TRUE";
then
	Rscript script/table_subset.R
fi
