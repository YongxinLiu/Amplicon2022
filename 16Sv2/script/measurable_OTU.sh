#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
input=result/otutab.txt
# 统计方法，默认edgeR基于负二项分布的检验，可选wilcoxon秩和检验，也叫‘Mann-Whitney’ test.
design=doc/design.txt
g1=groupID
g1_list=''
output=result/
execute=TRUE
abundance_threshold=0.01
# 矩阵自身标准化
normalization=TRUE
# 数据单位校正，默认为1不变，如万*0.01为百分比，100为百万
unit=1

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    measurable_OTU.sh
Version:     1.0
Date:        2020/1/14
Author:      Yong-Xin Liu
Email:       metagenome@126.com
Website:     https://blog.csdn.net/woodcorpse
Description: Filter OTU table by group
Notes:       Input OTU table mustbe in raw reads counts, and output reads count and norm table
-------------------------------------------------------------------------------
Copyright:   2020 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Jingying Zhang, Yong-Xin Liu, et al. 
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
Nature Biotechnology. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4
-------------------------------------------------------------------------------
Version 1.0 2020/1/14
Filter OTU table by group

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
1. OTUs table filtered in raw data
2. OTUs table filtered in percentage data

OPTIONS:
	-d design for each samples, default doc/design.txt
	-e execuate Rscript, default TRUE
	-i OTU table in reads counts, default result/otutab.txt
	-o output director, default result/
	-A group name
	-B group selected list, empty will not select
	-C group name2, alternative select for abundance
	-? show help of script

Example:
measurable_OTU.sh -i ${input} -t ${abundance_threshold} -d ${design} \
    -A genotype -B '${g1_list}' \
    -o result/

EOF
}


# 参数解析 Analysis parameter
while getopts "c:d:e:h:i:m:n:o:p:q:s:t:w:A:B:C:F:O:N:U:" OPTION
do
	case $OPTION in
		d)
			design=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		i)
			input=$OPTARG
			;;
		n)
			number=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		q)
			FDR=$OPTARG
			;;
		t)
			abundance_threshold=$OPTARG
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
		U)
			unit=$OPTARG
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

if [ ${g2} = ""]; then
	g2=${g1}
fi

# 建立脚本目录
mkdir -p script

# 开始写R统计绘图脚本
cat <<END >script/measurable_OTU.R
#!/usr/bin/env Rscript
# 
# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#    Jingying Zhang, Yong-Xin Liu, et al. 
#    NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#    Nature Biotechnology. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：高通量测序reads counts值筛选
# Functions: Calculate pvalue and FDR for each OTUs by edgeR or wilcon
# Main steps: 
# - Reads data matrix and design
# - Filter by groups
# - Save result table in result/otutab_measure.txt and otutab_measure_norm.txt

# 清空工作环境 Clean enviroment object
rm(list=ls())


# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 3. 读取输入文件

# 读取实验设计
design = read.table("${design}", header=T, row.names=1, sep="\t") # , comment.char=""
# 统一改实验列为group
design\$group = design\$${g1}

# 按实验组筛选 Select by manual set group
if ($select1){
	design = subset(design, group %in% c(${g1_list}))
# 调置组排序 Set group order
	design\$group  = factor(design\$group, levels=c(${g1_list}))
}

# 读取OTU表
otutab = read.table(paste("${input}", sep=""), header=T, row.names=1, quote = "", sep="\t", comment.char="") 
# Show features (include OTU/ASV/Taxonomy) numbers
print("Total features number")
print(dim(otutab)[1])

# 实验设计与输入文件交叉筛选
idx = rownames(design) %in% colnames(otutab)
design = design[idx, , drop = F]
otutab = otutab[,rownames(design)]

# 按丰度值按组中位数筛选OTU
# 标准化为百分比例，并转置
if (${normalization}){
	norm = t(otutab)/colSums(otutab,na=T)*100
}else{
	# 非标准化时为，默认抽样10000，除以100标准为百分比
	norm=t(otutab)*${unit}
}
# 检查样本标准化后是否为100
# rowSums(norm)

# 筛选组，按组求中位数
# grp = design[, "${g2}", drop=F] # 需要按第二条件筛选时使用
grp = design[, "${g2}", drop=F]
# 按行名合并
mat_t2 = merge(grp, norm, by="row.names")
mat_t2 = mat_t2[,-1]
# 按组求中位数，中位数筛选更有效去除异常值
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=median) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean\$group
colnames(mat_mean_final) = geno
# 按丰度按组中位数筛选
filtered = mat_mean_final[apply(mat_mean_final,1,max) >= ${abundance_threshold}, ] # select OTU at least one sample > 0.1%
otutab = otutab[rownames(filtered),]

# 按均值输出和保存相对丰度，汇总才为真实丰度
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean\$group
colnames(mat_mean_final) = geno

# 保存筛选后的OTU表
mat_mean_high = mat_mean_final[rownames(filtered),]
write.table(paste("OTUID\t",sep=""), file=paste("${output}", "otutab_measurable.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(otutab, file=paste("${output}", "otutab_measurable.txt",sep=""), append = T, quote = F, sep = '\t', row.names = T, col.names = T))

# 保存相对丰度
norm = as.data.frame(t(norm))
norm = norm[rownames(filtered),]

write.table(paste("OTUID\t",sep=""), file=paste("${output}", "otutab_measurable_norm.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(round(norm,5), file=paste("${output}", "otutab_measurable_norm.txt",sep=""), append = T, quote = F, sep = '\t', row.names = T, col.names = T))

print("Selected high abundance OTUs number")
print(dim(mat_mean_high)[1])

print(colSums(mat_mean_high))

print(mean(colSums(mat_mean_high)))

END


# 执行脚本，脚本运行目录即工作目录(与脚本位置无关)
if test "${execute}" == "TRUE";
then
	Rscript script/measurable_OTU.R
fi
