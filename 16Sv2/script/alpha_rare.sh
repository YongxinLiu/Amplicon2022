#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
input=result/alpha/rare.txt
design=doc/design.txt
g1=groupID
g1_list=''
output=result/alpha/
width=8 
height=5
text_size=7
execute=TRUE

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    alpha_rare.sh
Version:     1.0
Date:        2018/4/5
Author:      Yong-Xin Liu
Email:       metagenome@126.com
Website:     https://blog.csdn.net/woodcorpse
Description: Based on usearch -alpha_div_rare, draw alpha rarefracation curve
Notes:       Draw line plot in samples and groups
-------------------------------------------------------------------------------
Copyright:   2018 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
https://doi.org/10.1007/s11427-018-9284-4
-------------------------------------------------------------------------------
Version 1.0 2018/4/5
Based on usearch -alpha_div_rare and design, draw rarefracation curve line plot in samples and groups

# All input and output should be in default directory, or give relative or absolute path by -i/-d

# Input files: design.txt, index.txt

# 1. 实验设计 doc/design.txt, SampleID and groupID column is needed
#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	groupID	genotype
GroupAr1	ACGCTCGACA	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT
GroupAr2	ATCAGACACG	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT

# 2. Alpha稀释指数 result/alpha/rare.txt, calculate by usearch10 -alpha_div_rare
richness  Gr1    Gr10   Gr11    Gr12     Ga1    Ga0     Ga11    Ga12
1       74.0    95.0    90.0    116.0   74.0    84.0    69.0    105.0 
2       117.0   161.0   144.0   191.0   120.0   130.0   113.0   182.0 

# Output file
1. Rarefraction curve: result/alpha/rare_samples/gooups.pdf/png

OPTIONS:
	-d design for each samples, default doc/design.txt
	-e execuate Rscript, default TRUE
	-h figure height, default 5
	-i alpha rare file, default result/alpha/rare.txt
	-o output director, default result/alpha/
	-s text size, default 7
	-w figure width, default 8
	-A group name
	-B group selected list, empty will not select
	-? show help of script

Example:
alpha_rare.sh -i ${input} -d ${design} -A ${g1} -B '${g1_list}' -o ${output} -w ${width} -h ${height}

EOF
}


# 参数解析 Analysis parameter
while getopts "d:e:h:i:m:o:s:w:A:B:" OPTION
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
			input=$OPTARG
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
cat <<END >script/alpha_rare.R
#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：Alpha稀释曲线
# Functions: Rarefraction curve show alpha richness among samples and groups
# Main steps: 
# - Reads data table include rarefraction curve and design
# - Draw rarefraction curve of samples
# - Calculate group mean and SE (Standard Error)
# - Draw rarefraction curve of groups with SE


# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list = c("digest")
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

# 读取usearch alpha文件
rare = read.table("${input}", header=T, row.names=1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table("${design}", header=T, row.names=1, sep="\t")
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

# 实验设计与输入文件交叉筛选
idx = rownames(sub_design) %in% colnames(rare)
sub_design=sub_design[idx,]
sub_rare=rare[,rownames(sub_design)]

sampFile = as.data.frame(sub_design\$group,row.names = row.names(sub_design))
colnames(sampFile)[1] = "group"
rare=sub_rare

# 4. 绘图

# 4.1 直接展示样品 Sample
# 默认步长为1，折线不平滑，改为4减少锯齿
rare =rare[(1:25)*4,]
rare\$x = rownames(rare) # 添加x轴列
rare_melt = melt(rare, id.vars=c("x")) # 转换为长表格
rare_melt\$x = factor(rare_melt\$x, levels=1:100) # 设置x轴顺序

rare_melt3 = merge(sampFile,rare_melt, by.x="row.names", by.y="variable")
rare_melt3\$variable=rare_melt3\$Row.names

# 按样品分组，按组上色
p = ggplot(rare_melt3, aes(x = x, y = value, group = variable, color = group )) + 
  geom_line()+xlab("Rarefraction Percentage")+ylab("Richness (Observed OTUs)")+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10)+ theme_classic()
p
ggsave(paste("${output}rare", "_samples.pdf", sep=""), p, width = $width, height = $height)
ggsave(paste("${output}rare", "_samples.png", sep=""), p, width = $width, height = $height)
print("Alpha rarefraction curve: samples curve done!!!")

# 4.1 求组均值+标准误的曲线
# 求各组均值
# 读取样本筛选后的变量
rare = sub_rare
# 默认步长为1，折线不平滑，改为4减少锯齿
rare =rare[(1:25)*4,]
# 转置rare表格与实验设计合并，并去除第一列样品名
mat_t = merge(sampFile, t(rare), by="row.names")[,-1]
# 按第一列合并求均值
mat_mean = aggregate(mat_t[,-1], by=mat_t[1], FUN=mean)
# 修正行名
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean\$group
colnames(mat_mean_final) = geno

rare=as.data.frame(round(mat_mean_final))
rare\$x = rownames(rare)
rare_melt = melt(rare, id.vars=c("x"))

# 求各组标准误
# 转置rare表格与实验设计合并，并去除第一列样品名
se = function(x) sd(x)/sqrt(length(x)) # function for Standard Error
mat_se = aggregate(mat_t[,-1], by=mat_t[1], FUN=se) # se 为什么全是NA
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

rare_se=as.data.frame(round(mat_se_final))
rare_se\$x = rownames(rare_se)
rare_se_melt = melt(rare_se, id.vars=c("x"))

# 添加标准误到均值中se列
rare_melt\$se=rare_se_melt\$value

# 添加levels顺序，否则
rare_melt\$x = factor(rare_melt\$x, levels=c(1:100))

p = ggplot(rare_melt, aes(x = x, y = value, group = variable, color = variable )) + 
  geom_line()+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.5) +
  xlab("Percentage")+ylab("Richness (Observed OTUs)")+theme_classic()+
#  theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10) 
p
ggsave(paste("${output}rare", "_groups.pdf", sep=""), p, width = $width, height = $height)
ggsave(paste("${output}rare", "_groups.png", sep=""), p, width = $width, height = $height)

print("Alpha rarefraction curve: groups mean curve done!!!")

END



# 执行脚本，脚本运行目录即工作目录(与脚本位置无关)
if test "${execute}" == "TRUE";
then
	Rscript script/alpha_rare.R
fi
