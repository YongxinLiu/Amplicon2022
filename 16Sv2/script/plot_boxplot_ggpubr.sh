#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e 

# 默认参数 Default parameter
input=result/alpha/index.txt
method='"chao1","richness","shannon_e"'
design=doc/design.txt
g1=groupID
g1_list=''
output=result/alpha/
width=8 
height=5
text_size=7
execute=TRUE
transposition=FALSE
normlizatation=FALSE

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    plot_boxplot_ggpubr.sh
Version:     1.0
Date:        2018/5/9
Author:      Yong-Xin Liu
Email:       metagenome@126.com
Website:     https://blog.csdn.net/woodcorpse
Description: Based on usearch alph_div, draw alpha boxplot and statistics 
Notes:       Statistic alpha by aov and TukeyHSD
-------------------------------------------------------------------------------
Copyright:   2018 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
https://doi.org/10.1007/s11427-018-9284-4
-------------------------------------------------------------------------------
Version 1.0 2018/4/5
Based on usearch alph_div and design, draw alpha boxplot and statistics.

# All input and output should be in default directory, or give relative or absolute path by -i/-d

# Input files: design.txt, index.txt

# 1. 实验设计 doc/design.txt, SampleID and groupID column is needed
#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	groupID	genotype
GroupAr1	ACGCTCGACA	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT
GroupAr2	ATCAGACACG	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	GroupA	WT

# 2. Alpha指数 result/alpha/index.txt, calculate by usearch10 -alpha_div
Sample	berger_parker	buzas_gibson	chao1	dominance	equitability	jost	jost1	reads	richness
GroupAr1	0.0749	0.00262 994.1	0.968	0.633	42.1	78.3	29941.0 986.0	0.444	0.032	4.36
GroupAr10	0.0635	0.00505 1206.1	0.981	0.707	77.7	151.0	29886.0 1205.0	0.294	0.0188	5.02

# Output file
1. Alpha diversity boxplot: result/alpha/${method}.pdf/png, default method include chao1 richness shannon_e 
2. Alpha diversity statisitcs: result/alpha/${method}.txt, aov statistics, TukeyHSD test, conf.level = 0.95

OPTIONS:
	-d design for each samples, default doc/design.txt
	-e execuate Rscript, default TRUE
	-h figure height, default 5
	-i alpha diversity index, default result/alpha/index.txt
	-m index method, default "chao1","richness","shannon_e"; alpha指数种类14种 "berger_parker","buzas_gibson","chao1","dominance","equitability","jost","jost1","reads","richness","robbins","simpson","shannon_e","shannon_2","shannon_10"
	-n normalization 是否自身标准化为百分比
	-t transposition 是否转换，默认行样品列属性，OTU表需转置
	-o output director, default result/alpha/
	-s text size, default 7
	-w figure width, default 8
	-A group name
	-B group selected list, empty will not select
	-? show help of script

Example:
plot_boxplot_ggpubr.sh -i ${input} -m '${method}' -d ${design} -A ${g1} -B '${g1_list}' -o ${output} -w ${width} -h ${height}

EOF
}


# 参数解析 Analysis parameter
while getopts "d:e:h:i:m:n:o:s:t:w:A:B:" OPTION
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
		m)
			method=$OPTARG
			;;
		n)
			normlizatation=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		s)
			text_size=$OPTARG
			;;
		t)
			transposition=$OPTARG
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
cat <<END >script/plot_boxplot_ggpubr.R
#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：Alpha多样性箱线图绘制及组间统计
# Functions: Boxplot show alpha richness among groups
# Main steps: 
# - Reads data table input.txt
# - Calculate pvalue and save in output.txt
# - Draw boxplot and save in output.pdf



# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","devtools","bindrcpp",
				"ggthemes","agricolae","dplyr")
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
alpha = read.table("${input}", header=T, row.names=1, sep="\t", comment.char="") 

# 转置数据矩阵
if ($transposition){
	alpha = as.data.frame(t(alpha))
}

# 标准化数据矩阵
if ($normlizatation){
	alpha = as.data.frame(alpha/rowSums(alpha,na=T) * 100) # normalization to 1000
}
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
idx = rownames(sub_design) %in% rownames(alpha)
sub_design=sub_design[idx,]
sub_alpha=alpha[rownames(sub_design),]

# 合并Alpha指数与实验设计 add design to alpha
index = cbind(sub_alpha, sub_design) 



# 4. 循环对每种指数统计和绘图
method = c(${method})
for(m in method){
	# 统计各组间差异
	model = aov(index[[m]] ~ group, data=index)
	# 计算Tukey显著性差异检验
	Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
	# 提取比较结果
	Tukey_HSD_table = as.data.frame(Tukey_HSD\$group) 
	# 保存一个制表符，解决存在行名时，列名无法对齐的问题
	write.table(paste(m, "\n\t", sep=""), file=paste("${output}",m,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
	# 保存统计结果，有waring正常
	suppressWarnings(write.table(Tukey_HSD_table, file=paste("${output}",m,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

	# LSD检验，添加差异组字母
	out = LSD.test(model,"group", p.adj="none") # alternative fdr
	stat = out\$groups
	# 分组结果添入Index
	index\$stat=stat[as.character(index\$group),]\$groups
	# 设置分组位置为各组y最大值+高的5%
	max=max(index[,c(m)])
	min=min(index[,c(m)])
	x = index[,c("group",m)]
	y = x %>% group_by(group) %>% summarise_(Max=paste('max(',m,')',sep=""))
	y=as.data.frame(y)
	rownames(y)=y\$group
	index\$y=y[as.character(index\$group),]\$Max + (max-min)*0.05
	
	# 输出原始数据，方便筛选 
#	write.table(paste("SampleID\t", sep=""), file=paste("${output}",m,"_raw.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
#	write.table(index[,c(m,"group")], file=paste("${output}",m,"_raw.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)
#
#	p = ggplot(index, aes(x=group, y=index[[m]], color=group)) +
#		geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
#		labs(x="Groups", y=paste(m, "%")) + theme_classic() +
#		geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
#		geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
#	if (length(unique(sub_design\$group))>5){
#		p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
#	}
#	p
#	ggsave(paste("${output}", m, ".pdf", sep=""), p, width = $width, height = $height)
#	ggsave(paste("${output}", m, ".png", sep=""), p, width = $width, height = $height)

	# 采用ggpubr添加统计的箱线图
	library(ggpubr)
	p = ggboxplot(index, x="group", y="index[[m]]", color="group", palette = "jco", add = "jitter") +stat_compare_means()
	# 添加p-value, 默认是Wilcoxon test
	p
	ggsave(paste("${output}", "p_", m, ".pdf", sep=""), p, width = $width, height = $height)
	# 提示工作完成
	print(paste("Output in ${output}", m, ".txt/pdf finished.", sep = ""))
}

END



# 执行脚本，脚本运行目录即工作目录(与脚本位置无关)
if test "${execute}" == "TRUE";
then
	Rscript script/plot_boxplot_ggpubr.R
fi
