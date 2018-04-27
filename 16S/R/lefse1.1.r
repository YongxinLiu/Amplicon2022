#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>



# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory 设置工作目录为 data (分析项目根目录)

# 1.1 程序功能描述和主要步骤

# 程序功能：基于OTU表、物种注释和实验设计生成LefSe界门、纲、目、科、属文件
# Functions: Boxplot show alpha richness among groups
# Main steps: 
# - Reads data table input.txt
# - Calculate pvalue and save in output.txt
# - Draw boxplot and save in output.pdf

# 程序使用示例
# USAGE
# Default
# # 显示脚本帮助
# Rscript ./script/alpha_boxplot.r -h
# # 默认画richness
# Rscript ./script/alpha_boxplot.r
# # 完整参数，输出文件名默认为alpha指数类型
# Rscript ./script/alpha_boxplot.r -i alpha/alpha.txt -t richness \
# -d doc/design.txt -n group \
# -o alpha/richness \
# -w 4 -e 2.5
# # 绘制chao1
# Rscript ./script/alpha_boxplot.r -t chao1 
options(warn = -1) # Turn off warning



# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="result/otu_table.txt",
                help="Input table to read; OTU表 [default %default]"),
    make_option(c("-t", "--tax"), type="character", default="result/rep_seqs_tax.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="doc/design_indtej.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="groupID",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/lefse.txt",
                help="output directory or prefix; 输出文件前缀, 有统计表txt和矢量图pdf [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste("alpha/",opts$type, sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("Taxonomy file is ", opts$tax,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
 print(paste("The output file prefix is ", opts$output, sep = ""))
}

library(dplyr)

# 1. 读取OTU表
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "")

# 2. 读取物种注释
tax = read.table(opts$tax, header=T, row.names= 1, sep="\t",comment.char = "") 

# OTU筛选：按OTU表筛选注释
idx = rownames(tax) %in% rownames(otutab)
tax=tax[idx,]

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t")
design$group = design[,opts$group]

# OTU筛选表和实验设计
idx = rownames(design) %in% colnames(otutab)
design=design[idx,]
otutab=otutab[,rownames(design)]

# 标准化，并筛选高丰度菌0.01%
norm = t(t(otutab)/colSums(otutab,na=T))*100
colSums(norm)
idx = rowMeans(norm) > 0.1
HA = norm[idx,]
colSums(HA)

# 数据筛选
tax = tax[rownames(HA),]

# 转换为等级|连接格式
tax$Phylum=paste(tax$Kindom,tax$Phylum,sep = "|")
tax$Class=paste(tax$Phylum,tax$Class,sep = "|")
tax$Order=paste(tax$Class,tax$Order,sep = "|")
tax$Family=paste(tax$Order,tax$Family,sep = "|")
tax$Genus=paste(tax$Family,tax$Genus,sep = "|")

# 按Kindom合并
grp <- tax[rownames(tax), "Kindom", drop=F]
merge=cbind(HA, grp)
HA_Kindom = merge %>% group_by(Kindom) %>% summarise_all(sum)
colnames(HA_Kindom)[1]="Class"

# 按Phylum合并
grp <- tax[rownames(tax), "Phylum", drop=F]
merge=cbind(HA, grp)
HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
colnames(HA_Phylum)[1]="Class"

# 按Class合并
grp <- tax[rownames(tax), "Class", drop=F]
merge=cbind(HA, grp)
HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
colnames(HA_Class)[1]="Class"

# 按Order合并
grp <- tax[rownames(tax), "Order", drop=F]
merge=cbind(HA, grp)
HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
colnames(HA_Order)[1]="Class"

# 按Family合并
grp <- tax[rownames(tax), "Family", drop=F]
merge=cbind(HA, grp)
HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
colnames(HA_Family)[1]="Class"

# 按Genus合并
grp <- tax[rownames(tax), "Genus", drop=F]
merge=cbind(HA, grp)
HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
colnames(HA_Genus)[1]="Class"

# 合并6个分类级
all = rbind(HA_Kindom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus)

# 写实验组的编号为首行
colnames(all)=c("Class",as.vector(design[colnames(all)[-1],]$group ))

write.table(all, file=paste(opts$output,sep=""),append = F, quote = F, row.names = F, col.names = T, sep="\t")


# 检验结果是否正确


