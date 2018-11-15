#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# 写脚本scripts/beta_pcoa_group.r，筛选每组的距离，选前3输出


# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory 设置工作目录为 data (分析项目根目录)



# 1. 程序功能描述和主要步骤

# 程序功能：Beta多样性主坐标轴分析及组间统计
# Functions: PCoA analysis of samples and groups comparing
# Main steps: 
# - Reads distance matrix input.txt
# - Calculate orrdinate by PCoA and show in scatter plot
# - Adonis calculate significant between groups distance and group inner distance

# 程序使用示例
# USAGE
# # 展示样品间距离分布，统计组间是否显著，也用于异常样品筛选
# 
# Rscript ./script/beta_pcoa.r -h
# 
# # 默认基于bray_curtis距离
# Rscript ./script/beta_pcoa.r
# 
# # 完整默认参数
# Rscript ./script/beta_pcoa.r -i beta/bray_curtis.txt -t bray_curtis \
# -d doc/design.txt -n group \
# -o beta/pcoa_bray_curtis \
# -w 4 -e 2.5 
# 
# # 基于unifrac距离
# Rscript ./script/beta_pcoa.r -t unifrac
options(warn = -1)



# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T) 
}
# 解析命令行
if (TRUE){
  option_list = list(
    make_option(c("-t", "--trans"), type="logical", default=FALSE,
                help="矩阵转置; 默认TRUE [default %default]"),   
    make_option(c("-i", "--input"), type="character", default="result/alpha/index.txt",
                help="Input beta distance; 距离矩阵,默认beta目录下与t同名，可指定 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="snp/cubic_1404_hmp2plink_maf0.02.fam",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=4,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=2.5,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 有txt和矢量图pdf [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$input==""){opts$input=paste("beta/",opts$type, ".txt", sep = "")}
  if (opts$output==""){opts$output=paste(opts$input, ".fam", sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The transposition ", opts$input,  sep = ""))
  print(paste("Type transposition is ", opts$trans,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


# 2. 依赖关系检查、安装和加载

# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","vegan","dplyr")
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

# 读取距离矩阵文件
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="", stringsAsFactors = F) 
# colSums(otutab)
# 转置为样品列
if (opts$trans){
  t_otutab = t(otutab)
}else{
  t_otutab=otutab
}



# 读取实验设计
design = read.table(opts$design, header=F, sep=" ", comment.char = "", stringsAsFactors = F) 
# 只取前5列
design = design[,1:5]
rownames(design) = design$V2

# 按行名合并fam和otutab，保留fam所有行
gemma = merge(design, t_otutab, by='row.names', all.x = T)[,-1]

# 重新排序
rownames(gemma) = gemma$V2
gemma = gemma[rownames(design),]

# 写入文件和表头
write.table(gemma, file=paste(opts$output,sep = ""), append = F, sep="\t", quote=F, row.names=F, col.names=F)
write.table(colnames(gemma), file=paste(opts$output,".header",sep = ""), append = F, sep="\t", quote=F, row.names=F, col.names=F)
