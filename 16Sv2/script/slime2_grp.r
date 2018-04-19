#!/usr/bin/env Rscript

# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>
# v1.4 2018-03-03

# Rstudio运行脚本，先设置工作目录：使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory 设置工作目录为21compare目录 或个人项目根目录

# 1.1 程序功能描述和主要步骤

# 程序功能：基于实验设计、OTU表和比较组，筛选；输出silme2需要的分组文件
# Functions: Calculate pvalue and FDR for each OTUs, and draw volcano plot, heatmap and manhattan plot
# Main steps: 
# - Reads data matrix and design
# - Calculate pvalue and FDR
# - Save result table in edgeR*.txt, figure in edgeR*.pdf
# 1.0 组间比较，输出差异分析表格
# 1.1 添加火山图，显示上、下调数量
# 1.2 如果存在otus.tax文件，结果文件自动添加物种注释
# 1.3 画热图，并注释分组、变化类型和门分类
# 1.4 添加曼哈顿图，门水平默认前10名着色，可自定义top10文件设置显示门类别


# 程序使用示例
# USAGE
# Default 所有参数均有默认值
# Rscript compare_edgeR.r
# 计算OE-WT
# Rscript compare_edgeR.r -c OE-WT
# 按site分组下Beijing-Hainan
# Rscript compare_edgeR.r -g site -c Beijing-Sanya
# 完整的必须参数
# Rscript compare_edgeR.r -g site -c Beijing-Sanya
#   -i otutab.txt \
#   -d design.txt \
#   -t otus.tax \
#   -d doc/design.txt \
#   -o otuput #  filename prefix for output directory name 
options(warn = -1)



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
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="temp/slime2data.txt",
                help="OTU table in counts; 原始OTU表counts值 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="doc/design.txt",
                help="Design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="groupID",
                help="Group name; 分组列名 [default %default]"),
    make_option(c("-c", "--compare"), type="character", default="IND-TEJ",
                help="Groups comparison; 组间比较 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output prefix; 结果前缀.txt表/pdf图 [default %default]")
)
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste("temp/slime2grp_",opts$compare, sep = "")}
  
  # 显示输入输出参数，用户确认是否正确
  print("Parameters are as follows. Please check it!")
  print(paste("The input data matrix file is ", opts$input,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("Group name is ", opts$group,  sep = ""))
  print(paste("Group compare is ", opts$compare,  sep = ""))
  print(paste("The output file is ", opts$output, sep = ""))
  print("",quote = F)
}



# 2. 依赖关系检查、安装和加载

# 2.1 安装CRAN来源常用包
# 依赖包列表：差异分析、绘图、热图、数据变换和开发者工具
package_list = c("limma","ggplot2","pheatmap","dplyr","devtools")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
# 基于reads counts值组间差异分析包
package_list = c("edgeR")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.3 安装Github常用包
# ggplot2美化
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
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t") 

# 提取样品组信息,默认为group可指定
design$group=design[,opts$group]

# 数据筛选，筛选两文件中共有
idx = rownames(design) %in% colnames(otutab) # match design with alpha
design = design[idx,]
otutab = otutab[,rownames(design)] 

# 按实验设计筛选
group_list=strsplit(opts$compare,'-')[[1]]
idx = design$group %in% group_list
sub_design=design[idx,]
#sub_otutab=as.matrix(otutab[,rownames(sub_design)])  
sub_design$sampleID=rownames(sub_design)

output = sub_design[,c("sampleID","group")]

write.table(output,file=paste0(opts$output,".txt",sep=""),append = F,quote = F,sep = '\t',row.names = F, col.names = F)

