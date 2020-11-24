#!/usr/bin/env Rscript
# 
# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. 
# NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. 
# Nature Biotechnology. 2019, 37: 676-684.
# https://doi.org/10.1038/s41587-019-0104-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录
rm(list = ls())



#---- 1. 准备工作 #----

#----1.1 程序功能描述和主要步骤#---- 

# 程序功能：特征表排序
# Functions: Sort feature table
# Main steps: 
# - Read feature table
# - Sort
# - Write



#----1.2 安装CRAN来源常用包#----

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("optparse", "reshape2","ggplot2","devtools","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#----1.3 解析命令行#----

# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/tax/sum_g.txt",
                help="Feature table [default %default]"),
    # make_option(c("-t", "--type"), type="numeric", default="4",
    #             help="Column of pvalue [default %default]"),
    # make_option(c("-w", "--width"), type="numeric", default=183,
    #             help="Width of figure[default %default]"),
    # make_option(c("-v", "--height"), type="numeric", default=59,
    #             help="Height of figure [default %default]"),
    make_option(c("-s", "--sort"), type="logical", default=T,
                help="default decrease, F as increase [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste(opts$input,".sort", sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("Decreasing sort is ", opts$sort,  sep = ""))
  print(paste("The output file file is ", opts$output, sep = ""))
}


#---- 2. 统计与绘图 #----

#----2.1 读取OTU表#----

otutab = read.table(paste(opts$input, sep=""), header=T, row.names=1, quote = "", sep="\t", comment.char="") 

#----2.2 排序#----
idx = order(rowSums(otutab), decreasing = opts$sort)
otutab = otutab[idx,]

#----2.3 输出结果#----
write.table(paste("#OTU\t",sep=""), file=paste(opts$output, sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(otutab, file=paste(opts$output, sep=""), append = T, quote = F, sep = '\t', row.names = T))

