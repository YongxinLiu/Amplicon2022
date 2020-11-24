#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#---- 1. 准备工作 #----

#----1.1 程序功能描述和主要步骤#---- 

# 程序功能：OTU按实验设计分组，计算菌值
# Functions: Calculate mean OTU abundance for tree
# Main steps: 
# - Reads data table input.txt


options(warn = -1) # Turn off warning


#----1.2 解析命令行#----

# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="pheno/alpha_richness",
                help="Feature table [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory and prefix [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=89,
                help="Figure width [default %default]"),
    make_option(c("-v", "--height"), type="numeric", default=59,
                help="Figure height [default %default]"))
  opts <- parse_args(OptionParser(option_list=option_list))

  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste(opts$input, ".pdf", sep = "")}

  # 显示输入输出确认是否正确
  print(paste("Input phenotype file: ", opts$input,  sep = ""))
  print(paste("Output histogram PDF: ", opts$output, sep = ""))
}

# 0. 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list <- c("dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#---- 1.3 读取表型 #----
pheno = read.table(opts$input, header=F, row.names= NULL, sep="\t", comment.char = "", stringsAsFactors = F)
# head(otutab)
colnames(pheno) = c("family","variety","phenotype")
pheno = na.omit(pheno)
pheno$log2 = log2(pheno$phenotype)


#----2. 统计分析#----
# 正态性检验 Shapiro-Wilk normality test，保存p-value
# https://blog.csdn.net/zzminer/article/details/8858469
library(ggpubr)
pheno$geno=as.factor("A")

#----2.1 正对照，正态随机树#----
set.seed(1)
pheno$rnorm = rnorm(dim(pheno)[1])
pvalue = shapiro.test(pheno$rnorm)
# W = 0.99495, p-value = 0.7587，正态随机数，结果为0.75是正态
p1 = gghistogram(pheno, x="rnorm", add = "mean", color = "geno", fill = "geno", rug = TRUE, palette = c("#00AFBB") )
xposition = min(pheno$rnorm)+ abs((min(pheno$rnorm) - mean(pheno$rnorm)))*0.2
p1 = p1 + annotate("text",x=xposition,y=8, label=signif(pvalue$p.value, digits = 3))

#----2.2 负对照，非正态分布的指数分布#----
pheno$nonNorm = 2^(1:dim(pheno)[1])
pvalue = shapiro.test(pheno$nonNorm)
# W = 0.99495, p-value = 0.7587，正态随机数，结果为0.75是正态
p2 = gghistogram(pheno, x="nonNorm", add = "mean", color = "geno", fill = "geno", rug = TRUE, palette = c("#00AFBB") )
xposition = min(pheno$nonNorm)+ abs((max(pheno$nonNorm) - mean(pheno$nonNorm)))*0.2
p2 = p2 + annotate("text",x=xposition,y=8, label=signif(pvalue$p.value, digits = 3))


#----2.3 表型原始值检验#----
sw =  shapiro.test(pheno$phenotype)
# print(sw)
# W = 0.97716, p-value = 0.002832
if (sw$p.value > 0.05){
  print(paste0("Phenotype raw data is normality"))
  flag = "normality"
}else{
  print(paste0("Phenotype raw data isn't normality"))
  flag = "non-normality"
}
pvalue = shapiro.test(pheno$phenotype)
# W = 0.99495, p-value = 0.7587，正态随机数，结果为0.75是正态
p3 = gghistogram(pheno, x="phenotype", add = "mean", color = "geno", fill = "geno", rug = TRUE, palette = c("#00AFBB") )
xposition = min(pheno$phenotype)+ abs((min(pheno$phenotype) - mean(pheno$phenotype)))*0.2
p3 = p3 + annotate("text",x=xposition,y=8, label=signif(pvalue$p.value, digits = 3))

pvalue_text = paste0("Phenotype Shapiro-Wilk normality test\nW\t",sw[[1]],"\np-value\t",sw[[2]],"\nConclusion\t",flag,"\n" )
write.table(pvalue_text, file=paste(opts$output, "test.txt", sep="_"), append = F, sep="\t", quote=F, row.names=F, col.names=F)

pvalue1=sw[[2]]

#----2.4 表型log2值检验#----
sw = shapiro.test(pheno$log2)
# print(sw)
# W = 0.95744, p-value = 1.353e-05, 
if (sw$p.value > 0.1){
  print(paste0("Phenotype log2(data) is normality"))
  flag = "normality"
}else{
  print(paste0("Phenotype raw data isn't normality"))
  flag = "non-normality"
}

pvalue2=sw[[2]]

pvalue = shapiro.test(pheno$log2)
# W = 0.99495, p-value = 0.7587，正态随机数，结果为0.75是正态
p4 = gghistogram(pheno, x="log2", add = "mean", color = "geno", fill = "geno", rug = TRUE, palette = c("#00AFBB") )
xposition = min(pheno$log2)+ abs((min(pheno$log2) - mean(pheno$log2)))*0.2
p4 = p4 + annotate("text",x=xposition,y=8, label=signif(pvalue$p.value, digits = 3))

write.table(paste0("Phenotype log2(data) Shapiro-Wilk normality test\nW\t",sw[[1]],"\np-value\t",sw[[2]],"\nConclusion\t",flag,"\n" ), file=paste(opts$output, "test.txt", sep="_"), append = T, sep="\t", quote=F, row.names=F, col.names=F)

# 拼图
suppressPackageStartupMessages(library(cowplot))
p0 = plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2)
ggsave(opts$output, p0, width = opts$width*2, height = opts$height * 2, units = "mm")

# 保存特征转换前后Pvalue
write.table(paste(basename(opts$input),pvalue1,pvalue2,sep="\t" ), file=paste(dirname(opts$input), "/test.txt", sep=""), append = T, sep="\t", quote=F, row.names=F, col.names=F)

