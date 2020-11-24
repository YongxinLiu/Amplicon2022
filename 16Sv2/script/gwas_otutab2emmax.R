#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：转换OTU表为emmax的fam格式
# Functions: format otu table to emmax
# Main steps: 
# - OTU table: otutab.txt
# - OTU list: otuid
# - fam header: emmax/snp.tfam


options(warn = -1) # Turn off warning


# 1.2 解析命令行
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
    make_option(c("-i", "--input"), type="character", default="result/otutab_norm.txt",
                help="Feature table [default %default]"),
    make_option(c("-l", "--list"), type="character", default="pheno/otuid",
                help="OTU list [default %default]"),
    make_option(c("-f", "--fam"), type="character", default="../emmax/snp.tfam",
                help="Emmax fam example [default %default]"),
    make_option(c("-o", "--output"), type="character", default="pheno/",
                help="output directory and prefix [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))

  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste0(opts$input)}
}
print(opts)

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

#---- 1. 读取表型 #----
pheno = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
pheno = as.data.frame(t(pheno))

#---- 2. 读取ID列表，没有则跳过 #----
if (file.exists(opts$list)){
  otuid = read.table(opts$list, header=F, row.names= NULL, sep="\t", comment.char = "", stringsAsFactors = F)
  idx = colnames(pheno) %in% otuid$V1
  pheno = pheno[, idx]
}

#---- 3. 读取表型文件 #----
fam = read.table(opts$fam, header=F, row.names= NULL, sep="\t", comment.char = "", stringsAsFactors = F)
if (dim(fam)[2]<3){
  fam = read.table(opts$fam, header=F, row.names= NULL, sep=" ", comment.char = "", stringsAsFactors = F)
}
fam = fam[,1:2]
rownames(fam) = fam$V2


#---- 4. 批量生成表型文件 #----
pheno = cbind(rownames(pheno), pheno)
colnames(pheno)[1] = "ID"

merged = merge(fam, pheno, by.x = 'V2', by.y = 'ID', all.x = T)

for (i in 3:dim(merged)[2]){
  # print(i)
  # i=3
  otu_fam = merged[, c(1,2,i)]
  write.table(otu_fam, paste0(opts$output, colnames(otu_fam)[3]),
              col.names=F, row.names=F, quote=F, sep = "\t")
  # log2转换
  otu_fam[,3] = log2(otu_fam[,3])
  write.table(otu_fam, paste0(opts$output, "log2.", colnames(otu_fam)[3]),
              col.names=F, row.names=F, quote=F, sep = "\t")  
}
