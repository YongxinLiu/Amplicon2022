#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>



# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory 设置工作目录为 data (分析项目根目录)


# 1.1 程序功能描述和主要步骤

# 程序功能：Alpha稀释曲线
# Functions: Rarefraction curve show alpha richness among samples and groups
# Main steps: 
# - Reads data table include rarefraction curve and design
# - Draw rarefraction curve of samples
# - Calculate group mean and SE (Standard Error)
# - Draw rarefraction curve of groups with SE

# 程序使用示例
# USAGE
# Default
# # 显示帮助
# Rscript ./script/alpha_rare.r -h
# # 默认绘制4x2.5英寸的样品和组均值图
# Rscript ./script/alpha_rare.r
# # 指定输入文件和实验设计，实验组列，图片长宽和输出文件前缀
# Rscript ./script/alpha_rare.r -i alpha/alpha_rare.txt \
# -d doc/design.txt -n group \
# -o alpha/rare_ \
# -w 4 -e 2.5 
options(warn = -1)


# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T) 
}
# 参数不可用g，WARNING: unknown gui 'group', using X11 gui优先级高于R参数
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="alpha/alpha_rare.txt",
                help="Input rarefraction file to read; Usearch稀释抽样文件 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="doc/design.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=4,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=2.5,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="alpha/rare_",
                help="output directory or prefix; 输出文件前缀, 有样品和组矢量图pdf [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


# 2. 依赖关系检查、安装和加载

# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list <- c("reshape2","ggplot2")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list <- c("digest")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.3 安装Github常用包
# 参数解析、数据变换、绘图和开发包安装
package_list <- c("kassambara/ggpubr")
for(p in package_list){
  q=unlist(strsplit(p,split = "/"))[2]
  if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install_github(p)
    suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}





# 3. 读取输入文件

# 读取usearch alpha rare文件
rare = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char="") 

# 提取样品组信息,默认为group可指定
# design$genotype中genotype替换opts$group，可以写为design[,opts$group]，或design[[opts$group]]
sampFile = as.data.frame(design[,opts$group],row.names = row.names(design))
colnames(sampFile)[1] = "group"

# 4. 统计与绘图

# 统计各组间差异

# 直接展示样品 Sample
# 默认步长为1，折线不平滑，改为4减少锯齿
rare =rare[(1:25)*4,]
rare$x = rownames(rare) # 添加x轴列
rare_melt = melt(rare, id.vars=c("x")) # 转换为长表格
rare_melt$x = factor(rare_melt$x, levels=1:100) # 设置x轴顺序

rare_melt3 = merge(sampFile,rare_melt, by.x="row.names", by.y="variable")
rare_melt3$variable=rare_melt3$Row.names

# 按样品分组，按组上色
p = ggplot(rare_melt3, aes(x = x, y = value, group = variable, color = group )) + 
  geom_line()+xlab("Rarefraction Percentage")+ylab("Richness (Observed OTUs)")+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10)+ theme_classic()
p
ggsave(paste(opts$output, "samples.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, "samples.png", sep=""), p, width = opts$width, height = opts$height)

print("Alpha rarefraction curve: samples curve done!!!")



# 求组均值+标准误的曲线
# 求各组均值
# 读取usearch rarefraction文件，上面己经修改，必须重新读入
rare = read.table(opts$input, header=T, row.names= 1, sep="\t") 
# 默认步长为1，折线不平滑，改为4减少锯齿
rare =rare[(1:25)*4,]
# 转置rare表格与实验设计合并，并去除第一列样品名
mat_t = merge(sampFile, t(rare), by="row.names")[,-1]
# 按第一列合并求均值
mat_mean = aggregate(mat_t[,-1], by=mat_t[1], FUN=mean)
# 修正行名
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

rare=as.data.frame(round(mat_mean_final))
rare$x = rownames(rare)
rare_melt = melt(rare, id.vars=c("x"))

# 求各组标准误
# 转置rare表格与实验设计合并，并去除第一列样品名
se = function(x) sd(x)/sqrt(length(x)) # function for Standard Error
mat_se = aggregate(mat_t[,-1], by=mat_t[1], FUN=se) # se 为什么全是NA
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

rare_se=as.data.frame(round(mat_se_final))
rare_se$x = rownames(rare_se)
rare_se_melt = melt(rare_se, id.vars=c("x"))

# 添加标准误到均值中se列
rare_melt$se=rare_se_melt$value

# 添加levels顺序，否则
rare_melt$x = factor(rare_melt$x, levels=c(1:100))

p = ggplot(rare_melt, aes(x = x, y = value, group = variable, color = variable )) + 
  geom_line()+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.5) +
  xlab("Percentage")+ylab("Richness (Observed OTUs)")+theme_classic()+
#  theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10) 
p
ggsave(paste(opts$output, "groups.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, "groups.png", sep=""), p, width = opts$width, height = opts$height)

print("Alpha rarefraction curve: groups mean curve done!!!")
