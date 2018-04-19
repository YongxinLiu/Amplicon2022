#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>



# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory 设置工作目录为 data (分析项目根目录)



# 1. 程序功能描述和主要步骤

# 程序功能：样品和组的堆叠柱状图
# Functions: Stackplot of each taxonomy level
# Main steps: 
# - Reads distance matrix input.txt
# - Calculate orrdinate by PCoA and show in scatter plot
# - Adonis calculate significant between groups distance and group inner distance

# 程序使用示例
# USAGE
# Default
# # 显示帮助，主要是参数说明
# Rscript ./script/tax_stackplot.r -h
# 
# # 默认按phylum和前8类展示, 4X2.5
# Rscript ./script/tax_stackplot.r 
# # Legend too long, main text overlap. Increase figure size.
# 
# # 按目前10，图片宽6 x 4
# Rscript ./script/tax_stackplot.r -t order \
# -b 10 -w 6 -e 4 #  filename prefix for output directory name 
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
    make_option(c("-t", "--type"), type="character", default="p",
                help="Taxonomy level p c o f g; 分类学级别, 门phylum 纲class 目order 科family 属genus [default %default]"),
    make_option(c("-i", "--input"), type="character", default="",
                help="Merged taxonomy file; 分类学合并结果 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="doc/design.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-b", "--number"), type="numeric", default=8,
                help="Number taxonomy for showing; 展示分类数量 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=4,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=2.5,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 通常会有统计表txt、矢量图pdf和位图png [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$input==""){opts$input=paste("tax/sum_",opts$type,".txt", sep = "")}
  if (opts$output==""){opts$output=paste("tax/stackplot_",opts$type,"",sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("Merged taxonomy file is ", opts$input,  sep = ""))
  print(paste("Taxonomy level is ", opts$type,  sep = ""))
  print(paste("Number taxonomy for showing is ", opts$number,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


# 2. 依赖关系检查、安装和加载

# 2.1 安装CRAN来源常用包

# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","vegan")
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

# 读取样品分类学文件
tax_sample = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char="") 

# 提取样品组信息,默认为genotype可指定
sampFile = data.frame(group=design[,opts$group],
                      sample=row.names(design), 
                      row.names = row.names(design))

# 数据筛选，筛选两文件中共有
idx = rownames(sampFile) %in% colnames(tax_sample) # match design with alpha
sampFile = sampFile[idx,]
tax_sample = tax_sample[,rownames(sampFile)] 
# 检查，每组是否标准化为100%
# colSums(tax_sample)



# 4. 统计与绘图

# 4.1 每个样品堆叠图 Stackplot for each samples

# 按丰度降序排序
mean_sort = tax_sample[(order(-rowSums(tax_sample))), ]
mean_sort = as.data.frame(mean_sort)

# 筛选前7类，其它归为other，可设置不同组数
other = colSums(mean_sort[opts$number:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(opts$number - 1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[opts$number] = c("Low abundance")
# 再次检验计算是否出错
# colSums(mean_sort)

# 保存变量备份，并输出至文件
merge_tax=mean_sort
write.table("\t", file=paste(opts$output,"_sample.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(merge_tax, file=paste(opts$output,"_sample.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)

# 添加分类学列
mean_sort$tax = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("tax")))
# data_all$tax  = factor(data_all$tax, levels=rownames(mean_sort))   # set taxonomy order
data_all = merge(data_all, sampFile, by.x="variable", by.y = "sample")

p = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=1)+ 
  scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ group, scales = "free_x", switch = "x") +  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Percentage (%)")+ theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
p

# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, "_sample.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, "_sample.png", sep=""), p, width = opts$width, height = opts$height)
print(paste(opts$output, "_sample.pdf/txt finished.", sep = ""))



# 4.2 按组均值绘制柱状图

# 按组合并求均值

# 转置样品名添加组名，并去除多余的两个样品列
mat_t = t(merge_tax)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,c(-1,-3)]

# 按组求均值，转置，再添加列名
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

# 保存变量备份，并输出至文件
mean_sort=as.data.frame(mat_mean_final)
write.table("\t", file=paste(opts$output,"_group.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(merge_tax, file=paste(opts$output,"_group.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)

# 数据转换长表格并绘图
mean_sort$tax = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("tax")))
# 设置分类学顺序，默认字母，可选丰度或手动
# data_all$tax  = factor(data_all$tax, levels=rownames(mean_sort))   

p = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+ theme_classic()
p

# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, "_group.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, "_group.png", sep=""), p, width = opts$width, height = opts$height)
print(paste(opts$output, "_group.pdf/txt finished.", sep = ""))


# 5. 保存图表

# 提示工作完成

print("Taxonomy stackplot in sample done!!!")
