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
rm(list=ls()) 

#----1.2 解析命令行#----

# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息
# 参考数据目录 /mnt/bai/yongxin/rice/integrate16s/v2OTU/LN
if (TRUE){
  # 359 / 675
  option_list <- list(
    make_option(c("-g", "--genotype"), type="character", default="../snpTtest/4m27568855.ped",
                help="Feature table [default %default]"),
    make_option(c("-p", "--phenotype"), type="character", default="pheno/Burkholderia",
                help="Feature table [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory and prefix [default %default]"),
    make_option(c("-m", "--method"), type="character", default="t.test",
                help="wilcox or t.test [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=89,
                help="Figure width [default %default]"),
    make_option(c("-v", "--height"), type="numeric", default=59,
                help="Figure height [default %default]"))
  opts <- parse_args(OptionParser(option_list=option_list))

  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste(opts$genotype, ".pdf", sep = "")}

  # 显示输入输出确认是否正确
  print(paste("Input genotype file: ", opts$genotype,  sep = ""))
  print(paste("Input phenotype file: ", opts$phenotype,  sep = ""))
  print(paste("Test method: ", opts$method,  sep = ""))
  print(paste("Output histogram PDF: ", opts$output, sep = ""))
}

# 0. 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list <- c("dplyr","ggplot2")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

main_theme = theme(panel.background=element_blank(), panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"), 
                   axis.line.y=element_line(size=.5, colour="black"), 
                   axis.ticks=element_line(color="black"), axis.text=element_text(color="black", size=7),
                   legend.position="right", legend.background=element_blank(), legend.key=element_blank(), 
                   legend.text= element_text(size=),text=element_text(family="sans", size=7),
                   plot.title = element_text(hjust = 0.5))
# 读入表型文件
phenotype <- read.table(opts$phenotype, header = FALSE, sep="\t")
pheno_id = basename(opts$phenotype)

# 读取SNP基因型 
snp.table <- read.table(opts$genotype, header = FALSE, sep=" ")

# 提取SNP ID
snpid = gsub(".ped","",basename(opts$genotype))

# 基因型和表型按第2列合并
mer <- merge(phenotype,snp.table,by.x = "V2",by.y = "V2",all=FALSE)
# 添加基因型列，为7，8列合并
mer[,snpid] <- paste(mer$V7,mer$V8,sep ="")
# 查看基因型的各类
unique(mer[,snpid])

# 排除不可知的00类型SNP，只筛选ID，丰度和基因型3列
mer <- mer[!mer[,snpid] %in% c("00"), c("V2","V3.x",snpid)]
colnames(mer) <- c("sample","abundance","snp")
# 删除缺失 
mer <- na.omit(mer)

# 统计每种基因型的数量
c <- as.data.frame(table(mer$snp))
# 计算最小等位基因频率
# maf <- round((min(c$Freq)/sum(c$Freq))*100, digits = 4)
# # 生成x轴整数序列，适合Log2转换的OTU RPM数据，为1-10的范围
# br <- seq(floor(min(mer$abundance)),ceiling(max(mer$  abundance)),by=1)
# p <- ggplot(mer,aes(x=abundance,fill=snp))+
#   geom_histogram(binwidth=0.5,color="black",breaks=br)+main_theme+
#   ggtitle(label = paste("Phenotype", basename(opts$phenotype),"distribution",sep=" "))+
#   scale_x_continuous(limits=c(floor(min(mer$abundance)),ceiling(max(mer$abundance))), breaks=br)+
#   labs(x=paste("Abundance (log2)\n",snpid,"MAF=",maf,sep = " "),y="Freq")
# p
# ggsave(opts$output, p, width = opts$width, height = opts$height, units = "mm")

# 对数2转换，与GWAS变换的一致
# mer$abundance = log2(mer$abundance)

# t检验，数据乘除10000统计结果不变
# mer$abundance = mer$abundance/12640*100


# 绘制组间差异比较
library(ggpubr)
comp_t = compare_means(abundance ~ snp, data = mer, method = opts$method)

# 结果太多，仅输出显著的
if (comp_t$p < 0.05){
  
write.table(opts$genotype, file=paste0(opts$output,opts$method, ".txt"), append = F, sep="\t", quote=F, eol = "", row.names=F, col.names=F)
suppressWarnings(write.table(comp_t, file=paste0(opts$output,opts$method,".txt"), append = T, sep="\t", quote=F, row.names=F, col.names=T))
# comp_m = compare_means(abundance ~ snp, data = mer, method = "wilcox")
# write.table(opts$genotype, file=paste0(opts$genotype,".wilcox.txt"), append = F, sep="\t", quote=F, eol = "", row.names=F, col.names=F)
# suppressWarnings(write.table(comp_m, file=paste0(opts$genotype,".wilcox.txt"), append = T, sep="\t", quote=F, row.names=F, col.names=T))

p = ggboxplot(mer, x = "snp", y = "abundance", color = "snp", palette = "jco", add = "jitter") +
  stat_compare_means(aes(label=..p.format..),method = opts$method, label.x = 1.5)+ # p.signif
  xlab(snpid)+ylab("Relative abundance(%)")+labs(title=pheno_id)+
  theme(legend.position = "none",text = element_text(family = "sans", size = 7))
p
ggsave(paste0(opts$output,opts$method,".pdf"), p, width = opts$width, height = opts$height, units = "mm")
write.table(mer, file=paste0(opts$output,opts$method, ".tsv"), append = F, sep="\t", quote=F, eol = "", row.names=F, col.names=T)
# p
# p = ggboxplot(mer, x = "snp", y = "abundance", color = "snp", palette = "jco", add = "jitter") + stat_compare_means(method = "wilcox")
# ggsave(paste0(opts$genotype,".wilcox.pdf"), p, width = opts$width, height = opts$height, units = "mm")
}

