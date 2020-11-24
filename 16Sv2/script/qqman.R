#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录
rm(list = ls())



# 1.1 程序功能描述和主要步骤

# 程序功能：基于gemma结果绘制manhattan plot
# Functions: GWAS, gemma result plot
# Main steps: 
# - Reads gemma pvaule
# - Draw manhattanplot and save in output.pdf



# 1.2 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("optparse", "reshape2","ggplot2","devtools","qqman","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 1.3 解析命令行

# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="emmax/PositiveControl.qqman",
                help="GWAS results include chr, snp, position and pvalue [default %default]"),
    make_option(c("-t", "--type"), type="numeric", default="4",
                help="Column of pvalue [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=183,
                help="Width of figure[default %default]"),
    make_option(c("-v", "--height"), type="numeric", default=59,
                help="Height of figure [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste(opts$input, sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("Column of pvalue is ", opts$type,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


# 2. 统计与绘图
pvalue = read.table(opts$input, header=T, check.names=F)

# sub_table = subset(otu_table,otu_table$p_score <=0.001,select = c(1,2,3,opts$type))

# colnames(sub_table) = c("CHR", "SNP", "BP", "P")

# 我的默认参数画法
# pdf(file=paste(opts$output,".pdf",sep=""), width=opts$width, height=opts$height) # output to PDF or screen
# manhattan(pvalue, cex = 0.6,  # main = "Manhattan Plot",  ylim = c(3, 10),
#           cex.axis = 0.9, col = c("blue4", "orange3"))
# dev.off()


# 志文的参数
# 手动设置染色体数量相同的颜色
color<-c("#87CEFA", "#0000FF", "#F08080", "#FF0000", "#90EE90", "#008000", "#DA70D6", "#800080","#FFE4B5","#FFA500", "#0000FF", "#F08080")
# 文件名中提取特征名
traitname<-gsub(".qqman", "", opts$input, perl=T)

# 输出PNG格式manhattan plot，文件名、宽、高、分辨率和长度单位
png(file=paste(opts$output,".png",sep=""),width=opts$width,height=opts$height, res=300, units="mm")
par(mar = c(5,7,3,3),mgp = c(3,1.5,0))
manhattan(pvalue, col=color, suggestiveline=F, cex=1, 
          cex.axis=1, genomewideline=-log10(1e-6), 
          main=traitname, xlab="", ylab="", cex.main=1)
mtext("Chromsome", side = 1, line = 4, cex = 1)
mtext("-lg(P)", side = 2, line = 5, cex = 1)
dev.off()

pdf(file=paste(opts$output,".pdf",sep=""),width=opts$width/25.4,height=opts$height/25.4, paper = "special", pointsize=8)
par(mar = c(5,7,3,3),mgp = c(3,1.5,0))
manhattan(pvalue, col=color, suggestiveline=F, cex=1, 
          cex.axis=1, genomewideline=-log10(1e-6), 
          main=traitname, xlab="", ylab="", cex.main=1)
mtext("Chromsome", side = 1, line = 4, cex = 1)
mtext("-lg(P)", side = 2, line = 5, cex = 1)
dev.off()


##QQplot
png(file=paste(opts$output,"QQplot_",".png",sep=""),width=opts$width/25.4,height=opts$width/25.4, res=300, units="in")
par(mar = c(5,7,3,3),mgp = c(4,1.5,0))
qq(as.numeric(pvalue$P),main=traitname,xlim=c(0,10),ylim = c(0,10),
   pch=18, col ="blue4", cex.axis=1.0, cex.lab=1.0, cex.main=1.0,las = 1)
dev.off()

pdf(file=paste(opts$output,"QQplot_",".pdf",sep=""),width=opts$width/25.4,height=opts$width/25.4)
par(mar = c(5,7,3,3),mgp = c(4,1.5,0))
qq(as.numeric(pvalue$P),main=traitname,xlim=c(0,10),ylim = c(0,10),
   pch=18, col ="blue4", cex.axis=1.0, cex.lab=1.0, cex.main=1.0, las = 1)
dev.off()
