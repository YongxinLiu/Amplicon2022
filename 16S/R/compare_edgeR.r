#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>
# v1.4 2018-03-03

# Rstudio运行脚本，先设置工作目录：使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory 设置工作目录为21compare目录 或个人项目根目录

# 1.1 程序功能描述和主要步骤

# 程序功能：高通量测序reads counts值的组间比较并可绘制火山图、热图和曼哈顿图
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
    make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
                help="OTU table in counts; 原始OTU表counts值 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="doc/design.txt",
                help="Design file; 实验设计文件 [default %default]"),
    make_option(c("-t", "--taxonomy"), type="character", default="tax/taxtab.txt",
                help="Taxonomy file; 物种注释 [default %default]"),
    make_option(c("-T", "--top10tax"), type="character", default="tax_phylum.top10",
                help="Top 10 phylum; 自定义门图例 [default %default]"),    
    make_option(c("-n", "--group"), type="character", default="group",
                help="Group name; 分组列名 [default %default]"),
    make_option(c("-c", "--compare"), type="character", default="A-B",
                help="Groups comparison; 组间比较 [default %default]"),
    make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
                help="Threshold of P-value, 显著性阈值 [default %default]"),
    make_option(c("-f", "--fdr"), type="numeric", default=0.2,
                help="Threshold of FDR, 假阳性率阈值 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output prefix; 结果前缀.txt表/pdf图 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=8,
                help="Figure width; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=5,
                help="Figure heidth; 图片高 [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste("compare/edgeR_",opts$compare, sep = "")}
  
  # 显示输入输出参数，用户确认是否正确
  print("Parameters are as follows. Please check it!")
  print(paste("The input data matrix file is ", opts$input,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The taxonomy file is ", opts$taxonomy,  sep = ""))
  print(paste("Top 10 phylum file is ", opts$top10tax,  sep = ""))
  print(paste("Group name is ", opts$group,  sep = ""))
  print(paste("Group compare is ", opts$compare,  sep = ""))
  print(paste("Threshold of P-value is ", opts$pvalue,  sep = ""))
  print(paste("Threshold of FDR is ", opts$fdr,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
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

# 读取OTU表
dat = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "") 

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char = "") 

# 将选定的分组列统一命名为group
design$group=design[,opts$group]



# 4. 统计

# edgeR函数
run_edgeR = function(dat,design,compare){
  # 筛选比较组
  group_list=strsplit(opts$compare,'-')[[1]]
  idx = design$group %in% group_list
  sub_design=design[idx,]
  sub_dat=as.matrix(dat[,rownames(sub_design)])  

  d = DGEList(counts=sub_dat,group=factor(sub_design$group))
  d = calcNormFactors(d)
  # check samples is in right groups
  d$samples 

  # design.mat = model.matrix(~factor(sub_design$group))
  design.mat = model.matrix(~ 0 + factor(sub_design$group))
  rownames(design.mat)=colnames(sub_dat)
  colnames(design.mat)=levels(factor(sub_design$group))
  DAO = estimateDisp(d,design.mat)
  fit = glmFit(DAO,design.mat)
  BvsA <- makeContrasts(contrasts = opts$compare, levels=design.mat)
  lrt = glmLRT(fit,contrast=BvsA)
  # lrt = glmLRT(fit,coef=2)
  
  nrDAO=as.data.frame(topTags(lrt, n=nrow(dat)))
  nrDAO=as.data.frame(nrDAO)
  head(nrDAO)

  nrDAO$logFC=round(nrDAO$logFC,3)
  nrDAO$logCPM=round(nrDAO$logCPM,3)
  nrDAO$level = ifelse(nrDAO$logFC>0 & nrDAO$PValue<opts$pvalue & nrDAO$FDR < opts$fdr, "Enriched",ifelse(nrDAO$logFC<0 & nrDAO$PValue<opts$pvalue & nrDAO$FDR < opts$fdr, "Depleted","NotSig"))
  nrDAO$level=factor(nrDAO$level,levels = c("Enriched","Depleted","NotSig"))
  
  # 如果存在物种注释，添加注释信息
  if (file.exists(opts$taxonomy)){
  tax = read.table(opts$taxonomy, header=T, row.names= 1, sep="\t", comment.char = "") 
  tax = tax[rownames(nrDAO),]
  nrDAO=cbind(nrDAO,tax)
  
  # Draw manhattan plot and color by phylum
  x=nrDAO
  x$otu=rownames(x)
  x$neglogp=-log10(x$PValue)
  # order taxonomy
  x = arrange(x, Phylum, Class, Order, Family, Genus,otu)
  rownames(x) = x$otu
  
  # Taxonomy top 10, other in low abundance
  x$tax=x$Phylum
  if (file.exists(opts$top10)){
    top10 = read.table(opts$top10tax)
    top10 = as.vector(top10$V1)
  }else{
    top10=sort(c("Acidobacteria","Actinobacteria","Bacteroidetes",
            "Chloroflexi","Firmicutes","Proteobacteria",
            "Verrucomicrobia","Planctomycetes","Unassigned"))
  }
  print(paste("The top",length(top10)[1],"phylum as legends.", sep=" "))
  print(top10)
  print("",quote = F)
  x$tax=as.vector(x$tax)
  # levels(x$tax)=c(unique(x$tax),"Low Abundance")
  if (length(unique(x$tax)) > length(top10)){
    x[!(x$tax %in% top10),]$tax = "Low Abundance"
  }
  x$otu=factor(x$otu,levels = x$otu)
  FDR = min(x$neglogp[x$level=="Enriched"])
  x$Level=x$level
  x$Phylum=x$tax
  p = ggplot(x, aes(x=otu, y=neglogp, color=Phylum, size=logCPM, shape=Level)) +
    geom_point(alpha=.7) + 
    geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
    scale_shape_manual(values=c(17, 25, 20))+
    scale_size(breaks=c(5, 10, 15)) +
    labs(x="OTUs", y="-log10(P)", title=paste(group_list[1], "vs", group_list[1], sep=" ")) +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="right")
  p
  ggsave(paste(opts$output, "_manhattan.pdf", sep=""), p, 
         width = opts$width*2, height = opts$height)
  }
  
  # Add MeanA and MeanB in percentage
  # normlization to percentage
  norm = t(t(sub_dat)/colSums(sub_dat,na=T))*100
  # check norm is right?
  colSums(norm)
  # calculate groupA mean
  A_list = subset(sub_design, group %in% group_list[1])
  A_norm = norm[, rownames(A_list)]
  A_mean = as.data.frame(rowMeans(A_norm))
  colnames(A_mean)=c("MeanA")
  # calculate groupB mean
  B_list = subset(sub_design, group %in% group_list[2])
  B_norm = norm[, rownames(B_list)]
  B_mean = as.data.frame(rowMeans(B_norm))
  colnames(B_mean)=c("MeanB")
  # merge and reorder
  Mean = round(cbind(A_mean, B_mean, A_norm, B_norm),3)
  Mean = Mean[rownames(nrDAO),]   
  output=cbind(nrDAO[,-3],Mean)

  # write all OTU for volcano plot and manhattan plot
  write.table("OTUID\t", file=paste(opts$output,"_all.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  write.table(output,file=paste0(opts$output,"_all.txt",sep=""),append = T,quote = F,sep = '\t',row.names = T)

  # 计算上、下调OTUs数量
  NoE= dim(output[output$level=="Enriched",])[1]
  NoD= dim(output[output$level=="Depleted",])[1]
  # 绘制火山图，
  p = ggplot(output, aes(x=logFC, y=logCPM, color=level)) + 
    geom_point() + xlim(-4, 4) + theme_classic()+
    scale_colour_manual(values=c("red","green","grey")) + 
    labs(x="log2(fold change)", y="log2(count per million)", 
         title=paste(group_list[1], "vs", group_list[2], sep=" "))+ 
    annotate("text",x=-3,y=15,label=paste(NoD,sep=""))+ 
    annotate("text",x=3,y=15,label=paste(NoE,sep=""))
  p
  ggsave(paste(opts$output, "_volcano.pdf", sep=""), p, 
         width = opts$width, height = opts$height)

  # 数据筛选，pvalue < 0.05，FDR < 0.2
  output=output[output$PValue < opts$pvalue,]
  output=output[output$FDR < opts$fdr,]
  # 保存筛选结果于sig.txt结尾文件中
  write.table("OTUID\t", file=paste(opts$output,"_sig.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  write.table(output,file=paste0(opts$output,"_sig.txt",sep=""),append = T,quote = F,sep = '\t',row.names = T)
  
  
  if (file.exists(opts$taxonomy)){
  # 绘差异OTUs有分组、物种门的热图pheatmap
    
  # 制作注释行变化类型和分类学门水平的数据框
  anno_row=data.frame(Level = output$level, 
                            Taxonomy=output$Phylum, 
                            row.names = rownames(output))
  # 制作注释列分组信息
  anno_col=data.frame(Group = sub_design$group,
                      row.names = rownames(sub_design))
  
  # 绘制热图
  pheatmap(norm[rownames(output),],
           scale = "row",
           cutree_rows=2,cutree_cols = 2,
           annotation_col = anno_col, annotation_row = anno_row,
           filename = paste(opts$output, "_heatmap.pdf", sep=""),
           width=opts$width, height=opts$height, 
           annotation_names_row= T,annotation_names_col=F,
           show_rownames=F,show_colnames=T,
           main = paste("Differential abundance OTUs of",group_list[1], "vs", group_list[2],sep=" "),
           fontsize=7,display_numbers=F)
}
}

# edgeR计算组间差异
run_edgeR(dat,design,opts$compare)

# 提示工作完成
print(paste("Output in files in ", opts$output,"*.txt/pdf. All works done!!!", sep = ""))
print("",quote = F)
