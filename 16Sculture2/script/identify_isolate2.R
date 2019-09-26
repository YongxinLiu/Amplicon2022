#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology. 2019, 37: 676-684. doi:10.1038/s41587-019-0104-4

# Version 2.0, Based on OTU table, and taxonomy, analysis OTU pure pole and select candidate well for select bacterial 2019-7-2
# v2, 基于OTU表和物种注释，筛选每个OTU的最纯及备选来源孔，以及每个孔的菌
# 190704，增加每个孔中第三列菌

# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","grid","scales","vegan","dplyr"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
library("ggplot2") # load related packages
library("grid")
library("scales")
library("vegan")
library("dplyr")

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))
# Public file 1. "result/otutab.txt"  raw reads count of each OTU in each sample
otu_table = read.delim("result/otutab.txt", row.names= 1,  header=T, sep="\t")

# Public file 2. "result/taxonomy_8.txt"  taxonomy for each OTU, tab seperated
taxonomy = read.delim("result/taxonomy_8.txt", row.names= 1,header=T, sep="\t")
# colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","evalue")
taxonomy$Full=paste(taxonomy$Phylum,taxonomy$Class,taxonomy$Order,taxonomy$Family,taxonomy$Genus,taxonomy$Species,sep = ";")
#taxonomy$full=paste(taxonomy$phylum,taxonomy$family,taxonomy$genus,sep = ";")

# 表格交叉筛选
idx = rownames(otu_table) %in% rownames(taxonomy)
otu_table = otu_table[idx,]
taxonomy = taxonomy[rownames(otu_table),]


# 按顺序抽取样品与OTU数据关系，样品稀释曲线

# 计算抽取两个样品时计算结果
#设置初始值
count=c(1:2)
# 提取前两列
temp = otu_table[,1:2]
# 筛选零的行
temp1=temp[rowSums(temp)>0,]
# 显示当前行数
count[2]=dim(temp1)[1]

# 计算抽取所有样品的结果，替换为NBT或SL更快的算法
# 样本数量，默认只算前2000个，不然太慢
sample=dim(otu_table)[2]
if (sample > 2000){
  sample=2000
}
# 循环从2至所有sample，sample为2000时间还可以接受，但到4000时就很久不动了。
for(i in 2:sample) {
  # 提取前i列
  temp = otu_table[,1:i]
  # 筛选非零的行
  temp1=temp[rowSums(temp)>0,]
  # 计算当前行数,即OTU数量
  count[i]=dim(temp1)[1]
}
rare_sample=as.data.frame(x=cbind(c(1:sample),count))
colnames(rare_sample)=c("Sample","OTU")
p = ggplot(rare_sample,aes(x=Sample,y=OTU))+geom_line()+main_theme # +geom_point()
#p
ggsave(file=paste("result/rarefraction_curve.pdf", sep=""), p, width = 5, height = 3)
# ggsave(file=paste("rarefraction_curve.png", sep=""), p, width = 5, height = 3)



# 样品中OTU种类与丰度统计
# 统计每中的前三个OTU
sample=dim(otu_table)[2]
sample_stat=as.data.frame(cbind(1:sample,colnames(otu_table)))
rownames(sample_stat)=sample_stat$V2

for(i in 1:sample) {
  temp=otu_table[order(otu_table[,i],decreasing = TRUE),c(i,1)]
  sample_stat[i,3]=temp[1,1]
  sample_stat[i,4]=temp[1,1]/sum(temp[,1])*100
  sample_stat[i,5]=rownames(temp[1,])
  sample_stat[i,6]=taxonomy[rownames(temp[1,]),]$Full
  sample_stat[i,7]=temp[2,1]
  sample_stat[i,8]=temp[2,1]/sum(temp[,1])*100
  sample_stat[i,9]=rownames(temp[2,])
  sample_stat[i,10]=taxonomy[rownames(temp[2,]),]$Full
  sample_stat[i,11]=temp[3,1]
  sample_stat[i,12]=temp[3,1]/sum(temp[,1])*100
  sample_stat[i,13]=rownames(temp[3,])
  sample_stat[i,14]=taxonomy[rownames(temp[3,]),]$Full}
# 选菌原则：同一OTU，纯度优先、数据量第二，最多选前三；输出OTU, purity-count-
sample_stat=sample_stat[,-(1:2)]
colnames(sample_stat)=c("Count1","Purity1","OTU1","Taxonomy1","Count2","Purity2","OTU2","Taxonomy2","Count3","Purity3","OTU3","Taxonomy3")
sample_stat$ID=rownames(sample_stat)
sample_stat = arrange(sample_stat, OTU1, desc(Purity1), desc(Count1))
rownames(sample_stat)=sample_stat$ID
# sample_stat=sample_stat[,-9]

write.table("Rownames\t", file=paste("result/culture_bacteria.xls",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(sample_stat, file=paste("result/culture_bacteria2.xls", sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T))



# 读取孔信息，并先菌
sample_stat = read.delim("result/culture_bacteria2.xls", row.names= 1,  header=T, sep="\t")

# 查看top1的OTU数量
print("Top1 OTU number:")
length(unique(sample_stat$OTU1))
top1=unique(sample_stat$OTU1)
# 查看top2中unique的数量
# top2=unique(sample_stat$OTU2)
# common = intersect(top1, top2)
# top2_uniq = setdiff(top2, common)


# 将竖表整理为横表选菌，先选top3，ID+纯度+丰度

## 目标菌建立列表，至少两列矩阵为数据框
otu_stat=as.data.frame(cbind(1:2,1:2))

sample=dim(sample_stat)[1]
otu="OTU_0"
j=0 # 每个otu的候选孔数量，默认3
m=0 # OTU编号，从0起，一般2-5百
for(i in 1:sample) {
#  i=1
  if (otu != as.character( sample_stat[i,"OTU1"])){
    otu=as.character( sample_stat[i,"OTU1"])
    j=0
    m=m+1
    otu_stat[m,1]=as.character(sample_stat[i,"OTU1"])
    otu_stat[m,2]=as.character(sample_stat[i,"Taxonomy1"])
  }
  if (j==0){
    #otu_stat[otu,3]=sample_stat[i,"ID"] # 用rownames索引可读到值，但是无法写入，可能是因子的问题，但用坐标可以写入
    otu_stat[m,3]=as.character(sample_stat[i,"ID"])
    otu_stat[m,4]=as.integer(sample_stat[i,"Purity1"])
    otu_stat[m,5]=as.integer(sample_stat[i,"Count1"])
  }
  if (j==1){
    otu_stat[m,6]=as.character(sample_stat[i,"ID"])
    otu_stat[m,7]=as.integer(sample_stat[i,"Purity1"])
    otu_stat[m,8]=as.integer(sample_stat[i,"Count1"])
  }  
  if (j==2){
    otu_stat[m,9]=as.character(sample_stat[i,"ID"])
    otu_stat[m,10]=as.integer(sample_stat[i,"Purity1"])
    otu_stat[m,11]=as.integer(sample_stat[i,"Count1"])
  }  
  if (j==3){
    otu_stat[m,12]=as.character(sample_stat[i,"ID"])
    otu_stat[m,13]=as.integer(sample_stat[i,"Purity1"])
    otu_stat[m,14]=as.integer(sample_stat[i,"Count1"])
  }  
  if (j==4){
    otu_stat[m,15]=as.character(sample_stat[i,"ID"])
    otu_stat[m,16]=as.integer(sample_stat[i,"Purity1"])
    otu_stat[m,17]=as.integer(sample_stat[i,"Count1"])
  }  
  j=j+1
}

## 删除OTU中非数字，并转换为数值
otu_stat$V1=as.integer(gsub("OTU_","",otu_stat$V1,perl=TRUE))
## 按OTU数字顺序排序
otu_stat = arrange(otu_stat, V1)
## 重命名行名
colnames(otu_stat)=c("ID","Taxonomy","Well1","Purity1","Count1","Well2","Purity2","Count2","Well3","Purity3","Count3","Well4","Purity4","Count4","Well5","Purity5","Count5") # 
## 写入文件
# write.table("Rownames\t", file=paste("result/culture_select.xls",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(otu_stat, file=paste("result/culture_select2.xls", sep=""), append = T, sep="\t", quote=F, row.names=F, col.names=T))

print(paste("Output in result/culture_select2.xls finished. 选菌候选表", sep = ""))
print(paste("Output in result/culture_bacteria2.xls finished. 每个孔的信息，包括丰度最高两种菌", sep = ""))
