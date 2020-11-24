# Version 1.0, 2017-10-12 analysis OTU and select candidate well for select bacterial
# Version 2.0, 2018-12-06 analysis OTU table and select candidate well, test rarefracation boxplot for saturation


## Parameter part:



# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","grid","scales","vegan","dplyr"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd("result") # set work directory
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
# Public file 1. "otu_table.txt"  raw reads count of each OTU in each sample
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

# Public file 2. "result/taxonomy_8.txt"  taxonomy for each OTU, tab sepepurityd
taxonomy = read.delim("taxonomy_8.txt", row.names= 1, header=T, sep="\t")
taxonomy$full=paste(taxonomy[,1],taxonomy[,2],taxonomy[,3],taxonomy[,4],taxonomy[,5],taxonomy[,6],taxonomy[,7],sep = ";")


# 饱和度曲线

## 按顺序抽取样品与OTU数据关系，样品稀释曲线

source("/mnt/bai/yongxin/github/Amplicon/16Sv2/script/stat_plot_functions.R")
# 输入OTU表，可选参数count为真阈值3，最小至最大梯度30，每个梯度重复次数30
index = sample_rare(otu_table, 1, 30, 50)
p = ggplot(index, aes(x=sample, y=richness, color=sample)) +
  geom_boxplot(show.legend = F)+ # (alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") 
  labs(x="Samples number", y=paste("Richness")) + theme_classic() + main_theme +
  # geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
  theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))
p
ggsave(paste("sample_rarefracation_boxplot", ".pdf", sep=""), p, width = 8, height = 5)


# 旧方法：计算抽取两个样品时计算结果
#设置初始值
# count=c(1:2)
# # 提取前两列
# temp = otu_table[,1:2]
# # 显示列维数
# dim(temp)[2]
# # 筛选零的行
# temp1=temp[rowSums(temp)>0,]
# # 显示当前行数
# count[2]=dim(temp1)[1]
# # 计算抽取所有样品的结果
# # 样本数量
# sample=dim(otu_table)[2]
# # 循环从2至所有sample，sample为2000时间还可以接受，但到4000时就很久不动了。
# for(i in 2:sample) {
#   # 提取前i列
#   temp = otu_table[,1:i]
#   # 筛选非零的行
#   temp1=temp[rowSums(temp)>0,]
#   # 计算当前行数,即OTU数量
#   count[i]=dim(temp1)[1]
# }
# rare_sample=as.data.frame(x=cbind(c(1:sample),count))
# colnames(rare_sample)=c("Sample","OTU")
# p = ggplot(rare_sample,aes(x=Sample,y=OTU))+geom_line()+main_theme # +geom_point()
# p
# ggsave(file=paste("rarefraction_curve.pdf", sep=""), p, width = 5, height = 3)
# ggsave(file=paste("rarefraction_curve.png", sep=""), p, width = 5, height = 3)



# 样品中OTU种类与丰度统计
# 统计每中的前三个OTU

sample_stat=as.data.frame(cbind(1:sample,colnames(otu_table)))
rownames(sample_stat)=sample_stat$V2

for(i in 1:sample) {
  temp=otu_table[order(otu_table[,i],decreasing = TRUE),c(i,1)]
  sample_stat[i,3]=temp[1,1]
  sample_stat[i,4]=temp[2,1]
  sample_stat[i,5]=temp[1,1]/sum(temp[,1])*100
  sample_stat[i,6]=temp[2,1]/sum(temp[,1])*100
  sample_stat[i,7]=rownames(temp[1,])
  sample_stat[i,8]=rownames(temp[2,])
  sample_stat[i,9]=taxonomy[rownames(temp[1,]),]$full
  sample_stat[i,10]=taxonomy[rownames(temp[2,]),]$full
}
# 选菌原则：同一OTU，纯度优先、数据量第二，最多选前三；输出OTU, purity-count-
sample_stat=sample_stat[,-(1:2)]
colnames(sample_stat)=c("count1","count2","purity1","purity2","otu1","otu2","tax1","tax2")
sample_stat$id=rownames(sample_stat)
sample_stat = arrange(sample_stat, otu1, desc(purity1), desc(count1))
rownames(sample_stat)=sample_stat$id
write.table(sample_stat, file=paste("culture_bacteria.xls", sep=""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)



# 读取孔信息，并先菌
sample_stat = read.delim("culture_bacteria.xls", row.names= 1,  header=T, sep="\t")

# 查看top1的OTU数量
print("top1 OTU number:")
length(unique(sample_stat$otu1))
top1=unique(sample_stat$otu1)
# 查看top2中unique的数量
top2=unique(sample_stat$otu2)
common = intersect(top1, top2)
top2_uniq = setdiff(top2, common)


# 将竖表整理为横表选菌，先选top3，ID+纯度+丰度

## 目标菌建立列表，至少两列矩阵为数据框
otu_stat=as.data.frame(cbind(1:2,1:2))

sample=dim(sample_stat)[1]
otu="otu0"
j=0 # 每个otu的候选孔数量，默认3
m=0 # OTU编号，从0起，一般2-5百
for(i in 1:sample) {
#  i=1
  if (otu != as.character( sample_stat[i,"otu1"])){
    otu=as.character( sample_stat[i,"otu1"])
    j=0
    m=m+1
    otu_stat[m,1]=as.character(sample_stat[i,"otu1"])
    otu_stat[m,2]=as.character(sample_stat[i,"tax1"])
  }
  if (j==0){
    #otu_stat[otu,3]=sample_stat[i,"id"] # 用rownames索引可读到值，但是无法写入，可能是因子的问题，但用坐标可以写入
    otu_stat[m,3]=as.character(sample_stat[i,"id"])
    otu_stat[m,4]=as.integer(sample_stat[i,"purity1"])
    otu_stat[m,5]=as.integer(sample_stat[i,"count1"])
  }
  if (j==1){
    otu_stat[m,6]=as.character(sample_stat[i,"id"])
    otu_stat[m,7]=as.integer(sample_stat[i,"purity1"])
    otu_stat[m,8]=as.integer(sample_stat[i,"count1"])
  }  
  if (j==2){
    otu_stat[m,9]=as.character(sample_stat[i,"id"])
    otu_stat[m,10]=as.integer(sample_stat[i,"purity1"])
    otu_stat[m,11]=as.integer(sample_stat[i,"count1"])
  }  
  if (j==3){
    otu_stat[m,12]=as.character(sample_stat[i,"id"])
    otu_stat[m,13]=as.integer(sample_stat[i,"purity1"])
    otu_stat[m,14]=as.integer(sample_stat[i,"count1"])
  }  
  if (j==4){
    otu_stat[m,15]=as.character(sample_stat[i,"id"])
    otu_stat[m,16]=as.integer(sample_stat[i,"purity1"])
    otu_stat[m,17]=as.integer(sample_stat[i,"count1"])
  }  
  j=j+1
}

## 删除OTU中非数字，并转换为数值
otu_stat$V1=as.integer(gsub("OTU_","",otu_stat$V1,perl=TRUE))
## 按OTU数字顺序排序
otu_stat = arrange(otu_stat, V1)
## 重命名行名
colnames(otu_stat)=c("ID","taxonomy","well1","purity1","count1","well2","purity2","count2","well3","purity3","count3","well4","purity4","count4","well5","purity5","count5") # 
## 写入文件
write.table(otu_stat, file=paste("culture_select.xls", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)

print(paste("Output in culture_select.xls finished. 选菌候选表", sep = ""))
print(paste("Output in culture_bacteria.xls finished. 每个孔的信息，包括丰度最高两种菌", sep = ""))

