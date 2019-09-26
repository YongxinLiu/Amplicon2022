# 筛选某个指定的OTU，保留存在表达，并计算其相对丰度，排序编号 

otutab = read.delim("result/otutab.txt", row.names= 1,  header=T, sep="\t")
totutab =  as.data.frame(t(otutab))
totutab$Sum = rowSums(totutab)

# 筛选目标和总和
otu61 = totutab[,c("Sum","OTU_61")]
idx = otu61$OTU_61>0
summary(idx)
otu61 = otu61[idx,]
otu61$Percentage = otu61$OTU_61/otu61$Sum*100
otu61=otu61[order(otu61$Percentage,decreasing = T),]

## 写入文件

write.table("SampleID\t", file=paste("result/otu61.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(otu61, file=paste("result/otu61.txt", sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T))


