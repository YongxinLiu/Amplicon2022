library(dplyr)

# 1. 读取OTU表
otutab = read.table("result/otutab.txt", header=T, row.names= 1, sep="\t", comment.char = "")

# 2. 读取物种注释
tax = read.table("tax/taxtab.txt", header=T, row.names= 1, sep="\t",comment.char = "") 

# 标准化，并筛选高丰度菌0.01%
norm = t(t(otutab)/colSums(otutab,na=T))*100
colSums(norm)
idx = rowMeans(norm) > 0.0001
HA = norm[idx,]
colSums(HA)

# 数据筛选
tax = tax[rownames(HA),]

# 转换为等级|连接格式
tax$Phylum=paste(tax$Kindom,tax$Phylum,sep = "|")
tax$Class=paste(tax$Phylum,tax$Class,sep = "|")
tax$Order=paste(tax$Class,tax$Order,sep = "|")
tax$Family=paste(tax$Order,tax$Family,sep = "|")
tax$Genus=paste(tax$Family,tax$Genus,sep = "|")

# 按Kindom合并
grp <- tax[rownames(tax), "Kindom", drop=F]
merge=cbind(HA, grp)
HA_Kindom = merge %>% group_by(Kindom) %>% summarise_all(sum)
colnames(HA_Kindom)[1]="Class"

# 按Phylum合并
grp <- tax[rownames(tax), "Phylum", drop=F]
merge=cbind(HA, grp)
HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
colnames(HA_Phylum)[1]="Class"

# 按Class合并
grp <- tax[rownames(tax), "Class", drop=F]
merge=cbind(HA, grp)
HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
colnames(HA_Class)[1]="Class"

# 按Order合并
grp <- tax[rownames(tax), "Order", drop=F]
merge=cbind(HA, grp)
HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
colnames(HA_Order)[1]="Class"

# 按Family合并
grp <- tax[rownames(tax), "Family", drop=F]
merge=cbind(HA, grp)
HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
colnames(HA_Family)[1]="Class"

# 按Genus合并
grp <- tax[rownames(tax), "Genus", drop=F]
merge=cbind(HA, grp)
HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
colnames(HA_Genus)[1]="Class"

# 合并6个分类级
all = rbind(HA_Kindom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus)

# 修改样品名为组名：删除结尾的数字
colnames(all) = gsub("\\d+$","",colnames(all),perl=TRUE)

write.table(all, file="tax/lefse.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
