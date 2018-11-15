# ###### manually set working directory
# setwd("/Users/Yang/Documents/R_learning/Projects/CC_manuscript/20150901_Supplementary_Figure_4/")
# getwd()
# 
# # select top 100 OTUs in viriable top100.
# 
# # 文件1：OTU，是否培养，分类学门、纲、目、科、属
# all = read.table("1.1_current_leaf_all.txt", sep="\t", header = TRUE,stringsAsFactors = F)
# head(all)
# # 文件2：前100个OTU的列表 
# OTU_top_100 = read.table("1.2_leaf_study_top100.txt", sep="\t", header = F,stringsAsFactors = F)
# head(OTU_top_100)
# length(OTU_top_100[,1])
# intersect(x = OTU_top_100[,1], y = all[,1])
# 
# top100 = data.frame(stringsAsFactors = F)
# for ( i in 1:length(all[,1])){
#     print(i)
#     if (all[i,1] %in% OTU_top_100[,1]){
#         top100 = rbind(top100,all[i,])    
#     }
#         
# }
# head(top100)
# # 筛选top100，用awk或row.names都可以快速筛选
# write.table (top100, file="2_leaf_current_study_top100.txt", sep="\t", col.names=T, row.names=F, quote=F)


top100 = read.table("result/taxonomy_8.txt", sep="\t", header = TRUE, stringsAsFactors = F)
# make dot seperated tree file
tree = data.frame(top100[,3:6], top100[,1], stringsAsFactors = F)
head(tree)
## clarify taxonomy
tree[,1] = paste("p__",tree[,1],sep = "")
tree[,2] = paste("c__",tree[,2],sep = "")
tree[,3] = paste("o__",tree[,3],sep = "")

write.table (tree, file="1_tree_plain.txt", sep=".", col.names=F, row.names=F, quote=F)

## make family annotation files
# 列出现在有门、纲、科
Phylum = unique(tree[,1]) 
#Phylum
Class = unique(tree[,2])
#Class
family = unique(tree[,4])
#family

# 筛选四大菌门中的科
family_pro = tree[tree[,1]=="p__Proteobacteria",4]
#family_pro
family_act = tree[tree[,1]=="p__Actinobacteria",4] 
#family_act
family_bac = tree[tree[,1]=="p__Bacteroidetes",4]
#family_bac
family_fir = tree[tree[,1]=="p__Firmicutes",4]
#family_fir 
#unique(family_fir)

# 对每个科进行注释
annotation_family = data.frame(stringsAsFactors = F)
for (element in family)
{
#    element
    annotation_family_ind = data.frame(stringsAsFactors = F)
    annotation_family_ind[1,1] = element
    annotation_family_ind[1,2] = "annotation"
    annotation_family_ind[1,3] = "*"
    # 设置文字旋转90度
    annotation_family_ind[2,1] = element
    annotation_family_ind[2,2] = "annotation_rotation"
    annotation_family_ind[2,3] = "90"
    # 设置背景色，四大门各指定一种色，其它为灰色
    annotation_family_ind[3,1] = element
    annotation_family_ind[3,2] = "annotation_background_color" 
    
    if (element %in% family_pro)
    {
        annotation_family_ind[3,3] = "#85F29B"
    } else if (element %in% family_act)
    {
        annotation_family_ind[3,3] = "#F58D8D"   
    } else if (element %in% family_fir)
    {
        annotation_family_ind[3,3] = "#F7C875"  
    } else if (element %in% family_bac)
    {
        annotation_family_ind[3,3] = "#91DBF6"   
    } else {
        annotation_family_ind[3,3] = "grey"   
    }
    annotation_family = rbind(annotation_family,annotation_family_ind)
    
    
#    annotation_family_ind
#    annotation_family
    
}
write.table(annotation_family, "2_annotation_family.txt", sep = "\t", quote = F,col.names = F,row.names = F, na="")



# make match anotation file
# 添加注释符号

annotation_match = data.frame(stringsAsFactors = F)

for (i in 1:length(top100[,1])){
    annotation_match_ind = data.frame(stringsAsFactors = F)
#    print(i)
    if (top100[i,2]=="YES"){
        annotation_match_ind[1,1] = top100[i,1]
        annotation_match_ind[1,2] = "ring_shape"
        annotation_match_ind[1,3] = 1
        annotation_match_ind[1,4] = "R"
        annotation_match = rbind(annotation_match, annotation_match_ind)   
    }
    
    
}

#length(annotation_match[,1])
write.table (annotation_match, file="3_annotation_match.txt", sep = "\t", quote = F,col.names = F,row.names = F, na="")



rm(list=ls())
