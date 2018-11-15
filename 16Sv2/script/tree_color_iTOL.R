#Read taxonomy_8 and color by 4 phylum
# 按result/taxonomy_8.txt文件中4大菌门着色
taxonomy = read.table("result/taxonomy_8.txt", sep="\t", header = TRUE, stringsAsFactors = F)

#Set color type
phy_color = data.frame(Phylum = c("Proteobacteria", "Actinobacteria", "Bacteroidetes", "Firmicutes", "Other"), Color = c("#85F29B", "#F58D8D", "#F7C875", "#91DBF6", "#AAAAAA" ))

# Merge add color, other set to grey
tax_color = merge(taxonomy, phy_color, by = "Phylum", all.x = T)
tax_color[is.na(tax_color$Color),]$Color = "#AAAAAA"

write.table (tax_color[,c("OTUID","Phylum","Color")], file="result/taxonomy_8_color.txt", sep="\t", col.names=F, row.names=F, quote=F)
