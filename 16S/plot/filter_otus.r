
# Script for filter one of  OTU more than 0.1% abundance

# Input file "otu_table.txt"
# Output file "otu_id_k1.txt"
setwd("result")

thre=0.001 # 1/1000

otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")
norm = t(t(otu_table)/colSums(otu_table,na=T))
#colSums(norm) # check
k1 = norm[apply(norm,1,max) > thre, ] # select OTU at least one sample > 0.1%
dim(k1)
write.table(rownames(k1), file="otu_id_k1.txt", quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\n")

