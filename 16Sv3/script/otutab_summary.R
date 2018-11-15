# otutab_summary

# read raw read counts of otu table
otutab = read.table("otutab.txt", header=T, row.names=1, sep="\t", comment.char = "")

# OTU number:	
write.table(paste0("OTU number:\t", dim(otutab)[1]), file=paste0("otutab_sum.txt"), 
            append = F, quote = F, eol = "\n", row.names = F, col.names = F)

# Sample number:	
write.table(paste0("Sample number:\t", dim(otutab)[2]), file=paste0("otutab_sum.txt"), 
            append = T, quote = F, row.names = F, col.names = F)

# Total read counts:	
sample.size = as.data.frame(colSums(otutab))
total.read.counts = colSums(sample.size)
write.table(paste0("Total read counts:\t", total.read.counts), file=paste0("otutab_sum.txt"), 
            append = T, quote = F, row.names = F, col.names = F)

# Summary counts:
sample.size.summary = summary(sample.size)
write.table(paste0(sample.size.summary), file=paste0("otutab_sum.txt"), 
            append = T, quote = F, row.names = F, col.names = F)

# plot sample size by ggplot2
library(ggplot2)
df = sample.size
colnames(df)[1] = "ReadCounts"
idx = order(df$ReadCounts, decreasing = T)
df$SampleID = factor(rownames(sample.size), levels = row.names(df[idx,,drop=F]))


p = ggplot(df, aes(x=SampleID, y=ReadCounts, color="red")) +
  geom_point() + coord_flip() + 
  geom_text(label=df$ReadCounts, color = "blue", size = 3) +
  theme_bw()
# p  
ggsave(paste("otutab_sum.pdf", sep=""), p, width = 5, height = dim(otutab)[2]*0.2)
ggsave(paste("otutab_sum.png", sep=""), p, width = 5, height = dim(otutab)[2]*0.2)

