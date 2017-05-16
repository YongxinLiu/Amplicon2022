# Install related packages
# Work well in R3.3.3
if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggtree","colorspace"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("result")
library("ggtree")
library("colorspace")

tree <- read.tree("tax_rep_seqs.tree")
tax <- read.table("tax_rep_seqs.tax",row.names=1)
colnames(tax) = c("kingdom","phylum","class","order")


groupInfo <- split(row.names(tax), tax$phylum) # OTU and phylum for group
tree <- groupOTU(tree, groupInfo)
pdf(file="ggtree_phylum.pdf", width=8, height=8)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+
  scale_color_manual(values=c(rainbow_hcl(length(unique(tax$phylum))+1)), breaks=1:length(unique(tax$phylum)), labels=levels(tax$phylum))+
  theme(legend.position = "right") +geom_tiplab2(size=3)
dev.off()
png(file="ggtree_phylum.png", width=8, height=8, units = "in", res = 300)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+
  scale_color_manual(values=c(rainbow_hcl(length(unique(tax$phylum))+1)), breaks=1:length(unique(tax$phylum)), labels=levels(tax$phylum))+
  theme(legend.position = "right") +geom_tiplab2(size=3)
dev.off()



groupInfo <- split(row.names(tax), tax$class) # OTU and class for group
tree <- groupOTU(tree, groupInfo)
pdf(file="ggtree_class.pdf", width=8, height=8)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+
  scale_color_manual(values=c(rainbow_hcl(length(unique(tax$class))+1)), breaks=1:length(unique(tax$class)), labels=levels(tax$class))+
  theme(legend.position = "right") +geom_tiplab2(size=3)
dev.off()
png(file="ggtree_class.png", width=8, height=8, units = "in", res = 300)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+
  scale_color_manual(values=c(rainbow_hcl(length(unique(tax$class))+1)), breaks=1:length(unique(tax$class)), labels=levels(tax$class))+
  theme(legend.position = "right") +geom_tiplab2(size=3)
dev.off()



groupInfo <- split(row.names(tax), tax$order) # OTU and order for group
tree <- groupOTU(tree, groupInfo)
pdf(file="ggtree_order.pdf", width=8, height=8)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+
  scale_color_manual(values=c(rainbow_hcl(length(unique(tax$order))+1)), breaks=1:length(unique(tax$order)), labels=levels(tax$order))+
  theme(legend.position = "right") +geom_tiplab2(size=3)
dev.off()
png(file="ggtree_order.png", width=8, height=8, units = "in", res = 300)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group))+
  scale_color_manual(values=c(rainbow_hcl(length(unique(tax$order))+1)), breaks=1:length(unique(tax$order)), labels=levels(tax$order))+
  theme(legend.position = "right") +geom_tiplab2(size=3)
dev.off()


