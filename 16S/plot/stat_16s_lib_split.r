# Install related packages
if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","reshape2"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("result")
library("ggplot2")
library("reshape2")

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

# Public file 1. "design.txt"  Design of experiment
design = read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/design.txt", header=T, row.names= 1, sep="\t") 

# Public file 2. "otu_table.txt"  raw reads count of each OTU in each sample
sample_count = read.table("T4_split.count", header=F, sep="\t") # dataframe
colnames(sample_count)=c("Sample","Count")

# Set group style, single or combine
if (FALSE){
	design$group=paste(design$genotype,sub_design$batch,sep = ".")
}else{
	design$group=design$genotype
}

sample_count_group = cbind(sample_count, design[match(sample_count$Sample, rownames(design)), ]) 

# stat: count identify, position: stack dodge fill
p = ggplot(sample_count_group, aes(x=Sample, y = Count, fill=group))+ 
  geom_bar(stat = "identity",position="dodge", width=0.7)+ 
  xlab("Library")+ylab("Pair Reads count")+main_theme+ theme(axis.text.x = element_text(angle = 90))+labs(title="Library T4")
# + theme(axis.text.x = element_text(size = 15, family = "myFont", color = "green", face = "bold", vjust = 0.5, hjust = 0.5, angle = 45))
p
ggsave("stat_lib_split_T4.pdf", p, width = 8, height = 5)
ggsave("stat_lib_split_T4.png", p, width = 8, height = 5)

