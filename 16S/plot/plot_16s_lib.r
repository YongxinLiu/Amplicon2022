# Install related packages
if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2","reshape2"))
}

## Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location, default is runing directory
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("result") # set work directory
library("ggplot2")
library("reshape2")

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(size=.5, colour="black"),
                    axis.line.y=element_line(size=.5, colour="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black", size=8),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    legend.text= element_text(size=8),
                    text=element_text(family="sans", size=8))

sum = read.table("qc.sum", header=T, sep="\t") # dataframe
data_all = as.data.frame(melt(sum, id.vars=c("Library")))
colnames(data_all)=c("Library","Process","value")

# stat: count identify, position: stack dodge fill
p = ggplot(data_all, aes(x=Library, y = value, fill = Process ))+ 
  geom_bar(stat = "identity",position="dodge", width=0.7)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("Library")+ylab("Pair Reads count")+main_theme
p
ggsave("stat_lib_qc_sum.pdf", p, width = 8, height = 2.5)
ggsave("stat_lib_qc_sum.png", p, width = 8, height = 2.5)


lib_len = read.table("length.txt", header=F, sep="\t") # dataframe
colnames(lib_len)=c("Library","Count","Length")

# stat: count identify, position: stack dodge fill
p = ggplot(lib_len, aes(x=Length, y = Count, color = Library ,shapes = Library))+  geom_line() +
  xlab("Length")+ylab("Pair Reads count")+main_theme
p
ggsave("stat_lib_length.pdf", p, width = 8, height = 2.5)
ggsave("stat_lib_length.png", p, width = 8, height = 2.5)




