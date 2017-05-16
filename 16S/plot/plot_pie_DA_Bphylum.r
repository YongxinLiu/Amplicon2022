if (FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2"))
}

# Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("result_k1-c")
library("ggplot2")
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   text=element_text(family="sans"))

# Read topN ordered phylum, and format to dataframe according to family/otu
tax0 = read.table("tax_phylum.topN", header=F, sep="\t") # dataframe
tax0$V1=gsub("[\\w;_]+p__","",tax0$V1,perl=TRUE) # format taxonomy
tax=tax0$V1 # dataframe to vector
data=as.data.frame(tax) # reture to dataframe
data$count=rep(0,length(tax)) # set all is 0
data$tax  = factor(data$tax, levels=data$tax)   # set taxonomy order

# function plot pie in DA family/otu, color by phylum
plot_pie_familyBphylum <- function(sampleV){
	sampleA <- as.vector(sampleV$sampA)
	sampleB <- as.vector(sampleV$sampB)
	SampAvsB=paste(sampleA,"vs", sampleB, sep="")
    # DA family numbers in pie, and group by phylum: Read DA taxonomy, regexp taxonomy phylum, rename other to Low Abundance, add count 1, subset to new datafream and sum
    da_tax = read.table(paste("otu_",SampAvsB,"_enriched.txt",sep=""), header=T, row.names= 1, sep="\t")
    if ("otu"=="family"){
        da_tax$tax=row.names(da_tax)
        da_tax$tax=gsub("[\\w;_]+p__","",da_tax$tax,perl=TRUE) # remove before phylum
    }
        da_tax$tax=gsub("p__","",da_tax$tax,perl=TRUE) # remove before phylum
    da_tax$tax=gsub(";[\\w;_]+","",da_tax$tax,perl=TRUE) # nremove after phylum
    da_tax$tax = ifelse(da_tax$tax %in% tax, da_tax$tax, "Low Abundance") # non top to low abundance
    da_tax$count=rep(1,length(da_tax$tax)) # add each family count is 1
    sub_tax=da_tax[,c("tax","count")] # fetch tax and count to new dataframe

    sub_tax=rbind.data.frame(sub_tax,data)
    mat_mean <- aggregate(sub_tax[,-1], by=sub_tax[1], FUN=sum) # mean
    rownames(mat_mean)=mat_mean$tax
    mat_mean=mat_mean[tax,]

    nums <- mat_mean$x
    df <- data.frame(type = tax, nums = nums)  
    df$type=factor(df$type, levels=tax)
    p <- ggplot(data = df, mapping = aes(x = 'Content', y = nums, fill = type)) + geom_bar(stat = 'identity', position = 'stack', width = 1)  
    if ("hollow" == "basic") {
        # Draw legend on right
        label <- paste(df$type, df$nums, sep = ' ')  
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + scale_fill_discrete(labels = label)  +main_theme
    }else if ("hollow" == "hollow"){
        # Draw total number in inner, size = 3 is 9 pt
        label=rep(sum(df$nums),length(df$nums))
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + theme(legend.position = "none") + geom_text(aes(x = 0, label = label), size = 3) + theme(panel.background=element_blank(),panel.grid=element_blank())
    }else if ("hollow" == "count"){
        # Draw total number in inner
        label=df$nums
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + theme(legend.position = "none") + geom_text(aes(y = df$nums/2 + c(0, cumsum(df$nums)[-length(df$nums)]), x = sum(df$nums)/150, label = label), size = 3) + theme(panel.background=element_blank(),panel.grid=element_blank())
    }
    ggsave(file=paste("pie_otu_", SampAvsB, "_enriched.pdf", sep=""), p, width=1.5, height=1.5, useDingbats=F) # 1/4 of half page
    ggsave(file=paste("pie_otu_", SampAvsB, "_enriched.png", sep=""), p, width=1.5, height=1.5)
    print(paste("pie_otu_", SampAvsB, "_enriched.pdf is finished!!!", sep=""))
    da_tax = read.table(paste("otu_",SampAvsB,"_depleted.txt",sep=""), header=T, row.names= 1, sep="\t")
    if ("otu"=="family"){
        da_tax$tax=row.names(da_tax)
        da_tax$tax=gsub("[\\w;_]+p__","",da_tax$tax,perl=TRUE) # remove before phylum
    }
        da_tax$tax=gsub("p__","",da_tax$tax,perl=TRUE) # remove before phylum
    da_tax$tax=gsub(";[\\w;_]+","",da_tax$tax,perl=TRUE) # nremove after phylum
    da_tax$tax = ifelse(da_tax$tax %in% tax, da_tax$tax, "Low Abundance") # non top to low abundance
    da_tax$count=rep(1,length(da_tax$tax)) # add each family count is 1
    sub_tax=da_tax[,c("tax","count")] # fetch tax and count to new dataframe

    sub_tax=rbind.data.frame(sub_tax,data)
    mat_mean <- aggregate(sub_tax[,-1], by=sub_tax[1], FUN=sum) # mean
    rownames(mat_mean)=mat_mean$tax
    mat_mean=mat_mean[tax,]

    nums <- mat_mean$x
    df <- data.frame(type = tax, nums = nums)  
    df$type=factor(df$type, levels=tax)
    p <- ggplot(data = df, mapping = aes(x = 'Content', y = nums, fill = type)) + geom_bar(stat = 'identity', position = 'stack', width = 1)  
    if ("hollow" == "basic") {
        # Draw legend on right
        label <- paste(df$type, df$nums, sep = ' ')  
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + scale_fill_discrete(labels = label)  +main_theme
    }else if ("hollow" == "hollow"){
        # Draw total number in inner, size = 3 is 9 pt
        label=rep(sum(df$nums),length(df$nums))
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + theme(legend.position = "none") + geom_text(aes(x = 0, label = label), size = 3) + theme(panel.background=element_blank(),panel.grid=element_blank())
    }else if ("hollow" == "count"){
        # Draw total number in inner
        label=df$nums
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + theme(legend.position = "none") + geom_text(aes(y = df$nums/2 + c(0, cumsum(df$nums)[-length(df$nums)]), x = sum(df$nums)/150, label = label), size = 3) + theme(panel.background=element_blank(),panel.grid=element_blank())
    }
    ggsave(file=paste("pie_otu_", SampAvsB, "_depleted.pdf", sep=""), p, width=1.5, height=1.5, useDingbats=F) # 1/4 of half page
    ggsave(file=paste("pie_otu_", SampAvsB, "_depleted.png", sep=""), p, width=1.5, height=1.5)
    print(paste("pie_otu_", SampAvsB, "_depleted.pdf is finished!!!", sep=""))
}

if ("" == "FALSE") {
	compare_data <- as.vector(unique(sub_design$group))
	len_compare_data <- length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare <- as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			plot_pie_familyBphylum(tmp_compare)
		}
	}
}else {
	compare_data <- read.table("/mnt/bai/yongxin/ath/jt.terpene.16S/batch3/doc/group_compare.txt", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) <- c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){
    plot_pie_familyBphylum(compare_data[i,])}
}	
