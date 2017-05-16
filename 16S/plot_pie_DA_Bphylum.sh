#!/bin/bash
set -e

# Default parameter
compare_group=compare_group.txt
execute='TRUE'
ist='FALSE'
output='result'
# analysis level, default family, alternative otu
level=family
# pie styles, include basic, count and hollow, default hollow
style=hollow
tax_topn=tax_phylum.topN

usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    plot_pie_DA_Bphylum.sh, depended on taxonomy/DAOTU_egr.sh, incoperate plot_pie_family/otuBphylum.sh
Revision:    1.1
Date:        2017/5/11
Author:      Yong-Xin Liu
Email:       woodcorpse@163.com
Website:     http://bailab.genetics.ac.cn/
Description: This script is used to draw different abudance(DA) family in pie, and color by phylum; Ref: Figure 1C of Lebeis, S.L., Paredes, S.H., Lundberg, D.S., Breakfield, N., Gehring, J., McDonald, M., Malfatti, S., Glavina del Rio, T., Jones, C.D., Tringe, S.G., and Dangl, J.L. (2015). Salicylic acid modulates colonization of the root microbiome by specific bacterial taxa. Science 349, 860-864. 
Notes:       Only reads DA family, and group count by phylum; ouput pie hollow or basic
-------------------------------------------------------------------------------
Copyright:   2017 (c) Yong-Xin Liu
License:     GPL

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

If any changes are made to this script, please mail me a copy of the changes
-------------------------------------------------------------------------------
Version 1.0 2017/4/10
The first one , draw different abundance family and otu, and color by phylum. Pie style include basic, count, hollow. 
Version 1.1 2017/5/11
Add png format figure, for html report

# Input files
It requires at least 3 input files: 

## 1. otu_DM1vsWT_depleted/enriched.txt file, must have first column contain full family taxonomy, reads file list from compare_group.txt
#family
logFC   logCPM  LR      PValue  sig     level   percentage      fold    DM1.B1.1        DM1.B1.2        DM1.B1.3        DM1.B1.4        DM1.B1.5     
k__Bacteria;p__Bacteroidetes;c__Sphingobacteriia;o__Sphingobacteriales;f__      -1.2961421084962        11.1977987930675        16.0882020000147     
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;0.650  -5.09047570276644       6.22701576077493        14.277489561742 0.000157740037771433 
k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae      -4.03792934660057       7.64059991303981        18.4231432650522     
#otu
logFC   logCPM  LR      PValue  sig     level   otu     neglogp tax     percentage      fold    DM1.B1.1        DM1.B1.2        DM1.B1.3        DM1.B1.4
OTU_90  -1.98865595370181       9.3598532707379 20.4016554822535        6.27754668555716e-06    -1      depleted    OTU_90  11.9785313089758        p__Bacteroidete
OTU_168 -1.13924287034095       9.15055745516429        11.6522166211753        0.000641260789202598    -1      depleted    OTU_168 7.35207433635363        p__Bact
OTU_620 -3.09910064926046       7.22009539060197        13.9504820291733        0.000187689523524157    -1      depleted    OTU_620 8.58072143093136        Other

# 2. tax_phylum.topN, full name of phylum taxonomy
k__Bacteria;p__Proteobacteria
k__Bacteria;p__Actinobacteria
k__Bacteria;p__Bacteroidetes
k__Bacteria;p__Firmicutes
Low Abundance

#3. compare_group.txt, must have SampleA and SampleB in one line and seperated by tab, manually design
OE	WT
KO	WT
OE	KO

# Output file
1. pie_family/otu_AvB_enriched/depleted.pdf/png

OPTIONS:
	-c pairs needed to do comparasion, compare_group.txt
	-e execuate Rscript, TRUE or FALSE
	-h help
	-i install package TRUE or FALSE, default FALSE
	-l level, default family, alternative otu
	-o output director, default result/
	-s style, default hollow, alternative basic, count
	-t tax_phylum.topN
EOF
}

# Analysis parameter
while getopts "c:e:h:i:l:o:s:t:" OPTION
do
	case $OPTION in
		c)
			compare_group=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		h)
			usage
			exit 1
			;;
        i)
			ist=$OPTARG
			;;
		l)
			level=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		s)
			style=$OPTARG
			;;
		t)
			tax_topn=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done




cat <<END >plot_pie_DA_Bphylum.r
if ($ist){
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("ggplot2"))
}

# Basic plotting stuff
# Set working enviroment in Rstudio, select Session - Set working directory - To source file location
rm(list=ls()) # clean enviroment object
setwd(system("pwd", intern = T))
setwd("${output}")
library("ggplot2")
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   text=element_text(family="sans"))

# Read topN ordered phylum, and format to dataframe according to family/otu
tax0 = read.table("${tax_topn}", header=F, sep="\t") # dataframe
tax0\$V1=gsub("[\\\\w;_]+p__","",tax0\$V1,perl=TRUE) # format taxonomy
tax=tax0\$V1 # dataframe to vector
data=as.data.frame(tax) # reture to dataframe
data\$count=rep(0,length(tax)) # set all is 0
data\$tax  = factor(data\$tax, levels=data\$tax)   # set taxonomy order

# function plot pie in DA family/otu, color by phylum
plot_pie_familyBphylum <- function(sampleV){
	sampleA <- as.vector(sampleV\$sampA)
	sampleB <- as.vector(sampleV\$sampB)
	SampAvsB=paste(sampleA,"vs", sampleB, sep="")
    # DA family numbers in pie, and group by phylum: Read DA taxonomy, regexp taxonomy phylum, rename other to Low Abundance, add count 1, subset to new datafream and sum
END


for abundance in enriched  depleted
	do
cat <<END >>plot_pie_DA_Bphylum.r
    da_tax = read.table(paste("${level}_",SampAvsB,"_${abundance}.txt",sep=""), header=T, row.names= 1, sep="\t")
    if ("${level}"=="family"){
        da_tax\$tax=row.names(da_tax)
        da_tax\$tax=gsub("[\\\\w;_]+p__","",da_tax\$tax,perl=TRUE) # remove before phylum
    }
        da_tax\$tax=gsub("p__","",da_tax\$tax,perl=TRUE) # remove before phylum
    da_tax\$tax=gsub(";[\\\\w;_]+","",da_tax\$tax,perl=TRUE) # nremove after phylum
    da_tax\$tax = ifelse(da_tax\$tax %in% tax, da_tax\$tax, "Low Abundance") # non top to low abundance
    da_tax\$count=rep(1,length(da_tax\$tax)) # add each family count is 1
    sub_tax=da_tax[,c("tax","count")] # fetch tax and count to new dataframe

    sub_tax=rbind.data.frame(sub_tax,data)
    mat_mean <- aggregate(sub_tax[,-1], by=sub_tax[1], FUN=sum) # mean
    rownames(mat_mean)=mat_mean\$tax
    mat_mean=mat_mean[tax,]

    nums <- mat_mean\$x
    df <- data.frame(type = tax, nums = nums)  
    df\$type=factor(df\$type, levels=tax)
    p <- ggplot(data = df, mapping = aes(x = 'Content', y = nums, fill = type)) + geom_bar(stat = 'identity', position = 'stack', width = 1)  
    if ("${style}" == "basic") {
        # Draw legend on right
        label <- paste(df\$type, df\$nums, sep = ' ')  
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + scale_fill_discrete(labels = label)  +main_theme
    }else if ("${style}" == "hollow"){
        # Draw total number in inner, size = 3 is 9 pt
        label=rep(sum(df\$nums),length(df\$nums))
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + theme(legend.position = "none") + geom_text(aes(x = 0, label = label), size = 3) + theme(panel.background=element_blank(),panel.grid=element_blank())
    }else if ("${style}" == "count"){
        # Draw total number in inner
        label=df\$nums
        p = p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + theme(legend.position = "none") + geom_text(aes(y = df\$nums/2 + c(0, cumsum(df\$nums)[-length(df\$nums)]), x = sum(df\$nums)/150, label = label), size = 3) + theme(panel.background=element_blank(),panel.grid=element_blank())
    }
    ggsave(file=paste("pie_${level}_", SampAvsB, "_${abundance}.pdf", sep=""), p, width=1.5, height=1.5, useDingbats=F) # 1/4 of half page
    ggsave(file=paste("pie_${level}_", SampAvsB, "_${abundance}.png", sep=""), p, width=1.5, height=1.5)
    print(paste("pie_${level}_", SampAvsB, "_${abundance}.pdf is finished!!!", sep=""))
END
done

cat <<END >>plot_pie_DA_Bphylum.r
}

if ("${compare_pair}" == "FALSE") {
	compare_data <- as.vector(unique(sub_design\$group))
	len_compare_data <- length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare <- as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			plot_pie_familyBphylum(tmp_compare)
		}
	}
}else {
	compare_data <- read.table("${compare_group}", sep="\t", check.names=F, quote='', com='')
	colnames(compare_data) <- c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){
    plot_pie_familyBphylum(compare_data[i,])}
}	
END

if test "${execute}" == "TRUE";
then
	Rscript plot_pie_DA_Bphylum.r
fi
