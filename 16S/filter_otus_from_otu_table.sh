#!/bin/bash
set -e

usage()
{
cat <<EOF >&2
Usage:
This script is used to filter OTU table by percentage in one sample, such 0.001
It requires at least 1 input files: otu_table and threshold
Date: 2017-3-28, version 1.0

OPTIONS:
	-t threshold of filter, default 0.001
	-o output director, default result/
EOF
}

# Default parameter
output='result'
tax_1_thre=0.001

# Analysis parameter
while getopts "t:o:" OPTION
do
	case $OPTION in
		t)
			tax_1_thre=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		h)
			usage
			exit 1
			;;
		?)
			usage
			exit 1
			;;
	esac
done




cat <<END >filter_otus.r

# Script for filter one of  OTU more than 0.1% abundance

# Input file "otu_table.txt"
# Output file "otu_id_k1.txt"
setwd("${output}")

thre=${tax_1_thre} # 1/1000

otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")
norm = t(t(otu_table)/colSums(otu_table,na=T))
#colSums(norm) # check
k1 = norm[apply(norm,1,max) > thre, ] # select OTU at least one sample > 0.1%
dim(k1)
write.table(rownames(k1), file="otu_id_k1.txt", quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\n")

END



Rscript filter_otus.r
