#!/bin/bash
set -e

### Default parameter
compare_group=compare_group.txt
design=design.txt
execute='TRUE'
ist='FALSE' # install package, default FALSE
output='result_k1-c' # default work directory, remove low abundance < .1% and p__Cyanobacteria,p__Chloroflexi
width=4
height=2.5
text_size=8
pvalue=0.05

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    rmd_16s.sh
Revision:    1.0
Date:        2017/5/14
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: This script is used to perform calculate different abundance OTUs, and statistic by edgeR.
Notes:       Visualize in volcano plot, manhattan plot and heatmap
-------------------------------------------------------------------------------
Copyright:   2017 (c) Yong-Xin Liu
License:     GPL
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
If any changes are made to this script, please mail me a copy of the changes
-------------------------------------------------------------------------------
Version 1.0 2017/5/14
The first one , write R markdown

# All input and output should be in -o directory, or give relative -o path, or absolute path
# Input files: design.txt, otu_table.txt, rep_seqs_tax.txt, compare_group.txt and otu_table_css.txt

# 1. design.txt, grouping samples, must have SampleID and group info, group1/2 can give to parameter g1 and g2, manually design
SampleID	BarcodeSequence	group	gene	batch	description
WT.1	TAGCTT	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down
WT.2	GGCTAC	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down
WT.3	CGCGCG	WT	ggps9.10	2	double mutant of ggps9-ggps10, cause A/B down

# 2. compare_group.txt, must have SampleA and SampleB in one line and seperated by tab, manually design
OE	WT
KO	WT
OE	KO


# Output file
1. Rmd: report file, can format to html/pdf

OPTIONS:
	-a css OTU table, default in otu_table_css.txt
	-b pvalue threshold, default=0.05
	-c pairs needed to do comparasion, compare_group.txt
	-d design for each samples, default design.txt
	-e execuate Rscript, TRUE or FALSE
	-f OTU table, default in otu_table.txt
	-g group order by group1 list, default FALSE, not allow both TRUE with -m merge group, when B exist is TRUE
	-h figure height, default 2.5
	-i install package TRUE or FALSE, default FALSE
	-m merge group 1 and group 2 as new group name, default FALSE, not allow both TRUE with -g group order
	-n show top N taxonomy number, default 5, recommend phylum is 5-10
	-o output director, default result_k1-c/
	-p default TRUE for compare group, set FALSE to loop each pair group
	-s text size, default 8
	-t default rep_seqs_tax.txt, taxonomy of all OTU
	-w figure width, default 4
	-A group1 name
	-B group1 selected list
	-C group2 name
	-D group2 selected list
	-S sample abundance style, default none, alternative css and percentage
	-h/? show help of script

Example:
	rmd_16s.sh -e FALSE
EOF
}


# Analysis parameter
while getopts "a:b:c:d:e:f:g:h:i:m:o:p:s:t:w:A:B:C:D:" OPTION
do
	case $OPTION in
		a)
			otu_table_css=$OPTARG
			;;
		b)
			pvalue=$OPTARG
			;;
		c)
			compare_group=$OPTARG
			;;
		d)
			design=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		f)
			norm_otu=$OPTARG
			;;
		g)
			group_order=$OPTARG
			;;
		h)
			height=$OPTARG
			;;
		i)
			ist=$OPTARG
			;;
		m)
			merge_group=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		p)
			pair_group=$OPTARG
			;;
		s)
			text_size=$OPTARG
			;;
		t)
			taxonomy=$OPTARG
			;;
		w)
			width=$OPTARG
			;;
		A)
			g1=$OPTARG
			;;
		B)
			g1_list=$OPTARG
			group_order=TRUE
			select1=TRUE
			;;
		C)
			g2=$OPTARG
			;;
		D)
			g2_list=$OPTARG
			select2=TRUE
			;;
		S)
			style=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done



cat <<END >index.Rmd
--- 
title: "细菌16S扩增子分析"
author:
- 客户单位：中国科学院遗传与发育生物学研究所白洋组姜婷
- 服务单位：中国科学院遗传与发育生物学研究所白洋组刘永鑫
- 联系方式：yxliu@genetics.ac.cn
- 项目编号：20170428-1
- 项目周期：2017-04-28 ~ 2017-05-31
- 官方网站：http://bailab.genetics.ac.cn/
date: '\`r Sys.Date()\`'
documentclass: article
link-citations: yes
biblio-style: apalike
---

\`\`\`{r setup, include=FALSE}
library(knitr)
output <- opts_knit\$get("rmarkdown.pandoc.to")
html = FALSE
latex = FALSE
opts_chunk\$set(echo = FALSE, out.width="100%", fig.align="center", fig.show="hold", warning=FALSE, message=FALSE)
if (output=="html") {
	html = TRUE
}
if (output=="latex") {
	opts_chunk\$set(out.width="95%", out.height='0.7\textheight', out.extra='keepaspectratio', fig.pos='H')
	latex = TRUE
}
knitr::opts_chunk\$set(cache=TRUE, autodep=TRUE)

mtime <- function(files){
  lapply(Sys.glob(files), function(x) file.info(x)\$mtime)
}

set.seed(718)
\`\`\`

\`\`\`{asis, echo=html}
# Bailab, SKLPG/CEPAMS, IGDB, CAS {-}
\`\`\`

\`\`\`{r cover, eval=html, out.width="99%"}
figs_1 = paste0("figure/slide", c("1", "2"),"_raw.jpg")
knitr::include_graphics(figs_1)
\`\`\`

END



cat <<END >01-aim.Rmd
# 课题目的 {#project_aim}

比较分析拟南芥野生型、萜类合成(TPS)基因过表达和突变体各样品组间微生物组的物种丰富度、群落结构、各分类学级别及OTU水平上相对丰度和调控网络的差异。揭示萜类化合物对细菌微生物组的影响，以及在调控微生物群落结果中的规律和意义。

END



cat <<END >02-design.Rmd
# 课题设计 {#project_design}

样品准备如 Table \@ref(tab:design) 所示。[design.txt](doc/design.txt)

\`\`\`{r design}
table_design <- read.table("doc/design.txt", sep="\t", header=T)
knitr::kable(table_design, caption="样品详细及分组信息总结。", booktabs=TRUE)
\`\`\`

END



cat <<END >03-scheme.Rmd
# 课题方案 {#project_scheme}

我的材料有那些菌？taxonomy tree, phylogenetic tree.  

实验组和对照组间是否存在不同？alpha diversity, beta diversity.  

具体有那些不同？Differentially abundance taxonomy and OTU.  

整个分析流程包含以下10部分内容：报告测序数据质控；测序数据过滤及各步骤统计、样品数据量和长度分布；Alpha多样性分析: Shannon entropy和observed OTU；Beta多样性分析: 采用bray curtis和weighted unifrac距离计算距离的主坐标轴分析(PCoA/MDS)；限制条件的PCoA分析(CPCoA/CCA/RDA); 分类树及进化树展示OTU物种信息及进化关系；各分类级别丰度分析：包括门、纲、目、科、属水平；差异OTU分析：包括火山图、热图、曼哈顿图展示差异OTU数量、丰度、变化样式及分类学信息；组间差异OTU比较，观察不同组间的分类学样式，以及共有或特有OTU；其它有待进一步分析的内容，如OTU调控网络构建等。


**1. 测序reads数量和质量评估；Quality control of sequencing reads**

Table: (\#tab:seq-quality-explanatioan-ch) 测序质量评估结果解读方法

-----------------------------------------------------------------------------------
评估内容                   结果解释 (图例中会标记对应评估内容为PASS、WARN和FAIL, 具体处理方式详见下面中英文解释)
-------------------------  --------------------------------------------------------------------------------
Per base quality           测序reads从5'到3'的碱基的质量值 (Q)。该值越大越代表对应碱基测序准确度越高。假设p为一个碱基测序错误的概率，则Q=-10 * log10(p). 质量值为10时，对应碱基出错的概率为10%；质量值为20时，对应碱基出错的概率为1%。通常来讲，3'端的碱基质量会低于5'端；另外5'端最初几个碱基也会出现较大的质量值波动。我们在后期处理时，会去除低质量的碱基以保证分析结果的准确性。

Adaptor content            判断测序reads中是否残留接头序列。存在接头序列和不存在接头序列都是合理的，取决于测序数据下机后是否进行了接头去除和去除的是否完整。若在分析时检测到接头序列存在，我们会首先去除接头，然后进行后续分析，以保证分析结果的准确性。

Per sequence GC content    测序reads的GC含量。正常的测序reads的GC含量符合正态分布模式 (形如图中蓝色的倒钟形线)。若在平滑的曲线上存在一个尖峰表示测序样品存在特定的序列污染如混入了引物二聚体。若GC含量分布曲线比较平坦则代表可能存在不同物种的序列污染。当这一指标异常时，可能导致后期的序列比对或拼接存在问题，需要引起注意。

Per base sequence content  测序reads的碱基偏好性。正常的测序结果中一个序列不同的碱基没有偏好性，图中的线应平行。Bisulfite测序中存在甲基化的C到T的转变，会导致这一评估结果异常。我们人工核验无误后，可以忽略软件对这一检测结果的评价。
-----------------------------------------------------------------------------------


Table: (\#tab:seq-quality-explanatioan-en) Explanation of quality control by fastqc.

-----------------------------------------------------------------------------------
Analysis                   Explanation
-------------------------  --------------------------------------------------------------------------------
Per base quality           The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

Per sequence GC content    Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example),  which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species.

Adaptor content            Any library where a reasonable proportion of the insert sizes are shorter than the read length will trigger this module. This doesn't indicate a problem as such - just that the sequences will need to be adapter trimmed before proceeding with any downstream analysis.

Per base sequence content  In a random library you would expect that there would be little to no difference between the different bases of a sequence run,  so the lines in this plot should run parallel with each other. The relative amount of each base should reflect the overall amount of these bases in your genome,  but in any case they should not be hugely imbalanced from each other.
-----------------------------------------------------------------------------------


(ref:scheme-read-fastqc) 测序Reads质量评估。HiSeq2500产出Clean reads左端(A)和右端(B)各250 bp数据质量评估，选取测序reads碱基质量分布判断建库或测序有无异常。双端数据raw和clean reads左端(C)和右端(D)接头及引物污染情况分布，接头去除干净与否、及有效数据比例评估。  Quality control of raw reads (Andrews, 2010)

\`\`\`{r scheme-read-fastqc, fig.cap="(ref:scheme-read-fastqc)"}
knitr::include_graphics("figure/fig1.fastqc.png")
\`\`\`

**2. 样品提取及过滤各步骤统计；Statistics of reds filter processes**

(ref:scheme-read-summary) 统计文库处理过程及样品可用数据量。(A) 柱状图展示各文库数据标准化筛选各步骤有效数据分布。主要包括数据低质量及污染序列过滤、双端合并、筛选扩增子并统一序列方向、按barcode拆分样品、去除5’引物序列、去除3’引物序列为下一步分析的高质量样本序列；(B). 柱状图展示各样品的数据量分布，最小值也大于2万，大部分在12万左右，完全符合实验设计要求；(C) 可用数据的长度分布，可以看到本实验扩增子长度范围集中在360-390 bp，主峰位于370-380 bp间。  Statistics of reads filter processes in libraries and data size of samples. (A) Bar plot showing reads count of each library in read filter process; (B) Bar plot showing reads counts of each sample; (C) Length distribution of amplicons (Caporaso et al., 2010, Edgar, 2013).

\`\`\`{r scheme-read-summary, fig.cap="(ref:scheme-read-summary)"}
knitr::include_graphics("figure/fig2.summary.png")
\`\`\`

**3. Alpha多样性分析；Alpha (α) diversity**

(ref:scheme-sample-alpha) Alpha多样性展示各组间微生物多样性，方法采用(A) Shannon index，包括样品的可操作分类单元(operational taxonomic unit, OTU)数量及种类丰度信息；(B) Observed OTUs index，只包括样品OTU种类信息。图中KO(knock out)代表基因敲除突变体，OE(overexpression)代表基因过表达株系，WT(wild-type)代表野生型。附表有各种间t-test方法统计的p-value水平。此外还可计算chao1和PD whole tree等方法下的多样性分析。[各Alpha多样性计算方法详细](http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html)  Within sample diversity (α-diversity) measurements among each genotype. (A) Shannon index, estimated species richness and evenness; (B) Observed OTUs index, only calculate species richness. These results indicate genotype not significantly change microbial diversity. The horizontal bars within boxes represent median. The tops and bottoms of boxes represent 75th and 25th quartiles, respectively. The upper and lower whiskers extend 1.5× the interquartile range from the upper edge and lower edge of the box, respectively. All outliers are plotted as individual points (Edwards et al., 2015).

\`\`\`{r scheme-sample-alpha, fig.cap="(ref:scheme-sample-alpha)"}
knitr::include_graphics("figure/fig3.alpha.png")
\`\`\`

**4. Beta多样性分析；Beta (β) diversity **

(ref:scheme-sample-beta) 采用主坐标轴分析展示第1/2坐标轴下各组间微生物组差异(dissimilarity)，距离计算方法采用(A) bray curtis; (B) weighted unifrac. 如图A中可以看到坐标轴1可以解释24.15%的变异，坐标轴2可以解释12.32%的变异，KO与WT较为相似；而OE在第一轴方向上明显与WT分开，表明其微生物组呈现明显变化；同时还发现KO1中存在三个样品存在明显异常。  Principal coordinate analysis (PCoA) using the (A) bray curtis metric and (B) weighted unifrac metric shows dissimilarity of microbial communities. The result indicates that the largest separation is between WT and OE (PCoA 1) and the second largest source of variation is between WT and KO (PCoA 2). (Edwards et al., 2015)

\`\`\`{r scheme-sample-beta, fig.cap="(ref:scheme-sample-beta)"}
knitr::include_graphics("figure/fig4.beta.png")
\`\`\`

**5. 限制条件下的主坐标轴分析；Constrained principal coordinate analysis**

(ref:scheme-sample-CPCoA) 以基因型为条件分析贡献率及组间差异；分析表明基因型可解释微生物组的22.7%的变异，且各基因型间均可明显分开，且KO和OE的重复又能很好聚在一起，表明不同基因对微生物组的群落结构有明显的调控作用，且不同突变体和过表达株系的位点和生物学重复间表现出良好的可重复性。  Constrained principal coordinate analysis on bacterial microbiota. Variation between samples in Bray-Curtis distances constrained by genotype (22.7% of the overall variance; p < 0.05) (Bulgarelli et al., 2015).

\`\`\`{r scheme-sample-CPCoA, fig.cap="(ref:scheme-sample-CPCoA)"}
knitr::include_graphics("figure/fig5.CPCoA.png")
\`\`\`

**6. 分类树及进化树展示OTU物种信息及进化关系；Taxonomy and phylogenetic tree of OTU**

(ref:scheme-sample-tree) 样品中高丰度(>0.5%)OTU的分类树和系统发生学分析。(A)分类树，其中OTU按分类学的科水平进行背景高亮着色，显示本研究中主要丰度的细菌科；(B)系统发生树，按门水平进行着色，结果表明细菌的物种注释信息与16S的序列发生树的进化关系高度一致。  Taxonomy and phylogenetic tress show high abundance OTU (>0.5%), and their family and phylum annotation of taxonomy (Asnicar et al., 2015, Yu et al., 2017). 

\`\`\`{r scheme-sample-tree, fig.cap="(ref:scheme-sample-tree)"}
knitr::include_graphics("figure/fig6.tree.png")
\`\`\`

**7. 分类学不同分类级别的丰度分析；Differentially abundance of bacterial in each taxonomy level**

(ref:scheme-sample-tax) 柱状图展示各类微生物组分类学门水平相对丰度。(A) 堆叠柱状图，X轴为各样品组，Y轴为各门类相对百分比，只列出了丰度大于0.1%的门，其它所有门归入Low Abundance类。(B). 条形图展示最高丰度的五大菌门平均丰度及标准误，我们可以观察到与WT相比，各基因型的Proteobacteria丰度降低，而Actinobacteria丰度升高。注: 分类学注释可从门、纲、目、科、属五个级别进行丰度可视化及差异统计分析。  Bar plot showing phyla abundances in each genotype. (A). Stack plot showing high abundance (>0.1%) phylum; (B). Bar plot showing top 5 phylum abundance and error bar in each genotype. All the KO and OE were show enriched in Actinobacteria and depleted in Proteobacteria. Note: Differentially abundance taxonomy can analyze under phylum, order, class, family and genus level (Bulgarelli et al., 2015, Lebeis et al., 2015).

\`\`\`{r scheme-sample-tax, fig.cap="(ref:scheme-sample-tax)"}
knitr::include_graphics("figure/fig7.tax.png")
\`\`\`

**8. 差异OTUs分析；Differentially abundance OTUs**

(ref:scheme-sample-otu) KO1基因型存在一些丰度显著上调或下调的OTU (P & FDR < 0.05, GLM likelihood rate test)。(A) 火山图展示KO与WT相比OTU的变化，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝代表显著上下调，图中数字代表显著差异OTU数量，形状代表OTU的门水平物种注释；（B）热图展示KO与WT显著差异OTU在每个样品中丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低；可以看到我们找到的差异OTU在每组样品中重复非常好，同时也发现了在beta diversity分析中发现的KO1中存在的两个异常样品应该为KO1.7, KO1.8, 需检查实验材料准备了取材步骤有无问题？或补弃样品重复（C）曼哈顿图展示OTU的变化情况及在各门水平中的分布，x轴为OTU按物种门水平物种注释字母排序，y轴为pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，颜色为门水平注释，大小为相对丰度，形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调。  KO1 are enriched and depleted for certain OTUs (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of abundance and fold change of OTUs; (B) Heatmap showing differentially abundance OTUs of KO1 compared WT; (C) Manhattan plot showing phylum pattern of differentially abundance OTUs. These results show Actinobacterial has more enriched OTUs (Bai et al., 2015, Edwards et al., 2015, Zgadzaj et al., 2016).

\`\`\`{r scheme-sample-otu, fig.cap="(ref:scheme-sample-otu)"}
knitr::include_graphics("figure/fig8.otu.png")
\`\`\`

**9. 组间差异OTU比较；Compare differentially abundance OTUs among groups**

(ref:scheme-sample-overlap) 比较组间差异OTU的分类学样式、共有或特有。(A) 饼形图展示各种差异OTU细菌门水平分类比例。中间数字为所有显著差异OTU的数目。可以看到KO1与KO2样式相似，OE1与OE2样式相似。且上调OTU较多为Actinobacteria，而下调OTU绝大多数为Proteobacteria。(B) 维恩图展示各基因型差异OTUs间的共有和特有数量。图中所显各基因型组间重复间大部分OTUs共有；而且还发现KO和OE还存在一些相似变化样式的OTUs。  Taxonomy, common and unique OTUs in each group. (A) Pie charts show phyla of bacterial OTUs identified as either enriched or depleted in each genotype compared with WT. The number of OTUs in each category is noted inside each donut. (B) Venn diagrams show common and unique OTUs in each group (Lebeis et al., 2015).

\`\`\`{r scheme-sample-overlap, fig.cap="(ref:scheme-sample-overlap)"}
knitr::include_graphics("figure/fig9.overlap.png")
\`\`\`

**10. 其它数据分析过程中发现的有意思的点，商讨后，有意义的深入分析；Other points and ideas for further discussion and analysis **

**参考文献Reference**

Andrews, S. (2010) FastQC: a quality control tool for high throughput sequence data.  
Asnicar, F., Weingart, G., Tickle, T.L., Huttenhower, C. and Segata, N. (2015) Compact graphical representation of phylogenetic data and metadata with GraPhlAn. PeerJ, 3, e1029.  
Bai, Y., Müller, D.B., Srinivas, G., Garrido-Oter, R., Potthoff, E., Rott, M., Dombrowski, N., Münch, P.C., Spaepen, S., Remus-Emsermann, M., Hüttel, B., McHardy, A.C., Vorholt, J.A. and Schulze-Lefert, P. (2015) Functional overlap of the Arabidopsis leaf and root microbiota. Nature, 528, 364-369.  
Bulgarelli, D., Garrido-Oter, R., Munch, P.C., Weiman, A., Droge, J., Pan, Y., McHardy, A.C. and Schulze-Lefert, P. (2015) Structure and function of the bacterial root microbiota in wild and domesticated barley. Cell Host Microbe, 17, 392-403.  
Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Pena, A.G., Goodrich, J.K., Gordon, J.I., Huttley, G.A., Kelley, S.T., Knights, D., Koenig, J.E., Ley, R.E., Lozupone, C.A., McDonald, D., Muegge, B.D., Pirrung, M., Reeder, J., Sevinsky, J.R., Turnbaugh, P.J., Walters, W.A., Widmann, J., Yatsunenko, T., Zaneveld, J. and Knight, R. (2010) QIIME allows analysis of high-throughput community sequencing data. Nat. Methods, 7, 335-336.  
Edgar, R.C. (2013) UPARSE: highly accurate OTU sequences from microbial amplicon reads. Nat Meth, 10, 996-998.  
Edwards, J., Johnson, C., Santos-Medellín, C., Lurie, E., Podishetty, N.K., Bhatnagar, S., Eisen, J.A. and Sundaresan, V. (2015) Structure, variation, and assembly of the root-associated microbiomes of rice. Proceedings of the National Academy of Sciences, 112, E911-E920.  
Lebeis, S.L., Paredes, S.H., Lundberg, D.S., Breakfield, N., Gehring, J., McDonald, M., Malfatti, S., Glavina del Rio, T., Jones, C.D., Tringe, S.G. and Dangl, J.L. (2015) Salicylic acid modulates colonization of the root microbiome by specific bacterial taxa. Science, 349, 860-864.  
Yu, G., Smith, D.K., Zhu, H., Guan, Y. and Lam, T.T.Y. (2017) Ggtree: an r package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution, 8, 28-36.  
Zgadzaj, R., Garrido-Oter, R., Jensen, D.B., Koprivova, A., Schulze-Lefert, P. and Radutoiu, S. (2016) Root nodule symbiosis in Lotus japonicus drives the establishment of distinctive rhizosphere, root, and nodule bacterial communities. Proc. Natl. Acad. Sci. USA, 113, E7996-e8005.

END



cat <<END >04-a-sequenceQuality.Rmd
# 测序质量总结 {#sequencing_quality_summary}

## 测序质量评估 {#sub-sequence-qc}

说明：16S扩增子测序数据主要来自HiSeq2500产出的双端各250 bp (PE250)数据，因为读长长且价格便宜(性价比高)。HiSeqX PE150和MiSeq PE300也比较常见，但PE150过短分辨率低，而PE300价格高且末端序列质量过低。此外454在之前研究较多且设备已经停产，PacBio读长长可直接测序16S全长1.5kb代表未来的趋势。测序公司通常会返回raw data和clean data两种数据，raw data为测序获得的原始数据，而clean data则为去除含有接头序列及测序不确定N比例较高的结果，通常直接采用clean data进行质量评估及后续分析。数据质量评估结果中测序reads碱基质量分布图，常用于判断建库或测序有无异常。序列重复情况分布，判断原始序列的DNA质量、重复序列比例及PCR扩增重复情况，如重复序列较高可能某些菌高丰度或PCR扩增导致，对低丰度菌的结果影响较大。

END

list=`cat doc/library.txt|tr "\n" " "`

for library in ${list}
	do
cat <<END >>04-a-sequenceQuality.Rmd
(ref:quality-fastqc-${library}) 测序Reads质量评估文库${library}。Clean reads左端(A)和右端(B)数据质量评估；clean reads左端(C)和右端(D)序列重复情况分布。Quality control of clean reads [HTML report of library ${library}_1](clean_data/${library}_1_fastqc.html)  [HTML report of library ${library}_2](clean_data/${library}_2_fastqc.html)

\`\`\`{r quality-fastqc-${library}, fig.cap="(ref:quality-fastqc-${library})", out.width="49%"}
figs_1 = paste0("clean_data/${library}_", c("1_fastqc/Images/per_base_quality", "2_fastqc/Images/per_base_quality", "1_fastqc/Images/duplication_levels", "2_fastqc/Images/duplication_levels"),".png")
knitr::include_graphics(figs_1)
\`\`\`

END
done

cat <<END >>04-a-sequenceQuality.Rmd
## 样品数据量分布 {#sub-sequence-count}


END

for library in ${list}
	do
cat <<END >>04-a-sequenceQuality.Rmd
(ref:quality-split-${library}) 文库${library}各样品按barcode拆分获得的高质量测序数据，按实验设计组着色。Distribution of sequenced reads of samples in library ${library}. Samples were colored by group information. 1 Million = 10^6^. [PDF](result/stat_lib_split_${library}.pdf)

\`\`\`{r quality-split-${library}, fig.cap="(ref:quality-split-${library})"}
knitr::include_graphics("result/stat_lib_split_${library}.png")
\`\`\`

END
done

cat <<END >>04-a-sequenceQuality.Rmd
## 测序数据预处理总结 {#sub-sequence-summary}

(ref:quality-sum) 测序文库数据量和长度分布。(A) 柱状图展示各文库数据标准化筛选各步骤有效数据分布。主要包括数据低质量及污染序列过滤、双端合并、筛选扩增子并统一序列方向、按barcode拆分样品、去除5’引物序列、去除3’引物序列为下一步分析的高质量样本序列；(B) 折线图展示各测序文库中序列的长度分析。Data size and length distribution of sequencing libraries. (A) Bar plot showing reads count of each library in read filter process. (B) Line plot showing reads length distribution of each library. [Sum PDF](result/stat_lib_qc_sum.pdf)  [Length PDF](result/stat_lib_length.pdf)

\`\`\`{r quality-sum, fig.cap="(ref:quality-sum)"}
knitr::include_graphics(c("result/stat_lib_qc_sum.png","result/stat_lib_length.png"))
\`\`\`
END



if test "${execute}" == "TRUE";
then
	Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
fi
