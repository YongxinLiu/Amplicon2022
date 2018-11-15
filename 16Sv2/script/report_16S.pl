#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq!
Usage:    rmd_16s.pl -s summary.txt -d design.txt -l library.txt -c group_compare.txt -v group_venn.txt -t group_tern.txt -g group1 -b output_dir -S elite_report -a abundance_threshold  
Function: Write 16S report in Rmarkdown
Command:  -d design.txt
          -g group1
          -D group2_list, for samples in different batch or condition
          -F group3_list, for samples in different batch or condition
          -l library.txt
          -v venn.txt
          -s summary.txt
          -c compare.txt
          -t tern.txt
          -a abundance threshold, default 0.005
          -b version for output directory
          -S TRUE or FLASE, whether report simplified report, default FALSE
          -h header line number, default 0
          -p pvalue threshold
          -q qvalue or FDR threshold
Author:   Liu Yong-Xin, yxliu\@genetics.ac.cn, QQ:42789409
Version:  v1.0
Update:   2018/5/2
Notes:    1.0 show alpha, beta and taxonomy
\n!
    )
}
# 	time rmd_16s.pl -s ${summary} -d ${design} -l ${library} -c ${compare} -v ${venn} -t ${tern} -g ${g1} -b ${version} -S ${elite_report} -a ${thre}

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'a:i:o:f:p:q:d:e:l:c:v:m:n:h:s:g:D:F:t:b:S:', \%opts );
$opts{a}=0.01 unless defined($opts{a});
$opts{f}=1.2 unless defined($opts{f});
$opts{p}=0.05 unless defined($opts{p});
$opts{q}=0.1 unless defined($opts{q});
$opts{h}=0 unless defined($opts{h});
$opts{o}="./" unless defined($opts{o}); # work directory
$opts{s}="doc/summary.txt" unless defined($opts{s});
$opts{d}="doc/design.txt" unless defined($opts{d});
$opts{l}="doc/library.txt" unless defined($opts{l});
$opts{c}="doc/compare.txt" unless defined($opts{c});
$opts{m}="edgeR" unless defined($opts{m});
$opts{v}="doc/venn.txt" unless defined($opts{v});
$opts{t}="doc/tern.txt" unless defined($opts{t});
$opts{e}="TRUE" unless defined($opts{e}); # report elite version
$opts{g}="groupID" unless defined($opts{g}); # default group column name
$opts{S}="FALSE" unless defined($opts{S}); # report elite version

&usage unless (exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));

my %list;
if (-e "$opts{s}") {
open LIST,"<$opts{s}";
while (<LIST>) {
	chomp;
	my @tmp=split/\t/;
	$list{$tmp[0]}=$tmp[1];
}
}
close LIST;

# Set project default parameter
$list{title}="Analysis report of 16s rDNA Sequencing" unless defined($list{title});
$list{client}="Bailab member" unless defined($list{client});
$list{partner}=$list{client} unless defined($list{partner});
$list{analyst}="Dr. Yong-Xin Liu, Bailab, IGDB, CAS" unless defined($list{analyst});
$list{contact}="yxliu\@genetics.ac.cn 010-64808722" unless defined($list{contact});
$list{project}="170901-1" unless defined($list{project});
$list{period}="2018-05-01~ 2018-12-31" unless defined($list{period});
$list{website}="http://bailab.genetics.ac.cn/" unless defined($list{website});
$list{wechat}="[植物微生物组](http://mp.weixin.qq.com/s/QkgNlzK_rpauKzSs2fUd7A)" unless defined($list{wechat});
$list{logo}="Bailab, SKLPG/CEPAMS, IGDB, CAS {-}" unless defined($list{logo});

# Prepare relative files
print "Clean report enviroment and prepare relative files\n";
`rm -fr *.Rmd _bookdown_files $opts{b}`;
`cp -f /mnt/bai/yongxin/github/Amplicon/16Sv2/rmd/* ./`;
`sed -i 's/html/$opts{b}/g' _bookdown.yml`;


open OUTPUT,">$opts{o}index.Rmd";
print OUTPUT qq!
--- 
title: "$list{title}"
author:
- Partner实验数据：$list{partner}
- Analyst生信分析：$list{analyst}
- Contact联系方式：$list{contact}
- Project项目编号：$list{project}
- Period项目周期：$list{period}
- Website官方网站：$list{website}
- Wechat公众号：$list{wechat}
date: '`r Sys.Date()`'
documentclass: article
bibliography: [16s.bib]
link-citations: yes
biblio-style: apalike
---

```{r setup, include=FALSE}
library(knitr)
output = opts_knit\$get("rmarkdown.pandoc.to")
html = FALSE
latex = FALSE
opts_chunk\$set(echo = FALSE, out.width="100%", fig.align="center", fig.show="hold", warning=FALSE, message=FALSE)
if (output=="html") {
	html = TRUE
}
if (output=="latex") {
	opts_chunk\$set(out.width="95%", out.height='0.7\\textheight', out.extra='keepaspectratio', fig.pos='H')
	latex = TRUE
}
knitr::opts_chunk\$set(cache=TRUE, autodep=TRUE)
mtime = function(files){
  lapply(Sys.glob(files), function(x) file.info(x)\$mtime)
}
set.seed(315)
```

```{asis, echo=html}
# $list{logo}
```

```{r cover, eval=html, out.width="99%"}
knitr::include_graphics("http://bailab.genetics.ac.cn/img/banner.png")
```
!;
close OUTPUT;



# 第2部分：a. 多样性分析

open OUTPUT,">$opts{o}2a-diversity.Rmd";
print OUTPUT qq!

# 多样性分析 Diversity {#result-diversity}

**多样性指数 Alpha diversity index**

Usearch计算各样品常用14种Alpha多样性计算方法结果见如 Table \\\@ref(tab:alpha) 所示。[TXT](result/alpha/index.txt)

```{r alpha}
table_alpha = read.table("result/alpha/index.txt", sep="\\t", header=T, row.names = 1)
table_alpha=head(table_alpha,n=3)
knitr::kable(table_alpha, caption="样品14种Alpha多样性结果前3行预览", booktabs=TRUE)
```

## Alpha箱线图 Boxplot {#diversity-alpha-boxplot}

(ref:div-alpha) 箱线图展示各样品及组的微生物组Alpha多样性，方法采用(A) Shannon index，包括样品的可操作分类单元(operational taxonomic unit, OTU)种类(richness)及丰度(evenness)信息；(B) Observed OTUs (Richness) index，只包括样品OTU种类信息。(C) Chao1 index,基于样品测序中单拷贝OTU(饱合情况)估算样品物种种类的方法; 各组间分组采用R语言agricolae包的LSD.test函数统计，两组上方无相同字母代表组间存在显著差异(pvalue < 0.05)。附文本有t-test方法统计各组间是否存在显著差异的p-value水平。 Within sample diversity (α-diversity) measurements among each genotype. (A) Shannon index, estimated species richness and evenness; (B) Observed OTUs / Richness index, only calculate species richness; (C) Chao1 index, calculate richness based on observed, singletons and doubletons.
 Shannon [PDF](result/alpha/shannon_e.pdf) [TXT](result/alpha/shannon_e.txt)  Richness [PDF](result/alpha/richness.pdf) [TXT](result/alpha/richness.txt)  Chao1 [PDF](result/alpha/chao1.pdf) [TXT](result/alpha/chao1.txt)

```{r div-alpha, fig.cap="(ref:div-alpha)", out.width="99%"}
figs_2 = paste0("result/alpha/", c("shannon_e", "richness","chao1"),".png")
knitr::include_graphics(figs_2)
```


## Alpha稀释曲线 Rarefracation curve {#diversity-alpha-rare}

(ref:div-rare) 折线图展示各样品(A)及组(B)的物种丰富度稀释曲线。所有样品抽样基于数据量为10000条的标准化OTUs表，再进行重抽样1% ~ 100%的100个梯度，观测OTUs的数量，采用ggplot2绘制拆线图展示。
Rarefraction curve of richness (Observed OTUs) in samples and groups. [TXT](result/alpha/rare.txt)  [Samples PDF](result/alpha/rare_samples.pdf)  [Groups PDF](result/alpha/rare_groups.pdf)  

```{r div-rare, fig.cap="(ref:div-rare)", out.width="99%"}
figs_2 = paste0("result/alpha/rare_", c("samples", "groups"),".png")
knitr::include_graphics(figs_2)
```

## 主坐标轴分析 PCoA {#diversity-pcoa}

(ref:div-beta) 主坐标轴分析(PCoA)展示第1/2坐标轴下各样品间微生物组差异(dissimilarity)，距离计算方法采用(A) Bray-Curtis; (B) Unifrac; (C) Unweighted Unifrac。采用Adonis统计各样品组间的显著性差异P值(见Stat链接)。
Principal coordinate analysis (PCoA) using the (A) bray curtis metric, (B) unweighted unifrac metric and (C) weighted unifrac metric shows dissimilarity of microbial communities. bray_curtis [PC1/2](result/beta/bray_curtis.pdf) [PC3/4](result/beta/bray_curtis_34.pdf) [Stat](result/beta/bray_curtis.stat) [Label](result/beta/bray_curtis_label.pdf) weighted_unifrac [PC1/2](result/beta/weighted_unifrac.pdf) [PC3/4](result/beta/weighted_unifrac_34.pdf) [Stat](result/beta/weighted_unifrac.stat) [Label](result/beta/weighted_unifrac_label.pdf) unweighted_unifrac [PC1/2](result/beta/unweighted_unifrac.pdf) [PC3/4](result/beta/unweighted_unifrac_34.pdf) [Stat](result/beta/unweighted_unifrac.stat) [Label](result/beta/unweighted_unifrac_label.pdf)

```{r div-beta, fig.cap="(ref:div-beta)", out.width="99%"}
figs_2 = paste0("result/beta/", c("bray_curtis", "weighted_unifrac"),".png")
knitr::include_graphics(figs_2)
```

!;


$file = "result/beta/cpcoa_bray.pdf";
if (-e $file) {
print OUTPUT qq!
	
## 限制性 PCoA

(ref:div-CPCoA) 以基因型为条件分析其贡献率和样品组间差异(筛选至少在一个组中位数存在OTU丰度 > 0.05% )。variance代表当前分组条件下各样品间差异所占的比重或贡献率，P值示基因型各组间是否存在显著差异，各样品间距离计算方法为默认为Bray-Curtis距离，可选Jaccard距离。
Constrained principal coordinate analysis on bacterial microbiota (Only OTUs abundance median of one gorup > 0.05% were kept). Variation between samples in Bray-Curtis distances constrained by groupID. (Bulgarelli et al., 2015).[Bray-Curtis PDF](result/beta/cpcoa_bray.pdf)  [Bray-Curtis PDF labels](result/beta/cpcoa_bray_label.pdf) [jaccard PDF](result/beta/cpcoa_jaccard.pdf)  [jaccard PDF labels](result/beta/cpcoa_jaccard_label.pdf)  


```{r div-CPCoA, fig.cap="(ref:div-CPCoA)", out.width="99%"}
knitr::include_graphics("result/beta/cpcoa_bray.png")
```

!;

}
close OUTPUT;


## 样品组比较
open DATABASE,"<$opts{c}";
while (<DATABASE>) {
	chomp;
	#print "$_\n";	
	push @group,$_;
}
close DATABASE;


@tax_en=qw#p c o f g#;
@tax_cn=qw#门 纲 目 科 属#;

open OUTPUT,">$opts{o}2a-taxonomy.Rmd";
print OUTPUT qq!

# 生物分类学 Taxonomy {#result-taxonomy}

!;

foreach $i (0..4) {

# 开启精简模式，只输出门、目、属水平差异，纲、科跳过
if ($opts{S} eq "TRUE") {
	next if $i==1;
	next if $i==3;
}

print OUTPUT qq!
## $tax_cn[$i]水平差异 {#result-taxonomy-$tax_en[$i]}

(ref:taxonomy-$tax_en[$i]) 柱状图展示各样品组微生物组分类学$tax_cn[$i]水平相对丰度。(A) 堆叠柱状图展示各样品相对丰度。(B) 堆叠柱状图展示各组平均相对丰度，X轴为各样品组，Y轴为各$tax_cn[$i]类相对百分比，只列出了丰度大于0.1%的$tax_cn[$i]，其它所有$tax_cn[$i]归入Low Abundance类。(C) 条形图展示最高丰度的五大菌$tax_cn[$i]平均丰度及标准误，我们可以观察各样品组$tax_cn[$i]水平上相关丰度的差异及组内生物学重复间的波动范围。[stack sample PDF](result/tax/sum_$tax_en[$i]_sample.pdf)  [stack group PDF](result/tax/sum_$tax_en[$i]_group.pdf)  [raw Data](result/tax/sum_$tax_en[$i].txt)

```{r taxonomy-$tax_en[$i], fig.cap="(ref:taxonomy-$tax_en[$i])", out.width="99%"}
figs_2 = paste0("result/tax/sum_$tax_en[$i]_", c("sample","group"),".png")
knitr::include_graphics(figs_2)
```

组间差异$tax_en[$i]水平概述

样品组间显著差异OTUs数量(组内相对丰度中位数 > $opts{a}%, Fold-Change > $opts{f}, P-value < $opts{p}, FDR < $opts{q}, 统计方法为Wilcoxon rank sum test)如 Table \\\@ref(tab:otu-sum-$tax_en[$i]) 所示。[TXT](result/compare_$tax_en[$i]/summary.txt)；

```{r otu-sum-$tax_en[$i]}
table_otu = read.table("result/compare_$tax_en[$i]/summary.txt", sep="\\t", header=T)
knitr::kable(table_otu, caption="各样品组间差异OTUs数量汇总", booktabs=TRUE)
```

!;
foreach (@group) {
	#print $_,"\n";
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
	#print $tmp[0],"\n",$tmp[1],"\n";
	$file = "result/compare_$tax_en[$i]/$tmp[0]-$tmp[1]_sig.txt";

if (-e $file) {
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

$tmp[0]与$tmp[1]相比显著差异的分类单元信息如 Table \\\@ref(tab:taxonomy-$tmp[0]vs$tmp[1]-$tax_en[$i]) 所示。[All TXT](result/compare_$tax_en[$i]/$tmp[0]-$tmp[1]_all.txt) [Significant TXT](result/compare_$tax_en[$i]/$tmp[0]-$tmp[1]_sig.txt)

```{r taxonomy-$tmp[0]vs$tmp[1]-$tax_en[$i]}
m = read.table("result/compare_$tax_en[$i]/$tmp[0]-$tmp[1]_sig.txt", sep="\\t", header=T, row.names = 1, comment.char = "")
e = m[m\$logFC>0,]
e = head(e[order(-e\$MeanA),],n=10)
d = m[m\$logFC<0,]
d = head(d[order(-d\$MeanB),],n=10)
m = rbind(e,d)
m=m[,c("MeanA","MeanB","logFC","PValue","FDR")]
knitr::kable(m, row.names=T, caption="$tmp[0]与$tmp[1]显著差异的$tax_cn[$i]", booktabs=TRUE)
```

!;
}else{
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

无显著差异丰度分类单元；No significantlly differentially abundance taxonomy.

!;
}
}

# 维恩图
if (-e "$opts{v}") {
open DATABASE,"<$opts{v}";
my @venn;
while (<DATABASE>) {
	chomp;
	push @venn,$_;
}
close DATABASE;


## 如果venn非空，读列表文件
if (@venn>=1) {
print OUTPUT qq!
### $tax_cn[$i]维恩图  {#result-$tax_en[$i]-venn}

!;

my $j=0;
foreach (@venn) {
	chomp;
	$j++;
	my @tmp=split/\t/; # sampleA and sampleB
	my $venn_list2;
	foreach $tmp (@tmp) {
		$venn_list2.=$tmp;
	}
	$tmp[2]="C" unless defined($tmp[2]);
	$tmp[3]="D" unless defined($tmp[3]);
	$tmp[4]="E" unless defined($tmp[4]);
	$venn_list="$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]";

# 判断是否有维恩图，再输出
$file = "result/compare_$tax_en[$i]/diff.list.venn$venn_list.pdf";
if (-e $file) {
print OUTPUT qq!

#### $venn_list

(ref:$tax_en[$i]-venn-$j) 维恩图展示各比较组差异$tax_cn[$i]水平共有和特有数量。Venn diagrams show common and unique OTUs in each group. [PDF](result/compare_$tax_en[$i]/diff.list.venn$venn_list.pdf) [TXT](result/compare_$tax_en[$i]/diff.list.venn$venn_list2.xls.xls)  


```{r $tax_en[$i]-venn-$j, fig.cap="(ref:$tax_en[$i]-venn-$j)", out.width="99%"}
figs_2 = paste0("result/compare_$tax_en[$i]/diff.list.venn", "$venn_list", ".png")
knitr::include_graphics(figs_2)
```

!;
}
}
}
}

}





#open OUTPUT,">$opts{o}2b-taxonomy.Rmd";
#print OUTPUT qq!
## 目差异分析 {#result-o}
#
### 差异目概述 {#result-o-sum}
#
#样品组间显著差异目数量(P < 0.001, FDR < 0.05，统计方法为$opts{m})如 Table \\\@ref(tab:o-sum) 所示。[TXT](result/compare_o/summary.txt)；
#
#```{r o-sum}
#table_o = read.table("result/compare_o/summary.txt", sep="\t", header=T)
#knitr::kable(table_o, caption="各样品组间差异目数量汇总", booktabs=TRUE)
#```
#
#!;
#my $i=0;
#foreach (@group) {
#	chomp;
#	my @tmp=split/\t/; # sampleA and sampleB
#	$i++;
#
#print OUTPUT qq!
#
### $tmp[0] vs $tmp[1]
#
#(ref:o-$i) $tmp[0] vs $tmp[1]组间相对丰度显著差异目 (Pvalue < 0.001 & FDR < 0.05, GLM likelihood rate test(edgeR), or wilcoxon rank test)。(A) 火山图展示两组比较目的变化，x轴为o差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝分别代表显著上下调，灰为没有显著变化的目；(B) 热图展示$tmp[0]与$tmp[1]显著差异o在每个样品中丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低，黄色代表中间水平；(C) 曼哈顿图展示o的变化情况及在各门水平中的分布，x轴为o按物种门水平物种注释字母排序，y轴为Pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表o，颜色为门水平注释，大小为相对丰度，形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调；(D) 曼哈顿图按目水平上色。
#$tmp[0] are enriched and depleted for certain 目 (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of abundance and fold change of 目; (B) Heatmap showing differentially abundance 目; (C) Manhattan plot showing phylum pattern of differentially abundance 目, colored by phylum; (D) Manhattan plot colored by order.
#[Heatmap](result/compare_o/$tmp[0]-$tmp[1]_heatmap.pdf)
#
#```{r o-$i, fig.cap="(ref:o-$i)", out.width="99%"}
#figs = paste0("result/compare_o/$tmp[0]-$tmp[1]_", c("heatmap"),".png")
#knitr::include_graphics(figs)
#```
#
#$tmp[0]与$tmp[1]相比显著差异的目。A_mean, B_mean分别为各组内相对丰度均值的百分比(%)，如 Table \\\@ref(tab:tab-o-$i) 所示。[All](result/compare_o/$tmp[0]-$tmp[1]_all.txt)  [Significant](result/compare_o/$tmp[0]-$tmp[1]_sig.txt)
#
#```{r tab-o-$i}
#m = read.table("result/compare_o/$tmp[0]-$tmp[1]_sig.txt", sep="\\t", header=T, row.names = 1)
#e = m[m\$logFC>0,]
#e = head(e[order(-e\$MeanA),],n=10)
#d = m[m\$logFC<0,]
#d = head(d[order(-d\$MeanB),],n=10)
#m = rbind(e,d)
#m=m[,c("MeanA","MeanB","logFC","PValue","FDR")]
#knitr::kable(m, row.names=T, caption="$tmp[0]与$tmp[1]显著差异的目；Significantlly different 目.", booktabs=TRUE)
#```
#
#!;
#}

open OUTPUT,">$opts{o}2c-compare.Rmd";
print OUTPUT qq!
# OTUs差异分析 {#result-otu}

## 差异OTUs概述 {#result-otu-sum}

样品组间显著差异OTUs数量(组内Abundance中位数 > $opts{a}%, Fold-Change > $opts{f}, P-value < $opts{p}, FDR < $opts{q}, 统计方法为$opts{m})如 Table \\\@ref(tab:otu-sum) 所示。[TXT](result/compare/summary.txt)；

```{r otu-sum}
table_otu = read.table("result/compare/summary.txt", sep="\\t", header=T)
knitr::kable(table_otu, caption="各样品组间差异OTUs数量汇总", booktabs=TRUE)
```

!;
my $i=0;
foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
	$i++;
	$file = "result/compare/$tmp[0]-$tmp[1]_sig.txt";

if (-e $file) {
print OUTPUT qq!

## $tmp[0] vs $tmp[1]

(ref:otu-$i) $tmp[0] vs $tmp[1]组间相对丰度显著差异OTUs (Pvalue < 0.001 & FDR < 0.05, GLM likelihood rate test(edgeR), or wilcoxon rank test)。(A) 火山图展示两组比较OTUs的变化，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝分别代表显著上下调，灰为没有显著变化的OTUs；(B) 热图展示$tmp[0]与$tmp[1]显著差异OTU在每个样品中丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低，黄色代表中间水平；(C) 曼哈顿图展示OTU的变化情况及在各门水平中的分布，x轴为OTU按物种门水平物种注释字母排序，y轴为Pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，颜色为门水平注释，大小为相对丰度，形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调；(D) 曼哈顿图按目水平上色。
$tmp[0] are enriched and depleted for certain OTUs (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of abundance and fold change of OTUs; (B) Heatmap showing differentially abundance OTUs; (C) Manhattan plot showing phylum pattern of differentially abundance OTUs, colored by phylum; (D) Manhattan plot colored by order.
[VolcanoPlot](result/compare/$tmp[0]-$tmp[1]_volcano.pdf)  [Heatmap](result/compare/$tmp[0]-$tmp[1]_heatmap.pdf) [Mamnahttan plot](result/compare/$tmp[0]-$tmp[1]_all.txt_man_pc.pdf)

```{r otu-$i, fig.cap="(ref:otu-$i)", out.width="99%"}
figs = paste0("result/compare/$tmp[0]-$tmp[1]_", c("volcano","heatmap","all.txt_man_pc"),".png")
knitr::include_graphics(figs)
```

$tmp[0]与$tmp[1]相比显著差异的OTUs。A_mean, B_mean分别为各组内相对丰度均值的百分比(%)，如 Table \\\@ref(tab:tab-otu-$i) 所示。[All](result/compare/$tmp[0]-$tmp[1]_all.txt)  [Significant](result/compare/$tmp[0]-$tmp[1]_sig.txt)

```{r tab-otu-$i}
m = read.table("result/compare/$tmp[0]-$tmp[1]_sig.txt", sep="\\t", header=T, row.names = 1)
e = m[m\$logFC>0,]
e = head(e[order(-e\$MeanA),],n=10)
d = m[m\$logFC<0,]
d = head(d[order(-d\$MeanB),],n=10)
m = rbind(e,d)
m=m[,c("MeanA","MeanB","logFC","PValue","FDR","Phylum","Class","Order","Family","Genus")]
knitr::kable(m, row.names=T, caption="$tmp[0]与$tmp[1]显著差异的OTUs；Significantlly different OTUs.", booktabs=TRUE)
```

!;
}else{
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

无显著差异丰度分类单元；No significantlly differentially abundance taxonomy.

!;
}
}


my @venn;
my @venn2;
## 如果venn.txt文件存在，读列表文件
if (-e "$opts{v}") {
open DATABASE,"<$opts{v}";
while (<DATABASE>) {
	chomp;
	push @venn,$_;
}
close DATABASE;


## 如果venn非空，读列表文件
if (@venn>=1) {
print OUTPUT qq!
## OTUs维恩图  {#result-otu-venn}

!;

my $j=0;
foreach (@venn) {
	chomp;
	$j++;
	my @tmp=split/\t/; # sampleA and sampleB
	my $venn_list2;
	foreach $tmp (@tmp) {
		$venn_list2.=$tmp;
	}
	$tmp[2]="C" unless defined($tmp[2]);
	$tmp[3]="D" unless defined($tmp[3]);
	$tmp[4]="E" unless defined($tmp[4]);
	$venn_list="$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]";

print OUTPUT qq!

### $venn_list

(ref:otu-venn-$j) 维恩图展示各比较组差异OTU的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [PDF](result/compare/diff.list.venn$venn_list.pdf) [TXT](result/compare/diff.list.venn$venn_list2.xls.xls)  


```{r otu-venn-$j, fig.cap="(ref:otu-venn-$j)", out.width="99%"}
figs_2 = paste0("result/compare/diff.list.venn", "$venn_list", ".png")
knitr::include_graphics(figs_2)
```

!;
}
}
}

open OUTPUT,">$opts{o}4a-culture.Rmd";

print OUTPUT qq!

# 个性分析

!;

$file = "result/41culture/otu.txt";
if (-e $file) {
print OUTPUT qq!
## 可培养OTUs

OTUs与培养菌序列相似度大于97%(可培养)统计信息如 Table \\\@ref(tab:culture-sum) 所示。[Summary table](result/41culture/summary.txt)

```{r culture-sum}
tab = read.table("result/41culture/summary.txt", sep="\\t", header=F)
knitr::kable(tab, caption="可培养统计", booktabs=TRUE)
```

OTUs与培养菌序列相似度列表如 Table \\\@ref(tab:culture-tab) 所示。[Summary table](result/41culture/otu.txt)

```{r culture-tab}
tab = read.table("result/41culture/otu.txt", sep="\\t", header=T)
knitr::kable(head(tab, n= 30), caption="TUs与培养菌序列相似度列表Top30", booktabs=TRUE)
```

!;
}

close OUTPUT;
#(ref:otu-graphlan) OTUs物种树，外圈颜色势图代表OTUs丰度，内圈黑块标注此OTUs是否存在97%以上相似的菌保，仅相对丰度0.1%的OTUs展示在本图中[PDF](filter__k1/graphlan.pdf)。图中展示可培养数量和相比分析比例，详见[TXT](filter__k1/culture.txt)。
#Taxonomy tree showing high abundance OTUs. Outer circle heatmap show relative abundance ot OTUs, and inner circle of black block show OTUs have cultured similar stocks.
#
#```{r otu-graphlan, fig.cap="(ref:otu-graphlan)", out.width="99%"}
#figs_1 = paste0("filter__k1/graphlan",".png")
#knitr::include_graphics(figs_1)
#```


open OUTPUT,">$opts{o}9-references.Rmd";
print OUTPUT qq!
`r if (knitr:::is_html_output()) '# References {-}'`

# 附录

## 分析参数和配置信息

- [实验设计](doc/design.txt)
- [分析参数](makefile)
- [比较组](doc/compare.txt)
- [韦恩图](doc/venn.txt)

## 分析结果重要文件下载

- OTU表 [Count TXT](result/otutab.txt) [Count biom](result/otutab.biom) [抽样1万标准化 TXT](result/otutab_norm.txt) [抽样1万标准化 biom](result/otutab_norm.biom) [Greengene Count txt](result/otutab_gg.txt)  
- 代表性序列 [fasta](result/otu.fa) [tree](result/otu.tree)
- 物种注释 [2列格式](result/taxonomy_2.txt) [8列格式](result/taxonomy_8.txt)

!;
close OUTPUT;

# 记录程序运行时间
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

# 备份分析参数文件
#$parameter=`readlink  makefile`;
#`cp $parameter $opts{b}/makefile.txt`;


# 编译Rmarkdown为HTML
if ($opts{e} eq "TRUE") {
	`Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"`;
}

# 添加htaccess方式加密访问
open ACCESS,">$opts{b}/.htaccess";
print ACCESS qq!
# AuthName must have, "" not allow blank, not need change
AuthName "Need user name and password"
AuthType Basic
AuthUserFile /mnt/bai/yongxin/bin/config/users
require valid-user
!;
`chmod +x $opts{b}/.htaccess`;

print "Result please visiting http://210.75.224.110/report/16Sv2/$opts{b}\n";
