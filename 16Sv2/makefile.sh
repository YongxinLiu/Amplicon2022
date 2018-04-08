SHELL:=/bin/bash

# 16S扩增子第二版流程配置文件，请根据项目具体情况修改 v2.0 2018/4/3
# Config file of 16S Amplicon pipeline version 2, please modify according to experiment design. v2.0 2018/4/3



# 1. 测序文库参数 Standard pipeline parameter

# 工作目录 working directory
# 修改wd为当前工作目录pwd
wd=`pwd`
# make init # 建立分析所需子目录
# 设置任务最大运行任务/线程数，超过CPU数量效率反而会降低
p=36
# 数据库
# Greengene 13 May database, fa for usearch format, udb for usearch index
usearch_gg=/mnt/bai/public/ref/gg_13_5_otus/97_otus_usearch.udb
## Silva 132 database, fa for usearch format, udb for usearch index
usearch_silva=/mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.udb
usearch_rdp=/mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb

# 1.1 lane_split 拆分下机数据为文库 Split lane into library
# 需要将lane文件放入seq目录，对应的index和文库放在doc/library.txt
# 注意：务必查看文库文件中Index具体格式，默认为#Index，其它情况需修改主流程源代码main_pipeine.sh
# 文库名
lane=lane

# 1.2 library_split 拆分文库为样品 Split library into sample
# 默认只拆分单左端barcode类型的样品:先匹配左端，再提取序列ID，再提取右端，最后改名，注意实验设计要严格规范无空格

# 1.3 sample_merge 双端序列合并 Merge pair-end reads
# 如果是pair-end reads是phred64，需要先使用fastp转换为33且关闭质控(质控影响序列长度)，再使用usearch10 mergepair

# 1.4 fq_trim 切除引物与条型码 Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
stripleft=29
stripright=18

# 1.5 fq_qc 质量控制fastq filter
# 默认错误率<0.01 keep reads error rates less than 1%
fastq_maxee_rate=0.01

# 1.6 fa_unqiue 序列去冗余 Remove redundancy
# 最小序列频率miniuniqusize默认为8，去除低丰度，增加计算速度，整lane的序列可更改为30，甚至100
minuniquesize=30

# 1.7 otu_pick 挑选OTU Pick OTUs
# 可选97% cluster_otus 和 unoise3 ，默认unoise3
otu_method=cluster_otus
# OTU日志文件，记录从挑选至过滤为最终的过程
otu_log=result/otu.log

# 1.8 chimera_ref 参考去嵌合 remove chimera
# 此处推荐使用大数据，如SILVA132，其它数据库如rdp_gold是错误的
# /mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.fasta # 99%非冗余1.1G, 内存使用5.8G, 8min, 16.2% chimeras
# /mnt/bai/public/ref/silva/SILVA_132_SSURef_tax_silva.fasta # 全部3.3G, 内存使用15.9G, 30min, 16.2% chimeras
# 使用非冗余99%的省内存和时间，但结果差不多
chimera_ref=${usearch_silva}
# 模式，嵌合体比例由大到小: high_confidence specific balanced sensitive sensitive
chimera_mode=balanced

# 1.9 host_rm 去宿主 remove host
# 去宿主方法选择 blast / sintax_gg / sintax_silva，推荐：sintax_silva
host_method=sintax_silva
# 方法1. blast宿主基因组(含叶绿体/线粒体)去除同源序列，如水稻微生物，需要提供水稻基因组；可调相似度和覆盖度的阈值(百分数)
host=/mnt/bai/public/ref/rice/msu7/all.con
host_similarity=90
host_coverage=90
# 方法2. 基于gg注释结果筛选, 去除叶绿体Chloroplast、线粒体mitochondria，默认为usearch_gg数据库
# 方法3. silva注释可识线粒体Mitochondria、叶绿体Chloroplast和真核生物Eukaryota(包括宿主、真菌、原生动物等)，默认为usearch_silva数据库

# 1.10 otutab_create 生成OTU表 Creat OTUs table
# 有 usearch10 和 vsearch 两个软件可选，默认usearch10，vsearch多线程会更快些
map_method=vsearch
map_identify=0.97

# 1.11 otutab_filter OTU表筛选 Filter OTU table
# OTU表筛选日志文件
log_otutable=result/otutab.log
# 按样本量筛选，默认5000，根据otu_stats结果调整
min_sample_size=5000
# 按矩阵中每个点count, freq筛选，低于阈值变为0
# 按OTU丰度和频率筛选，如OTU测序量至少8次，相对丰度百万分之一(建议挑选序列去冗余部分调高阈值更合理)
min_otu_size=8
# 按频率筛选，推荐十万分之一0.00001，范围千一至百分一0.001 - 0.000001之间
min_otu_freq=0.000001
# 抽样标准化的值，推荐最小10000，根据统计结果选择筛选后最小值或可保留大部分样品的值
sample_size=30000

# 1.12 tax_assign 物种注释 Assign taxonomy
# 物种注释推荐使用小而准的数据库，如rdp trainset 16(由Robert整理)
# 可选gg, silva，分别从官网下载并shell调整格式
sintax_db=${usearch_silva}
# 分类准确度阈值，默认0.8，注释太少最小可改0.5，发现有明显错误可最高上升为0.95
sintax_cutoff=0.8

# 1.13 tax_sum 物种注释统计 Taxonomy summary
# 按门、纲、目、科、属水平分类汇总，结果位于：result/tax/sum_*.txt
# 输出可读的taxonomy与OTU对应文件，有2列和8列版，result/taxonomy_2/8.txt

# 1.14 tree_make 多序列比对和进化树 Multiply alignment and make_phylogeny

# 1.15 alpha_calc Alpha多样性指数计算 Calculate alpha diversity index
# alpha指数计算结果为 result/alpha/index.txt
# 稀释梯度抽样方法 richness (observed OTUs)-method fast / with_replacement / without_replacement , 结果位于 result/alpha/rare.txt
rare_method=without_replacement

# 1.16 beta_calc Beta多样性进化树和距离矩阵计算 Beta diversity tree and distance matrix
# 距离矩阵计算方法，34种可选： abund_jaccard, binary_chisq, binary_chord, binary_euclidean, binary_hamming, binary_jaccard, binary_lennon, binary_ochiai, binary_otu_gain, binary_pearson, binary_sorensen_dice, bray_curtis, bray_curtis_faith, bray_curtis_magurran, canberra, chisq, chord, euclidean, gower, hellinger, kulczynski, manhattan, morisita_horn, pearson, soergel, spearman_approx, specprof, unifrac, unifrac_g, unifrac_g_full_tree, unweighted_unifrac, unweighted_unifrac_full_tree, weighted_normalized_unifrac, weighted_unifrac
# 默认使用4种
dis_method=bray_curtis,binary_jaccard,weighted_unifrac,unweighted_unifrac

# 1.17 otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
# 比对方法和相似度同1.10 mapping
otutab_gg=/mnt/bai/public/ref/gg_13_5_otus/rep_set/97_otus.fasta



# 2. 统计绘图

# 绘图通用参数
# 实验设计文件位置，全局，其它图默认调此变量，也可单独修改；并选择表中的组列和具体分组
design=${wd}/doc/design.txt
g1=groupID
# tail -n+2 doc/design.txt|cut -f 5|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
# 绘图使用的实验组，顺序即图中显示顺序；为空时使用所有组和默认顺序
g1_list='"Col","ThasKO2","ThahKO","ThadKO","ACT2KO"'

# 组间比较列表
compare=${wd}/doc/compare.txt
# 组间共有、特有比较列表
venn=${wd}/doc/venn.txt
# 图片长宽，按nature全面版、半版页面设置
# 图片长宽和字体大小，7组以下用默认，7组以上改为8x5或更大； figure size, recommend 4x2.5, 5x3(default), 8x5, 16x10, text_size 6, 7(default), 8
width=8
height=5
text_size=7

# 图中显示legend, 如taxonomy的数量，5，8(default)，10
legend_number=10
# 差异统计按丰度过滤 abundance filter，如丰度按万分之五过滤，减少计算量，提高OTU的FDR值
abundance_thre=0.0005
# 差异比较方法，默认是edgeR的 lrt ，可选 wilcoxon 秩和检验
compare_method="lrt"
# 显著性P值过滤 threshold of P-value，可选0.05, 0.01, 0.001。采用FDR校正，此参数意义不大，即使0.001也没有FDR < 0.2过滤严格
pvalue=0.01
# 统计检验方式fdr
fdr=0.05
# 差异倍数logFC常用1.5, 2, 4倍，对应0.585, 1, 2；菌丰度变化倍数不明显，还可用1.3和1.7倍对应0.379和0.766
logFC=0.379

# 统计绘图和网页报告版本控制
sub="v1"
doc=doc/${sub}
# 报告输出目录
version=report_AC_${sub}


# 2.1 alpha_boxplot Alpha多样性指数箱线图 Alpha index in boxplot
# alpha箱线图绘制参数
ab_input=${wd}/result/alpha/index.txt
# alpha指数种类14种 head -n1 result/alpha/index.txt | tr '\t' '\n'|tail -n+2|awk '{print "\""$1"\""}'|tr "\n" ","
# "berger_parker","buzas_gibson","chao1","dominance","equitability","jost","jost1","reads","richness","robbins","simpson","shannon_e","shannon_2","shannon_10"
# 默认使用chao1, richness和shannon_e三种，可修改ab_type添加或减少
ab_method='"chao1","dominance","richness","simpson","shannon_e"'
ab_design=${design}
ab_group_name=${g1}
ab_group_list=${g1_list}
ab_output=${wd}/result/alpha/
ab_width=${width}
ab_height=${height}

# 2.2 alpha_rare Alpha丰富度稀释曲线 Alpha rarefracation curve
ar_input=${wd}/result/alpha/rare.txt
ar_design=${design}
ar_group_name=${g1}
ar_group_list=${g1_list}
ar_output=${wd}/result/alpha/
ar_width=${width}
ar_height=${height}

# 2.3 beta_pcoa 主坐标轴分析距离矩阵 PCoA of distance matrix
bp_input=${wd}/result/beta/
bp_method='"binary_jaccard","bray_curtis","unweighted_unifrac","weighted_unifrac"'
bp_design=${design}
bp_group_name=${g1}
bp_group_list=${g1_list}
bp_output=${wd}/result/beta/
bp_width=${width}
bp_height=${height}
# 散点图是否按组添加置信椭圆，TRUE添加，FALSE不添加，默认T
bp_ellipse=TRUE
# 实验比较组，可用默认，也可设置单独文件，没有则不计算
bp_compare=${compare}

# 2.4 beta_cpcoa 限制性主坐标轴分析: OTU表基于bray距离和CCA  CCA of bray distance matrix
# 输入的OTU表，可原始count，也可以标准化的结果，无差异
bc_input=${wd}/result/otutab.txt
# Method from vegdist() of vegan: "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis"
bc_method='"bray","jaccard"'
bc_design=${design}
bc_group_name=${g1}
bc_group_list=${g1_list}
bc_output=${wd}/result/beta/
bc_width=${width}
bc_height=${height}
# 散点图是否按组添加置信椭圆，TRUE添加，FALSE不添加，默认T
bc_ellipse=TRUE

# 2.5 tax_stackplot 样品和组分类学各级别的堆叠柱状图 Stackplot showing taxonomy in each level
ts_input=${wd}/result/tax/sum_
ts_level='"p","c","o","f","g"'
ts_design=${design}
ts_group_name=${g1}
ts_group_list=${g1_list}
ts_output=${ts_input}
ts_width=${width}
ts_height=${height}
# 显列图例的数量，推荐6，8，10，默认10
ts_number=${legend_number}
# 设置图例的顺序，默认FALSE按分类单元字母顺序排列，TRUE则按丰度由到大小排列
ts_order=FALSE

# 2.6 DA_compare 组间差异比较
Dc_input=${wd}/result/otutab.txt
# 差异比较方法edgeR or wilcox，默认edgeR
Dc_method='edgeR'
Dc_design=${design}
Dc_group_name=${g1}
Dc_group_list=${g1_list}
Dc_output=${wd}/result/compare/
# 设置差异OTUs顺序，默认按P值排序，TRUE为按丰度排序
Dc_order=FALSE













## 文库列表 library list
## ll clean_data/*.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|sort -u|tr "\n" " "
## tail -n+2 doc/library.txt |cut -f 1|tr "\n" " "
#list=L1 L2 L3 L4 L5 L6 L7 L8
#
## 文库建库方法
## 类型，单左侧barcode选择barcode_single_end，单右端和双端均选择barcode_paired_stitched。barcode_paired_stitched for barcode in end and both
#lib_type=barcode_single_end
## 正向barcode长度 forword barcode length
#bc1=10
## 反向barcode长度 reverse barcode length
#bc2=0
## Barcode类型，值为前两个barcode长度加和 barcode type, usually length equal barcode 1 plus barcode 2
#bt=10
## 质量值类型，分33或64两种；determine_phred-score.pl clean_data/L1_1.fq.gz
#phred=64
## OTUs聚类方式，主要包括cluster_otus，unoise3
#cluster=unoise3
#
## 2. 分组信息 Group
#
## 主要分组列名 Primary design group column
#g1=groupID
## 主要分组筛选，不筛选可为空
## 获取组信息 tail -n+2 doc/design.txt|cut -f 5|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
## 主要分组总表：
## "b1Col","b1ThasKO","b2ACT2KO","b2ACT3GK","b2ACT3KO","b2ALDHKO","b2BS","b2Col","b2ThasKO","b2ThasOE","b3ACT1KD","b3ACT1OE1","b3ACT1OE2","b3ACT2CR","b3ACT2KO","b3ACT2OE1","b3ACT2OE3","b3ACT2ThahDK","b3ACT3GK","b3ACT3KO","b3ACT3OE1","b3ACT3OE6","b3ACT3OE8","b3ALDHKO","b3BS","b3Col","b3ColAC","b3ThadKO","b3ThadOE1","b3ThadOE2","b3ThahKO","b3ThahOE1","b3ThahOE2","b3ThasJTOE3","b3ThasKO1","b3ThasKO2","b3ThasOE1","b3ThasOE3","b3ThasOE6","b3ThasThahDO","b3TL1KO","b3TL1OE"
## 筛除OE "b1Col","b1ThasKO","b2ACT2KO","b2ACT3KO","b2ALDHKO","b2BS","b2Col","b2ThasKO","b3ACT1KD","b3ACT2CR","b3ACT2KO","b3ACT2ThahDK","b3ACT3GK","b3ACT3KO","b3ALDHKO","b3BS","b3Col","b3ThadKO","b3ThahKO","b3ThasKO1","b3ThasKO2","b3TL1KO"
## batch3 root:  "b3ACT1KD","b3ACT2CR","b3ACT2KO","b3ACT2ThahDK","b3ACT3GK","b3ACT3KO","b3ALDHKO","b3Col","b3ThadKO","b3ThahKO","b3ThasKO1","b3ThasKO2","b3TL1KO"
## subgroup:
## "b3Col","b3ThahKO","b3ThasKO1","b3ThasKO2"
## "b3Col","b3ACT1KD","b3ACT2CR","b3ACT2KO"
## "b3Col","b3ACT2ThahDK","b3ACT3GK","b3ACT3KO"
## "b3Col","b3ALDHKO","b3ThadKO","b3TL1KO"
## 筛选4组6genes："b3Col","b3ThasKO1","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO","b3ACT2CR"
## 筛选4组4genes："b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2CR"
## 筛选4组4genes修改ACT2CR为KO："b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO"
#g1_list='"b3Col","b3ThasKO2","b3ThahKO","b3ThadKO","b3ACT2KO"'
## 次要分组列名 Secondary design group column，没有先真batch
#g2=genotype
## 次要分组筛选，不筛选可为空
## 次要分组总表：tail -n+2 doc/design.txt|cut -f 6 |sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
## "ACT1","ACT2","ACT2Thah","ACT3","ALDH","BS","Col","ColAC","Thad","Thah","Thas","ThasThah","TL1","WT"
#g2_list='"ACT1","ACT2","ACT2Thah","ACT3","ALDH","BS","Col","ColAC","Thad","Thah","Thas","ThasThah","TL1","WT"'
#
## 第三分组，可按此分组分别画Constrained PCoA，本示例是在不同土壤类型下画品种间差异
## 第三分组总表：tail -n+2 doc/design.txt|cut -f 7 |sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
#g3=type
## "BS","CR","DK","DO","GK","JTOE","KD","KO","OE","WT"
#g3_list=''
#
## 合并主要和次要分组 default FALSE, if merge_group=TRUE, must group_order=FALSE
#merge_group=FALSE
## 图例是否按主要组筛选样式排序，default TRUE ，但必须与merge_group保持相反
#group_order=TRUE
## 成对比较，TRUE为默认比对group_compare.txt，而FALSE则自动两两比对
#pair_compare=TRUE
## 仅批次按形状显示，默认按分组形状 Only Shape batch, default FALSE
#batch=FALSE
#
#
#
## 3. 报告参数
## 报告输入信息目录 
#sub="b3_3"
#doc=doc/${sub}
## 报告输出目录
#version=AC_b3_all_${sub}_lrt
## 报告输出是否精简版 report elite report, if FALSE report all figure and table, TRUE report frequently used figure
#elite_report=FALSE
## 图片长宽和字体大小，7组以下用默认，7组以上改为8x5或更大； figure size, recommend 4x2.5, 5x3(default), 8x5, 16x10, text_size 6, 7(default), 8
#width=8
#height=5
#text_size=7
## 图中显示taxonomy的数量，5，8(default)，10
#tax_number=10
## 按丰度和分类单元过滤OTU OTU taxonomy and abundance filter parameter
## 丰度按万分之五过滤 # threshold of filter low abundance OTU，OTU太多计算过慢可改为千一
#thre=0.0005
## 物种目前只去除叶绿体(k__Bacteria     p__Cyanobacteria        c__Chloroplast )和线粒体(k__Bacteria     p__Proteobacteria       c__Alphaproteobacteria  o__Rickettsiales        f__mitochondria)。不过滤p_xxx # filter some phylum	p__Actinobacteria,p__Bacteroidetes,p__Firmicutes,p__Proteobacteria p__Cyanobacteria,p__Chloroflexi c__Chloroplast,f__mitochondria
#taxonomy=c__Chloroplast,f__mitochondria
## 显著性P值过滤 # threshold of filter differentially abundance OTU，结果多可选0.01, 0.001
#pvalue=0.05
## 统计检验方式fdr, in edgeR have fdr or nonw
#adjust_method="fdr"
#fdr=0.2
## logFC常用1.5, 2, 4倍，对应0.585, 1, 2；此外还有1.3和1.7倍对应0.379和0.766
#logFC=0.379
#ellipse=TRUE
## 差异比较方法，默认是edgeR的 lrt ，可选 wilcoxon 秩和检验
#compare_method="lrt"
#
#
## 4. 不常用参数
### 输入输出目录文件 Input and output directory and files
### 可变配置文件目录，包括6个文本文件，主要个性group_*.txt来设置比较组、维恩图和三元图；可在doc在建子目录，复制并编辑，修改此处目录
#seq=clean_data
#summary=${wd}/doc/summary.txt
#library=${wd}/doc/library.txt
#design=${wd}/doc/design.txt
#compare=${wd}/${doc}/group_compare.txt
#venn=${wd}/${doc}/group_venn.txt
#tern=${wd}/${doc}/group_tern.txt
#temp=temp
#result=result
### 过滤OTU表结果目录 result based on filter OTU table
#result_f=result_k1-c
#
### 日志文件，记录数据量整体过滤和OTU过滤 log file for basic statistics
#log_reads=result/log_reads.txt
#log_otus=result/log_otus.txt
#log_usearch=result/log_usearch.txt
#
### 过滤序列质量>19为99%准确度 base quality, accurate > 99%; 29 means 99.9%
#quality=19
### 过滤N的数字，默认0，在引物中有N时最多3
#N=3
### 16S primers F799 and R1192 
## 5` primer used for 16S
#primer5=AACMGGATTAGATACCCKG
## 3` primer used for 16S, must reverse compliment
#primer3=GGAAGGTGGGGATGACGT 
## 引物匹配错误率，建议0.15，可调到最高0.25
#er=0.25
#
### 保留扩增子的最小长度，细菌799-1192用300，真菌ITS1用220 # min length, recommend 300 for bacterial 16S and 220 for ITS
#min_len=300
### 最小样本量 # sample min count, filter samples less than thre_count
#thre_count=5000
### 用于聚类的序列最低8，可选样品重复数量(replications numbers)、3xrep，甚至 1/1M # min count of unique reads, reconmend 1/1000000?
#minuniquesize=15
### 比对序列的相似度阈值，默认0.97 # similarity of cluster OTU
#sim=0.97
### 最大使用计算机线程数，主要给clustero多序列比对使用 # threads number used: 32
#p=32
### 用于筛选绘制圈图ggtree和igraphlan的OTU # filter OTU percentage > 0.5% for draw taxonomy and phylogenetic tree, 0.1% about 150 OTU is too much to show
#tax_per=0.005
### OTU物种注释的方法 # rdp, blast, rtax, mothur, uclust, sortmerna , default=uclust, recommend rdp is better
#method=rdp
### Alpha多样性分析的抽样数据 # alpha rarefaction count, recoomend 10000, at least 5000
#rarefaction=10000
### OTU表达丰度样式，默认为百分比percentage，可选css, rpm # add css or percentage mean add normlized value, default no sample data
#otu_stat_style=percentage
#
## 数据库 database；目前RDP虽然全，但是出现分类明显不准确且缺少Firmicute问题？暂时用gg13.8
### silva 128 99%, 492M uchime2_ref 中建议不去，如果去用最大的数据库，替换原29M rdp为新492M silva 128 99%
#rdp=/mnt/bai/public/ref/silva/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/99/99_otus_16S.fasta
### 绿色基因细菌16S数据库多序列比对文件，用于建立多序列比对和进化树 97% 763M, 99% 1.5G, SILVA 128 99% 13G, RDP 11.5仅细菌比对文件有78G，过G计算会很困难
#gg_align=/mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta
#### RDP 11.5 16S细菌和古菌序列: 注释比例高，但出现分类时间过长、进化树聚类不一致、缺少Firmicute等问题？
##gg_seq=/mnt/bai/public/ref/rdp/Bacteria_Archaea_seq.97
#### RDP 11.5 16S细菌和古菌物种注释信息
##gg_tax=/mnt/bai/public/ref/rdp/Bacteria_Archaea_tax.txt
#
#### RDP数据库用于去除嵌合体 rdp gold database, for remove chimera
##rdp=/mnt/bai/public/ref/rdp_gold.fa
#### 绿色基因细菌16S数据库多序列比对文件，用于建立多序列比对和进化树  greengene bacterial 16S database
##gg_align=/mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta 
### 绿色基因细菌16S数据库 greengene bacterial 16S database，虽然旧、不完整，但快、准。
#gg_seq=/mnt/bai/public/ref/gg_13_8_otus/rep_set/97_otus.fasta
### 绿色基因细菌16S数据库物种注释信息 greengene bacterial 16S database
#gg_tax=/mnt/bai/public/ref/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
#
#
## culture_graphlan
## 筛选指定组样品并鉴定培养比例,且样品也要对应
#type=""
#filter=filter_${type}_k1
#thre2=0.0005
#otu_table=${wd}/${result_f}/otu_table.txt
#cluture_db=/mnt/bai/yongxin/culture/ath/result/${type}culture_select.fa
#
#
## ## 不同列计算PCoA, compartment(第一轴分开), genotype(1/2/3/4分不开), site(第三轴分开), day
#group_color=day
#time_course=TRUE



# 载入16S扩增子主流程
# 放在makefile最好，保证所有参数己加载 
# Loading amplicon 16S main pipeline. Don't change the following file unless necessary.
include /mnt/bai/yongxin/github/Amplicon/16Sv2/main_pipeline.sh


