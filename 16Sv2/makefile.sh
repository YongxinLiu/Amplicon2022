#SHELL:=/bin/bash
#
## 16S扩增子流程配置文件，请根据项目具体情况修改 v1.4 2017/12/23
## Config file of 16S Amplicon pipeline, please modify according to experiment design. v1.4 2017/12/23
#
#
#
## 1. 测序文库参数 Standard pipeline parameter
#
## 工作目录 working directory
## pwd # 修改wd为新工作目录
#wd=/mnt/bai/yongxin/ath/jt.HuangAC/batch3all
## make init # 建立流程所需目录
#
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
#
#
#
## 载入16S扩增子主流程，修改请保证向前兼容，最后只添加新分枝流程
## 放在makefile最好，保证所有参数己加载 
## Loading amplicon 16S main pipeline. Don't change the following file unless necessary.
#include /mnt/bai/yongxin/ref/amplicon/16s/makefile_16s
#
#
