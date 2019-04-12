# Project config file

help:
	# 16S Amplicon pipeline config file for bacterial identify
	# Version 1.0 2017/9/26
	# 流程使用方法
	# 1. 初始化：复制makefile和makefile.cfg至工作目录，设置初始目录(wd=pwd)，建立相关目录：make init
	# 2. 原始数据：上传原始数据*.fq.gz至clean_data/；原始数据文件名重命名为L1/2/3_1/2.fq.gz(可选)；
	#	批量删除文件名中共有部分可使用rename，例如：rename 's/-AN-RS-//g;s/_HW5WYBCXY_L1//g;s/clean\.//g' *.gz 
	# 3. 实验设计：将16s_design_template.xlsx中每个表认真填好，将前6个表格保存为同名.txt文件于doc/目录；
	# write_mappingfile_culture.pl -o doc/R2A3.txt -s Medicago -p 42 # 按板数量和物种生成mapping file
	# cp doc/R2A3.txt doc/L1.txt
	# cat <(head -n1 doc/L1.txt|sed 's/#//g') doc/L* |grep -v '#' > doc/design.txt
	#	批量改文库名：make rename
	#	合并mapping file为实验设计design：cat <(head -n1 doc/L1.txt|sed 's/#//g') doc/L* |grep -v '#' > doc/design.txt
	#	获取文库列表给list变量：ll clean_data/*.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|sort -u|tr "\n" " "
	#	添写g1/g2实验的主要和次要分组列表，如groupID, genotype, compartment, soiltype，没有g2可使用默认batch
	#	获取实验组信息总表并手动筛选给g1_list变量：tail -n+2 doc/design.txt|cut -f 5 |sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
	#	仔细检查下面样品和参数信息，按下列make系列命令运行
	# English help
	# 1. Set work directory, Upload clean data to directory clean_data/, copy makefile* and make init;
	# 2. Upload mapping file according to each library and design to doc/
	# cat <(head -n1 doc/L1.txt|sed 's/#//g') doc/L* |grep -v '#' > doc/design.txt # merge mappingfile to design
	# 3. Manually modify each library name or using rename regexp batch rename
	# rename 's/yoyo\-//g;s/_HM223BCXY_L2//g;s/clean\.//g' clean_data/*.gz
	# ll clean_data/*.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|sort -u
	# lib=`ls clean_data/*.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|uniq|tr "\n" " "` # get library name list
	# 4. Manually set variable in the following part
	# make rename # creat directory
	# make stat # stat each library
	# make otu_table # make otu table
	# make diversity # alpha diversity
	# make rmd # write report 


# 标准流程参数 Standard pipeline parameter
## 工作目录 working directory
wd=/mnt/bai/yongxin/culture/wheat
## 文库列表 library list
list=L10 L5 L6 L7 L8 L9
## 文库建库方法
### 文库建库方式，单左侧barcode选择barcode_single_end，单右端和双端均选择barcode_paired_stitched # barcode_paired_stitched for barcode in end
lib_type=barcode_paired_stitched
### 正向barcode长度 forword barcode length, according to experiment design
bc1=10
### 反向barcode长度 forword barcode length, according to experiment design
bc2=6
### Barcode类型，值为前两者相加和 barcode type, usually length equal barcode 1 add barcode 2
bt=16
#let bt=bc1+bc2
## 主要分组列名 Primary design group column
g1=plate
## 次要分组列名 Secondary design group column
g2=batch
## 主要分组筛选，如不筛选可注释此行；
### 获取组信息，根据design中列调整cut后面列位置：tail -n+2 doc/design.txt|cut -f 5 |sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
### 获取的总表用于筛选：
## tail -n+2 doc/design.txt|cut -f 5 |sort|uniq|grep -P 'R$' | awk '{print "\""$1"\""}'|tr "\n" ","
g1_list='"RnfpR","soil"'
## 次要分组筛选，如不筛选可注释此行；方法同主要分组tail -n+2 doc/design.txt|cut -f 10 |sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
###"1","2","3"
g2_list='"3"'
## 合并主要和次要分组 default FALSE, if merge_group=TRUE, must group_order=FALSE
merge_group=FALSE
## 图例是否按主要组筛选样式排序，default TRUE ，但必须与merge_group保持相反
group_order=TRUE
## 成对比较，TRUE为默认比对group_compare.txt，而FALSE则自动两两比对
pair_compare=TRUE
## Shape batch, default FALSE
batch=FALSE

## 报告参数
### 报告输入目录
doc=doc
### 报告输出目录
version=test
### 报告输出版本，是否精简版 report elite report, default FALSE, report all figure and table, TRUE report frequently used figure
elite_report=TRUE
## 图片长宽和字体大小，7组以下用默认，7组以上改为8x5或更大； figure size, recommend 4x2.5, 5x3(default), 8x5, 16x10, text_size 6, 7(default), 8
width=40
height=5
text_size=1
### 图中显示taxonomy的数量，5，8(default)，10
tax_number=8
## 按丰度和分类单元过滤OTU OTU taxonomy and abundance filter parameter
### 丰度按万分之五过滤 # threshold of filter low abundance OTU
thre=0.0005
### 物种目前去除蓝细菌和绿细菌门 # filter some phylum	p__Actinobacteria,p__Bacteroidetes,p__Firmicutes,p__Proteobacteria
taxonomy=p__Cyanobacteria,p__Chloroflexi
### 显著性P值过滤 # threshold of filter differentially abundance OTU
pvalue=0.05



# 实验流程不常用参数
## 输入输出目录文件 Input and output directory and files
## 可变配置文件目录，包括6个文本文件，主要个性group_*.txt来设置比较组、维恩图和三元图；可在doc在建子目录，复制并编辑，修改此处目录
seq=clean_data
summary=${wd}/${doc}/summary.txt
library=${wd}/doc/library.txt
design=${wd}/doc/design.txt
compare=${wd}/${doc}/group_compare.txt
venn=${wd}/${doc}/group_venn.txt
tern=${wd}/${doc}/group_tern.txt
temp=temp
result=result
## 过滤OTU表结果目录 result based on filter OTU table
result_f=result_k1-c

## 日志文件，记录数据量整体过滤和OTU过滤 log file for basic statistics
log=result/readme.log
## 过滤序列质量>19为99%准确度 base quality, accurate > 99%; 29 means 99.9%
quality=19
## 16S primers F799 and R1192 
primer5=AACMGGATTAGATACCCKG # 5` primer used for 16S
primer3=GGAAGGTGGGGATGACGT # 3` primer used for 16S, must reverse compliment
## 保留扩增子的最小长度，细菌799-1192用300，真菌ITS1用220 # min length, recommend 300 for bacterial 16S and 220 for ITS
min_len=300
## 最小样本量 # sample min count, filter samples less than thre_count
thre_count=20
## 用于聚类的序列最低丰度，最小2，单样品推荐4，混样推荐8
minuniquesize=8
## 聚类序列的相似度阈值，默认0.99 # similarity of database
sim=0.99
## 最大使用计算机线程数，主要给clustero多序列比对使用 # threads number used: 32
p=32
## 用于筛选绘制圈图ggtree和igraphlan的OTU # filter OTU percentage > 0.5% for draw taxonomy and phylogenetic tree, 0.1% about 150 OTU is too much to show
tax_per=0.005
## OTU物种注释的方法 # rdp, blast, rtax, mothur, uclust, sortmerna , default=uclust, recommend rdp is better
method=rdp
## Alpha多样性分析的抽样数据 # alpha rarefaction count, recoomend 10000, at least 5000
rarefaction=1000
## OTU表达丰度样式，默认为百分比percentage，可选css, rpm # add css or percentage mean add normlized value, default no sample data
otu_stat_style=percentage

# 数据库 database
## silva 128 99%, 492M uchime2_ref 中建议不去，如果去用最大的数据库，替换原29M rdp为新492M silva 128 99%
rdp=/mnt/bai/public/ref/silva/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/99/99_otus_16S.fasta
## 绿色基因细菌16S数据库多序列比对文件，用于建立多序列比对和进化树 97% 763M, 99% 1.5G, SILVA 128 99% 13G, RDP 11.5仅细菌比对文件有78G，过G计算会很困难
gg_align=/mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta
## RDP 11.5 16S细菌和古菌序列
gg_seq=/mnt/bai/public/ref/rdp/Bacteria_Archaea_seq.97
## RDP 11.5 16S细菌和古菌物种注释信息
gg_tax=/mnt/bai/public/ref/rdp/Bacteria_Archaea_tax.txt


# graphlan
## 自然样品结果位置
nature=/mnt/bai/yongxin/wheat/NP

include /mnt/bai/yongxin/github/Amplicon/16Sculture/pipeline.md
