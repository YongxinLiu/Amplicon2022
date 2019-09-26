SHELL:=/bin/bash

	# 16S扩增子第二版流程配置文件，请根据项目具体情况修改 v2.0 2018/4/3
	# Config file of 16S Amplicon pipeline version 2, please modify according to experiment design. v2.0 2018/4/3


# 1. 标准流程参数 Parameters of standard pipeline

	## 工作目录 Working directory
	# 修改wd为当前工作目录pwd
	wd=`pwd`
	# make init # 建立分析所需子目录

	# 设置任务最大运行任务/线程数，超过CPU数量效率反而会降低
	p=32
	
	# 数据库
	# Greengene 13 May database, fa for usearch format, udb for usearch index
	usearch_gg=/mnt/bai/public/ref/gg_13_5_otus/97_otus_usearch.udb
	## Silva 132 database, fa for usearch format, udb for usearch index
	usearch_silva=/mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.udb
	usearch_rdp=/mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb

## 1.1 实验设计检查 Validate mapping file

## 1.2 文库双端合并 Merge clean reads

## 1.3 提取Barcode

	# 文库建库方法
	# 文库建库方式，单左侧barcode选择barcode_single_end，单右端和双端均选择barcode_paired_stitched, barcode_paired_stitched for barcode in end
	lib_type=barcode_paired_stitched
	# 正向barcode长度 forword barcode length
	bc1=10
	# 反向barcode长度 forword barcode length
	bc2=6

## 1.4  拆分文库为样品 Split library into sample
	# Barcode类型，值为前两者相加和 barcode type, usually length equal barcode 1 add barcode 2
	
	bt=16
	## 过滤序列质量>19为99%准确度 base quality, accurate > 99%; 29 means 99.9%
	quality=19
	# 统计 split_libraries_stat
	# Split library into sample
	# 默认只拆分单左端barcode类型的样品:先匹配左端，再提取序列ID，再提取右端，最后改名，注意实验设计要严格规范无空格
	# 每个文库中数量统计为result/split/L?.txt，结果有文本、PDF和PNG见result/split目录
	sample_split=result/sample_split.log

## 1.5. fq_trim 切除引物和标签

	# primer V5 19bp in left, and primer V7 18bp in right
	stripleft=19
	stripright=18

## 1.6. fa_unqiue 序列去冗余

	# Remove redundancy
	# 最小序列频率默认为8，去除低丰度，增加计算速度，整lane的序列推荐1/1M，即上一步最后一行的数据量
	minuniquesize=8

## 1.7. **otu_pick 挑选OTU**

	# Pick OTUs
	# 可选97% cluster_otus 和 unoise3 ，默认unoise3
	otu_method=unoise3
	# OTU日志文件，记录从挑选至过滤为最终的过程
	otu_log=result/otu.log

## 1.8. chimera_ref 参考去嵌合

	# Remove chimeras
	# 此处推荐使用大数据，如SILVA132，其它数据库如rdp_gold是错误的
	# /mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.fasta # 99%非冗余1.1G, 内存使用5.8G, 8min, 16.2% chimeras
	# /mnt/bai/public/ref/silva/SILVA_132_SSURef_tax_silva.fasta # 全部3.3G, 内存使用15.9G, 30min, 16.2% chimeras
	# 使用非冗余99%的省内存和时间，但结果差不多
	chimera_ref=${usearch_silva}
	# 模式，嵌合体比例由大到小: high_confidence specific balanced sensitive sensitive, 可选none关闭
	chimera_mode=none

## 1.9. host_rm 去宿主

	# Remove host original sequences
	# 去宿主方法选择 blast / sintax_gg / sintax_silva / sintax_silva_its / sintax_unite / none，推荐：sintax_silva
	host_method=sintax_silva
	# 方法1. blast宿主基因组(含叶绿体/线粒体)去除同源序列，如水稻微生物，需要提供水稻基因组；可调相似度和覆盖度的阈值(百分数)
	host=/mnt/bai/public/ref/rice/msu7/all.con
	host_similarity=90
	host_coverage=90
	# 方法2. 基于gg注释结果筛选, 去除叶绿体Chloroplast、线粒体mitochondria，默认为usearch_gg数据库
	# 方法3. silva注释可识线粒体Mitochondria、叶绿体Chloroplast和真核生物Eukaryota(包括宿主、真菌、原生动物等)，默认为usearch_silva数据库


## 1.10. otutab_create 生成OTU表

	# Creat OTUs table
	# 有 usearch10 和 vsearch 两个软件可选，默认 usearch10 ，vsearch多线程会更快些
	map_method=usearch10
	map_identify=0.97

## 1.11. otutab_filter OTU表筛选

	# Filter OTU table
	# OTU表筛选日志文件
	log_otutable=result/otutab.log
	# 按样本量筛选，默认30，根据otu_stats结果调整
	min_sample_size=30
	# 按矩阵中每个点count, freq筛选，低于阈值变为0
	# 按OTU丰度和频率筛选，如OTU测序量至少8次，相对丰度百万分之一(建议挑选序列去冗余部分调高阈值更合理)
	min_otu_size=8
	# 按频率筛选，推荐十万分之一0.00001，范围千一至百分一0.001 - 0.000001之间
	min_otu_freq=0.000001
	# 抽样标准化的值，推荐最小1000，根据统计结果选择筛选后最小值或可保留大部分样品的值
	sample_size=1000

## 1.12. tax_assign 物种注释

	# Assign taxonomy
	# 物种注释推荐使用小而准的数据库，如rdp trainset 16(由Robert整理)
	# 可选gg, silva, rdp分别从官网下载并shell调整格式，gg较准但旧，silva全但不准，rdp少而准，比较通用
	sintax_db=${usearch_rdp}
	# 分类准确度阈值，默认0.8，注释太少最小可改0.5，发现有明显错误可最高上升为0.95，改为零为最大化显示物种注释
	sintax_cutoff=0.6

## 1.13. tax_sum 物种注释统计

	# Taxonomy summary
	# 按门、纲、目、科、属水平分类汇总，结果位于：result/tax/sum_*.txt
	# 输出可读的taxonomy与OTU对应文件，有2列和8列版，result/taxonomy_2/8.txt

## 1.14. tree_make 多序列比对和进化树

	# Multiply alignment and make_phylogeny



include /mnt/bai/yongxin/github/Amplicon/16Sculture2/pipeline.md
