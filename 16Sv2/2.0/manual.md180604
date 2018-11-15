<!-- TOC -->

- [1. 准备工作 Preparation](#1-准备工作-preparation)
    - [1.1. 准备流程配置文件](#11-准备流程配置文件)
    - [1.2. 初始化工作区](#12-初始化工作区)
    - [1.3. 准备原始数据](#13-准备原始数据)
- [2. 处理序列 Processing sequencing data](#2-处理序列-processing-sequencing-data)
    - [2.1. 按实验设计拆分lane为文库](#21-按实验设计拆分lane为文库)
    - [2.2. 按实验设计拆分文库为样品](#22-按实验设计拆分文库为样品)
    - [2.3. 样品双端合并、重命名、合并为单一文件](#23-样品双端合并重命名合并为单一文件)
    - [2.4. 切除引物与标签](#24-切除引物与标签)
    - [2.5. 质量控制](#25-质量控制)
    - [2.6. 序列去冗余](#26-序列去冗余)
    - [2.7. 挑选OTU](#27-挑选otu)
    - [2.8. 有参去嵌合体](#28-有参去嵌合体)
    - [2.9. 去除宿主](#29-去除宿主)
    - [2.10. 生成OTU表](#210-生成otu表)
    - [2.11. 过滤样本和OTUs](#211-过滤样本和otus)
    - [2.12. 物种注释](#212-物种注释)
    - [2.13. 物种统计](#213-物种统计)
    - [2.14. 多序列比对和进化树](#214-多序列比对和进化树)
    - [2.15. Alpha多样性指数计算](#215-alpha多样性指数计算)
    - [2.16. Beta多样性距离矩阵计算](#216-beta多样性距离矩阵计算)
    - [2.17. 有参考构建OTU表](#217-有参考构建otu表)
- [3. 统计绘图 Statistics and plot](#3-统计绘图-statistics-and-plot)
    - [3.1. Alpha多样性指数箱线图](#31-alpha多样性指数箱线图)
    - [3.2. Alpha丰富度稀释曲线](#32-alpha丰富度稀释曲线)
    - [3.3. 主坐标轴分析距离矩阵](#33-主坐标轴分析距离矩阵)
    - [3.4. 限制性主坐标轴分析](#34-限制性主坐标轴分析)
    - [3.5. 样品和组各级分类学堆叠柱状图](#35-样品和组各级分类学堆叠柱状图)
    - [3.6. 组间差异比较](#36-组间差异比较)

<!-- /TOC -->

# 1. 准备工作 Preparation

## 1.1. 准备流程配置文件

    # Prepare config file of pipeline
	# 设置工作目录并准备流程文件
	wd=maize/180604breed
	cd # 切换至家目录
	mkdir -p $wd # 创建工作目录
	mkdir -p github/Work/$wd # 创建脚本备份目录
	cp /mnt/bai/yongxin/github/Amplicon/16Sv2/parameter.md github/Work/$wd # 复制参数文件至同步备份目录
	ln -s `pwd`/github/Work/$wd/parameter.md $wd/makefile # 建立工作目录软链
	cp /mnt/bai/yongxin/github/Amplicon/16Sv2/manual.md github/Work/$wd # 复制帮助和实验记录文件至同步备份目录
	ln -s `pwd`/github/Work/$wd/manual.md $wd/manual.sh # 建立工作目录软链
	cd $wd # 进入工作目录

## 1.2. 准备实验设计

	# Initialize the working directory
	make init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt
	# 按library中第二列index准备测序文库
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/180528.lane11/Clean/CWHPEPI00001823/seq/"$2"_1.fq seq/"$1"_1.fq");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/180528.lane11/Clean/CWHPEPI00001823/seq/"$2"_2.fq seq/"$1"_2.fq");}' <(tail -n+2 doc/library.txt )

	# 标准多文库实验设计拆分，保存模板中design页为doc/design_raw.txt
	split_design.pl -i doc/design_raw.txt
	# 从其它处复制实验设计
	cp ~/ath/jt.HuangAC/batch3/doc/L*.txt doc/
	# 删除多余空格，windows换行符等
	sed -i 's/ //g;s/\r/\n/' doc/*.txt 
	head -n3 doc/L01.txt
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L01.txt | sed 's/#//g') <(cat doc/L* |grep -v '#') > doc/design.txt
	# 检查是否相等
	wc -l doc/design.txt
	cut -f 1 doc/design.txt|sort|uniq|wc -l

## 1.3. 准备原始数据

	# 拆lane和质量转换归为原始seq目录中处理
	# Prepare raw data
	#ln ~/seq/180210.lane9.ath3T/Clean/CWHPEPI00001683/lane_* ./
	#cp ~/ath/jt.HuangAC/batch3/doc/library.txt doc/
	
	# 检查数据质量，转换为33
	#determine_phred-score.pl seq/lane_1.fq.gz
	# 如果为64，改原始数据为33
	rename 's/lane/lane_33/' seq/lane_*
	# 关闭质量控制，主要目的是格式转换64至33，不然usearch无法合并
	#time fastp -i seq/lane_64_1.fq.gz -I seq/lane_64_2.fq.gz \
	#	-o seq/lane_1.fq.gz -O seq/lane_2.fq.gz -6 -A -G -Q -L -w 9
	# 1lane 80GB, 2 threads, 102min


# 2. 处理序列 Processing sequencing data

## 2.1. 按实验设计拆分lane为文库

	# Split lane into libraries
	# lane文件一般为seq/lane_1/2.fq.gz
	# lane文库信息doc/library.txt：至少包括编号、Index和样品数量三列和标题
	# head -n3 doc/library.txt
	#LibraryID	IndexRC	Samples
	#L1	CTCAGA	60
	
	# 按library.txt拆分lane为library
	# make lane_split


## 2.2. 按实验设计拆分文库为样品


	# 拆分样品
	head -n3 doc/L01.txt
	# 按L1/2/3...txt拆分library为samples
	make library_split
	make library_split_stat

## 2.3. 样品双端合并、重命名、合并为单一文件

	# Merge paired reads, renames and merge all samples
	# 样品双端合并、重命名、合并为单一文件, 注意fastq为33格式，64位采用fastp转换
	make sample_merge
	make sample_merge_stat


## 2.4. 切除引物与标签

	# Cut primers and lables
	# 切除左端标签和引物，右端 引物
	# Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
	make fq_trim


## 2.5. 质量控制

	# Quality control
	# 过滤序列中预期累计错误率>1%的序列
	make fq_qc


## 2.6. 序列去冗余

	# Remove redundancy, get unique reads
	make fa_unqiue


## 2.7. 挑选OTU

	# Pick OTUs
	# unoise3速度比cluster_otus慢上百倍
	make otu_pick


## 2.8. 有参去嵌合体

	# Remove chimiras by silva database
	# 基于SILVA数据库去除
	make chimera_ref


## 2.9. 去除宿主

	# Remove host
	# 根据SILVA注释去除线粒体、叶绿体、真核生物18S和未知序列(非rRNA)
	make host_rm


## 2.10. 生成OTU表
	
	# Create OTUs table
	# 默认使用vsearch更快10倍，可选usearch10，线程不可超48
	make otutab_create


## 2.11. 过滤样本和OTUs

	# OTU table filter samples and OTU
	# 推荐过滤低测序量<5000的样本，筛选大于1RPM的OTU
	make otutab_filter 


## 2.12. 物种注释

	# Assign taxonomy
	# 默认使用RDP trainset快而准，GG太旧，Silva太慢
	# 推荐阈值为0.6保证注释更完整
	make tax_assign


## 2.13. 物种统计
	
	# Taxonomy summary
	# 必须所有物种有注释，否则有可能报错
	make tax_sum


## 2.14. 多序列比对和进化树
	
	# Multiply alignment and make_phylogeny
	# usearch10/culsterO结果不同可能影响多样性分析(usearch unifrac结果更可信)
	# 进化树，用于树图和多样性分析
	make tree_make

## 2.15. Alpha多样性指数计算
	
	# Calculate alpha diversity index
	# alpha指数计算结果为 result/alpha/index.txt
	# 稀释梯度结果位于 result/alpha/rare.txt
	make alpha_calc

## 2.16. Beta多样性距离矩阵计算
	
	# Beta diversity tree and distance matrix
	# 最好用usearch，结果unifrac分类更好；clustero+fastree结果PCoA较差
	make beta_calc
	# ---Fatal error--- ../calcdistmxu.cpp(32) assert failed: QueryUniqueWordCount > 0 致信作者; 改用qiime1

## 2.17. 有参考构建OTU表

	# Reference based OTU table
	# otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
	make otutab_gg



# 3. 统计绘图 Statistics and plot

## 3.1. Alpha多样性指数箱线图
	
	# Alpha index in boxplot
	make alpha_boxplot

## 3.2. Alpha丰富度稀释曲线
	
	# Alpha rarefracation curve
	make alpha_rare

## 3.3. 主坐标轴分析距离矩阵
	
	# PCoA of distance matrix
	make beta_pcoa

## 3.4. 限制性主坐标轴分析

	# Constrained PCoA / CCA of bray distance matrix
	# OTU表基于bray距离和CCA，至少3个组 
	make beta_cpcoa

## 3.5. 样品和组各级分类学堆叠柱状图

	# Stackplot showing taxonomy in each level
	make tax_stackplot

## 3.6. 组间差异比较 
	
	# Group compareing by edgeR or wilcox
	# 可选负二项分布，或wilcoxon秩和检验
	make DA_compare
