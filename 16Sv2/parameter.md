SHELL:=/bin/bash

	# 16S扩增子第二版流程配置文件，请根据项目具体情况修改 v2.0 2018/4/3
	# Config file of 16S Amplicon pipeline version 2, please modify according to experiment design. v2.0 2018/4/3


# 1. 标准流程参数 Parameters of standard pipeline

	## 工作目录 Working directory
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

## 1.1. lane_split 拆分下机数据为文库

	# Split lane into library
	# 需要将lane文件放入seq目录，对应的index和文库放在doc/library.txt
	# 注意：务必查看文库文件中Index具体格式，默认为#Index，其它情况需修改主流程源代码main_pipeine.sh
	# 文库名
	lane=lane
	lib_log=result/library.log

## 1.2. library_split_stat 拆分文库为样品

	# Split library into sample
	# 默认只拆分单左端barcode类型的样品:先匹配左端，再提取序列ID，再提取右端，最后改名，注意实验设计要严格规范无空格
	# 每个文库中数量统计为result/split/L?.txt
	sample_split=result/sample_split.log

## 1.3. sample_merge_stat 双端序列合并

	# Merge pair-end reads
	# 如果是pair-end reads是phred64，需要先使用fastp转换为33且关闭质控(质控影响序列长度)，再使用usearch10 mergepair
	sample_merge=result/sample_merge.log

## 1.4. fq_trim 切除引物和标签

	# Cut barcode 10bp + primer V5 19bp in left, and primer V7 18bp in right
	stripleft=29
	stripright=18

## 1.5. fq_qc 质量控制
	
	# fastq filter
	# 默认错误率<0.01 keep reads error rates less than 1%
	fastq_maxee_rate=0.01

## 1.6. fa_unqiue 序列去冗余

	# Remove redundancy
	# 最小序列频率miniuniqusize默认为8，去除低丰度，增加计算速度，整lane的序列可更改为30，甚至100
	minuniquesize=100

## 1.7. otu_pick 挑选OTU

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
	# 模式，嵌合体比例由大到小: high_confidence specific balanced sensitive sensitive
	chimera_mode=balanced

## 1.9. host_rm 去宿主

	# Remove host original sequences
	# 去宿主方法选择 blast / sintax_gg / sintax_silva，推荐：sintax_silva
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
	# 按样本量筛选，默认5000，根据otu_stats结果调整
	min_sample_size=5000
	# 按矩阵中每个点count, freq筛选，低于阈值变为0
	# 按OTU丰度和频率筛选，如OTU测序量至少8次，相对丰度百万分之一(建议挑选序列去冗余部分调高阈值更合理)
	min_otu_size=8
	# 按频率筛选，推荐十万分之一0.00001，范围千一至百分一0.001 - 0.000001之间
	min_otu_freq=0.000001
	# 抽样标准化的值，推荐最小10000，根据统计结果选择筛选后最小值或可保留大部分样品的值
	sample_size=30000

## 1.12. tax_assign 物种注释

	# Assign taxonomy
	# 物种注释推荐使用小而准的数据库，如rdp trainset 16(由Robert整理)
	# 可选gg, silva，分别从官网下载并shell调整格式
	sintax_db=${usearch_silva}
	# 分类准确度阈值，默认0.8，注释太少最小可改0.5，发现有明显错误可最高上升为0.95
	sintax_cutoff=0.6

## 1.13. tax_sum 物种注释统计

	# Taxonomy summary
	# 按门、纲、目、科、属水平分类汇总，结果位于：result/tax/sum_*.txt
	# 输出可读的taxonomy与OTU对应文件，有2列和8列版，result/taxonomy_2/8.txt

## 1.14. tree_make 多序列比对和进化树

	# Multiply alignment and make_phylogeny

## 1.15. alpha_calc Alpha多样性指数

	# Calculate alpha diversity index
	# alpha指数计算结果为 result/alpha/index.txt
	# 稀释梯度抽样方法 richness (observed OTUs)-method fast / with_replacement / without_replacement , 结果位于 result/alpha/rare.txt
	rare_method=without_replacement

## 1.16. beta_calc Beta多样性距离矩阵

	# Beta diversity tree and distance matrix
	# 距离矩阵计算方法，34种可选： abund_jaccard, binary_chisq, binary_chord, binary_euclidean, binary_hamming, binary_jaccard, binary_lennon, binary_ochiai, binary_otu_gain, binary_pearson, binary_sorensen_dice, bray_curtis, bray_curtis_faith, bray_curtis_magurran, canberra, chisq, chord, euclidean, gower, hellinger, kulczynski, manhattan, morisita_horn, pearson, soergel, spearman_approx, specprof, unifrac, unifrac_g, unifrac_g_full_tree, unweighted_unifrac, unweighted_unifrac_full_tree, weighted_normalized_unifrac, weighted_unifrac
	# 默认使用4种
	dis_method=bray_curtis,binary_jaccard,weighted_unifrac,unweighted_unifrac

## 1.17. otutab_ref 有参比对生成OTU表

	# 如Greengenes，可用于picurst, bugbase分析
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
	# 差异统计按丰度过滤 abundance filter，如丰度按万分之一过滤，减少计算量，提高OTU的FDR值，可选十万/百万之一
	abundance_thre=0.0001
	# 差异比较方法，默认是 edgeR ，可选 wilcox 秩和检验
	compare_method="edgeR"
	# 显著性P值过滤 threshold of P-value，可选0.05, 0.01, 0.001。采用FDR校正，此参数意义不大，即使0.001也没有FDR < 0.2过滤严格
	pvalue=0.01
	# 统计检验方式fdr
	FDR=0.05
	# 差异变化倍数常用1.5, 2, 4倍，对应logFC为0.585, 1, 2；菌丰度变化倍数不明显，还可用1.3和1.7倍对应0.379和0.766
	FC=1.3

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
	Dc_compare=${compare}
	Dc_method=${compare_method}
	Dc_pvalue=${pvalue}
	Dc_FDR=${FDR}
	Dc_FC=${FC}
	Dc_thre=${abundance_thre}
	Dc_design=${design}
	Dc_group_name=${g1}
	Dc_group_list=${g1_list}
	Dc_output=${wd}/result/compare/

	# 2.7 plot_volcano 基于差异OTU表绘制火山图
	pv_input=${wd}/result/tax/sum_
	pv_design=${design}
	pv_output=${pv_input}
	pv_width=5
	pv_height=7
	# 显列图例的数量，推荐6，8，10，默认10
	pv_number=${legend_number}
	# 设置图例的顺序，默认FALSE按分类单元字母顺序排列，TRUE则按丰度由到大小排列
	pv_order=FALSE



include /mnt/bai/yongxin/github/Amplicon/16Sv2/main_pipeline.sh
