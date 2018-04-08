# 0 准备工作

# 准备流程配置文件
cp ~/github/Amplicon/16Sv2/makefile.sh makefile
cp ~/github/Amplicon/16Sv2/manual.sh ./

# 建立初始目录
make init



# 1 原始序列预处理

# 1.1 准备lane测序文件及文库实验设计
ln ~/seq/180210.lane9.ath3T/Clean/CWHPEPI00001683/lane_* ./
cp ~/ath/jt.HuangAC/batch3/doc/library.txt doc/
head -n3 doc/library.txt
#LibraryID	IndexRC	Samples
#L1	CTCAGA	60
#L2	GACGAC	60
# 按library.txt拆分lane为library
make lane_split
## Split lane运行时间太长，改用头1000000万行在temp中测试
#cd temp
#mkdir -p seq doc temp result
#zcat ../seq/lane_1.fq.gz | head -n 1000000 | gzip > seq/lane_1.fq.gz
#zcat ../seq/lane_2.fq.gz | head -n 1000000 | gzip > seq/lane_2.fq.gz
#cp ../doc/library.txt doc/
#time parallel -j ${p} --no-notice "zcat seq/lane_1.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > seq/{1}_1.fq" ::: `tail -n+2 doc/library.txt | cut -f 2`
#time parallel -j ${p} --no-notice "zcat seq/lane_2.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$' > seq/{1}_2.fq" ::: `tail -n+2 doc/library.txt | cut -f 2`
## 修改索引为文库编号 rename index to library ID
##awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$2"_1.fq seq/"$1"_1.fq")}' <(tail -n+2 doc/library.txt) # 左端
##awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$2"_2.fq seq/"$1"_2.fq")}' <(tail -n+2 doc/library.txt) # 右端
##awk 'BEGIN{OFS=FS="\t"}{system("rename 's/"$2"/"$1"/' seq/"$2"_?.fq")}' <(tail -n+2 doc/library.txt) # 调用rename无效，还是用mv好
#awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$2"_1.fq seq/"$1"_1.fq");system("mv seq/"$2"_2.fq seq/"$1"_2.fq")}' <(tail -n+2 doc/library.txt) # 左右两次system

# 1.2 准备每个文库的实验设计
# 标准多文库实验设计拆分
split_design.pl -i doc/design_raw.txt
# 从其它处复制实验设计
cp ~/ath/jt.HuangAC/batch3/doc/L?.txt doc/
sed -i 's/ //g;s/\r/\n/' doc/*.txt # 删除多余空格，windows换行符等
head -n3 doc/L1.txt
##SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	groupID	genotype	type	lines	replicate	compartment	batch	specieDescription
#ThasOE1r1	ACGCTCGACA	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	ThasOE1	Thas	OE	1	1	root	3	Aarabidopsis	ThasOE1r1
#ThasOE1r2	ATCAGACACG	AACMGGATTAGATACCCKG	ACGTCATCCCCACCTTCC	ThasOE1	Thas	OE	1	2	root	3	Aarabidopsis	ThasOE1r2
# 按L1/2/3...txt拆分library为samples
make library_split
## 测试本功能的源代码：初始化变量测试拆样
#library=L1
#mkdir -p seq/sample
#barcode=ACGCTCGACA
#grep -A 2 -B 1 -P "^${barcode}" seq/${library}_1.fq|grep -v -P '^--$' > seq/sample/${barcode}_1.fq
## 并行匹配左端
#parallel -j 12 --no-notice "grep -A 2 -B 1 -P "^{1}" seq/${library}_1.fq|grep -v -P '^--$' > seq/sample/{1}_1.fq" ::: `tail -n+2 doc/${library}.txt | cut -f 2`
## 提取序列名
#awk 'NR%4==1' seq/sample/${barcode}_1.fq |sed 's/^@//;s/1$/2/' > temp/${barcode}.id
## 筛选右端序列
#usearch10 -fastx_getseqs seq/${library}_2.fq -labels temp/${barcode}.id -fastqout seq/sample/${barcode}_2.fq
## 并行提取ID，按ID提取右端，再改名
#parallel "awk 'NR%4==1' seq/sample/{1}_1.fq |sed 's/^@//;s/1$/2/' > temp/{1}.id" ::: `tail -n+2 doc/${library}.txt | cut -f 2`
#parallel "usearch10 -fastx_getseqs seq/${library}_2.fq -labels temp/{1}.id -fastqout seq/sample/{1}_2.fq" ::: `tail -n+2 doc/${library}.txt | cut -f 2`
#awk 'BEGIN{OFS=FS="\t"}{system("mv seq/sample/"$2"_1.fq seq/sample/"$1"_1.fq");system("mv seq/sample/"$2"_2.fq seq/sample/"$1"_2.fq")}' <(tail -n+2 doc/${library}.txt)

# 1.3 实验设计样品逐个双端合并、重命名、合并为单一文件
# 依据各文库L*.txt文件生成实验设计
cat <(head -n1 doc/L1.txt) <(cat doc/L* |grep -v '#') > doc/design.txt
# 运行本节
make sample_merge
## 本节流程测试过程
#usearch10 -fastq_mergepairs seq/sample/Colr1_1.fq -reverse seq/sample/Colr1_2.fq -fastqout temp/Colr1.merged.fq -relabel Colr1.
## CharToIntQual('e') Phred score 68 out of range 0..41 # 64质量不支持，需转换
## fastp质量转换：-6转换质量64至33，-A禁用切接头，-G禁用切除polyG，-Q禁用质控，-L禁用长度过滤
#fastp -i seq/sample/Colr1_1.fq -I seq/sample/Colr1_2.fq -o temp/Colr1_1.fq -O temp/Colr1_2.fq -6 -A -G -Q -L

# 1.4 切除引物与条型码 Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
make fq_trim
# 05:49 37Mb    100.0% Processing, 145 (0.0%) too short

# 1.5 质量控制fastq filter
make fq_qc
# 单端：15:43 1.7Gb   100.0% Filtering, 96.1% passed
# 双端合并：11:59 1.7Gb   100.0% Filtering, 99.9% passed 合并后质量上升明显，99.9%质量合格

# 1.6 序列去冗余 Remove redundancy
make fa_unqiue
# 09:29 83.6Gb Min size 1, median 1, max 2695385, avg 12.65
# 117685 uniques written, 6748488 clusters size < 30 discarded (79.0%)
#05:21 50.1Gb 54745090 seqs, 9256698 uniques, 7887911 singletons (85.2%)
#05:21 50.1Gb Min size 1, median 1, max 2178111, avg 5.91
#75652 uniques written, 7299707 clusters size < 30 discarded (78.9%)

# 1.7 挑选OTU Pick OTUs
make otu_pick
# unoise3: 02:01:31 520Mb   100.0% 31802 good, 14758 chimeras 速度比cluster_otus慢上百倍
# cluster_otus: 01:09 76Mb    100.0% 4285 OTUs, 5051 chimeras

# 1.8 去嵌合数据比较
make chimera_ref
# 02:22 6.9Gb   100.0% Chimeras 689/4285 (16.1%), in db 850 (19.8%), not matched 2746 (64.1%)
## 测试环境: rice timecourse
#cd ~/rice/timecourse
## silva132: Nr99,1.1M | 07:50 5.8Gb   100.0% Chimeras 3984/24553 (16.2%), in db 4652 (18.9%), not matched 15917 (64.8%)
#chimera_ref=/mnt/bai/public/ref/silva/SILVA_132_SSURef_Nr99_tax_silva.fasta
#usearch10 -uchime2_ref temp/otus.fa -db ${chimera_ref} -strand plus -mode balanced -chimeras temp/otus_chimeras.fa
## silva132: All, 3.3M | 30:36 15.9Gb  100.0% Chimeras 3974/24553 (16.2%), in db 4991 (20.3%), not matched 15588 (63.5%)
#chimera_ref=/mnt/bai/public/ref/silva/SILVA_132_SSURef_tax_silva.fasta
#usearch10 -uchime2_ref temp/otus.fa -db ${chimera_ref} -strand plus -mode balanced -chimeras temp/otus_chimeras.fa
## rdp16s: All, 21M | 300:58 831Mb   100.0% Chimeras 4371/24553 (17.8%), in db 757 (3.1%), not matched 19425 (79.1%)
#chimera_ref=/mnt/bai/public/ref/rdp/rdp_16s_v16_sp.fa
#usearch10 -uchime2_ref temp/otus.fa -db ${chimera_ref} -strand plus -mode balanced -chimeras temp/otus_chimeras.fa

# 1.9 去除宿主 remove host
make host_rm
# -sintax: 00:46 2.5Gb   100.0% Processing
# 测试比对水稻基因组去宿主，只去年了43个>97%相似，且覆盖度>90%；选择是否是真宿主的阈值不好确定
#blastn -query temp/otus.fa -db ${host} -out temp/otus_no_chimeras.blastn \
#	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' \
#	-num_alignments 1 -evalue 1 -num_threads ${p}
#	awk '$3>97 && $13>90' temp/otus_no_chimeras.blastn |cut -f 1 > $temp/otus_no_chimeras.id
## 可以找到高丰度的菌
## 正对照
#grep '_mitochondria' result/rep_seqs_tax.txt
#grep 'Chloroplast' result/rep_seqs_tax.txt
#grep -P -v 'mitochondria|Chloroplast' temp/otus_no_chimeras.tax | cut -f 1
# 去宿主没有采用gg_13_5的注释结果去除的全面

# 1.10 生成OTU表 Creat OTUs table
make otutab_create
# usearch10 : 05:21:05 1.8Gb   100.0% Searching stripped.fq, 79.6% matched 43622567 / 54819798 mapped to OTUs (79.6%)
# vsearch: Matching query sequences: 43647730 of 54745090 (79.73%) real	24核 54m; 48核 31min; 96核 27min

# 1.11 OTU表样本和OTU过滤 OTU table filter
make otutab_filter 
# 推荐过滤低测序量<5000的样本，筛选至少大于百万分之一丰度的OTU

# 1.12 物种注释 Assign taxonomy
make tax_assign
# RDP train set 16, 21M, 24s, 1.9Gb # 发现OTU_1原来gg注释为线粒体，而rdp中注释不出来，在线rdp也注释不出来。所有之前用gg注释一次并移除宿主
# 三种数据比较
# rdp trainset 16 21Mb: 00:03 106Mb  14470 names, tax levels min 3, avg 6.7, max 7; WARNING: 2 taxonomy nodes have >1 parent; 00:06 1.9Gb
# gg_13_5 97% 151Mb: 2995 names, tax levels min 1, avg 5.0, max 7; WARNING: 1 taxonomy nodes have >1 parent; 00:46 2.5Gb
# silva132 99% 1.1G: 114457 names, tax levels min 2, avg 6.8, max 7; WARNING: 2324 taxonomy nodes have >1 parent; 06:26 6.7Gb; index 04:44 9.0Gb

# 1.13 物种统计 Taxonomy summary
make tax_sum
# ---Fatal error---Line 487 has 3 fields (need 4)
#sed -n '487p' temp/otu.fa.tax|cat -A # 此行缺少4列，以制表符结尾
# 解决方法：1删除这些行；2用第一列补充# grep '\t\n' 

# 1.14 tree_make 多序列比对和进化树 Multiply alignment and make_phylogeny
make tree_make

# 1.15 Alpha多样性指数计算 Calculate alpha diversity index
make alpha_calc
# alpha指数计算结果为 result/alpha/index.txt
# 稀释梯度抽样方法 richness (observed OTUs)-method fast / with_replacement / without_replacement , 结果位于 result/alpha/rare.txt
# 03:30 12Mb

# 1.16 beta_calc Beta多样性进化树和距离矩阵计算 Beta diversity tree and distance matrix
make beta_calc
# ---Fatal error--- ../calcdistmxu.cpp(32) assert failed: QueryUniqueWordCount > 0 致信作者; 改用qiime1

# 1.17 otutab_gg 有参比对，如Greengenes，可用于picurst, bugbase分析
make otutab_gg



# 2. 统计绘图

# 2.1 alpha_boxplot Alpha多样性指数箱线图 Alpha index in boxplot
make alpha_boxplot

# 2.2 alpha_rare Alpha丰富度稀释曲线 Alpha rarefracation curve
make alpha_rare

# 2.3 beta_pcoa 主坐标轴分析距离矩阵 PCoA of distance matrix
make beta_pcoa

# 2.4 beta_cpcoa 限制性主坐标轴分析: OTU表基于bray距离和CCA  CCA of bray distance matrix
make beta_cpcoa

# 2.5 tax_stackplot 样品和组分类学各级别的堆叠柱状图 Stackplot showing taxonomy in each level
make tax_stackplot