## 快速分析 Quick Start(所需文件准备好)

    # 清理零字节文件，用于从头重新分析项目清空makefile点位文件
	find . -name "*" -type f -size 0c | xargs -n 1 rm -f 
    rm -r temp result
    # 建立程序必须目录
	make 10init

    # 1.1 实验设计检查
	make 11validate_mapping
    # 1.2 文库双端合并
	make 12library_merge
    # 1.3 提取Barcode
	make 13extract_barcodes
    # 1.4  拆分文库为样品并质控
	make 14split_libraries
	make 14split_libraries_stat
    # 1.5 切除引物
	make 15fq_trim
    # 1.6 序列去冗余
	make 16fa_unqiue
    # 1.7 挑选OTU 
	make 17otu_pick
    # 1.8 基于参考序列去嵌合
	make 18chimera_ref
    # 1.9 去除宿主 remove host
	make 19host_rm

    # 2.1 生成OTU表
	make 21otutab_create
    # 2.2 OTU表筛选 Filter OTU table
	make 22otutab_filter
    # 2.3 OTU表抽样标准化
	make 23otutab_norm
    # 2.4 物种注释 Assign taxonomy
	make 24tax_assign
    # 2.5 物种分类汇总 Taxonomy summary
	make 25tax_sum
    # 2.6 多序列比对和进化树
	make 26tree_make
    # 2.7 筛选菌identify bac
	make 27identify_isolate


# 1. 处理序列 Processing sequences

	# 0. 准备工作 Preparation

	## 0.1 准备流程配置文件

	# 创建环境代码见~/github/Work/initial_project.sh
	# 设置工作目录
	wd=culture/medicago/190626
	## 准备实验设计
	cd ~/$wd
	# Initialize the working directory
	make 10init

	# 保存模板中basic页中3. 测序文库列表library为doc/library.txt，参考SeqLibraryList.xlsx，与/mnt/bai/yongxin/seq/amplicon目录中文件对应
	# 按library中第二列index准备测序文库，如果压缩要添加.gz，并用gunzip解压
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
	awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
    # 检查数据链接，全红为错误，绿色为正常
    ll seq/*
	# 如果压缩文件，要强制解压链接
	gunzip -f seq/*.gz


	## 写mappingfile, s为物种，p为板数；多个library需要多个mappingfile
	# 可指定品种、部位和培养基类型
	# 单个文库
	write_mappingfile_culture2.pl -o doc/L1.txt -s medicago -L L1 -v A17 -c Rhizosphere -m R2A -B 1 -p 1
	# 批量相同属性文库
	for i in `seq 1 10`; do write_mappingfile_culture2.pl -o doc/L${i}.txt -s medicago -L L${i} -v A17 -c Rhizosphere -m R2A -B 1 -p 48; done
	# 按Library信息批量生成
	# awk '{if(NR>1){system("echo "$1" "$6" "$7" "$8)}}' doc/library.txt
	awk '{if(NR>1){system("write_mappingfile_culture2.pl -o doc/"$1".txt -s medicago -L "$1" -v "$6" -c "$7" -m "$8" -B "$9" -p "$10)}}' doc/library.txt
	# L9, L10手动修改个性化数据，在Excel中手动修改 


	# 删除多余空格，windows换行符等(MAC用户禁用)
	sed -i 's/ //g;s/\r//' doc/*.txt
	# 查看数据格式、列名，了解和检查实验设计
	head -n3 doc/L1.txt
	# 依据各文库L*.txt文件生成实验设计
	cat <(head -n1 doc/L1.txt | sed 's/#//g') <(cat doc/L* |grep -v '#'|grep -v -P '^SampleID\t') > doc/design.txt
	# 检查是否相等
	wc -l doc/design.txt
	cut -f 1 doc/design.txt|sort|uniq|wc -l
	# 查看冗余的列(仅上方不等时使用)
	cut -f 1 doc/design.txt|sort|uniq -c| less


## 1.2. 按实验设计拆分文库为样品


	# 拆分样品
	head -n3 doc/L1.txt
	# 按L1/2/3...txt拆分library为samples
	# 输入为seq/L*.fq，输出为seq/sample/*.fq
	make library_split
	make library_split_stat
	# 统计结果见result/split有txt/pdf/png，推荐看png方便快速查看每张位图
	# 查看样本量排序
	sort -k2,2n result/sample_split.log|less

## 1.3. 样品双端合并、重命名、合并为单一文件

	# Merge paired reads, renames and merge all samples
	# 样品双端合并、重命名、合并为单一文件, 注意fastq为33格式，64位采用fastp转换
	# 输入为seq/sample/*.fq，输出为seq/all.fq
	make sample_merge
	make sample_merge_stat
	# result/sample_merge.log中有每个样本合并后的序列数量


## 1.4. 切除引物与标签

	# Cut primers and lables
	# 切除左端标签和引物，右端 引物
	# Cut barcode 10bp + V5 19bp in left， and V7 18bp in right
	# 输入为seq/all.fq，输出为temp/stripped.fq
	make fq_trim


## 1.5. 质量控制

	# Quality control
	# 过滤序列中预期累计错误率>1%的序列
	# 输入为temp/stripped.fq，输出为temp/filtered.fa
	make fq_qc



    # (第一阶段结束，获得纯净扩增子序列temp/filtered.fa，可提供此文件从下面开始)


## 1.6. 序列去冗余

	# Remove redundancy, get unique reads
	# 输入为temp/filtered.fa，输出为temp/uniques.fa
	make fa_unqiue


## 1.7. 挑选OTU

	# Pick OTUs
	# unoise3速度比cluster_otus慢上百倍，更精细但结果也更多
	# 输入为temp/uniques.fa，输出为temp/Zotus.fa
	make otu_pick


## 1.8. 有参去嵌合体

	# Remove chimiras by silva database
	# 基于SILVA数据库去除
	make chimera_ref


## 1.9. 去除宿主

	# Remove host
	# 根据SILVA注释去除线粒体、叶绿体、真核生物18S和未知序列(非rRNA)
	make host_rm


    # (第二阶段结束，获得OTU代表序列result/otu.fa，可提供此文件和测序数据temp/filtered.fa从下方起始)


## 1.10. 生成OTU表
	
	# Create OTUs table
	# 默认使用vsearch更快10倍，可选usearch10，线程不可超48
	make otutab_create


## 1.11. 过滤样本和OTUs

	# OTU table filter samples and OTU
	# 推荐过滤低测序量<5000的样本，筛选大于1RPM的OTU
	make otutab_filter 


## 1.12. 物种注释

	# Assign taxonomy
	# 默认使用RDP trainset快而准，GG太旧，Silva太慢
	# 推荐阈值为0.6保证注释更完整
	make tax_assign


## 1.13. 物种统计
	
	# Taxonomy summary
	# 必须所有物种有注释，否则有可能报错
	make tax_sum


## 1.14. 多序列比对和进化树
	
	# Multiply alignment and make_phylogeny
	# usearch10/culsterO结果不同可能影响多样性分析(usearch unifrac结果更可信)
	# 进化树，用于树图和多样性分析
	make tree_make


# 与实验菌建立联系

## 注释分菌孔信息
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$8"\t"$9"\t"$10} NR>FNR {print $0,a[$1]}' doc/design.txt result/culture_bacteria.xls > result/culture_bacteria_anno.xls

## 筛选纯菌并建索引
	tail -n+2 result/culture_select.xls | cut -f 1,2|sed 's/^/OTU_/g;s/;/\t/g'|less>result/culture_select.tax
	filter_fasta.py -f result/otu.fa -o result/culture_select.fa -s result/culture_select.tax
	sed -i 's/OTU/COTU/' result/culture_select.fa
	makeblastdb -dbtype nucl -in result/culture_select.fa
	makeblastdb -dbtype nucl -in result/otu.fa


## 2019/7/3 与最新10个菌库比较，找OTU_2 Bacillales
    cwd=culture10
    mkdir -p ${cwd}
    blastn -query ~/medicago/AMF2/result/otu.fa -db result/culture_select.fa -out ${cwd}/otu_culture.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 # 输出13列为coverage
    # 查看关注菌对应的编号 
    less -S ${cwd}/otu_culture.blastn 
    blastn -query ~/medicago/AMF2/result/otu.fa -db result/otu.fa -out ${cwd}/otu_culture_all.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 # 输出13列为coverage

    # 进一步筛选第三位出现的菌是否存在有OTU_61，修改程序名，输出文件名+2，每个孔有Top3的结果
    Rscript ~/github/Amplicon/16Sculture2/script/identify_isolate2.R
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$8"\t"$9"\t"$10} NR>FNR {print $0,a[$1]}' doc/design.txt result/culture_bacteria2.xls > result/culture_bacteria2_anno.xls
    # 检查了第三列中OTU_61，数据只只有1-5条，且丰度在0.5%以下

    ## 筛选某个指定OTU
    Rscript ~/github/Amplicon/16Sculture2/script/select_OTU.R
    # 添加孔实验设计注释
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$8"\t"$9"\t"$10} NR>FNR {print $0,a[$1]}' doc/design.txt result/otu61.txt > result/otu61.anno.txt
    # 添加孔详细
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print $0,a[$1]}' result/culture_bacteria2_anno.xls result/otu61.txt > result/otu61.anno.txt
    # 添加每个COTU与OTU1的相似度
    echo 'OTU_2' > ${cwd}/otu2.id
	filter_fasta.py -f ~/medicago/AMF2/result/otu.fa -o ${cwd}/bacillus.fa -s ${cwd}/otu2.id
    makeblastdb -dbtype nucl -in ${cwd}/bacillus.fa
    blastn -query result/otu.fa -db ${cwd}/bacillus.fa -out ${cwd}/cotu_otu2.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 # 输出13列为coverage
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$3} NR>FNR {print $0,a[$8]}' ${cwd}/cotu_otu2.blastn result/otu61.anno.txt > result/otu61.anno2.txt

    # 与其它物种菌库比对
    # 拟南芥
    blastn -query ${cwd}/bacillus.fa -db ~/culture/ath/result/Rootculture_select.fa -out ${cwd}/ath_otu2.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9
    less  ${cwd}/ath_otu2.blastn # 有COTU_317 100%匹配
    # 与菌保比对
    blastn -query ${cwd}/bacillus.fa -db ~/culture/ath/wet/root_leaf.fa -out ${cwd}/ath_otu2_stock.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9
    less  ${cwd}/ath_otu2_stock.blastn # 有COTU_317 100%匹配

    # 水稻
    blastn -query ${cwd}/bacillus.fa -db ~/culture/rice/result/culture_select.fasta -out ${cwd}/rice_otu2.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9
    less  ${cwd}/rice_otu2.blastn # 有rice_311 100%匹配
    # 直接比对水稻菌保
    blastn -query ${cwd}/bacillus.fa -db ~/culture/rice/190626/sanger_stock/16S_rice_culture_collection.fasta -out ${cwd}/rice_otu2.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9
    less  ${cwd}/rice_otu2.blastn # 与3菌株100%一致

    # 小麦
    blastn -query ${cwd}/bacillus.fa -db ~/culture/wheat/result/culture_select.fa -out ${cwd}/wheat_otu2.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9
    less  ${cwd}/wheat_otu2.blastn # 有rice_311 100%匹配，OTU_128


## 2019/10/30 与最新10个菌库比较，苜蓿更新为AMF3，找OTU_3 Gammaproteobacteria	Pseudomonadales	Pseudomonadaceae	Pseudomonas下降，但OTU_1 Gammaproteobacteria	Pseudomonadales	Pseudomonadaceae	Rhizobacter上升
    cwd=culture10pseudo
    mkdir -p ${cwd}
    blastn -query ~/medicago/AMF3/result/otu.fa -db result/culture_select.fa -out ${cwd}/otu_culture.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 # 输出13列为coverage
    # 查看关注菌对应的编号，OTU——3最接近的是COTU_111 100%和COTU_543 99.47%,2个错配
    less -S ${cwd}/otu_culture.blastn 
    # 查看和总表中，还有更好的吗，发现OTU_111存在100%
    blastn -query ~/medicago/AMF3/result/otu.fa -db result/otu.fa -out ${cwd}/otu_culture_all.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 10 -evalue 1 -num_threads 9 # 输出13列为coverage
    less result/otu.fa ${cwd}/otu_culture_all.blastn
    # 进一步筛选第三位出现的菌是否存在有OTU_111和OTU_543，修改程序名，输出文件名+2，每个孔有Top3的结果
    Rscript ~/github/Amplicon/16Sculture2/script/identify_isolate2.R
    awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$8"\t"$9"\t"$10} NR>FNR {print $0,a[$1]}' doc/design.txt result/culture_bacteria2.xls > result/culture_bacteria2_anno.xls


## Sanger鉴定 2019/10/22

    mkdir -p wet && cd wet
    mkdir -p temp
    # 本次使用拼接的全长

    # 第一批Sanger
    rm -r batch1
    # cp -r ~/medicago/AMF/wet/stock_sanger/list_detail.txt batch1.txt # 有修改
    cp -r ~/medicago/AMF/wet/stock_sanger/contig batch1/
    # 基本统计
    ls batch1/ | wc -l # 307条序列
    grep 'GGATTAGATACCC' batch1/* | wc -l # 正向173
    grep 'ACGTCATCCCCACCTTCC' batch1/* | wc -l # 反向115

    # 序列存在正反，正向移动到目录，反向调整
    rm -r 190130
    mkdir -p 190130
    for i in `grep 'GGATTAGATACCC' batch1/*| cut -f 1 -d ':'`; do
        mv $i 190130/
    done
    # 筛选反向互补
    rm temp/* 
    for i in `grep 'ACGTCATCCCCACCTTCC' batch1/*| cut -f 1 -d ':'`; do
        mv $i temp/
    done
    for i in `ls temp/`; do
    	revseq -sequence temp/$i -outseq 190130/$i -tag N
    done
    # 处理其它引物不匹配序列，直接截取，长度正常如下：
    # 5个正向的 MP02C3 MP02C7 MP02C8 MP02C8 MP03D12
    mv batch1/MP02* 190130/
    mv batch1/MP03* 190130/
    mv batch1/MP05* 190130/
    # 取反向
    for i in `ls batch1/`; do
    	revseq -sequence batch1/$i -outseq 190130/$i -tag N
    done
    # 序列改名
    rename 's/SEQ/seq/' 190130/*.txt
    rename 's/27F.//' 190130/*.txt
    awk 'BEGIN{OFS=FS="\t"}{system("mv 190130/"$2".seq.txt 190130/"$1".seq")}' batch1.txt
    # 合并为fa文件 seq merge and format to fasta
    sed -i '/>/d' 190130/*
    rm -rf 190130/MP*
    format_seq2fasta.pl -i "190130/*" -o 190130.fa 
    # 提取V5-V7区域
	# 检查测序区域，以 180614.fa 为列 V5: AACM GGATTAGATACCC KG V7: ACGTCATCCCCACCTTCC ，匹配正向引物成功
	cutadapt -g AACMGGATTAGATACCCKG -e 0.2 190130.fa -o temp/cut5.fa --discard-untrimmed
	# 反向时先切反向引物，1193R ACGTCATCCCCACCTTCC，正常为1492R GGTTACCTTGTTACGACTT找不到，可选 1391R GACGGGCGGTGTGTRCA F TGYACACACCGCCCGTC
	cutadapt -a GGAAGGTGGGGATGACGT -e 0.2 temp/cut5.fa -o 190130V5-7.fa
	# 长度分布，370-385，多数为379
	grep -v '>' 190130V5-7.fa|awk '{print length($0)}'|sort -n|uniq -c
	# 冗余度
	grep -c '>' 190130V5-7.fa # 304
	grep -v '>' 190130V5-7.fa | sort|uniq -c|wc -l # 153


    # 第二批菌 2019/10/22
    cp -r /mnt/bai/jingmei/Bacteria/* ./batch2
    # 基本统计
    ls batch2/ | wc -l # 366条序列
    grep 'GGATTAGATACCC' batch2/* | wc -l # 正向352
    grep 'ACGTCATCCCCACCTTCC' batch2/* | wc -l # 反向0
    # 批量改名
    ## 生成改名列表
    paste <(ls batch2/ | sed 's/-2/b/;s/-3/c/' | cut -f 2 -d '-' | cut -f 1 -d '.' | sed 's/^/M/;s/MM/M/') <(ls batch2/ | sed 's/.seq//') > filelist.txt
    mkdir -p 191022
    awk 'BEGIN{OFS=FS="\t"}{system("ln batch2/"$2".seq 191022/"$1".seq")}' filelist.txt
    # 文件内有多余的序列名，删除
    sed -i '/>/d' 191022/*.seq
    # 合并为fa文件 seq merge and format to fasta
    format_seq2fasta.pl -i "191022/*.seq" -o 191022.fa 
	# 检查测序区域，以 180614.fa 为列 V5: AACM GGATTAGATACCC KG V7: ACGTCATCCCCACCTTCC，匹配正向引物成功
	cutadapt -g AACMGGATTAGATACCCKG -e 0.2 191022.fa -o temp/cut5.fa
	# 反向时先切反向引物，1193R ACGTCATCCCCACCTTCC，正常为1492R GGTTACCTTGTTACGACTT找不到，可选 1391R GACGGGCGGTGTGTRCA F TGYACACACCGCCCGTC
	cutadapt -a GGAAGGTGGGGATGACGT -e 0.2 temp/cut5.fa -o 191022V5-7.fa
	# 长度分布，370-385，多数为379
	grep -v '>' 191022V5-7.fa|awk '{print length($0)}'|sort -n|uniq -c
	# 冗余度
	grep -c '>' 191022V5-7.fa # 366
	grep -v '>' 191022V5-7.fa | sort|uniq -c|wc -l # 168

    # 第三批菌 2020/5/14
    # 使用拼接好的序列，删除其它seq序列
    rm batch3/*.seq
    sed -i '/>/d' batch3/*.txt
    rename 's/txt/seq/' batch3/*.txt
    # 基本统计
    ls batch3/ | wc -l # 12条序列
    grep 'GGATTAGATACCC' batch3/* | wc -l # 正向12
    grep 'ACGTCATCCCCACCTTCC' batch3/* | wc -l # 反向0
    # 批量改名
    ## 生成改名列表
    paste <(ls batch3/ | sed 's/-2/b/;s/-3/c/' | cut -f 2 -d '-' | cut -f 1 -d '.' | sed 's/^/M/;s/MM/M/') <(ls batch3/ | sed 's/.seq//') > filelist.txt
    mkdir -p 200514
    awk 'BEGIN{OFS=FS="\t"}{system("ln batch3/"$2".seq 200514/"$1".seq")}' filelist.txt
    # 文件内有多余的序列名，删除
    sed -i '/>/d' 200514/*.seq
    sed -i 's/\r//' 200514/*.seq
    # 合并为fa文件 seq merge and format to fasta
    format_seq2fasta.pl -i "200514/*.seq" -o 200514.fa 
	# 检查测序区域，以 180614.fa 为列 V5: AACM GGATTAGATACCC KG V7: ACGTCATCCCCACCTTCC，匹配正向引物成功
	cutadapt -g AACMGGATTAGATACCCKG -e 0.2 200514.fa -o temp/cut5.fa
	# 反向时先切反向引物，1193R ACGTCATCCCCACCTTCC，正常为1492R GGTTACCTTGTTACGACTT找不到，可选 1391R GACGGGCGGTGTGTRCA F TGYACACACCGCCCGTC
	cutadapt -a GGAAGGTGGGGATGACGT -e 0.2 temp/cut5.fa -o 200514V5-7.fa
	# 长度分布，370-385，多数为379
	grep -v '>' 200514V5-7.fa|awk '{print length($0)}'|sort -n|uniq -c
	# 冗余度
	grep -c '>' 200514V5-7.fa # 12
	grep -v '>' 200514V5-7.fa | sort|uniq -c|wc -l # 5

    # 多批合并分析
    cat *V5-7.fa > merge.fa
    grep -c '>' merge.fa # 682
    makeblastdb -dbtype nucl -in  merge.fa

    # 菌保分析
	extract.pl merge.fa > temp/merge_uniq.fa
    cut -f 1 -d ',' temp/merge_uniq.fa > merge_uniq.fa
    grep -c '>' merge_uniq.fa # 243个V5-V7非冗余序列
    # fasta转换
	awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\t":$0 }' merge_uniq.fa | sed 's/>//' > temp/fasta1line.txt
    # 添加物种注释
    usearch10 -sintax merge_uniq.fa \
	-db /mnt/bai/public/ref/rdp/rdp_16s_v16_sp.udb -sintax_cutoff 0.6 -strand both \
	-tabbedout temp/otu.fa.tax -threads 32
    mkdir -p result
    cut -f 1,4 temp/otu.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > result/taxonomy_2.txt
    awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="Unclassified";a["p"]="Unclassified";a["c"]="Unclassified";a["o"]="Unclassified";a["f"]="Unclassified";a["g"]="Unclassified";a["s"]="Unclassified"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' result/taxonomy_2.txt | sed '1 i #OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' | sed '1 s/#OTU ID/OTUID/' > result/taxonomy_8.txt
    # 制作来源文件，添加多个来源 b[$1].="," ;print $0,b[$1]
    # awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$4} NR>FNR {split($1,x,";");for(i in x){print "1"i"2"a[i]}}' merge_metadata.txt result/taxonomy_8.txt 
    # ./detect_seq_distribution.pl 志文缩写脚本，也需要修改
    ./detect_seq_distribution.pl -i temp/merge_uniq.fa -d merge_metadata.txt -o result/variety.txt
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3} NR>FNR {print $0,a[$1]}' result/variety.txt result/taxonomy_8.txt | sed 's/A17/Yes/;s/NA/NO/;s/R108/Yes/;' | sed '1 s/$/Yes\tYes/' > result/taxonomy_9.txt


