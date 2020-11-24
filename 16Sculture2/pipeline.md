SHELL:=/bin/bash

	# 16S扩增子分析流程第二版 —— 培养菌鉴定 16S Amplicon pipeline version 2 —— Culture identify

	# 帮助文档
	
	# 流程所需的软件、脚本及数据库版本
	# Help: Pipeline dependency version of softwares, scripts and databases

version: 

	# 2019/6/26 2.1 Standard 16S anlysis report
	# 
	# 软件 Softwares
	# QIIME 1.9.1
	# fastqc -v # 0.11.5 测序数据质量评估
	# parallel # 并行/多线程任务管理
	# fastp -v # version 0.12.5 fastq文件质控、质量类型转换
	# usearch10 --help # v10.0.240_i86linux64 扩增子分析软件
	# clustalo --version # 1.2.1 多序列对齐
	# filter_alignment.py # from qiime1.9.1 筛选比对序列及区域
	# make_phylogeny.py # from qiime1.9.1 默认调用fastree建树
	# biom --version # 2.1.5， OTU表格式转换
	# 
	# 脚本 Scripts
	# alpha_boxplot.sh # 1.0 基于usearch alpha_div绘制箱线图
	#
	# 数据库 Databases
	# greengene13_5 # reference, and/or taxonomy database 
	# rdp train set 16 # taxonomy database
	# Silva132 # chimera reference, and/or taxonomy database 



# 1. 标准流程 Standard pipeline

	# 建立程序必须目录 Create work directory

10init:

	touch $@
	mkdir -p seq doc temp result script


## 1.1 实验设计检查 Validate mapping file

11validate_mapping: 

	touch $@
	mkdir -p doc/validate
	parallel -j ${p} "validate_mapping_file.py -o doc/validate/ -m doc/{1}.txt" \
		::: `tail -n+2 doc/library.txt | cut -f 1`


## 1.2 文库双端合并 Merge clean reads

12library_merge: 11validate_mapping

	touch $@
	# 双端序列合并
	mkdir -p temp/merge
	parallel -j ${p} "usearch10 -fastq_mergepairs seq/{1}_1.fq -reverse seq/{1}_2.fq \
		-fastqout temp/merge/{1}.fq -relabel {1}. -threads 1" ::: `tail -n+2 doc/library.txt | cut -f 1`


## 1.3 提取Barcode

13extract_barcodes: 12library_merge

	touch $@
	parallel -j ${p} "extract_barcodes.py -f temp/merge/{1}.fq -m doc/{1}.txt \
		-o temp/{1}_barcode -c barcode_paired_stitched --bc1_len 10 --bc2_len 6 -a --rev_comp_bc2" \
		::: `tail -n+2 doc/library.txt | cut -f 1`


## 1.4  拆分文库为样品并质控 Split library into sample and quality control

14split_libraries: 13extract_barcodes
	touch $@
	parallel -j ${p} "split_libraries_fastq.py -i temp/{}_barcode/reads.fastq -b temp/{}_barcode/barcodes.fastq \
		-m doc/{1}.txt -o temp/{}_split/ \
		-q 19 --max_barcode_errors 0 --barcode_type 16 --phred_offset 33" \
		::: `tail -n+2 doc/library.txt | cut -f 1`

14split_libraries_stat: 14split_libraries
	touch $@
	mkdir -p result/split
	parallel -j ${p} "tail -n+16 temp/{}_split/split_library_log.txt|head -n-4 > result/split/{}.txt" \
		::: `tail -n+2 doc/library.txt | cut -f 1`
	# 绘图：高5，长40，字体为1，才能看清近5000个样本编号
	parallel -j 1 "stat_16s_lib_split2.sh -o result/split/ -A plate -d `pwd`/doc/design.txt -l {} -h 5 -w 40 -s 1" \
		::: `tail -n+2 doc/library.txt | cut -f 1`
	# 转换qiime结果为usearch格式
	rm -rf temp/qc.fa
	parallel -j 1 "cut -f 1 -d ' ' temp/{}_split/seqs.fna | sed 's/_/./' >>  temp/qc.fa" \
		::: `tail -n+2 doc/library.txt | cut -f 1`


## 1.5 切除引物 Cut primers and quality filter

	# Cut barcode 10bp + V5 19bp in left and V7 18bp in right

15fq_trim: 14split_libraries_stat

	touch $@
	usearch10 -fastx_truncate temp/qc.fa \
		-stripleft ${stripleft} -stripright ${stripright} \
		-fastaout temp/filtered.fa 


## 1.6 序列去冗余 Remove redundancy
	# miniuniqusize为8，去除低丰度，增加计算速度

16fa_unqiue: 15fq_trim

	touch $@
	vsearch --derep_fulllength temp/filtered.fa \
		--relabel Uni --minuniquesize ${minuniquesize} --sizeout \
		--output temp/uniques.fa 
	echo -ne 'Unique reads\t' > ${otu_log}
	grep -c '>' temp/uniques.fa >> ${otu_log}
	cat ${otu_log}


## 1.7 挑选OTU Pick OTUs

	# 可选97% cluster_otus，或100% unoise3，默认unoise3，不支持多线程
	# ifeq条件必须顶格，否则报错

17otu_pick: 16fa_unqiue

	touch $@
	echo -e "OTU method\t${otu_method}" >> ${otu_log}
ifeq (${otu_method}, unoise3)
	# 类似100%聚类，只去除嵌合体和扩增及错误，保留所有高丰度序列
	usearch10 -unoise3 temp/uniques.fa -zotus temp/Zotus.fa -minsize ${minuniquesize}
else ifeq (${otu_method}, cluster_otus)
	# cluster_otus无法修改聚类参数，想使用不同聚类相似度，使用cluster_smallmem命令
	usearch10 -cluster_otus temp/uniques.fa -otus temp/Zotus.fa
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in unoise3 or cluster_otus") 
endif
	awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' temp/Zotus.fa > temp/otus.fa
	echo -ne 'OTU number\t' >> ${otu_log}
	grep -c '>' temp/otus.fa >> ${otu_log}
	cat ${otu_log}


## 1.8 基于参考序列去嵌合 Remove chemira by silva

	# 可选，推荐使用最新SILVA大数据库 https://www.arb-silva.de/
	# 官方模式推荐sensitive错，推荐balanced, high_confidence

18chimera_ref: 17otu_pick

	touch $@
ifeq (${chimera_mode}, none)
	cp temp/otus.fa temp/otus_no_chimeras.fa
	echo -ne 'not remove chimeras\t' >> ${otu_log}
else
	usearch10 -uchime2_ref temp/otus.fa \
		-db ${chimera_ref} -strand plus -mode ${chimera_mode} \
		-chimeras temp/otus_chimeras.fa -threads ${p}
	# 获得非嵌合体序列ID
	cat temp/otus.fa temp/otus_chimeras.fa| grep '>' | sort | uniq -u | sed 's/>//' > temp/no_chimeras.id
	# 筛选非嵌合体
	usearch10 -fastx_getseqs temp/otus.fa -labels temp/no_chimeras.id -fastaout temp/otus_no_chimeras.fa
	echo -ne 'no chimeras\t' >> ${otu_log}
endif
	grep -c '>' temp/otus_no_chimeras.fa >> ${otu_log}
	cat ${otu_log}


## 1.9 去除宿主 remove host

19host_rm: 18chimera_ref

	touch $@
ifeq (${host_method}, blast)
	# 方法1. 基于宿主基因组(含叶绿体/线粒体)比对
	blastn -query temp/otus_no_chimeras.fa -db ${host} -out temp/otus_no_chimeras.blastn \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' \
	-num_alignments 1 -evalue 1 -num_threads ${p}
	awk '$$3>${host_similarity} && $$13>${host_coverage}' temp/otus_no_chimeras.blastn | cut -f 1 | sort | uniq > temp/otus_host.id
	cat <(grep '>' temp/otus_no_chimeras.fa|sed 's/>//') temp/otus_host.id | sort | uniq -u > temp/otus_no_host.id
else ifeq (${host_method}, sintax_silva)
	# 方法2. 基于silva132的usearch注释结果筛选，排除线粒体、叶绿体、真核和非16S。仍有结果RDP会注释Chloroplast
	usearch10 -sintax temp/otus_no_chimeras.fa \
		-db ${usearch_silva} -sintax_cutoff ${sintax_cutoff} -strand both \
		-tabbedout temp/otus_no_chimeras.tax -threads ${p}
	grep -P -v 'Mitochondria|Chloroplast|Eukaryota|\t$$' temp/otus_no_chimeras.tax | cut -f 1 > temp/otus_no_host.id
else ifeq (${host_method}, sintax_unite)
	# 方法3. 基于unite筛选真菌
	# usearch使用unite注释ITS筛选真菌
	usearch10 -sintax temp/otus_no_chimeras.fa \
		-db ${usearch_unite} -sintax_cutoff ${sintax_cutoff} -strand both \
		-tabbedout temp/otus_no_chimeras.tax -threads ${p}
#	grep -P -v 'Fungi|\t$$' temp/otus_no_chimeras.tax | cut -f 1 > temp/otus_no_fungi.id
#	cat <(grep '>' temp/otus_no_chimeras.fa|sed 's/>//') temp/otus_host.id temp/otus_no_fungi.id | sort | uniq -u > temp/otus_no_host.id
	grep 'd:Fungi,p' temp/otus_no_chimeras.tax | cut -f 1 > temp/otus_no_host.id
else ifeq (${host_method}, none)
	# 方法4. 不过滤宿主和质体
	cut -f 1 temp/otus_no_chimeras.tax > temp/otus_no_host.id
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in blast, usearch_gg or usearch_silva") 
endif
	# 最终筛选结果复制入result
	usearch10 -fastx_getseqs temp/otus_no_chimeras.fa -labels temp/otus_no_host.id -fastaout temp/otus_no_host.fa
	echo -ne 'no host\t' >> ${otu_log}
	grep -c '>' temp/otus_no_host.fa >> ${otu_log}
	cat ${otu_log}
	cp temp/otus_no_host.fa result/otu.fa


## 2.1 生成OTU表 Creat OTUs table

21otutab_create: 19host_rm

	touch $@
ifeq (${map_method}, usearch10)
	# 方法1：usearch10, temp/stripped.fq建议新流程使用，与旧体系只有filtered.fa兼容
	usearch10 -otutab temp/filtered.fa -otus result/otu.fa -id ${map_identify} \
		-otutabout temp/otutab.txt -threads ${p}
else ifeq (${map_method}, vsearch)
	# 方法2：vsearch，比usearch10更快
	time vsearch --usearch_global temp/filtered.fa --db result/otu.fa --id ${map_identify} \
		--otutabout temp/otutab.txt --threads ${p}
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in usearch10 or vsearch") 
endif


## 2.2 OTU表筛选 Filter OTU table

22otutab_filter: 21otutab_create

	touch $@
	# 统计OTU表的基本信息
	usearch10 -otutab_stats temp/otutab.txt -output temp/otutab.txt.stat
	echo -ne "\nOTU table summary\nType\tThreshold\tReads\tSamples\tOTUs\nTotal\t-\t" > ${log_otutable}
	head -n3 temp/otutab.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	biom convert -i temp/otutab.txt -o temp/otutab.biom --table-type="OTU table" --to-json
	biom summarize-table -i temp/otutab.biom > temp/otutab.biom.sum
	# 按样本测序量筛选：通常低于5000的样本会删除
	usearch10 -otutab_trim temp/otutab.txt -min_sample_size ${min_sample_size} -output temp/otutab_trim1.txt
	usearch10 -otutab_stats temp/otutab_trim1.txt -output temp/otutab_trim1.txt.stat
	echo -ne "SampleSize\t${min_sample_size}\t" >> ${log_otutable}
	head -n3 temp/otutab_trim1.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	# 按OTU测序量筛选：通常低于8的OTU会删除
	usearch10 -otutab_trim temp/otutab_trim1.txt -min_otu_size ${min_otu_size} -output temp/otutab_trim2.txt
	usearch10 -otutab_stats temp/otutab_trim2.txt -output temp/otutab_trim2.txt.stat
	echo -ne "OtuSize\t${min_otu_size}\t" >> ${log_otutable}
	head -n3 temp/otutab_trim2.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	# 按OTU相对丰度筛选：通常低于1 RPM 的OTU会删除
	usearch10 -otutab_trim temp/otutab_trim2.txt -min_otu_freq ${min_otu_freq} -output temp/otutab_trim3.txt
	usearch10 -otutab_stats temp/otutab_trim3.txt -output temp/otutab_trim3.txt.stat
	echo -ne "OtuFreq\t${min_otu_freq}\t" >> ${log_otutable}
	head -n3 temp/otutab_trim3.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	# 复制最终版OTU表到结果目录
	cp temp/otutab_trim3.txt result/otutab.txt
	# 转换为biom格式
	biom convert -i result/otutab.txt -o result/otutab.biom --table-type="OTU table" --to-json
	# 统计OTU表
	biom summarize-table -i result/otutab.biom > result/otutab.biom.sum
	head -n 30 result/otutab.biom.sum


## 2.3 OTU表抽样标准化

23otutab_norm: 22otutab_filter

	touch $@
	# 依据最小样本量，设置OTU表标准化的阈值，如我们看到最小样品数据量为3.1万，可以抽样至3万
	usearch11 -otutab_rare result/otutab.txt -sample_size ${sample_size} -output result/otutab_norm.txt 
	usearch10 -otutab_stats result/otutab_norm.txt -output result/otutab_norm.txt.stat
	echo -ne "OtuNorm\t${sample_size}\t" >> ${log_otutable}
	head -n3 result/otutab_norm.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	cat ${log_otutable}
	# 转换为biom格式
	biom convert -i result/otutab_norm.txt -o result/otutab_norm.biom --table-type="OTU table" --to-json


## 2.4 物种注释 Assign taxonomy

24tax_assign: 23otutab_norm

	touch $@
	usearch10 -sintax result/otu.fa \
		-db ${sintax_db} -sintax_cutoff ${sintax_cutoff} -strand both \
		-tabbedout temp/otu.fa.tax -threads ${p}


## 2.5 物种分类汇总 Taxonomy summary

25tax_sum: 24tax_assign

	touch $@
	mkdir -p result/tax
	# 未分类的添加末注释标记，否则汇总时报错
	sed -i 's/\t$$/\td:Unassigned/' temp/otu.fa.tax
	# 按门、纲、目、科、属水平分类汇总
	# 默认用otutab_norm.txt，用otutab.txt是不是更好呢？
	for i in p c o f g;do \
		usearch10 -sintax_summary temp/otu.fa.tax -otutabin result/otutab.txt -rank $${i} \
			-output result/tax/sum_$${i}.txt; \
	done
	# 删除Taxonomy中异常字符如() " - /
	sed -i 's/(//g;s/)//g;s/\"//g;s/\/Chloroplast//g;s/\-/_/g;s/\//_/' result/tax/sum_*.txt
	# 格式化物种注释：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
	cut -f 1,4 temp/otu.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > result/taxonomy_2.txt
	# 生成物种表格：注意OTU中会有末知为空白，补齐分类未知新物种为Unassigned
	awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned"; split($$2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $$1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' result/taxonomy_2.txt | sed '1 i #OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > result/taxonomy_8.txt
	# 去除#号和空格，会引起读取表格分列错误
	sed -i 's/#//g;s/ //g' result/taxonomy_8.txt
	# 添加物种注释
	biom add-metadata -i result/otutab.biom --observation-metadata-fp result/taxonomy_2.txt -o result/otutab_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
	# 添加物种注释
	biom add-metadata -i result/otutab_norm.biom --observation-metadata-fp result/taxonomy_2.txt -o result/otutab_norm_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
	# 制作门+变形菌纲混合比例文件
	cat <(grep -v 'Proteobacteria' result/tax/sum_p.txt) <(grep 'proteobacteria' result/tax/sum_c.txt) > result/tax/sum_pc.txt


## 2.6 多序列比对和进化树 Multiply alignment and make_phylogeny

26tree_make: 25tax_sum

	touch $@
	# clustalo+qiime1.9.1
	clustalo -i result/otu.fa -o temp/otu_align.fa --seqtype=DNA --full --force --threads=${p}
	make_phylogeny.py -i temp/otu_align.fa -o result/otu.tree


## 2.7 筛选菌identify bac

27identify_isolate: 26tree_make

	touch $@
	Rscript /mnt/bai/yongxin/github/Amplicon/16Sculture2/script/identify_isolate.r


## 9.9 清理中间文件

	# 完成大数据分析后，有太多临时文件，需要清楚节约空间
clean:

	pigz seq/L*.fq
	rm -r temp/
