## 16S扩增子主流程，修改请保证向前兼容，最好只添加新分枝流程
## Amplicon 16S main pipeline. Modify need meet the previous running good, you would better only add new command
#
#
#
#version: 
#	# 帮助文档：流程所需的软件、脚本及数据库版本
#	# 16S Amplicon pipeline v1.4 2017/12/23
#	# 2017/6/28 1.0 Standard 16S reports, include alpha/beta diversity, differentially abundance (DA) taxonomy in each level, DA OTU, and high-level analysis include pie, venn and ternary plot.
#	# 2017/9/25 1.1 add format_taxonomy2full.pl to degub uncomplete taxonomy; usearch8 cluster_otus must add -ids 0.97, because default is 0.87
#	# 2017/10/13 1.2 change cluster_otus to unoise3, greengene13.8 to rdp11.5 97%
#	# 2017/12/1 1.3 Add multipleQC, 
#	# 2017/12/23 1.4 Add remove old file of merge library, split, cutadapt for debug error for exist file; add usearch rarefraction & plot, split main makefile to ref/amplicon/16s
#	#
#	# Operation System: Ubuntu 16.04.3 LTS
#	#
#	# # Software
#	# dos2unix -V # 6.0.4 (2013-12-30) 去除文本中的windows特有换行符^M
#	# fastqc -v # 0.11.5 测序数据质量控制
#	# multiqc --version # version 1.3.dev0 测序报告汇总
#	# qiime --version # 1.9.1 扩增子测序分析流程
#	#	join_paired_ends.py 双端序列合并为单端
#	#	extract_barcodes.py 提取样品序列单端/两端的barcode 
#	#	split_libraries_fastq.py 按barcode对文库进行样品拆分
#	#	align_seqs.py 多序列比对
#	#	filter_fasta.py 过滤OTU代表序列fasta
#	#	filter_samples_from_otu_table.py 按样本测序量筛选
#	#	filter_alignment.py 过滤多序列比对结果中75%相似序列及保守区
#	#	make_phylogeny.py 调用fasttree进行建树
#	#	single_rarefaction.py 对OTU表进行单次等量抽样
#	#	multiple_rarefactions.py 对OTU表进行多次等量抽样
#	#	alpha_diversity.py 计算多种方法的alpha多样性
#	#	collate_alpha.py 将多次抽样的alpha结果合并
#	#	make_rarefaction_plots.py 绘制稀释曲线
#	#	normalize_table.py OTU表进行CSS标准化
#	# cutadapt --version # 1.9.1 剪切扩增子双端的引物序列
#	# usearch8 --help # 8.0.1517_i86linux64 扩增子分析软件
#	# usearch10 --help # v10.0.240_i86linux64 扩增子分析软件
#	#	derep_fulllength 序列去冗余
#	#	cluster_otus 聚类OTU
#	#	uchime_ref 去除嵌合体
#	#	usearch_global 生成OTU表
#	#	uc2otutab.py 转换uc至OTU表
#	# biom --version # 2.1.5 biom文件OTU表统计和格式转换
#	# clustalo --version # 1.2.1
#	#
#	# # Self script
#	# stat_16s_lib_split.sh # 1.0 对每个文库拆分样品的结果进行统计和绘图
#	# plot_16s_lib.sh #1.0 对文库处理流程中的数据量和序列长度分布统计和绘图
#	# format_taxonomy2lefse.pl # 转换taxonomy至lefse格式
#	#
#	# # Database
#	# All database list on the makefile.config
#
## 0.1 建立程序必须目录 Create work directory
#init:
#	touch $@
#	mkdir -p ${seq} 
#	mkdir -p ${doc}
#	mkdir -p ${temp}
#	mkdir -p ${result}
#	mkdir -p ${result_f}
#
## 0.2 按doc/library.txt改文库名，方便批量检索 Batch rename seq/*fq file according to doc/design, make samples ID more meaningful
#rename: init
#	rename_fq.pl -i ${library} -s ./${seq}/ -t ./${seq}/ -m mv -h 1
#	dos2unix ${doc}/*
#
#
#
## 1. 并行文库质控拆分 Bacth parallel split each lib
#
## 1.1 质控 Quality control
#${lib}.qc:
#	touch $@
#	fastqc --threads 2 --quiet ${seq}/${lib}_*.fq.gz --extract
#
## 1.2 双端合并 Merge clean reads
#${lib}.merge: ${lib}.qc 
#	touch $@
#	# Method1. qiime script - join_paired_ends.py
#	rm -fr ${temp}/${lib}_join # Delete old output
#	join_paired_ends.py -f ${seq}/${lib}_1.fq.gz -r ${seq}/${lib}_2.fq.gz -m fastq-join -o ${temp}/${lib}_join
#	# Method2. Usearch 10 x64
##	zcat ${seq}/${lib}_1.fq.gz > ${seq}/${lib}_1.fq
##	zcat ${seq}/${lib}_2.fq.gz > ${seq}/${lib}_2.fq
##	mkdir -p temp/${lib}_join
##	usearch10 -fastq_mergepairs ${seq}/${lib}_1.fq -reverse ${seq}/${lib}_2.fq -fastqout ${temp}/${lib}_join/fastqjoin.join.fastq -relabel ${lib}
#	rm -f ${temp}/${lib}.fq # Delete old hardlink
#	ln ${temp}/${lib}_join/fastqjoin.join.fastq ${temp}/${lib}.fq
#	fastqc -quiet ${temp}/${lib}.fq
#
## 1.3 提取barcode Extract barcodes
#${lib}.extract_barcodes: ${lib}.merge
#	touch $@
#	validate_mapping_file.py -o temp/ -m doc/${lib}.txt
#	rm -fr ${temp}/${lib}_barcode
#	extract_barcodes.py -f ${temp}/${lib}_join/fastqjoin.join.fastq \
#		-m doc/${lib}.txt \
#		-o ${temp}/${lib}_barcode \
#		-c ${lib_type} --bc1_len ${bc1} --bc2_len ${bc2} -a --rev_comp_bc2
#
## 1.4 拆分样品 Split library (华大小麦剖面数据报incorrect value for phred_offset，但两种类型33/64均使用也无法通过)
#${lib}.split: ${lib}.extract_barcodes
#	touch $@
#	# Method1. split_libraries_fastq.py
#	rm -fr ${temp}/${lib}_split
#	split_libraries_fastq.py -i ${temp}/${lib}_barcode/reads.fastq \
#		-b ${temp}/${lib}_barcode/barcodes.fastq \
#		-m doc/${lib}.txt \
#		-o ${temp}/${lib}_split/ \
#		-q ${quality} --max_bad_run_length 3 --min_per_read_length_fraction 0.75 \
#		--max_barcode_errors 0 --barcode_type ${bt} --phred_offset=${phred} --sequence_max_n ${N}
#	tail -n+16 ${temp}/${lib}_split/split_library_log.txt|head -n-4>${result}/${lib}_split.count
#	# Method2. Usearch
#	# 准备barcode样品对应文件 samples barcode in fasta
##	tail -n+2 doc/${lib}.txt | cut -f 1-2 | sed 's/^/>/;s/\t/\n/' > ${temp}/${lib}.barcode.fa 
##	usearch10 -fastx_demux ${temp}/${lib}_barcode/reads.fastq -index ${temp}/${lib}_barcode/barcodes.fastq \
##	-barcodes ${temp}/${lib}.barcode.fa -fastqout temp/${lib}_demux.fq >> ${log_usearch} 2>&1 # label samples
##	grep 'sample' temp/${lib}_demux.fq | cut -f 2 -d '=' | cut -f 1 -d ';' | sort | uniq -c | awk '{print $$2"\t"$$1}' > ${result}/${lib}_split.count
#
## 1.5 切除接头引物序列 Cut adaptors
#${lib}.cutadapt: ${lib}.split
#	touch $@
#	rm -f ${temp}/${lib}_P5.fa ${temp}/${lib}_P53.fa
#	cutadapt -g ${primer5} -e ${er} --discard-untrimmed ${temp}/${lib}_split/seqs.fna -o ${temp}/${lib}_P5.fa
#
#${lib}.cutadapt3: ${lib}.cutadapt
#	touch $@
#	cutadapt -a ${primer3} -e ${er} --discard-untrimmed -m ${min_len} ${temp}/${lib}_P5.fa -o ${temp}/${lib}_P53.fa
#
## 1.6 统计每步 Statistics 1-5 each process
#${lib}.stat: ${lib}.cutadapt3
#	touch $@
#	stat_16s_lib_split.sh -o ${result} -A ${g1} -C ${g2} -d ${design} -l ${lib} -m ${merge_group}
#	echo 'Merged clean reads:' > ${lib}.stat
#	grep -c -P '^\+$$' ${temp}/${lib}_join/fastqjoin.join.fastq >> ${lib}.stat
#	echo 'Oriented reads:' >> ${lib}.stat
#	grep -c -P '^\+$$' ${temp}/${lib}_barcode/reads.fastq >> ${lib}.stat
#	echo 'Splitted lib reads:' >> ${lib}.stat
#	grep -c '>' ${temp}/${lib}_split/seqs.fna >> ${lib}.stat
#	echo 'Remove 5` primer:' >> ${lib}.stat
#	grep -c '>' ${temp}/${lib}_P5.fa >> ${lib}.stat
#	echo "Remove 3\` primer and length less than ${min_len} nt:" >> ${lib}.stat
#	grep -c '>' ${temp}/${lib}_P53.fa >> ${lib}.stat
#	grep -v '>' ${temp}/${lib}_P53.fa|awk '{print length($$0)}'|sort -n|uniq -c|sed 's/^ *//g;s/ /\t/g;s/^/${lib}\t/g' > ${temp}/length_${lib}.txt
#	rm -f stat_16s_lib*
#
## 1.7 Test unit for new function
#${lib}.test:
#	ln ${temp}/${lib}_join/fastqjoin.join.fastq ${temp}/${lib}.fq
#	fastqc -quiet ${temp}/${lib}.fq
#
## Batch run 1 - 6 step for each library
#qc merge extract_barcodes split cutadapt cutadapt3 stat test:
#	${foreach var, ${list}, make ${var}.$@ lib=${var} &}
#
#
#test_cmd:
#	touch $@
#	grep 'sample' temp/L10_demux.fq | cut -f 2 -d '=' | cut -f 1 -d ';' | sort | uniq -c | awk '{print $$2"\t"$$1}' > ${result}/L10_split.count
#
#
#
## 2. Make OTUs and calculate diversity
#
## 2.1 总结各文库双端合并前后质量评估 multiple QC
## 启动python3环境，make中无法启动，修改multiqc文件中首行指向python3即可，--pdf输出PDF同时html失去交互能力
#multiqc:
#	touch $@
#	# Summary pair-end fastq
#	multiqc clean_data/ -f -o clean_data/ --pdf
#	# Summary merged fastq
#	multiqc temp/ -f -o temp/ --pdf
#
## 2.2 合并所有样品 Merge all libraries, and format to usearch
#merge_library: 
#	touch $@
#	# Show filter read process: ${result}/bar_qc_sum.pdf/png
#	plot_16s_lib.sh -o ${result} -l "${list}"
#	cat ${temp}/*_P53.fa | sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$$/;/g' > ${temp}/seqs_usearch.fa
#	echo -ne $(date)"\nTotal_reads\t" > ${log_reads}
#	grep -c '>' ${temp}/seqs_usearch.fa >> ${log_reads}
#	cat ${log_reads}
#
## 2.3 序列去冗余 dereplication
#derep: merge_library
#	touch $@
##	usearch8 -derep_fulllength ${temp}/seqs_usearch.fa \
##		-fastaout ${temp}/seqs_unique.fa \
##		-minuniquesize ${minuniquesize} -sizeout >> ${log_usearch} 2>&1
#	usearch10 -fastx_uniques ${temp}/seqs_usearch.fa -sizeout -relabel U \
#		-fastaout ${temp}/seqs_unique.fa -minuniquesize ${minuniquesize} # >> ${log_usearch} 2>&1
#	echo -ne $(date)"\nUnique reads\t" > ${log_otus}
#	grep -c '>' ${temp}/seqs_unique.fa >> ${log_otus}
#	cat ${log_otus}
#
## 2.4 鉴定生物特征序列 unoise
#unoise: derep
#	touch $@
#	usearch10 -unoise3 ${temp}/seqs_unique.fa -zotus ${temp}/Zotus.fa -minsize ${minuniquesize}
#	awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' ${temp}/Zotus.fa > ${temp}/otus.fa
##	usearch8 -cluster_otus ${temp}/seqs_unique.fa \
##		-otus ${temp}/otus.fa \
##		-sizein -sizeout  -id 0.97 #  -uparseout ${temp}/otus.up -id ${sim} >>${log_usearch} 2>&1
#	echo -ne 'Cluster_OTU\t' >> ${log_otus}
#	grep -c '>' ${temp}/otus.fa >> ${log_otus}
#	cat ${log_otus}
#	ln -f ${temp}/otus.fa ${temp}/otus_rdp_align.fa
#
## 2.5 可先是unoise还是cluster_otu进行聚类
## 注意下面条件语句不可以缩进，ifeq else条件必须顶头才能运行，命令必须制表符经缩进
#cluster_otu: derep
#	touch $@
#ifeq (${cluster}, unoise3)
#	echo ${cluster}
#	usearch10 -unoise3 ${temp}/seqs_unique.fa -zotus ${temp}/Zotus.fa -minsize ${minuniquesize}
#	awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' ${temp}/Zotus.fa > ${temp}/otus.fa
#else
#	usearch10 -cluster_otus temp/seqs_unique.fa -otus ${temp}/otus.fa >> ${log_usearch} 2>&1
#endif
#	echo -ne 'Cluster_OTU\t' >> ${log_otus}
#	grep -c '>' ${temp}/otus.fa >> ${log_otus}
#	cat ${log_otus}
#	ln -f ${temp}/otus.fa ${temp}/otus_rdp_align.fa
#
#
## 设置初始值不存在时，变为unoise：无效，仍然还是原来的值
#all: 
#ifneq (${cluster}, cluster_otus)
#	${cluster}=unoise3
#	echo ${cluster}
#endif
#
## 6. (可选)基于参考库去嵌合-容易有假阳性 Remove chimeras by silva database, make false negative, option
#rm_chimeras: unoise
#	touch $@
##	usearch8 -uchime_ref ${temp}/otus.fa  \
##		-nonchimeras ${temp}/otus_rdp.fa \
##		-uchimeout ${temp}/otus_rdp.uchime -db ${rdp} -strand plus >> ${log_usearch} 2>&1
#	usearch10 -uchime2_ref ${temp}/otus.fa -db ${rdp} -chimeras ${temp}/otus_chimeras.fa -strand plus -mode balanced -threads ${p} # sensitive 23.5% chimeras, OTU_1 chimeras; balanced 8% Chimeras, OTU_1 still; high_confidence 1.5% Chimeras, OTU_1 still
#	cat ${temp}/otus.fa ${temp}/otus_chimeras.fa | grep '>' | sed 's/>//g'| sort | uniq -u > ${temp}/otus_nonchimeras.id # output non-chemias ID
#	usearch10 -fastx_getseqs ${temp}/otus.fa -labels ${temp}/otus_nonchimeras.id -fastaout ${temp}/otus_rdp.fa
#	echo -ne 'No_chimeras\t' >> ${log_otus}
#	grep -c '>' ${temp}/otus_rdp.fa >> ${log_otus}
#	cat ${log_otus}
#	align_seqs.py -i ${temp}/otus_rdp.fa -t ${gg_align} -o ${temp}/aligned/
#	grep '>' ${temp}/aligned/otus_rdp_aligned.fasta|cut -f 1 -d ' '|sed 's/>//g' > ${temp}/aligned/otus_rdp_aligned.id
#	filter_fasta.py -f ${temp}/otus_rdp.fa -o ${temp}/otus_rdp_align.fa -s ${temp}/aligned/otus_rdp_aligned.id
#	# fasta_subtraction.pl -i ${temp}/otus_rdp.fa -d ${temp}/aligned/otus_rdp_failures.fasta -o ${temp}/otus_rdp_align.fa
#	echo -ne 'No_bac\t' >> ${log_otus}
#	grep '>' -c ${temp}/otus_rdp_align.fa >> ${log_otus}
#	cat ${log_otus}
#
## 7. Generate representitive sequences and OTU table, remove low abundance samples
#otu_table: cluster_otu
#	touch $@
#	awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' ${temp}/otus_rdp_align.fa > ${result}/rep_seqs.fa
##	usearch8 -usearch_global ${temp}/seqs_usearch.fa -db ${result}/rep_seqs.fa -uc ${temp}/otu_table.uc -strand plus -id ${sim} >> ${log_usearch} 2>&1
##	uc2otutab.py ${temp}/otu_table.uc > ${temp}/otu_table_raw.txt 
#	# 比对时会丢弃一部分OTUs usearch10 mapping will discard some OTUs
##	usearch10 -otutab ${temp}/seqs_usearch.fa -otus ${result}/rep_seqs.fa -otutabout ${temp}/otu_table_raw.txt -threads ${p}
#	vsearch --usearch_global temp/seqs_usearch.fa --db ${result}/rep_seqs.fa --id 0.97 --otutabout ${temp}/otu_table_raw.txt --threads ${p}
#	biom convert -i ${temp}/otu_table_raw.txt -o ${temp}/otu_table_raw.biom --table-type="OTU table" --to-json
#	echo "Summary of otu_table_raw is in ${temp}/otu_table_raw.sum"
#	biom summarize-table -i ${temp}/otu_table_raw.biom >> ${temp}/otu_table_raw.sum
#	# 过滤数据量过小的样本，默认排除reads数小于5k的样品
#	filter_samples_from_otu_table.py -i ${temp}/otu_table_raw.biom -o ${result}/otu_table.biom -n ${thre_count}
#	# 统计过滤样品后的OTU表为终表
#	echo "Summary of otu_table_raw is in ${result}/otu_table.sum"
#	biom summarize-table -i ${result}/otu_table.biom > ${result}/otu_table.sum
#	biom convert -i ${result}/otu_table.biom -o ${result}/otu_table.txt --table-type="OTU table" --to-tsv
#	sed -i '/# Const/d' ${result}/otu_table.txt
#	usearch10 -otutab_stats ${result}/otu_table.txt -output ${result}/otu_table.report # Summary OTUs table
#	usearch10 -otutab_norm ${result}/otu_table.txt -sample_size ${rarefaction} -output ${result}/otu_table_norm.txt # normlize by subsample to 10000
#	# 删除表头，和首行的#及多余字符，方便R读取
##	sed -i 's/#OTU //g' ${result}/otu_table.txt
#
## 8. Taxonomy assignment
#assign_tax: otu_table
#	touch $@
#	assign_taxonomy.py -i ${result}/rep_seqs.fa -r ${gg_seq} -t ${gg_tax} -m ${method} -o result --rdp_max_memory 900000
#	# 每个OTU均编号，对结果按分类级合并影响很大，全改为UnAssigned
##	format_taxonomy2full.pl -i ${result}/rep_seqs_tax_assignments.txt -o ${result}/rep_seqs_tax_assignments.txt.full
#	# 全改为Unassigned
#	format_taxonomy2full_unassigned.pl -i ${result}/rep_seqs_tax_assignments.txt -o ${result}/rep_seqs_tax_assignments.txt.full
##	awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";;a["s"]="Unassigned";\
##	  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
##	  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"];}' \
##	  ${result}/rep_seqs_tax_assignments.txt > ${result}/rep_seqs_tax_assignments.txt.full
#	# debug rdp
#	sed -i 's/p__Bacteria/p__Firmicutes/g' ${result}/rep_seqs_tax_assignments.txt.full 
#	sed 's/;/\t/g;s/ //g' ${result}/rep_seqs_tax_assignments.txt.full > ${result}/rep_seqs_tax.txt # format for R read
##	mv ${result}/rep_seqs_tax_assignments.log ${temp}/rep_seqs_tax_assignments.log
#	# biom添加的taxonomy注释使用最新修正的taxonomy(.full)前后保持一致
#	biom add-metadata -i ${result}/otu_table.biom --observation-metadata-fp ${result}/rep_seqs_tax_assignments.txt.full -o ${result}/otu_table_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy # add taxonomy to biom
#	biom convert -i ${result}/otu_table_tax.biom -o ${result}/otu_table_tax.txt --to-tsv --header-key taxonomy
#	summarize_taxa.py -i ${result}/otu_table_tax.biom -o ${result}/sum_taxa # summary each level percentage
#	rm ${result}/sum_taxa/*.biom
#	sed -i '/# Const/d;s/#OTU //g' ${result}/sum_taxa/* # format for R read
#
## 9. Phylogeny tree
#tree: assign_tax
#	touch $@
#	clustalo -i ${result}/rep_seqs.fa -o ${temp}/rep_seqs_align.fa --seqtype=DNA --full --force --threads=${p}
#	filter_alignment.py -i ${temp}/rep_seqs_align.fa -o ${temp}/  # rep_seqs_align_pfiltered.fa, only very short conserved region saved
#	make_phylogeny.py -i ${temp}/rep_seqs_align_pfiltered.fasta -o ${result}/rep_seqs.tree # generate tree by FastTree
#
## 10. Alpha diversity
#alpha: tree
#	touch $@
#	# rarefaction=`head -n 7 ${result}/otu_table.sum|tail -n 1|cut -f 3 -d ' '|cut -f 1 -d '.'`
#	single_rarefaction.py -i ${result}/otu_table.biom -o ${temp}/otu_table_rare.biom -d ${rarefaction}
#	biom convert -i ${temp}/otu_table_rare.biom -o ${temp}/otu_table_rare.txt --table-type="OTU table" --to-tsv
#	sed -i '/# Const/d;s/#OTU //g' ${temp}/otu_table_rare.txt
#	alpha_diversity.py -i ${temp}/otu_table_rare.biom -o ${result}/alpha.txt -t ${result}/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree
#	# Calculate all alpha diversity, details in http://www.drive5.com/usearch/manual/alpha_metrics.html
#	usearch10 -alpha_div ${result}/otu_table_norm.txt -output ${result}/alpha.usearch 
#	# Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
#	usearch10 -alpha_div_rare ${result}/otu_table_norm.txt -output ${result}/alpha_rare.txt  -method without_replacement # 取1%-100%的序列中OTUs数量
#
## 10.1 Alpha diversity rarefactions
#alpha_qiime_rare: alpha
#	touch $@
#	#multiple_rarefactions for rarefraction curve
#	multiple_rarefactions.py -i ${result}/otu_table.biom -m 2000 -x 50000 -s 2000 -n 25 -o ${result}/a_rare/ # min, max, step
#	alpha_diversity.py -i ${result}/a_rare/ -m shannon,chao1,observed_otus,PD_whole_tree -o ${result}/a_div/ -t ${result}/rep_seqs.tree 
#	collate_alpha.py -i ${result}/a_div/ -o ${result}/a_collated/
#	#make_rarefaction_plots.py -i ${result}/a_collated/ -m doc/${library}.txt -o ${result}/a_rare_plots # error, but good in virtualbox and docker
#	cat <(echo -n "#") doc/design.txt > doc/design_rare.txt # make mapping file
#	docker run --rm -v `pwd`:/home --name=qiime yoshikiv/basespace-qiime-191-dev make_rarefaction_plots.py -i home/${result}/a_collated/ -m home/doc/design_rare.txt -o home/result
#
## 11. Beta diversity
#beta: alpha
#	touch $@
#	normalize_table.py -i ${result}/otu_table.biom -o ${temp}/otu_table_css.biom -a CSS
#	biom convert -i ${temp}/otu_table_css.biom -o ${result}/otu_table_css.txt --table-type="OTU table" --to-tsv
#	sed -i '/# Const/d;s/#OTU //g' ${result}/otu_table_css.txt
#	beta_diversity.py -i ${temp}/otu_table_css.biom -o ${result}/beta/ -t ${result}/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
#	sed -i 's/^\t//g' ${result}/beta/*
#
## 12. Taxonomy tree - GraPhlAn
#graphlan: beta
#	touch $@
#	filter_otus_from_otu_table.py --min_count_fraction ${tax_per} -i ${result}/otu_table.biom -o ${temp}/tax_otu_table.biom
#	filter_fasta.py -f ${result}/rep_seqs.fa -o ${temp}/tax_rep_seqs.fa -b ${temp}/tax_otu_table.biom 
#	echo "Number of OTU abundance > ${tax_per} :" >> ${log_otus}
#	grep -c '>' ${temp}/tax_rep_seqs.fa >> ${log_otus}
#	grep '>' ${temp}/tax_rep_seqs.fa|sed 's/>//g' > ${temp}/tax_rep_seqs.id
#	awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$$1]=$$0} NR>FNR {print a[$$1]}' ${result}/rep_seqs_tax_assignments.txt.full ${temp}/tax_rep_seqs.id|cut -f 2-3|grep 's__'|sed 's/; */\|/g' > ${temp}/tax_full_anno.txt 
#	echo "Number of OTU abundance > ${tax_per} with fully annotation :" >> ${log_otus}
#	wc -l ${temp}/tax_full_anno.txt >> ${log_otus}
#	echo "Number of OTU abundance > ${tax_per} with fully annotation unique:" >> ${log_otus}
#	sort ${temp}/tax_full_anno.txt|cut -f 1|uniq|wc -l >> ${log_otus}
#	format_taxonomy2lefse.pl -i ${temp}/tax_full_anno.txt -o ${temp}/tax_lefse.txt 
#	## order
#	export2graphlan.py -i ${temp}/tax_lefse.txt --tree ${temp}/tax_order.tree --annotation ${temp}/tax_order.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 4 --min_clade_size 1 --min_font_size 5
#	graphlan_annotate.py --annot ${temp}/tax_order.annot ${temp}/tax_order.tree ${temp}/tax_order.xml
#	sed -i 's/ref="A:1">o  /ref="A:1">/g' ${temp}/tax_order.xml
#	graphlan.py --dpi 300 ${temp}/tax_order.xml ${result}/tax_order.pdf --external_legends
#	graphlan.py --dpi 300 ${temp}/tax_order.xml ${result}/tax_order.png --external_legends
#	mv ${result}/tax_order_legend.* ${temp}/ 
#	## family
#	export2graphlan.py -i ${temp}/tax_lefse.txt --tree ${temp}/tax_family.tree --annotation ${temp}/tax_family.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 5 --min_clade_size 1 --min_font_size 4
#	graphlan_annotate.py --annot ${temp}/tax_family.annot ${temp}/tax_family.tree ${temp}/tax_family.xml
#	sed -i 's/ref="A:1">f  /ref="A:1">/g' ${temp}/tax_family.xml
#	graphlan.py --dpi 300 ${temp}/tax_family.xml ${result}/tax_family.pdf --external_legends
#	graphlan.py --dpi 300 ${temp}/tax_family.xml ${result}/tax_family.png --external_legends
#	mv ${result}/tax_family_legend.* ${temp}/ 
#	## genus
#	export2graphlan.py -i ${temp}/tax_lefse.txt --tree ${temp}/tax_genus.tree --annotation ${temp}/tax_genus.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 6 --min_clade_size 1 --min_font_size 3
#	graphlan_annotate.py --annot ${temp}/tax_genus.annot ${temp}/tax_genus.tree ${temp}/tax_genus.xml
#	sed -i 's/ref="A:1">g  /ref="A:1">/g' ${temp}/tax_genus.xml
#	graphlan.py --dpi 300 ${temp}/tax_genus.xml ${result}/tax_genus.pdf --external_legends
#	graphlan.py --dpi 300 ${temp}/tax_genus.xml ${result}/tax_genus.png --external_legends
#	mv ${result}/tax_genus_legend.* ${temp}/ 
#
## 13. Phylogenetic tree - ggtree
#ggtree: graphlan
#	touch $@
#	clustalo -i ${temp}/tax_rep_seqs.fa -o ${temp}/tax_rep_seqs_clus.fa --seqtype=DNA --full --force --threads=$p
#	make_phylogeny.py -i ${temp}/tax_rep_seqs_clus.fa -o ${temp}/tax_rep_seqs.tree
#	sed "s/'//g" ${temp}/tax_rep_seqs.tree > ${result}/tax_rep_seqs.tree # remove '
#	grep '>' ${temp}/tax_rep_seqs_clus.fa|sed 's/>//g' > ${temp}/tax_rep_seqs_clus.id
#	awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$$1]=$$0} NR>FNR {print a[$$1]}' ${result}/rep_seqs_tax_assignments.txt.full ${temp}/tax_rep_seqs_clus.id|sed 's/;/\t/g'|cut -f 1-5 |sed 's/p__/p_/g;s/c__/c_/g;s/o__/o_/g' > ${result}/tax_rep_seqs.tax
#	ggtree.sh -e TRUE -b ${pvalue} -d ${design} -m ${merge_group} -g ${group_order} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${result} -h ${height} -w ${width} -s ${text_size} -S percentage
#
#
## 14. Visuallize diversity, draw alpha, beta and Constrain PCoA
#diversity: ggtree
#	touch $@
#	diversity.sh -d ${design} -m ${merge_group} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -E ${g3} -F ${g3_list} -o ${result} -g ${group_order} -h ${height} -w ${width} -s ${text_size} -T ${batch} -L ${ellipse}
#	alpha_rare_usearch.sh -d ${design} -m ${merge_group} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${result} -g ${group_order} -h ${height} -w ${width} -s ${text_size}
#
#
### 15. Visuallize taxonomy, draw barplot+error bar, stack plot, first using qimme + limma; need update to count and edgeR
##taxonomy: diversity
##	touch $@
##	taxonomy_egr.sh -d ${design} -m ${merge_group} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${result} -g ${group_order} -n ${tax_number} -h ${height} -w ${width} -s ${text_size}
##
### 16. Visuallize DEOTU, draw volcano, manhattan, heatmap, venn
##DAOTU: taxonomy
##	touch $@
##	DAOTU_egr.sh -d ${design} -m ${merge_group} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${result}
#
#
#
### Filter OTU by abundance and reanalyze, such as 0.001%
## 1. Subset OTU table by abundance and taxonomy
#filter:
#	touch $@
#	rm -r ${result_f}
#	mkdir -p ${result_f}
#	# 按丰度筛选：建议添加排除宿主序列 genome+chloroplast+mitochondria、指定菌门，再进行丰度筛选
#	filter_otus_from_otu_table.sh -t ${thre} -o ${result} -d ${design} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} 
#	filter_otus_from_otu_table.py -i ${result}/otu_table_tax.biom -o ${temp}/k1.biom --otu_ids_to_exclude_fp ${result}/otu_id_k1.txt --negate_ids_to_exclude
#	echo 'Summary of otu_table_k1, one of sample OTU > 0.1%:' >> ${log_otus}
#	biom summarize-table -i ${temp}/k1.biom  > ${temp}/k1.biom.sum
#	echo -ne "${thre}\t" >> ${log_otus}
#	cat <(head -n2 ${temp}/k1.biom.sum|tail -n1|cut -f 3 -d ' ') >> ${log_otus}
#	cat ${log_otus}
##	# 去除宿主 remove host ITS OTU
##	blastn -query ${result}/rep_seqs.fa -db ${host_its} -out ${temp}/rep_seqs.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 # 输出13列为coverage
##	awk '$$3>97 && $$13>90' ${temp}/rep_seqs.blastn|cut -f 1 > ${result}/otu_id_host.txt
##	filter_otus_from_otu_table.py -i ${temp}/k1.biom -o ${temp}/k1_host.biom --otu_ids_to_exclude_fp ${result}/otu_id_host.txt
#	# -p是保留某些菌，-n是去除某类菌
#	filter_taxa_from_otu_table.py -i ${temp}/k1.biom -o ${result_f}/otu_table.biom -n ${taxonomy}
#	biom summarize-table -i ${result_f}/otu_table.biom > ${result_f}/otu_table.sum
#	# 统计去除指定菌门后的OTU数量
#	echo -ne 'Remove'${taxonomy}"\t" >> ${log_otus}
#	cat <(head -n2 ${result_f}/otu_table.sum|tail -n1|cut -f 3 -d ' ') >> ${log_otus}
#	cat ${log_otus}
#	filter_fasta.py -f ${result}/rep_seqs.fa -o ${result_f}/rep_seqs.fa -b ${result_f}/otu_table.biom
#	ln -f ${result_f}/otu_table.biom ${result_f}/otu_table_tax.biom
#	summarize_taxa.py -i ${result_f}/otu_table_tax.biom -o ${result_f}/sum_taxa
#	biom convert -i ${result_f}/otu_table_tax.biom -o ${result_f}/otu_table_tax.txt --to-tsv --header-key taxonomy
#	sed -i '/# Const/d;s/#OTU //g' ${result_f}/otu_table_tax.txt
#	rm ${result_f}/sum_taxa/*.biom
#	sed -i '/# Const/d;s/#OTU //g' ${result_f}/sum_taxa/*
#	biom convert -i ${result_f}/otu_table.biom -o ${result_f}/otu_table.txt --table-type="OTU table" --to-tsv
#	sed -i '/# Const/d;s/#OTU //' ${result_f}/otu_table.txt 
#	cut -f 1 ${result_f}/otu_table.txt | tail -n+2 > ${temp}/k1_t.id
#	awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$$1]=$$0} NR>FNR {print a[$$1]}' ${result}/rep_seqs_tax.txt ${temp}/k1_t.id > ${result_f}/rep_seqs_tax.txt
#
## 2. Re-analyze new OTU table
#rediv: filter
#	touch $@
#	clustalo -i ${result_f}/rep_seqs.fa -o ${temp}/rep_seqs_align.fa --seqtype=DNA --full --force --threads=${p}
#	filter_alignment.py -i ${temp}/rep_seqs_align.fa -o ${temp}/
#	make_phylogeny.py -i ${temp}/rep_seqs_align_pfiltered.fasta -o ${result_f}/rep_seqs.tree
#	single_rarefaction.py -i ${result_f}/otu_table.biom -o ${temp}/otu_table_rare.biom -d ${rarefaction}
#	# 小宁用Public安装了MPI软件导致alpha_diversity依赖关系变化而报错，卸载后恢复
#	alpha_diversity.py -i ${temp}/otu_table_rare.biom -o ${result_f}/alpha.txt -t ${result_f}/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree
#	#改用docker代替，但报错缺少hdf5无法处理biom格式
#	#docker run --rm -v `pwd`:/home --name=qiime yoshikiv/basespace-qiime-191-dev alpha_diversity.py -i home/${temp}/otu_table_rare.biom -o home/${result_f}/alpha.txt -t home/${result_f}/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree
#	normalize_table.py -i ${result_f}/otu_table.biom -o ${temp}/otu_table_css.biom -a CSS
#	biom convert -i ${temp}/otu_table_css.biom -o ${result_f}/otu_table_css.txt --table-type="OTU table" --to-tsv
#	sed -i '/# Const/d;s/#OTU //g' ${result_f}/otu_table_css.txt
#	beta_diversity.py -i ${temp}/otu_table_css.biom -o ${result_f}/beta/ -t ${result_f}/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
#	sed -i 's/^\t//g' ${result_f}/beta/*
#
## 3. redraw all figure
#
## Draw diversity related figures, include alpha boxplot and rarefraction, PCA, PCoA, LDA, CCA
#draw_div: rediv
#	touch $@
#	diversity.sh -d ${design} -m ${merge_group} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -E ${g3} -F ${g3_list} -o ${result_f} -g ${group_order} -h ${height} -w ${width} -s ${text_size} -T ${batch} -L ${ellipse}
#
## Draw taxonomy related figure
#draw_tax: draw_div
#	touch $@
#	taxonomy_egr.sh -b ${pvalue} -d ${design} -m ${merge_group} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${result_f} -g ${group_order} -n ${tax_number} -h ${height} -w ${width} -s ${text_size} -F ${adjust_method} -L ${logFC}
#	taxonomy_phylumpro.sh -b ${pvalue} -d ${design} -m ${merge_group} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${result_f} -g ${group_order} -n ${tax_number} -h ${height} -w ${width} -s ${text_size} -F ${adjust_method} -L ${logFC}
#	plot_pie_DA_Bphylum.sh -c ${compare} -l family -o ${result_f}
#	batch_venn.pl -i ${venn} -d ${result_f}/order.txt
#	rm -f ${result_f}/order.txt*venn*.xls.xls
#	batch2.pl -i "${result_f}/order.txt*venn*.xls" -d ${result_f}/database_order.txt -o ${result_f}/ -p vennNumAnno.pl
#	batch_venn.pl -i ${venn} -d ${result_f}/family.txt
#	rm -f ${result_f}/family.txt*venn*.xls.xls
#	batch2.pl -i "${result_f}/family.txt*venn*.xls" -d ${result_f}/database_family.txt -o ${result_f}/ -p vennNumAnno.pl
#
#draw_otu: draw_tax
#	touch $@
#	DAOTU_egr.sh -b ${pvalue} -d ${design} -m ${merge_group} -g ${group_order} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${result_f} -h ${height} -w ${width} -s ${text_size} -F ${fdr} -L ${logFC} -M ${compare_method}
#
#draw_venn: draw_otu
#	touch $@
#	plot_pie_DA_Bphylum.sh -c ${compare} -l otu -o ${result_f}
#	batch_venn.pl -i ${venn} -d ${result_f}/otu.txt
#	rm -f ${result_f}/otu.txt*venn*.xls.xls
#	batch2.pl -i '${result_f}/otu.txt*venn*.xls' -d ${result_f}/database.txt -o ${result_f}/ -p vennNumAnno.pl
#
#draw_ter: draw_venn
#	touch $@
#	ternary_16s.sh -d ${design} -m ${merge_group} -t ${tern} -A ${g1} -C ${g2} -D ${g2_list} -o ${result_f} 
#
### 4. Write RMD report of all result
#rmd: draw_ter
#	time rmd_16s.pl -s ${summary} -d ${design} -l ${library} -c ${compare} -v ${venn} -t ${tern} -g ${g1} -D ${g2_list} -F ${g3_list} -b ${version} -S ${elite_report} -a ${thre}
#	ln ${wd}/${version}/ /var/www/html/report/16s/${version} -sf
#	rm -f ${version}/${version}
#
## 补充实验阶段功能
#
### 选择不同列计算PCoA, compartment, genotype, site, day
#draw_div_pcoa: rediv
##	touch $@
#	diversity_pcoa.sh -d ${design} -m ${merge_group} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -E ${g3} -F ${g3_list} -o ${result_f} -g ${group_order} -h ${height} -w ${width} -s ${text_size} -T ${batch} -L ${ellipse} -G ${group_color} -M ${time_course}
#
### 单独添加phylum+pro冲击图
#draw_tax_pro: rediv
#	taxonomy_phylumpro.sh -b ${pvalue} -d ${design} -m ${merge_group} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${result_f} -g ${group_order} -n ${tax_number} -h ${height} -w ${width} -s ${text_size} -F ${fdr} -L ${logFC}
#
### 高级定制 High-Level
#
## 3.1 可培养比例分析 Culture
#culture_graphlan: 
#	# 建立完整OTU比对数据库，用于OTU可培养注释，报告会根据是否有此数来添加差异OTU的注释
#	blastn -query result/rep_seqs.fa -db ${cluture_db} -out result/rep_seqs.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 # 输出13列为coverage
#	awk '$$3*$$13>=9700' result/rep_seqs.blastn|cut -f 1-3 > result/rep_seqs.97
#	OTU_culture_anno.pl -i result/rep_seqs.fa -d result/rep_seqs.97 -o result/otu_cultured.txt
#	# 筛选根际土、根的k1 OTU,并在相应库中匹配培养比例；
#	mkdir -p ${filter}
#	echo -ne "Total OTUs:\t" > ${filter}/culture.txt
#	grep '>' -c result/rep_seqs.fa >> ${filter}/culture.txt
#	echo -ne "OTUs similar to stocks:\t" >> ${filter}/culture.txt
#	wc -l result/rep_seqs.blastn >> ${filter}/culture.txt
#	echo -ne "OTUs similar > 97% in stocks:\t" >> ${filter}/culture.txt
#	wc -l result/rep_seqs.97 >> ${filter}/culture.txt
#	filter_otus_from_otu_table.sh -t ${thre2} -o ${filter} -f ${otu_table} -d ${design} -F 'TRUE' -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list}
#	filter_fasta.py -f ${result}/rep_seqs.fa -o ${filter}/rep_seqs.fa.top -s ${filter}/otu_table_ha.id
#	echo -ne "Nature > ${thre2} OTUs:\t" >> ${filter}/culture.txt
#	grep -c '>' ${filter}/rep_seqs.fa.top >> ${filter}/culture.txt
#	# 分析这些OTU中可培养的比例
#	blastn -query ${filter}/rep_seqs.fa.top -db ${cluture_db} -out ${filter}/rep_seqs.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 # 输出13列为coverage
#	awk '$$3*$$13>=9700' ${filter}/rep_seqs.blastn|cut -f 1 > ${filter}/otu_cultured.txt
#	echo -ne "Stocked_OTUs:\t" >> ${filter}/culture.txt
#	grep -c 'OTU' ${filter}/otu_cultured.txt >> ${filter}/culture.txt
#	echo -ne "Nature_HA_abundance:\t" >> ${filter}/culture.txt
#	awk '{a=a+$$2} END {print a}' ${filter}/otu_table_ha.mean >> ${filter}/culture.txt # total is 0.835
#	echo -ne "Stocked_abundance:\t" >> ${filter}/culture.txt
#	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$$1]="culture"} NR>FNR {print $$0,a[$$1]}' ${filter}/otu_cultured.txt ${filter}/otu_table_ha.mean |grep 'culture'|awk '{a=a+$$2} END {print a}' >> ${filter}/culture.txt 
#	# 绘制graphlan
#	graphlan_culture.pl -i ${filter}/otu_table_ha.id -d ${filter}/otu_cultured.txt -t result/rep_seqs_tax_assignments.txt.full -o 0_ha_otu_culture.txt
#	Rscript /mnt/bai/yongxin/bin/graphlan_culture.R # 生成1树, 2科注释, 3培养注释文件
#	sed 's/\t/\tring_alpha\t3\t/g' ${filter}/otu_table_ha.zscore > ${filter}/abundance_heat.txt # 柱状用log2，热图用zscore
#	cat /mnt/bai/yongxin/culture/rice/graphlan/global.cfg 2_annotation_family.txt /mnt/bai/yongxin/culture/rice/graphlan/ring1.cfg 3_annotation_match.txt /mnt/bai/yongxin/culture/rice/graphlan/abundance_heat.cfg ${filter}/abundance_heat.txt > ${filter}/5_annotation.txt
#	graphlan_annotate.py --annot ${filter}/5_annotation.txt 1_tree_plain.txt ${filter}/graphlan.xml
#	graphlan.py ${filter}/graphlan.xml ${filter}/graphlan.pdf --size 5
#	graphlan.py ${filter}/graphlan.xml ${filter}/graphlan.png --size 5
#	cat ${filter}/culture.txt
#
#
## 3.2以GG为参考的OTU表预测宏基因组KEGG基因分类
#picrust_calc: 
#	touch $@
#	mkdir -p function
#	# 以gg13.5比对
#	#usearch8 -usearch_global temp/seqs_usearch.fa -db /mnt/bai/public/ref/gg_13_5_otus/rep_set/97_otus.fasta -uc function/gg135_otu_table.uc -strand plus -id 0.97 -threads 96
#	vsearch --usearch_global temp/seqs_usearch.fa --db /mnt/bai/public/ref/gg_13_5_otus/rep_set/97_otus.fasta --id 0.97 --otutabout function/gg135_otu_table.txt --threads ${p}
#	#uc2otutab.py function/gg135_otu_table.uc > function/gg135_otu_table.txt 
#	biom convert -i function/gg135_otu_table.txt -o function/gg135_otu_table.biom --table-type="OTU table" --to-json
#	biom summarize-table -i function/gg135_otu_table.biom > function/gg135_otu_table.sum
#	# 校正拷贝数
#	normalize_by_copy_number.py -i function/gg135_otu_table.biom -o function/normalized_otus.biom -c ~/test/example_PE250/picrust/picrust/data/16S_13_5_precalculated.tab.gz
#	# 预测宏基因组KO表
#	predict_metagenomes.py -i function/normalized_otus.biom -o function/metagenome_predictions.biom -c ~/test/example_PE250/picrust/picrust/data/ko_13_5_precalculated.tab.gz
#	# 转换结果为txt查看
#	biom convert -i function/metagenome_predictions.biom -o function/metagenome_predictions.txt --table-type="OTU table" --to-tsv
#	sed -i '/# Const/d;s/#OTU //g' function/metagenome_predictions.txt 	# 调整为标准表格
#	#less -S function/metagenome_predictions.txt # 预览KO表，和OTU表类似，只是OTU变为了KO
#	#wc -l function/metagenome_predictions.txt # 统计KO数据6909个
#	# 按功能级别分类汇总, -c指输出类型，有KEGG_Pathways, COG_Category, RFAM三种，-l是级别，分4级
#	categorize_by_function.py -i function/metagenome_predictions.biom -c KEGG_Pathways -l 3 -o function/metagenome_predictions.L3.biom
#	# 转换结果为txt查看
#	biom convert -i function/metagenome_predictions.L3.biom -o function/metagenome_predictions.L3.txt --table-type="OTU table" --to-tsv
#	sed -i '/# Const/d;s/#OTU //g' function/metagenome_predictions.L3.txt # 调整为标准表格
#
#picrust_draw: picrust_calc
#	#touch $@
#	# 差异分析与统计绘图：热图
#	DAKO_egr.sh -b ${pvalue} -d ${design} -m ${merge_group} -g ${group_order} -c ${compare} -p ${pair_compare} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o function/ -h ${height} -w ${width} -s ${text_size} -a ${adjust_method} -F ${lfc} -f ${picrust}
#
#
#
#
## 新功能测试区
#
#
## OTU表过滤
## makefile命令不能与目录同名
#filter_otu:
#	# 0. 设置输出文件位置，删除旧文件夹并新建
#	rm -rf ${of}
#	mkdir -p ${of}
#	# 1. 过滤宿主线粒体和叶绿体
#	# a. 直接过滤chloroplast+mitochondri名字的OTU
#	# -p是保留某些菌，-n是去除某类菌，通常基于Greengene13.5为c__Chloroplast,f__mitochondria
#	filter_taxa_from_otu_table.py -i result/otu_table_tax.biom -o temp/otu_table_tax_rmMt.biom -n ${taxonomy}
#	biom summarize-table -i temp/otu_table_tax_rmMt.biom > temp/otu_table_tax_rmMt.sum
#	head -n20 temp/otu_table_tax_rmMt.sum
#	biom convert -i temp/otu_table_tax_rmMt.biom -o ${of}/1otu_table_tax_rmMt.txt --to-tsv --header-key taxonomy
#	sed -i '/# Const/d' ${of}/1otu_table_tax_rmMt.txt
#	# b. 与宿主基因组比对过滤全部
#	# 2. 分级别统计
#	summarize_taxa.py -i temp/otu_table_tax_rmMt.biom -o ${of}/sum_taxa
#	rm ${of}/sum_taxa/*.biom
#	sed -i '/# Const/d' ${of}/sum_taxa/*
#
#filter_otu1:
#	# 3. 按样品和组筛选与合并
#	# 按样品名筛选，可用QIIME脚本
##	filter_otus_from_otu_table.py -i ${result}/otu_table_tax.biom -o ${temp}/k1.biom --otu_ids_to_exclude_fp ${result}/otu_id_k1.txt --negate_ids_to_exclude
#	# 按实验设计列1/2列筛选
#	filter_otus_by_design.sh -f ${of}/1otu_table_tax_rmMt.txt -d ${design} -A ${g1} -B ${g1_list} -C ${g2} -D ${g2_list} -o ${of}/2otu_table_filter.txt
#
#
#
## Usearch结合RDP种水平注释 usearch sintax + rdp 16s v16 sp database: 
#usearch_tax:
#	mkdir -p sintax
#	# rdp_16s_v16_sp.fa结果会报错：WARNING: 2 taxonomy nodes have >1 parent
#	usearch10 -sintax result/rep_seqs.fa -db /mnt/bai/public/ref/rdp/rdp_16s_v16_sp.fa -tabbedout sintax/sintax.txt -strand both -sintax_cutoff 0.8
#	# Line 1532 has 3 fields (need 4),连界都没过阈值，缺少4列，补全为Id:Unassigned
#	sed -i 's/\t$/\td:Unassigned/' sintax/sintax.txt
#	for i in p c o f g;do
#	usearch10 -sintax_summary sintax/sintax.txt -rank ${i} \
#	-otutabin otu_filter/5otu_table_sort_LN.txt \
#	-output sintax/sum_${i}.txt
#	done
#	sed -i 's/(//g;s/)//g;s/\"//g' sintax/sum_*.txt
#
