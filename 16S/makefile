SHELL:=/bin/bash
include makefile.config

## Bacth parallel split each library
# 1. Quality control
$(library).qc:
	touch $@
	fastqc -t 9 clean_data/$(library)_1.fq.gz
	fastqc -t 9 clean_data/$(library)_2.fq.gz
	unzip -o clean_data/$(library)_1_fastqc.zip -d clean_data/
	unzip -o clean_data/$(library)_2_fastqc.zip -d clean_data/
	# -o is overwrite exist file, -d is output directory

# 2. Merge clean reads
$(library).merge: $(library).qc 
	touch $@
	join_paired_ends.py -f clean_data/$(library)_1.fq.gz -r clean_data/$(library)_2.fq.gz -m fastq-join -o temp/$(library)_join

# 3. Split library
$(library).split: $(library).merge
	touch $@
	extract_barcodes.py -f temp/$(library)_join/fastqjoin.join.fastq \
		-m doc/$(library).mappingfile.txt \
		-o temp/$(library)_barcode \
		-c barcode_paired_stitched --bc1_len $(bc1) --bc2_len $(bc2) -a --rev_comp_bc2
	split_libraries_fastq.py -i temp/$(library)_barcode/reads.fastq \
		-b temp/$(library)_barcode/barcodes.fastq \
		-m doc/$(library).mappingfile.txt \
		-o temp/$(library)_split/ \
		-q $(quality) --max_bad_run_length 3 --min_per_read_length_fraction 0.75 --max_barcode_errors 0 --barcode_type $(bt)
	#cat temp/$(library)_split/split_library_log.txt >> $(library).stat

# 3.1 Split library stat
# Check samples reads is good? result/bar_split.pdf/png
$(library).split.stat: $(library).split
	touch $@
	mkdir -p result
	tail -n+16 temp/$(library)_split/split_library_log.txt|head -n-4>result/$(library)_split.count
	stat_16s_lib_split.sh -o $(result) -A $(g1) -C $(g2) -d $(design) -l $(library)
	
# 4. Remove adaptor
$(library).cutp: $(library).split.stat
	touch $@
	cutadapt -g $(primer5) -e 0.15 --discard-untrimmed temp/$(library)_split/seqs.fna -o temp/$(library)_P5.fa
	cutadapt -a $(primer3) -e 0.15 --discard-untrimmed -m $(min_len) temp/$(library)_P5.fa -o temp/$(library)_P53.fa

# 4.1 Statistics 1-4 each process
$(library).stat: $(library).cutp
	touch $@
	echo 'Merged clean reads:' > $(library).stat
	grep -c -P '^\+$$' temp/$(library)_join/fastqjoin.join.fastq >> $(library).stat
	echo 'Oriented reads:' >> $(library).stat
	grep -c -P '^\+$$' temp/$(library)_barcode/reads.fastq >> $(library).stat
	echo 'Splitted library reads:' >> $(library).stat
	grep -c '>' temp/$(library)_split/seqs.fna >> $(library).stat
	echo 'Remove 5` primer:' >> $(library).stat
	grep -c '>' temp/$(library)_P5.fa >> $(library).stat
	echo "Remove 3\` primer and length less than $(min_len) nt:" >> $(library).stat
	grep -c '>' temp/$(library)_P53.fa >> $(library).stat
	grep -v '>' temp/$(library)_P53.fa|awk '{print length($$0)}'|sort -n|uniq -c|sed 's/^ *//g;s/ /\t/g;s/^/$(library)\t/g' > temp/length_$(library).txt

# Batch run 1, 2, 3, 4
qc merge split split.stat cutp stat:
	$(foreach var, $(list), make $(var).$@ library=$(var) &)


## Main pipeline: Merge library, cluster otu, ......
# Statistics library
# Show filter read process: result/bar_qc_sum.pdf/png
stat_lib: init
	touch $@
	/mnt/bai/yongxin/bin/stat_16s_lib.sh
	/mnt/bai/yongxin/bin/plot_16s_lib.sh -o $(result)

# Merge all libraries, and format to usearch
merge_library: stat_lib
	touch $@
	cat temp/*_P53.fa | sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$$/;/g' > temp/seqs_usearch.fa
	echo -e $(date)"\nFinished splitting all libraries.\nTotal reads of merge each library :" > $(log)
	grep -c '>' temp/seqs_usearch.fa >> $(log)


# 5. Cluster OTU by Usearch
cluster_otus: merge_library
	touch $@
	usearch8 -derep_fulllength temp/seqs_usearch.fa \
		-fastaout temp/seqs_unique.fa \
		-minuniquesize $(minuniquesize) -sizeout >> $(log) 2>&1
	echo 'Unique reads:' >> $(log)
	grep -c '>' temp/seqs_unique.fa >> $(log)
	usearch8 -cluster_otus temp/seqs_unique.fa \
		-otus temp/otus.fa \
		-uparseout temp/otus.up -sizein -sizeout >> $(log) 2>&1
	echo 'Cluster OTU:' >> $(log)
	grep -c '>' temp/otus.fa >> $(log)

# 6. Remove chimeras by rdp_gold database
rm_chimeras: cluster_otus
	touch $@
	usearch8 -uchime_ref temp/otus.fa  \
		-nonchimeras temp/otus_rdp.fa \
		-uchimeout temp/otus_rdp.uchime -db $(rdp) -strand plus >> $(log) 2>&1
	echo 'Remove chimeras by rdp_gold database:' >> $(log)
	grep -c '>' temp/otus_rdp.fa >> $(log)
	align_seqs.py -i temp/otus_rdp.fa -t $(gg_align) -o temp/aligned/
	grep '>' temp/aligned/otus_rdp_aligned.fasta|cut -f 1 -d ' '|sed 's/>//g' > temp/aligned/otus_rdp_aligned.id
	filter_fasta.py -f temp/otus_rdp.fa -o temp/otus_rdp_align.fa -s temp/aligned/otus_rdp_aligned.id
	# fasta_subtraction.pl -i temp/otus_rdp.fa -d temp/aligned/otus_rdp_failures.fasta -o temp/otus_rdp_align.fa
	echo 'Remove non-bac seq by align_seqs.py:' >> $(log)
	grep '>' -c temp/otus_rdp_align.fa >> $(log)

# 7. Generate representitive sequences and OTU table, remove low abundance samples
otu_table: rm_chimeras
	touch $@
	awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' temp/otus_rdp_align.fa > result/rep_seqs.fa
	usearch8 -usearch_global temp/seqs_usearch.fa -db result/rep_seqs.fa -uc temp/otu_table.uc -strand plus -id $(sim) >> $(log) 2>&1
	uc2otutab.py temp/otu_table.uc > temp/otu_table_raw.txt 
	biom convert -i temp/otu_table_raw.txt -o temp/otu_table_raw.biom --table-type="OTU table" --to-json
	echo 'Summary of otu_table_raw:' >> $(log)
	biom summarize-table -i temp/otu_table_raw.biom >> $(log)
	filter_samples_from_otu_table.py -i temp/otu_table_raw.biom -o result/otu_table.biom -n $(thre_count)
	echo 'Summary of otu_table:' >> $(log)
	biom summarize-table -i result/otu_table.biom >> $(log)
	biom summarize-table -i result/otu_table.biom > result/otu_table.sum
	biom convert -i result/otu_table.biom -o result/otu_table.txt --table-type="OTU table" --to-tsv
	sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table.txt

# 8. Taxonomy assignment
assign_tax: otu_table
	touch $@
	assign_taxonomy.py -i result/rep_seqs.fa -r $(gg_seq) -t $(gg_tax) -m $(method) -o result
	sed 's/;/\t/g;s/ //g' result/rep_seqs_tax_assignments.txt > result/rep_seqs_tax.txt # format for R read
	mv result/rep_seqs_tax_assignments.log temp/rep_seqs_tax_assignments.log
	biom add-metadata -i result/otu_table.biom --observation-metadata-fp result/rep_seqs_tax_assignments.txt -o result/otu_table_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy # add taxonomy to biom
	biom convert -i result/otu_table_tax.biom -o result/otu_table_tax.txt --to-tsv --header-key taxonomy
	summarize_taxa.py -i result/otu_table_tax.biom -o result/sum_taxa # summary each level percentage
	rm result/sum_taxa/*.biom
	sed -i '/# Const/d;s/#OTU ID.//g' result/sum_taxa/* # format for R read

# 9. Phylogeny tree
make_tree: assign_tax
	touch $@
	clustalo -i result/rep_seqs.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=${p}
	filter_alignment.py -i temp/rep_seqs_align.fa -o temp/  # rep_seqs_align_pfiltered.fa, only very short conserved region saved
	make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fasta -o result/rep_seqs.tree # generate tree by FastTree

# 10. Alpha diversity
alpha: make_tree
	touch $@
	# rarefaction=`head -n 7 result/otu_table.sum|tail -n 1|cut -f 3 -d ' '|cut -f 1 -d '.'`
	single_rarefaction.py -i result/otu_table.biom -o temp/otu_table_rare.biom -d $(rarefaction)
	alpha_diversity.py -i temp/otu_table_rare.biom -o result/alpha.txt -t result/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree

# 11. Beta diversity
beta: alpha
	touch $@
	normalize_table.py -i result/otu_table.biom -o temp/otu_table_css.biom -a CSS
	biom convert -i temp/otu_table_css.biom -o result/otu_table_css.txt --table-type="OTU table" --to-tsv
	sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_css.txt
	beta_diversity.py -i temp/otu_table_css.biom -o result/beta/ -t result/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
	sed -i 's/^\t//g' result/beta/*

# 12. Taxonomy tree - GraPhlAn
graphlan: beta
	touch $@
	filter_otus_from_otu_table.py --min_count_fraction $(tax_per) -i result/otu_table.biom -o temp/tax_otu_table.biom
	filter_fasta.py -f result/rep_seqs.fa -o temp/tax_rep_seqs.fa -b temp/tax_otu_table.biom 
	echo "Number of OTU abundance > $(tax_per) :" >> $(log)
	grep -c '>' temp/tax_rep_seqs.fa >> $(log)
	grep '>' temp/tax_rep_seqs.fa|sed 's/>//g' > temp/tax_rep_seqs.id
	awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$$1]=$$0} NR>FNR {print a[$$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs.id|cut -f 2-3|grep 's__'|sed 's/; */\|/g' > temp/tax_full_anno.txt 
	echo "Number of OTU abundance > $(tax_per) with fully annotation :" >> $(log)
	wc -l temp/tax_full_anno.txt >> $(log)
	echo "Number of OTU abundance > $(tax_per) with fully annotation unique:" >> $(log)
	sort temp/tax_full_anno.txt|cut -f 1|uniq|wc -l >> $(log)
	format_taxonomy2lefse.pl -i temp/tax_full_anno.txt -o temp/tax_lefse.txt 
	## order
	export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_order.tree --annotation temp/tax_order.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 4 --min_clade_size 1 --min_font_size 5
	graphlan_annotate.py --annot temp/tax_order.annot temp/tax_order.tree temp/tax_order.xml
	sed -i 's/ref="A:1">o  /ref="A:1">/g' temp/tax_order.xml
	graphlan.py --dpi 300 temp/tax_order.xml result/tax_order.pdf --external_legends
	graphlan.py --dpi 300 temp/tax_order.xml result/tax_order.png --external_legends
	mv result/tax_order_legend.* temp/ 
	## family
	export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_family.tree --annotation temp/tax_family.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 5 --min_clade_size 1 --min_font_size 4
	graphlan_annotate.py --annot temp/tax_family.annot temp/tax_family.tree temp/tax_family.xml
	sed -i 's/ref="A:1">f  /ref="A:1">/g' temp/tax_family.xml
	graphlan.py --dpi 300 temp/tax_family.xml result/tax_family.pdf --external_legends
	graphlan.py --dpi 300 temp/tax_family.xml result/tax_family.png --external_legends
	mv result/tax_family_legend.* temp/ 
	## genus
	export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_genus.tree --annotation temp/tax_genus.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 6 --min_clade_size 1 --min_font_size 3
	graphlan_annotate.py --annot temp/tax_genus.annot temp/tax_genus.tree temp/tax_genus.xml
	sed -i 's/ref="A:1">g  /ref="A:1">/g' temp/tax_genus.xml
	graphlan.py --dpi 300 temp/tax_genus.xml result/tax_genus.pdf --external_legends
	graphlan.py --dpi 300 temp/tax_genus.xml result/tax_genus.png --external_legends
	mv result/tax_genus_legend.* temp/ 

# 13. Phylogenetic tree - ggtree
ggtree: graphlan
	touch $@
	clustalo -i temp/tax_rep_seqs.fa -o temp/tax_rep_seqs_clus.fa --seqtype=DNA --full --force --threads=$p
	make_phylogeny.py -i temp/tax_rep_seqs_clus.fa -o temp/tax_rep_seqs.tree
	sed "s/'//g" temp/tax_rep_seqs.tree > result/tax_rep_seqs.tree # remove '
	grep '>' temp/tax_rep_seqs_clus.fa|sed 's/>//g' > temp/tax_rep_seqs_clus.id
	awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$$1]=$$0} NR>FNR {print a[$$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs_clus.id|sed 's/;/\t/g'|cut -f 1-5 |sed 's/p__//g;s/c__//g;s/o__//g' > result/tax_rep_seqs.tax
	ggtree.sh -e FALSE # need R3.3.3 on windows, server not work well


# 14. Visuallize diversity, draw alpha, beta and Constrain PCoA
diversity: ggtree
	touch $@
	diversity.sh -d $(design) -m $(merge_group) -c $(compare_group) -p $(pair_compare) -A $(g1) -B $(g1_list) -C $(g2) -D $(g2_list) -o $(result) -g $(group_order) -h $(height) -w $(width) -s $(text_size)

# 15. Visuallize taxonomy, draw barplot+error bar, stack plot, first using qimme + limma; need update to count and edgeR
taxonomy: diversity
	touch $@
	taxonomy_egr.sh -d $(design) -m $(merge_group) -c $(compare_group) -p $(pair_compare) -A $(g1) -B $(g1_list) -C $(g2) -D $(g2_list) -o $(result) -g $(group_order) -n $(tax_number) -h $(height) -w $(width) -s $(text_size)

# 16. Visuallize DEOTU, draw volcano, manhattan, heatmap, venn
DAOTU: taxonomy
	touch $@
	DAOTU_egr.sh -d $(design) -m $(merge_group) -c $(compare_group) -p $(pair_compare) -A $(g1) -B $(g1_list) -C $(g2) -D $(g2_list) -o $(result)
	


## Filter OTU by abundance and reanalyze, such as 0.001%
# 1. Subset OTU table by abundance and taxonomy
filter:
	touch $@
	mkdir -p $(result_f)
	filter_otus_from_otu_table.sh -t $(thre) -o $(result)
	filter_otus_from_otu_table.py -i $(result)/otu_table_tax.biom -o temp/k1.biom --otu_ids_to_exclude_fp $(result)/otu_id_k1.txt --negate_ids_to_exclude
	echo 'Summary of otu_table_k1, one of sample OTU > 0.1%:' >> $(log)
	biom summarize-table -i temp/k1.biom  >> $(log)
	filter_taxa_from_otu_table.py -i temp/k1.biom -o $(result_f)/otu_table.biom -n $(taxonomy)
	echo 'Summary of otu_table_k1 remove:'$(taxonomy) >> $(log)
	biom summarize-table -i $(result_f)/otu_table.biom >> $(log)
	biom summarize-table -i $(result_f)/otu_table.biom > $(result_f)/otu_table.sum
	filter_fasta.py -f $(result)/rep_seqs.fa -o $(result_f)/rep_seqs.fa -b $(result_f)/otu_table.biom
	ln -f $(result_f)/otu_table.biom $(result_f)/otu_table_tax.biom
	summarize_taxa.py -i $(result_f)/otu_table_tax.biom -o $(result_f)/sum_taxa
	rm $(result_f)/sum_taxa/*.biom
	sed -i '/# Const/d;s/#OTU ID.//g' $(result_f)/sum_taxa/*
	biom convert -i $(result_f)/otu_table.biom -o $(result_f)/otu_table.txt --table-type="OTU table" --to-tsv
	sed -i '/# Const/d;s/#OTU ID.//' $(result_f)/otu_table.txt 
	cut -f 1 $(result_f)/otu_table.txt | tail -n+2 > temp/k1_t.id
	awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$$1]=$$0} NR>FNR {print a[$$1]}' $(result)/rep_seqs_tax.txt temp/k1_t.id > $(result_f)/rep_seqs_tax.txt

# 2. Re-analyze new OTU table
rediv: filter
	touch $@
	clustalo -i $(result_f)/rep_seqs.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=$(p)
	filter_alignment.py -i temp/rep_seqs_align.fa -o temp/
	make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fasta -o $(result_f)/rep_seqs.tree
	single_rarefaction.py -i $(result_f)/otu_table.biom -o temp/otu_table_rare.biom -d $(rarefaction)
	alpha_diversity.py -i temp/otu_table_rare.biom -o $(result_f)/alpha.txt -t $(result_f)/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree
	normalize_table.py -i $(result_f)/otu_table.biom -o temp/otu_table_css.biom -a CSS
	biom convert -i temp/otu_table_css.biom -o $(result_f)/otu_table_css.txt --table-type="OTU table" --to-tsv
	sed -i '/# Const/d;s/#OTU //g;s/ID.//g' $(result_f)/otu_table_css.txt
	beta_diversity.py -i temp/otu_table_css.biom -o $(result_f)/beta/ -t $(result_f)/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
	sed -i 's/^\t//g' $(result_f)/beta/*

# 3. redraw all figure
draw_div: rediv
	touch $@
	diversity.sh -d $(design) -m $(merge_group) -c $(compare_group) -p $(pair_compare) -A $(g1) -B $(g1_list) -C $(g2) -D $(g2_list) -o $(result_f) -g $(group_order) -h $(height) -w $(width) -s $(text_size)

draw_tax: draw_div
	touch $@
	taxonomy_egr.sh -d $(design) -m $(merge_group) -c $(compare_group) -p $(pair_compare) -A $(g1) -B $(g1_list) -C $(g2) -D $(g2_list) -o $(result_f) -g $(group_order) -n $(tax_number) -h $(height) -w $(width) -s $(text_size)
	plot_pie_DA_Bphylum.sh -c $(compare_group) -l family -o $(result_f)
	batch_venn.pl -i $(venn) -d $(result_f)/family.txt

draw_otu: draw_tax
	touch $@
	DAOTU_egr.sh -d $(design) -m $(merge_group) -c $(compare_group) -p $(pair_compare) -A $(g1) -B $(g1_list) -C $(g2) -D $(g2_list) -o $(result_f) -h $(height) -w $(width) -s $(text_size)
	plot_pie_DA_Bphylum.sh -c $(compare_group) -l otu -o $(result_f)
	batch_venn.pl -i $(venn) -d $(result_f)/otu.txt
	
