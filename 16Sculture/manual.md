# Data prepare
cd ~/culture/wheat1903
## 复制所有相关脚本
ln ~/github/Amplicon/16Sculture/parameter.md makefile
ln ~/github/Amplicon/16Sculture/manual.md manual.sh

pwd # 查看当前工作目录并修改变量wd
make init

## 准备原始数据并重命名，准备doc/library.txt，存文库L1/2与目前名的对应关系，按列表获得数据
awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_1.fq.gz seq/"$1"_1.fq.gz");}' <(tail -n+2 doc/library.txt )
awk 'BEGIN{OFS=FS="\t"}{system("ln -s /mnt/bai/yongxin/seq/amplicon/"$2"_2.fq.gz seq/"$1"_2.fq.gz");}' <(tail -n+2 doc/library.txt )
# 显示库名列表
ls seq/*.gz|sed 's/\t/\n/g'|cut -f 2 -d '/'|cut -f 1 -d '_'|uniq|tr '\n' ' '

## 写mappingfile, s为物种，p为板数；多个library需要多个mappingfile
# 有规律的文章用程序编写，没规律的手动填写
#for RPM in `seq 5 10`; do
#	write_mappingfile_culture.pl -o doc/L${RPM}.txt -s rice -p 48 -L L${RPM}
#done
# 其它手动保存，如L1
cat <(head -n1 doc/L1.txt|sed 's/#//g') doc/L* |grep -v '#' > doc/design.txt
# 检查是否有重名ID
wc -l doc/design.txt
cut -f 1 doc/design.txt|sort|uniq |wc -l
# 如果数变小，检查重复名ID
cut -f 1 doc/design.txt|sort|uniq -d

# 多文库并行分析
make qc # 序列质量评估
make merge # 双端序列合并
make split # 按barcode拆分样品
make split.stat # 拆分样品统计绘图
make cutp # 切除双端引物
make stat # 统计各步信息

# 合并数据后单线程分析
make merge_library # 合并所有文库数据，并统计表格
make derep # 
make unoise # OTU聚类，采用unoise方法
make otu_table # 生成OTU表，采用otutab方法
make assign_tax
make make_tree
make identify_isolate
# 数据质量量统计，和分析查看
firefox clean_data/L*_2*.html # 查看所有分析报告右端，直接看低质量即可
firefox result/stat_lib_split*.pdf # 查看所有数据分布，更建议下载本地放大查看

# 检查每一个库中reads>30的样本数据，发现过30的样品就有6k+，OTU表怎么只有不到5k；检查merge_library处合并的seqs_usearch.fa并没有各库总合大，重新运行
for RPM in `seq 1 6`; do
	tail -n+16 temp/L${RPM}_split/split_library_log.txt|head -n-3|awk '$2>30'|wc -l
done

# 统计OTU注释的物种信息
for RPM in `seq 2 8`; do
cut -f ${RPM} result/rep_seqs_tax.txt|sort|uniq|wc -l # 一共有多少种
cut -f ${RPM} result/rep_seqs_tax.txt|sort|uniq -c|awk '{print $2"\t"$1}'>result/count_tax_${RPM}.txt # 每种的组成
done

# 统计样品中的正负对照：A11为+，A12为-
grep -P 'A11' result/otu_table_raw.sum
grep -P 'A12' result/otu_table_raw.sum

for j in A B C D E F G H; do
	echo -ne "${j}\n"
for i in `seq 1 12`; do
	echo -ne "${i}\t"
	grep "${j}${i}\:" result/otu_table_raw.sum|wc -l 
done
done

# 对OTU表进行每孔注释，再选OTU，程序为identify_isolate.r
grep 'A1[12]' result/culture_bacteria.xls # 三个孔纯度低
grep 'A1[12]' result/culture_select.xls 
grep -c '>' result/rep_seqs.fa # 统计Zotu数量 
wc -l result/*.xls # 统计样品大于30的wells和候选选OTU数量

# 按分组分别找菌，并筛选可培养菌
temp=temp
result=result
i=5
for lib in RHTANP RHTANp RHTAnP BDNP BDNp BDnP; do
ln doc/L${i}.txt doc/${lib}.txt
filter_samples_from_otu_table.py -i ${result}/otu_table.biom -o ${result}/${lib}_otu_table.biom --sample_id_fp doc/${lib}.txt
biom convert -i ${result}/${lib}_otu_table.biom -o ${result}/${lib}_otu_table.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU //g' ${result}/${lib}_otu_table.txt
identify_isolate.sh -f ${lib}_otu_table.txt -o ${lib}
tail -n+2 result/${lib}culture_select.xls | cut -f 1,2|sed 's/^/OTU_/g;s/;/\t/g'|less>result/${lib}culture_select.tax
filter_fasta.py -f result/rep_seqs.fa -o result/${lib}culture_select.fa -s result/${lib}culture_select.tax
makeblastdb -dbtype nucl -in result/${lib}culture_select.fa
((i=$i+1))
done
# for total 
lib=""
tail -n+2 result/${lib}culture_select.xls | cut -f 1,2|sed 's/^/OTU_/g;s/;/\t/g'|less>result/${lib}culture_select.tax
filter_fasta.py -f result/rep_seqs.fa -o result/${lib}culture_select.fa -s result/${lib}culture_select.tax
makeblastdb -dbtype nucl -in result/${lib}culture_select.fa


# 6种条件下合并，并添加丰度
# 合6批库的表，每库选两个
cd ~/culture/wheat/result
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} NR>FNR {print $1,a[$1]}' <(cut -f 1,3-8 1BDNPculture_select.xls|sed '1 s/well/1BDNP/g') 0culture_select.xls | sed 's/\t$/\t\t\t\t\t\t/g' > culture_select1.xls
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} NR>FNR {print $0,a[$1]}' <(cut -f 1,3-8 2BDNpculture_select.xls|sed '1 s/well/2BDNp/g') culture_select1.xls | sed 's/\t$/\t\t\t\t\t\t/g' > culture_select2.xls
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} NR>FNR {print $0,a[$1]}' <(cut -f 1,3-8 3BDnPculture_select.xls|sed '1 s/well/3BDnP/g') culture_select2.xls | sed 's/\t$/\t\t\t\t\t\t/g' > culture_select3.xls
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} NR>FNR {print $0,a[$1]}' <(cut -f 1,3-8 4RHTANPculture_select.xls|sed '1 s/well/4RHTANP/g') culture_select3.xls | sed 's/\t$/\t\t\t\t\t\t/g' > culture_select4.xls
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} NR>FNR {print $0,a[$1]}' <(cut -f 1,3-8 5RHTANpculture_select.xls|sed '1 s/well/5RHTA/g') culture_select4.xls | sed 's/\t$/\t\t\t\t\t\t/g' > culture_select5.xls
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} NR>FNR {print $0,a[$1]}' <(cut -f 1,3-8 6RHTAnPculture_select.xls|sed '1 s/well/6RHTA/g') culture_select5.xls | sed 's/\t$/\t\t\t\t\t\t/g' > culture_select6.xls
wd=/mnt/bai/yongxin/wheat/NP
cluture_db=/mnt/bai/yongxin/culture/wheat/result/culture_select.fa
cp $wd/result/otu_table.txt ./
filter_otus_from_otu_table.sh -t 0.0001 -o ./ # 过滤万分之一的 946
# 比对建立联系
blastn -query ${wd}/result/rep_seqs.fa -db ${cluture_db} -out rep_seqs.blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -num_alignments 1 -evalue 1 -num_threads 9 # 输出13列为coverage
# 筛选97%相似菌
awk '$3*$13>=9700' rep_seqs.blastn|cut -f 1-3 |sed 's/\tOTU_/\t/g' > otu_cultured.txt
# 添加自然样品丰度，存在一个菌对应多个OTU，按丰度从小到大排序，读取为保留最大的
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {printf ($0"\t%.4f\n",(a[$1]*100))}' otu_table.txt.mean otu_cultured.txt | sort -k4,4n > otu_cultured.info
# 添加相似度和丰度
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$2]=$3"\t"$4} NR>FNR {print $0,a[$1]}' otu_cultured.info culture_select6.xls | sed '1 s/$/similarity\tabundance/' > culture_select7.xls
# 添加物种注释
awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$1]=$2} NR>FNR {print $0,a[$1]}' 0culture_select.xls culture_select7.xls > culture_select8.xls

