# Summary split and mapping 
mkdir -p result/qc
#touch temp/library.stat
cp /mnt/bai/yongxin/ref/amplicon/qc.sum result/qc.sum
list=`ls clean_data/*.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|uniq|tr "\n" " "` # library name list
for library in ${list}
    do
    echo ${library}
    echo "${library}"> temp/${library}.count # library ID
#    unzip clean_data/${library}_1_fastqc.zip
#    cp clean_data/${library}_1_fastqc.html result/qc/
#    cp ${library}_1_fastqc/Images/per_base_quality.png result/qc/${library}_1_quality.png
#    cp ${library}_1_fastqc/Images/duplication_levels.png result/qc/${library}_1_duplication.png
#    unzip clean_data/${library}_2_fastqc.zip
#    cp clean_data/${library}_2_fastqc.html result/qc/
#    cp ${library}_2_fastqc/Images/per_base_quality.png result/qc/${library}_2_quality.png
#    cp ${library}_2_fastqc/Images/duplication_levels.png result/qc/${library}_2_duplication.png
    grep 'Total Sequences' clean_data/${library}_1_fastqc/fastqc_data.txt|cut -f 2 >> temp/${library}.count
#    rm -r ${library}_1_fastqc 
#    rm -r ${library}_2_fastqc
    grep -P '^\d+' ${library}.stat >> temp/${library}.count
    #paste temp/library.stat temp/${library}.count > temp/library.stat
    temp=`cat temp/${library}.count|tr "\n" "\t" ` # library count list
    echo ${temp}>>result/qc.sum
    done
sed -i 's/ /\t/g' result/qc.sum
cat result/qc.sum
cat temp/length_* > result/length.txt
