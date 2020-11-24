#!/usr/bin/perl -w
# 加载时间管理，参数管理，文件名和路径处理的基础包，无须安装
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Scripts usage and about.
# 程序的帮助文档，良好的描述是程序重用和共享的基础，也是程序升级和更新的前提
###############################################################################
# 参考 /mnt/zhou/zhiwen/maize_16s_lanesplit/ASV/v2/JT/result/GWAS/otu_count/emmax_result/1/test_gwas_pipline/eQTL_merge.pl
sub usage {
    die(
        qq!
Usage:    gwas_eQTL_merge.pl -i emmax/OTU_12.snp_cluster -o emmax/OTU_12.QTL_merge -d emmax/OTU_12.ld
Function: Merge QTL based on LD 根据LD值，进一步合并QTL
Command:  -i snp_merge (Must)
          -o QTL_merge (Must)
          -d LD input
Author:   Liu Yong-Xin, liuyongxin_bio\@163.com, QQ:42789409
Version:  v1.0
Update:   2019/12/31
Notes:    
\n!
    )
}

###############################################################################
#命令行参数据的定义和获取，记录程序初始时间，设置参数默认值
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "LD file is $opts{d}\n" if defined($opts{d});


###############################################################################
#Main text.
###############################################################################
# 正文部分，读取输入文件，列出输入和输入文件的三行作为示例，方便编程处理数据
open INPUT,"<$opts{i}";
#input	chr	start	end	range	snpID	pvalueMin	count
#emmax/OTU_12.sig.ps     1       37862036        37862069        33      1m37862036      3.0197e-07      2
#emmax/OTU_12.sig.ps     2       11184799        11292399        107600  2m11226790      1.5769e-07      55
#emmax/OTU_12.sig.ps     8       496298  534009  37711   8m496298        2.8287e-08      13

open OUTPUT,">$opts{o}";
#input	chr	start	end	range	snpID	pvalueMin	count
#OTU_106 chr9    39722701        39722722        21      chr9.s_39722722 4.0029e-07      2


# 统计输入文件QTL行数
my $line = `wc -l $opts{i}`;
if ($line =~ /^(\d+)/){
	$line = $1;
}

# 当存在两个以上QTL时，才读取LD文件
if ($line > 2){
# 输入文件plink根据SNP列表计算LD的R2，只保留大于0.2的值
# CHR_A         BP_A             SNP_A  CHR_B         BP_B             SNP_B           R2 
#     1     10368725   chr1.s_10368725      1     10414058   chr1.s_10414058     0.824355 
#     1    173241870  chr1.s_173241870      1    173362007  chr1.s_173362007     0.809039 
#     1    173241870  chr1.s_173241870      1    173816122  chr1.s_173816122     0.744552 

# 保存LD文件结果于hash散列中
my %hash;

my $lineLD = `wc -l $opts{d}`;
if ($lineLD =~ /^(\d+)/){
	$lineLD = $1;
}

if ($lineLD > 1){
open LD,"<$opts{d}";
readline LD;
while(<LD>){
	chomp;
	$_ =~ s/\s+/\t/g;
	$_ =~ s/^\t//;
	my ($snp1,$snp2) = (split "\t",$_)[2,5];
	$hash{"$snp1#$snp2"} = 1;
}
close LD;
}

my $FL = readline INPUT;
chomp($FL);
my ($CHR,$START,$END,$ID,$PV,$NUM) = (split "\t",$FL)[1,2,3,5,6,7];
my $SNP = $ID;
while(<INPUT>){
	chomp;
	my ($chr,$start,$end,$id,$pv,$num) = (split "\t",$_)[1,2,3,5,6,7];
	if (exists $hash{"$ID#$id"}){
		$ID = $id;
		if ($PV > $pv){
			$PV = $pv;
			$SNP = $id;
		}
		$NUM = $NUM + $num;
		$END = $end;
	}else{
		my $range = $END - $START;
		print OUTPUT "$opts{i}\t$CHR\t$START\t$END\t$range\t$SNP\t$PV\t$NUM\n";
		$ID = $id;
		$CHR = $chr;
		$START = $start;
		$END = $end;
		$SNP = $id;
		$PV = $pv;
		$NUM = $num;
		}
	}
	close INPUT;
	close OUTPUT;
}elsif ($line == 2){
	# 只有两行的LD时
	`mv $opts{i} $opts{o}`;
}

###############################################################################
#Record the program running time!
# 输出程序运行时间
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

