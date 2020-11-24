#!/usr/bin/perl -w
# 加载时间管理，参数管理，文件名和路径处理的基础包，无须安装
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Scripts usage and about.
# 程序的帮助文档，良好的描述是程序重用和共享的基础，也是程序升级和更新的前提
###############################################################################
# 参考 /mnt/zhou/zhiwen/maize_16s_lanesplit/ASV/v2/JT/result/GWAS/otu_count/emmax_result/1/test_gwas_pipline/SNP_merge.pl
sub usage {
    die(
        qq!
Usage:    gwas_snp_merge.pl -i emmax_result_sig.ps -o cluster -d 20000
Function: Identify significantlly SNP distance < 20k as cluster 鉴定20K内的显著SNP作为QTL
Command:  -i emmax result < 1e-6 SNP (Must)
          -o output file name (Must)
          -d distance, default 20000
          -h header line number, default 0
Author:   Liu Yong-Xin, liuyongxin_bio\@163.com, QQ:42789409
Version:  v1.0
Update:   2019/12/30
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
print "Distance file is $opts{d}\n" if defined($opts{d});
# 调置参数的初始值，可以添加更多参数的默认值
$opts{h}=1 unless defined($opts{h});
$opts{d}=20000 unless defined($opts{h});

###############################################################################
#读入的数据或注释文件，用于与输入文件比较或注释(可选)，提供三种方式
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
# 1. 散列结构数据库，要求数据文件有唯一ID并且无顺序要求
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}
# 2. 数组结构数据库，无唯一ID，但有顺序要求
#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;
# 3. 批量数据文件，读取一批有相似结构的文件
#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
# 正文部分，读取输入文件，列出输入和输入文件的三行作为示例，方便编程处理数据
open INPUT,"<$opts{i}";
# 输入文件为emmax的标准输出结果，4列格式，最重要的是第1列SNP ID和最后一列Pvalue
#SNP	beta	beta SE	Pvalue
#1m36655774      1.1076  0.19174 3.195e-08
#1m36659693      0.87607 0.16428 2.8148e-07
#3m27929713      0.93875 0.18157 6.0614e-07
open OUTPUT,">$opts{o}";
#chrm    snppos          ref     mat_gtyp        pat_gtyp        c_gtyp  phase   mat_all pat_all cA      cC      cG      cT      winning SymCls  SymPval BindingSite     cnv
#1       4648    C       A       C       M       PHASED  C       A       0       11      0       0       M       Asym    0.0009765625    -1      0.902113

# 按取首行，设置起始点
my $FL = readline INPUT;
# 提取1，4列ID和P值
my ($ID1,$PV1) = (split "\t",$FL)[0,3];
# 拆分ID为染色体和位置，不同文件分隔符可能不同，需修改，如水稻中为m，玉米中为点.
my ($CHR,$P1) = split /\./,$ID1;
chomp($PV1);
$P1 =~ s/s_//gi;
my ($start,$end);
$start = $P1;
$end = $P1;
my $n = 1;
my $range = 0;
while(<INPUT>){
	chomp;
	my ($ID,$PV) = (split "\t",$_)[0,3];
	my ($chr,$pos) = split /\./,$ID;
	$pos =~ s/s_//gi;
	if($chr eq $CHR){
		if ($pos - $P1 <= $opts{d}){
			$n++;
			$P1 = $pos;
			$end = $pos;
			if($PV < $PV1){
				$PV1 = $PV;
				$ID1 = $ID;
			}
		}else{
			$range = $end - $start;
			print OUTPUT "$ARGV[0]\t$CHR\t$start\t$end\t$range\t$ID1\t$PV1\t$n\n";
			$P1 = $pos;
			$n = 1;
			$start = $pos;
			$end = $pos;
			$PV1 = $PV;
			$ID1 = $ID;
		}
	}else{
		$range = $end - $start;
		print OUTPUT "$ARGV[0]\t$CHR\t$start\t$end\t$range\t$ID1\t$PV1\t$n\n";
		$P1 = $pos;
		$CHR = $chr;
		$n = 1;
		$start = $pos;
		$end = $pos;
		$PV1 = $PV;
		$ID1 = $ID;
	}
}
$range = $end - $start;
print OUT "$ARGV[0]\t$CHR\t$start\t$end\t$range\t$ID1\t$PV1\t$n\n";
close INPUT;
close OUTPUT;

###############################################################################
#Record the program running time!
# 输出程序运行时间
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

