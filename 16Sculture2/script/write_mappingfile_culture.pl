#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

# 编写分菌鉴定mappingfile的脚本
# 通常每个文库有48个板index，每板96孔
# 基本信息为#SampleID       BarcodeSequence LinkerPrimerSequence    ReversePrimer   barcodeF        barcodeR        plate
# 可选信息有variety	compartment	medium    batch   species Description

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq!
Usage:    write_mappingfile_culture.pl -i genome_chr.len (-c chr -s start -e end) -o output_file.bed -b bin_size -l slid_step
Function: write a mapping file for identifying culture microbe
Command:  -i barcodeF file for hole
          -b barcodeR file for plate
          -p plate number, default=48
          -o output file name (Must), such as L1, L2
          -F forward, primer, default AACMGGATTAGATACCCKG, 16S 799F
          -R reverse, primer, default ACGTCATCCCCACCTTCC, 16S 1193R
          -c compartment type, default root, also can rhizosphere, soil, leaf
          -L library, default L1
          -v variety, default A50
          -m medium, default TSB
          -s species, default rice
          -d description, default WildType
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v2.0
Update:   2019/6/26
Notes:    Change hole 1-96 to A1-H12, Add variety and medium between compartment column
\n!
    )
}

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:b:p:F:R:c:s:B:L:P:v:m:', \%opts );
&usage unless ( exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{i}="/mnt/bai/yongxin/github/Amplicon/16Sculture2/doc/barcodeF96.txt" unless defined($opts{i});
$opts{b}="/mnt/bai/yongxin/github/Amplicon/16Sculture2/doc/barcodeR48.txt" unless defined($opts{b});
$opts{p}=48 unless defined($opts{p});
$opts{o}="L1.txt" unless defined($opts{o});
$opts{F}="AACMGGATTAGATACCCKG" unless defined($opts{F});
$opts{R}="ACGTCATCCCCACCTTCC" unless defined($opts{R});
$opts{c}="root" unless defined($opts{c});
$opts{B}="1" unless defined($opts{B});
$opts{s}="rice" unless defined($opts{s});
$opts{d}="WildType" unless defined($opts{d});
$opts{L}="L1" unless defined($opts{L});

open INPUT,"<$opts{i}";
my @barcodeF; #database in array
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	push @barcodeF,$tmp[0];
}
close INPUT;

open DATABASE,"<$opts{b}";
my @barcodeR; #database in array
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @barcodeR,$tmp[0];
}
close DATABASE;

my @hole;
@line=("A","B","C","D","E","F","G","H");
my $i=1;
foreach  $l(@line) {
	foreach  $r(1..12) {
          $hole[$i]="$l$r";
#          print $hole[$i],"\n";
          $i++;
	}
}

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
print OUTPUT "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tbarcodeF\tbarcodeR\tplate\tvariety\tcompartment\tmedium\tbatch\tspecies\tDescription\n";
foreach $i(1..$opts{p}) {
	foreach $j(1..96) {
          $F=$j-1;
          $R=$i-1;
          printf OUTPUT "$opts{L}P%02d$hole[$j]\t$barcodeF[$F]$barcodeR[$R]\t$opts{F}\t$opts{R}\t$barcodeF[$F]\t$barcodeR[$R]\t$opts{L}P%02d\t$opts{v}\t$opts{c}\t$opts{m}\t$opts{B}\t$opts{s}\t$opts{d}\n",$i,$i;
	}
}
close OUTPUT;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

