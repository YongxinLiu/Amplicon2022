#!/usr/bin/env perl
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
#use List::MoreUtils qw(none any all);
sub usage {
    die(
        qq!
Usage:    merge_id_with_same_sequence.pl -i inpute_file -o output_file 
Function: Merge id with the same sequence to one id
Command:  -i inpute file name (Must)
          -o output file name (Must)
	  -d database file name
Author:   Zhai Zhi-wen, arvin16\@126.com, QQ:786049393
Version:  v1.0
Update:   2018/10/23
Notes:    
\n!
    )
}
my %opts;
getopts('i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
# 调置参数的初始值，可以添加更多参数的默认值
$opts{h}=1 unless defined($opts{h});

open DATABASE,"<$opts{d}";
# 1. 散列结构数据库，要求数据文件有唯一ID并且无顺序要求
my %database; #database in hash
while (<DATABASE>) {
    chomp;
    my @tmp=split/\t/;
    $database{$tmp[5]} = $tmp[6];	
}
close DATABASE;

open (IN,"<$opts{i}");
open (OT,">$opts{o}");
while (<IN>){
        chomp;s/\r//;
	my %hash;
	next unless $_ =~ /^>/;
	my $merge_id = $_;
	$merge_id =~ s/^\>//;
       	my @temp = split(/\,/,$merge_id);
		foreach my $id (@temp){
			if ($database{$id} eq "A50"){$hash{A50} += 1}else{$hash{A50} += 0}; 
			if ($database{$id} eq "IR24"){$hash{IR24} += 1}else{$hash{IR24} += 0};
		}
		foreach my $tkey (keys %hash){
			if ($hash{$tkey} ge 1){$hash{$tkey} = "$tkey"}else{$hash{$tkey} = NA}
		}
	print OT "$temp[0]\t$hash{A50}\t$hash{IR24}\n";
	} 
close IN;
close OT;
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
