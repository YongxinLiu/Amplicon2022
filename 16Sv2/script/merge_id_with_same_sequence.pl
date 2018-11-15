#!/usr/bin/env perl
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
sub usage {
    die(
        qq!
Usage:    merge_id_with_same_sequence.pl -i inpute_file -o output_file 
Function: Merge id with the same sequence to one id
Command:  -i inpute file name (Must)
          -o output file name (Must)
Author:   Zhai Zhi-wen, arvin16\@126.com, QQ:786049393
Version:  v1.0
Update:   2018/10/23
Notes:    
\n!
    )
}
my %opts;
getopts('i:o:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
# 调置参数的初始值，可以添加更多参数的默认值
$opts{h}=1 unless defined($opts{h});

open (IN,"<$opts{i}");
open (OT,">$opts{o}");
my %hash;
my @a;
while (<IN>){
        chomp;s/\r//;
       	if ($_ =~ /^>/ ){
    	$id = $_ ;
	$id =~ s/^>//
	}else{
	my $seq = $_;
	push @{$hash{$seq}},$id
	}
	}
foreach (keys %hash){
	my $idd = join(',', @{$hash{$_}});
	print OT "\>$idd\n$_\n";
	}
close IN;
close OT;
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
