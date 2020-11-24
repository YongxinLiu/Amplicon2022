#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
##database in hash
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}
##database in array
#my (@tmp1,@tmp2);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#>Cluster 0
#0       1473nt, >1467... *
#>Cluster 1
#0       1457nt, >1335-2... at +/99%
#1       1460nt, >1990... *

open OUTPUT,">$opts{o}";
my %seq;
my %name;
my %rep;
while (<INPUT>) {
	if (/^>/){
		chomp;
		$id=$_;
		$id=~s/>Cluster //;
		next;
	}else{
		$seq{$id}++;
		my @tmp=split/\s+/;
		$tmp[2]=~s/>//;
		$tmp[2]=~s/\.//g;
		#print $tmp[2],",";
		$name{$id}.=$tmp[2].",";
		if(/\*$/){
		$rep{$id}=$tmp[2];
		}
	}
}
foreach  (sort { $a <=> $b } keys %seq) {
#cluster	representitive	number	ID
#0       1467    1       1467,
#1       1990    2       1335-2,1990,
#2       2318    2       1862,2318,
	print OUTPUT "$_\t$rep{$_}\t$seq{$_}\t$name{$_}\n";
}
close OUTPUT;


###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq/
Usage:    cd-hit_clstr_count_ID.pl -i input_file -o output_file
Function: format fasta into sequence in 1 line format
Command:  -i cd-hit clstr (Must)
          -o output (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2020-07-22
Notes:    
\n/
    )
}
