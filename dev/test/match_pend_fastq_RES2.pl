#!/usr/bin/perl
use strict;
use warnings;

####################################################################################################
#
#	Author: K. Wyatt McMahon, Ph.D.
#	Date: December 19, 2011
#
#   modified by Robert E Settlage, Ph.D on May 30, 2012
#	match_pend_fastq.pl - makes sure that the same sequences are in each file in paired fastq files
#
#	use: perl /groups/DAC/useful_perl/match_pend_fastq_RES2.pl <pe1.fastq>
#
#   output is fastq_base.paired_matched.tag.fastq
#####################################################################################################


#if(scalar(@ARGV)<1){
#	die "use: perl match_pend_fastq.pl <pe1.fastq>\n";
#}

my $read_indicator1_len = 2;
my $read_indicator2_len = 2;
my $search ="";

my $file1_NAME=$ARGV[0];
my $file2_NAME=$file1_NAME;
$file2_NAME=~s/_R1_/_R2_/;

my $zippedFlag = 0;
if ($file1_NAME =~ m/\.gz$/) {
	`gzip -d -v $file1_NAME`;
	$file1_NAME=~s/\.gz$//;
	`gzip -d -v $file2_NAME`;
	$file2_NAME=~s/\.gz$//;
	$zippedFlag = 1;
}

print "new files are:\n";
print $file1_NAME."\n";
print $file2_NAME."\n";

open(ONE, $file1_NAME)|| die "Couldn't open $file1_NAME\n";
open(TWO, $file2_NAME)|| die "Couldn't open $file2_NAME\n";

my (%onehash, %twohash);
my $line=<ONE>;

if ($line=~/^@/) {
	$search=substr $line, 0, 8;  ##this is getting a longer tag to find header lines
	if ($line=~/:N:0:/) {
		my @temp_array=split(/ /, $line);
		print "read tag is: ".$temp_array[1]."\n";
		$read_indicator1_len = length($temp_array[1]);
		$read_indicator2_len = length($temp_array[1]);
	}
}else{
	print "something is funky with the file, check it out\n";
	print "line 1 is $line\n";
	print "exiting\n";
	exit;
}

while($line){
	if($line=~/$search/){
		my $line_len = length($line);
		$line=substr $line, 0, $line_len - $read_indicator1_len;
		$onehash{$line}=1;
	}
	$line=<ONE>;
}

close ONE;
my $base = $file2_NAME;
$base=~ s/fastq$//;
my $out=$base . "paired_matched.fastq";
print "parsing $file2_NAME to get $base to output to $out\n";
open(OUT, ">$out") || die "Couldn't open $out\n";

my $line2=<TWO>;
my $count2=0;
while($line2){
	$count2++;
	my $line_len2 = length($line2);
	my $lookup2=substr $line2, 0, $line_len2 - $read_indicator2_len;
	#print $lookup2;
	if(exists($onehash{$lookup2})){
		$twohash{$lookup2}=1;
		print OUT $line2;	
		my $seq=<TWO>;
		print OUT $seq;
		my $title2=<TWO>;
		print OUT "+\n";
		my $scores=<TWO>;
		print OUT $scores;
		}
	$line2=<TWO>;
}
close TWO;

open(IN, $ARGV[0])|| die "couldn't open $ARGV[0]\n";


my $line3=<IN>;
my $base2 = $ARGV[0];
$base2=~ s/fastq$//;
my $out2=$base2 . "paired_matched.fastq";
print "parsing $ARGV[0] to get $base2 to output to $out2\n";
open (OUT2, ">$out2") || die "Couldn't open $out2\n";

while($line3){
	my $line_len3 = length($line3);
	my $lookup3 = substr $line3, 0, $line_len3 - $read_indicator1_len;
	if(exists($twohash{$lookup3})){
		print OUT2 "$line3";
		my $seq2=<IN>;
		print OUT2 $seq2;
		my $title3=<IN>;
		print OUT2 "+\n";
		my $scores2=<IN>;
		print OUT2 $scores2;
	}
	$line3=<IN>;
}


if ($zippedFlag eq 1) {
	`gzip -v $file1_NAME`;
	`gzip -v $file2_NAME`;
	`gzip -v $out`;
	`gzip -v $out2`;
}

exit;