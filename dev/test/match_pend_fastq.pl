#!/usr/bin/perl
use strict;
use warnings;

####################################################################################################
#
#	Author: K. Wyatt McMahon, Ph.D.
#
#
#	Date: December 19, 2011
#
#
#	match_pend_fastq.pl - makes sure that the same sequences are in each file in paired fastq files
#
#	use: perl match_pend_fastq.pl <pe1.fastq> <pe2.fastq> [searchline]
#
#
#####################################################################################################


if(scalar(@ARGV)<2){
	die "use: perl match_pend_fastq.pl <pe1.fastq> <pe2.fastq>\n";
}

my $search='@M00317';

if(length($ARGV[2])>0){
	$search=$ARGV[2];
	chomp $search;
}

open(ONE, $ARGV[0])|| die "Couldn't open $ARGV[0]\n";
open(TWO, $ARGV[1])|| die "Couldn't open $ARGV[1]\n";

my (%onehash, %twohash);
my $line=<ONE>;

while($line){
	
	if($line=~/$search/){
		$onehash{$line}=1;
	}
	$line=<ONE>;
}

close ONE;
my ($base, $fastq)=split(/\./, $ARGV[1]);
my $out="./matched/" . $base . ".matched.fastq";
open(OUT, ">$out");

my $line2=<TWO>;

my $count2=0;
while($line2){
	$count2++;
	my $lookup=$line2;
	$lookup=~s/\/2/\/1/;
	print $lookup;
	if(exists($onehash{$lookup})){

		$twohash{$lookup}=1;
		print OUT $line2;	
		my $seq=<TWO>;
		print OUT $seq;
		my $title2=<TWO>;
		print OUT $title2;
		my $scores=<TWO>;
		print OUT $scores;
		}
	$line2=<TWO>;
}
close TWO;

open(IN, $ARGV[0])|| die "couldn't open $ARGV[0]\n";


my $line3=<IN>;

my ($base2, $fastq2)=split(/\./, $ARGV[0]);
my $out2="./matched/" . $base2 . ".matched.fastq";
open (OUT2, ">$out2");

while($line3){
	
	if(exists($twohash{$line3})){
		print OUT2 "$line3";
		my $seq2=<IN>;
		print OUT2 $seq2;
		my $title3=<IN>;
		print OUT2 $title3;
		my $scores2=<IN>;
		print OUT2 $scores2;
	}
	$line3=<IN>;

}
