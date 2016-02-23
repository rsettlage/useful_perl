#!/usr/bin/perl -w
use strict;


#########################################################################
#
#
#           Author: K. Wyatt McMahon, Ph.D.
#           Date:  October 29, 2012
#
#           use: perl get_reads_stats.pl <filename><tag>
#
#########################################################################

if(scalar(@ARGV!=2)){

    die "use: perl get_reads_stats.pl <filename><tag>\n";
   }
open (IN, $ARGV[0]);
chomp $ARGV[0];
my $outfile=$ARGV[0] . "_out";

my $tag=$ARGV[1];
chomp $tag;
open (OUT, ">$outfile");

my $line=<IN>;
my $count=0;
my $length=0;
my $qscore=0;
my $q30count=0;
my $q20count=0;

print "working on $ARGV[0]\n";
while($line){

    if($line=~/$tag/){
	$count++;
	my $seq=<IN>;
	$length=$length + length($seq);
	my $plus=<IN>;
	my $qline=<IN>;

	my @q=split('', $qline);
	my $q2=0;
	foreach my $a(@q){
	    $q2=$q2+(ord($a)-33);
	}
	my $meanqscore=$q2/length($qline);

	if($meanqscore>20){
	    $q20count++;
	}

	if($meanqscore>30){
	    $q30count++;
	}

	$qscore=$qscore+$meanqscore;

	
    }
    $line=<IN>;
}

my $avglen=$length/$count;
my $avgqscore=$qscore/$count;

my $percgt20=($q20count/$count)*100;
my $percgt30=($q30count/$count)*100;

print OUT "$ARGV[0]\t$count\t$avglen\t$avgqscore\t$percgt20\t$percgt30\n";
close OUT;
close IN;
