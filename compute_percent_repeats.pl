#!/usr/bin/perl
# Lauren McIver
# 02/10/2014

# This script will compute the percent of a input sequence that is made up of repeats.
# This script uses TRF which was developed at Boston University. Please cite the paper at: http://tandem.bu.edu/trf/trf.html

use strict;
use warnings;

# These are the parameters we use for our microsatellite papers.
# The numbers in the CMMD should be the same as that found in OUTPUT_EXT.
use constant CMMD=>"./trf404.linux64.exe FILE 2 5 5 80 10 14 6 -h";
use constant OUTPUT_EXT=>".2.5.5.80.10.14.6.dat";

#NOTE: input file should be in fasta format (like example temp.fa).
(my $file)=@ARGV;

# run trf on input file
my $cmmd=CMMD;
$cmmd=~s/FILE/$file/;

print "\n\nExecuting $cmmd\n";
my $ret=`$cmmd`;
print "$ret\n";

my $trf_output_file=$file.OUTPUT_EXT;

# count the number of bases in input file that are part of repeats
open(INFILE,"< $trf_output_file") or die("Could not open file $trf_output_file\n");

my $total_repeat_bases=0;
while(my $line=<INFILE>)
{
	chomp($line);
	if($line=~/^[0-9]/)
	{
		my @tokens=split(" ",$line);
		my $seq=pop(@tokens);
		$total_repeat_bases=$total_repeat_bases+length($seq);
	}
}
close(INFILE);

# count the number of bases in the input file
open(INFILE,"< $file") or die("Could not open file $file\n");

my $total_bases=0;
while(my $line=<INFILE>)
{
	chomp($line);
	if($line=~/^[A|T|G|C]/)
	{
		$total_bases=$total_bases+length($line);
	}
}

my $percent=($total_repeat_bases/$total_bases)*100;
print "\n\nFasta file $file is composed of $percent % repeats.\n";
print "Repeat sequences for $file are located in $trf_output_file\n";
