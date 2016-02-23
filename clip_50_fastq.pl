#!/usr/bin/perl -w
use strict;

print "clipping $ARGV[0]\n";

open (IN, $ARGV[0]);

my $outfile=$ARGV[0] . ".minus50.fastq";
open (OUT, ">$outfile");

my $line=<IN>;

while($line){

    chomp $line;
    if($line=~/^\@HWI/){

	my $title=$line;
	my $seq=<IN>;
	my $plus=<IN>;
	my $qline=<IN>;
	$seq=substr($seq,0,51);
	$qline=substr($seq,0,51);

	print OUT "$title\n$seq\n$plus$qline\n";

    }
        $line=<IN>;
}

close IN;
close OUT;
