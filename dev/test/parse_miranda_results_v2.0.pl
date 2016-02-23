#!/usr/bin/perl -w
use strict;

my $file=$ARGV[0];
open(IN, $file)|| die "couldn't open $file\n";
my $outfile=$file . ".2.out";
open(OUT, ">$outfile")|| die "Couldn't open $outfile\n";

my $line=<IN>;
while($line){

    if($line=~/>>/){

	print OUT $line;
	$line=<IN>;
    }else{
	$line=<IN>;
    }
}
