#!/usr/bin/perl -w
use strict;

my $file=$ARGV[0];
my $out=$file . ".out.fa";

open(IN, $file);
open(OUT, ">$out")|| die "Couldn't open $out\n";

my $line=<IN>;

while($line){

    if($line=~/>/){

	my $name=$line;
	my $seq="";
	$line=<IN>;

	while($line!~/>/){
	    chomp $line;
	    $line=~s/\r//g;
	    $seq.=$line;
	    $line=<IN>;
	}
	print OUT "$name$seq\n";
#	$line=<IN>;
    }else{
	$line=<IN>;
    }
}
