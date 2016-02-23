#!/usr/bin/perl -w
use strict;

my ($infile, $outfile)=@ARGV;
chomp ($infile);

open (IN, $infile)|| die "Can't open $infile\n";
open (OUT, ">$outfile")|| die "Can't open $outfile\n";

my $line=<IN>;
my $contigname="";


while($line){
    if ($line ne "" && $line=~/Seq name/){
	chomp $line;
	my @line=split(/:/, $line);
	$contigname=$line[1];
	$contigname=~s/\r//g;
	$contigname=~s/ //g;
	$line=<IN>;
    }elsif($line ne "" && $line=~/CDS/){
	my @data=split(/\s+/, $line);
	if($data[4]=~/CDS/){
	print OUT "$contigname\tfgenesh\t$data[4]\t$data[5]\t$data[7]\t$data[8]\t$data[2]\t1\tgene_id=$contigname" . "_" . "$data[1]\n";
	}
	$line=<IN>;
    }else{
    	$line=<IN>;
    }
}
