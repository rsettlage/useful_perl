#!/usr/bin/perl -w
use strict;

#####################################################################
#
#
#  convertfilenames.pl - This script takes as arguments an infile and an outfile name.  
#            It renames all unmerged sequences (those with the name contig or s_) as contigs
#            and it renames all merged sequences as "original."
#
#
######################################################################


my $file=$ARGV[0];
chomp $file;
open (IN, $file) || die "Couldn't open $file\n";

my $outfile=$ARGV[1];
chomp $outfile;
open (OUT, ">$outfile") || die "Couldn't open $outfile\n";

my $line=<IN>;

while (<IN>){

    if ($line !~ />/){
	print OUT $line;
	$line=<IN>;
    }else{

	if ($line=~/merged/){

	    my ($nada, $name)=split(/>/, $line);
	    $name = ">original_" . $name;
	    print OUT $name;
	    $line=<IN>;
	}elsif($line=~/contig/){
	    my ($nada1, $name1)=split(/>/, $line);
	    $name1=">contig_" . $name1;
	    print OUT $name1;
	    $line=<IN>;
	}elsif ($line=~/s_/){
	    my ($nada2, $name2)=split(/>/, $line);
	    $name2=">contig_" . $name2;
	    print OUT $name2;
	    $line=<IN>;
	}else{
	    print "*************Missed one     $line*********************************\n";
	    $line=<IN>;
	}
    }
}
