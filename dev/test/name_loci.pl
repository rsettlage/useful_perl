#!/usr/bin/perl -w
use strict;

open(BLAST, $ARGV[0]);
#open(QUERY, $ARGV[1]);
#open(SUBJECT,$ARGV[2]);

my $outfile=$ARGV[0] . "_named";

open (OUT, ">$outfile");

# First, we will determine, for each unique locus, the total number of unique hits (evalue <-10) there were.

#Then we will determine whether there is enough data to name this locus this gene name.

my %loci2;
my %loci;

my $line=<BLAST>;

while($line){
my @data=split(/\t/, $line);
my $qu=$data[0];
my @qudata=split('_', $qu);
my $locus=$qudata[0] . "_" . $qudata[1] . "_" . $qudata[2];


if(!exists($loci{$locus})){
    $loci{$locus}=join("\t",$data[1],@data[1..11]);
}else{
#    print "$locus\t$data[1]\n";
}

my $locus2=$qudata[0] . "_" . $qudata[1] . "_" . $qudata[2];

unless(exists($loci2{$locus2})){
    $loci2{$locus2}=1;}

$line=<BLAST>;
}

my @array=keys(%loci2);
my @array2=keys(%loci);
print "$ARGV[0]\t$#array2\t$#array\n";

while(my($k, $v)=each(%loci)){
    print OUT "$k\t$v\n";
}
