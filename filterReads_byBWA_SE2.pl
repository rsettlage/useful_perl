#!/usr/bin/perl
use strict;
use warnings;

####################################################################################################
#
# Robert E Settlage, Ph.D June 13, 2012
#	filterReads_byBWA_SE.pl  filters reads using SAM file -- right now a hit is considered good enough to toss
#
#	use: perl filterReads_byBWA_SE.pl <reads1.fastq> <samfile>
#   note: sam file needs to be readfile.sam
#	this version does not use the sam header to validate a read.
#####################################################################################################


if(scalar(@ARGV)<1){
	die "use: perl filterReads_byBWA_SE.pl <reads1.fastq> \n";
}

my $fastqfile = $ARGV[0];
my $samfile = $ARGV[1];
my $zipped=0;

##make sure file is not zipped
my $file_EXT=$fastqfile;
$file_EXT=~ s/\.gz$// ;

if ($file_EXT ne $fastqfile) {
	`gzip -d $fastqfile` ;
	$fastqfile=$file_EXT;
	$samfile = $fastqfile.".sam";
	$zipped=1;
}

open(fastq_FILE1, $fastqfile)|| die "Couldn't open $fastqfile\n";
open(sam_FILE1, $samfile)|| die "Couldn't open $samfile\n";


my $base=$fastqfile;
##$base=~ s/fastq$//; ##my ($base, $fastq=split(/\./, $fastqfile);  ##this is dangerous, should just do a s/fastq//
my $out_hit=$base .$ARGV[1]. ".sam_hit.fastq";
my $out_nohit=$base .$ARGV[1]. ".sam_nohit.fastq";
open(OUT_HIT, ">$out_hit") || die "Couldn't open $out_hit\n";
open(OUT_NOHIT, ">$out_nohit") || die "Couldn't open $out_nohit\n";

print "reading sam file $samfile\n";

##get list of reads that get a hit
my %sam1hash;
my %valid_reads_hash;
my @temp1;
my @temp2;
my $line1=<sam_FILE1>;
while ($line1=~/^@/) {
	$line1=<sam_FILE1>;
}
while ($line1) {
	@temp1 = split(/\s+/, $line1);
	$valid_reads_hash{$temp1[0]} = 1;
	$line1=<sam_FILE1>;
}
close sam_FILE1;

##now scan through fastq and dump data to hit or nohit
print "scanning read file $fastqfile for hits\n";
my $line2=<fastq_FILE1>;
while($line2){
	@temp1 = split(/ /, $line2);
	$temp1[0] =~ s/^\@//;
	if(exists($valid_reads_hash{$temp1[0]})){
		print OUT_HIT $line2;	
		my $seq=<fastq_FILE1>;
		print OUT_HIT $seq;
		my $title2=<fastq_FILE1>;
		print OUT_HIT "+\n";
		my $scores=<fastq_FILE1>;
		print OUT_HIT $scores;
	}else{
		print OUT_NOHIT $line2;	
		my $seq=<fastq_FILE1>;
		print OUT_NOHIT $seq;
		my $title2=<fastq_FILE1>;
		print OUT_NOHIT "+\n";
		my $scores=<fastq_FILE1>;
		print OUT_NOHIT $scores;
	}
	$line2=<fastq_FILE1>;
}
close fastq_FILE1;
close OUT_HIT;
close OUT_NOHIT;

##dont forget to rezip
if ($zipped eq 1) {
	##`gzip $fastqfile` ;
}

exit;