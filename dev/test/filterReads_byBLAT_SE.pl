#!/usr/bin/perl
use strict;
use warnings;

####################################################################################################
#
# Robert E Settlage, Ph.D 17 Sept 2012
#	filterReads_byBLAT_SE.pl  filters reads using BLAT (NOT BLAT) table file -- right now a hit is considered good enough
#
#	use: perl /groups/DAC/useful_perl/filterReads_byBLAT_SE.pl <reads1.fastq>
#   note: BLAT file needs to be reads1.fastq.bln.txt
#
#####################################################################################################


if(scalar(@ARGV)<1){
	die "use: perl filterReads_byBLAT_SE.pl <reads1.fastq> \n";
}

my $fastqfile = $ARGV[0];
my $BLATfile = $fastqfile .".bln.txt";   #####<----------------------change this

open(fastq_FILE1, $fastqfile)|| die "Couldn't open $fastqfile\n";
open(BLAT_FILE1, $BLATfile)|| die "Couldn't open $BLATfile\n";

my $base=$fastqfile;
##$base=~ s/fastq//; ##my ($base, $fastq=split(/\./, $fastqfile);  ##this is dangerous, should just do a s/fastq//
my $out_hit=$base . ".BLAT_hit.fastq";  #####<----------------------change this
my $out_nohit=$base . ".BLAT_nohit.fastq";  #####<----------------------change this
open(OUT_HIT, ">$out_hit") || die "Couldn't open $out_hit\n";
open(OUT_NOHIT, ">$out_nohit") || die "Couldn't open $out_nohit\n";

print "reading BLAT file $BLATfile\n";

##get list of reads that get a hit
my %BLAT1hash;
my %valid_reads_hash;
my @temp1;
my @temp2;
my $line1=<BLAT_FILE1>;
while ($line1) {
	@temp1 = split(/\s+/, $line1);
	##print "--$temp1[0]--\n";
	##if (exists $BLAT1hash{$temp1[0]}) {
		$valid_reads_hash{$temp1[0]} = 1;
	##}
	$line1=<BLAT_FILE1>;
}
close BLAT_FILE1;
##now scan through fastq and dump data to hit or nohit
print "scanning read file $fastqfile for hits\n";
my $line2=<fastq_FILE1>;
while($line2){
	@temp1 = split(/ /, $line2);
	$temp1[0] =~ s/^\@//;
	##print "--$temp1[0]--\n";
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

exit;