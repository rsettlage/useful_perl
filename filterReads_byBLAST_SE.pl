#!/usr/bin/perl
use strict;
use warnings;

####################################################################################################
#
# Robert E Settlage, Ph.D June 13, 2012
#	filterReads_byBLAST_SE.pl  filters reads using BLAST table file -- right now a hit is considered good enough
#
#	use: perl /groups/DAC/useful_perl/filterReads_byBLAST_SE.pl <reads1.fastq>
#   note: blast file needs to be reads1.fastq.bln.txt
#
#####################################################################################################


if(scalar(@ARGV)<1){
	die "use: perl filterReads_byBLAST_SE.pl <reads1.fastq> \n";
}


my $fastqfile = $ARGV[0];
my $zippedFlag = 0;

if ($fastqfile =~ m/\.gz$/) {
	##`gzip -d -v $fastqfile`;
	##$fastqfile=~s/\.gz$//;
	$zippedFlag = 1;
	open(fastq_FILE1, "gunzip -c $fastqfile |") || die "Couldn't open $fastqfile\n";
	}else{
	open(fastq_FILE1, $fastqfile) || die "Couldn't open $fastqfile\n";
}

my $blastfile = $fastqfile .".HF_GoI.bln.txt.nodup";   #####<----------------------change this

#open(fastq_FILE1, $fastqfile)|| die "Couldn't open $fastqfile\n";
open(blast_FILE1, $blastfile)|| die "Couldn't open $blastfile\n";

my $base=$fastqfile;
##$base=~ s/fastq$//; ##my ($base, $fastq=split(/\./, $fastqfile);  ##this is dangerous, should just do a s/fastq//
my $out_hit=$base . ".GOF_hit.fastq";  #####<----------------------change this
my $out_nohit=$base . ".GOF_nohit.fastq";  #####<----------------------change this
open(OUT_HIT, ">$out_hit") || die "Couldn't open $out_hit\n";
open(OUT_NOHIT, ">$out_nohit") || die "Couldn't open $out_nohit\n";

print "reading blast file $blastfile\n";

##get list of reads that get a hit
my %blast1hash;
my %valid_reads_hash;
my @temp1;
my @temp2;
my $line1=<blast_FILE1>;
while ($line1) {
	@temp1 = split(/\s+/, $line1);
	##print "--$temp1[0]--\n";
	##if (exists $blast1hash{$temp1[0]}) {
		$valid_reads_hash{$temp1[0]} = 1;
	##}
	$line1=<blast_FILE1>;
}
close blast_FILE1;
##now scan through fastq and dump data to hit or nohit
print "scanning read file $fastqfile for hits\n";
my $line2=<fastq_FILE1>;
while($line2){
	chomp $line2;
	@temp1 = split(/ /, $line2);
	$temp1[0] =~ s/^\@//;
	##print "--$temp1[0]--\n";
	if(exists($valid_reads_hash{$temp1[0]})){
		print OUT_HIT $line2."\n";	
		my $seq=<fastq_FILE1>;
		print OUT_HIT $seq;
		my $title2=<fastq_FILE1>;
		print OUT_HIT "+\n";
		my $scores=<fastq_FILE1>;
		print OUT_HIT $scores;
	}else{
		print OUT_NOHIT $line2."\n";	
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

if ($zippedFlag eq 1) {
	`gzip -v $out_hit`;
	`gzip -v $out_nohit`;
}

exit;
