#!/usr/bin/perl
use strict;
use warnings;

####################################################################################################
#
# Robert E Settlage, Ph.D June 13, 2012
#	read_continuity_map.pl  creates a map of what bases are connected -- 
#   will tally the evidence, ie 5 reads show two bases are connected
#   will only look to the right, ie is base in current position connected to the base to the right
#
#	use: perl read_continuity_map.pl <sam.file> <reference.fasta>
#
#####################################################################################################


if(scalar(@ARGV)<2){
	die "use: perl read_continuity_map.pl <reads1.fastq> <reference.fasta>\n";
}

my $samfile = $ARGV[0];
my $fastafile = $ARGV[1];

open(SAM_FILE, $samfile)|| die "Couldn't open $samfile\n";
open(fasta_FILE, $fastafile)|| die "Couldn't open $fastafile\n";

my $base=$fastqfile;
##$base=~ s/fastq//; ##my ($base, $fastq=split(/\./, $fastqfile);  ##this is dangerous, should just do a s/fastq//
my $out_hit=$base . ".blast_hit.fastq";  #####<----------------------change this
my $out_nohit=$base . ".blast_nohit.fastq";  #####<----------------------change this
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