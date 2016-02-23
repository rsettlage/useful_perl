#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 15 S
# fasta_fastq_read_parser
#
#given a fasta or fastq file, program creates a .read1 and .read2 file,
#makes no assumptions on which (fasta or fastq) the input file is
#
#usage perl fasta_fastq_read_parser.pl <fasta|fastq>
#example perl fasta_fastq_read_parser.pl F4PCEZU01pairedtrimmedpaired.fa
#
#####################################################

open (FASTQ, $ARGV[0])|| die "Couldn't open $ARGV[0]\n";

my $tfastqoutREAD1 = $ARGV[0].".READ1";
my $tfastqoutREAD2 = $ARGV[0].".READ2";
open (TFASTQ_1, ">$tfastqoutREAD1") || die "Couldn't open $tfastqoutREAD1\n";
open (TFASTQ_2, ">$tfastqoutREAD2") || die "Couldn't open $tfastqoutREAD2\n";

my $line=<FASTQ>;
$read_FLAG = 1;

while($line) {
    chomp($line);

	$read = substr($line, -2, 1);
	#print "\n\n$read\n";
	
	if($read eq "1"){
		$read_FLAG = 1;
	}elsif($read eq "2"){
		$read_FLAG = 2;
	}
	
	if($read_FLAG eq "1"){
		print TFASTQ_1 "$line\n";
	}else{
		print TFASTQ_2 "$line\n";
    }
    $line=<FASTQ>;
}

close TFASTQ_1;
close TFASTQ_2;