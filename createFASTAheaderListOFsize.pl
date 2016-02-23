#!/usr/bin/perl

use warnings;

#usage
#perl createFASTAheaderListOFsize.pl test.fa
#

$fasta = $ARGV[0];
$fasta_out = $fasta.".list";

print "filtering $fasta\n";
open (FASTA, "$fasta") || die "Can't open $fasta $!\n";
open (FILTERED, ">$fasta_out") || die "Can't open $fasta_out $!\n";

	$sequence = "";	
	$first_line = 0;
	while (<FASTA>) {
		chomp;
		$test_header = substr($_, 0, 1);
		# print "\n$test_header\n";
		if ($test_header eq ">") {
			if ($first_line eq "1") {
				$seq_len = length($sequence); 
				printf FILTERED "$header\t$seq_len\n";
				print "$header\t$seq_len\n";
			}
			else {
				$first_line = 1;
			}
			$header = $_;
			$sequence = "";
		}
		else {
			$sequence =  $sequence.$_;
		}	
	}
	$seq_len = length($sequence); 
	printf FILTERED "$header\t$seq_len\n";
	print "$header\t$seq_len\n";

close FASTA;
close FILTERED;

exit;
