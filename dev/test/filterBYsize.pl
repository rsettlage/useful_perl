#!/usr/bin/perl

use warnings;

#usage
#perl filterBYsize.pl test.fa 500 
#

$min_length = $ARGV[1];
$fasta = $ARGV[0];
$fasta_out = $fasta.".filtered.fa";
$filtered_out = $fasta.".filtered_short.fa";

print "filtering $fasta\n";
open (FASTA, "$fasta") || die "Can't open $fasta $!\n";
open (FILTERED, ">$fasta_out") || die "Can't open $fasta_out $!\n";
open (FILTERED_OUT, ">$filtered_out") || die "Can't open $filtered_out $!\n";

	$sequence = "";	
	$first_line = 0;
	while (<FASTA>) {
		chomp;
		$test_header = substr($_, 0, 1);
		# print "\n$test_header\n";
		if ($test_header eq ">") {
			if ($first_line eq "1") {
				do_filter();
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
close FASTA;
close FILTERED;
close FILTERED_OUT;

exit;

sub do_filter {
	$seq_len = length($sequence); 
        # print "\n$header\n";
        # print "$sequence\n";
        # print "$seq_len";
        if ($seq_len > $min_length) {
                   printf FILTERED "$header\n";
                   printf FILTERED "$sequence\n";
                   }
        else {
                   printf FILTERED_OUT "$header\n";
                   printf FILTERED_OUT "$sequence\n";
        }
        # print "\ncontinue?\n";
        # $answer=<STDIN>;
        # chomp $answer;
        # if (($answer eq "no") || ($answer eq "n")) {
        #           print "exiting due to user input\n";
        #           exit;
        # }
}

