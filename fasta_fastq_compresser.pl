#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 28 Feb 2012
# fastq/fasta compresser
#
#given a fasta or fastq file, program creates a .read1 and .read2 file,
#
#usage perl fasta_fastq_compresser.pl <switch for fasta|fastq> <file> <common header>
#example perl fasta_fastq_compresser.pl fastq /groups/DAC/bobs_useful_stuff/perl_scripts/test.fq R0174892:8:
#
#####################################################

$file_type = $ARGV[0];
$file_NAME = $ARGV[1];
$common_HEADER = "@".$ARGV[2];

print "using parameters $file_type $file_NAME $common_HEADER for analysis\n";

open (INFILE, "$file_NAME")|| die "Couldn't open $file_NAME\n";

my $file_OUT = $ARGV[0].".compressed";

open (OUTFILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";

if ($file_type eq 'fastq') {
		###lets check if the line count is mod 4 and then assume things are good, if remainder is 2, lets notify the user and allow then to change thier mind, otherwise exit
	$count = `wc -l < $file_NAME`;
	$x = $count % 4;
	if ($x eq 0) {
		print "\n looks like a proper fastq file with $count lines moving on!!\n\n\n";
		print OUTFILE "$common_HEADER\n";
	}elsif($x eq 2){
		print "\n looks like a fasta file, what would you like to do...exiting for now\n\n\n";
	}else{
		print "not sure what this file is, exiting\n\n\n";
		exit;
	}
	while(my $line=<INFILE>) {
		####simplify the header
			chomp($line);
			$header = $line;
			$header =~ s/$common_HEADER//;
			print "header $line changed to $header\n";
			print OUTFILE "$header\n";
		####change the sequence line
			$line=<INFILE>;
			chomp($line);
			$sequence = $line;
			$sequence =~ tr/AGCT/0123/;
			$length_sequence = length ($sequence);
			print "sequence $line changed to $sequence\n";
			$x = $length_sequence % 8;
			print OUTFILE "$sequence\n";
			
		####dont use the optional header
			$line=<INFILE>;
			chomp($line);
			print "optional header $line is removed\n";
		####change the quality score line, assuming true phred
			$line=<INFILE>;
			chomp($line);
			$quality = $line;
			$quality =~ tr[\x20-\x7e][000000000044444444448888888888C];
			print OUTFILE "$quality\n";
			print "quality $line is changed to $quality\n\n";

			$sequence = substr($sequence, 0, 8);
			$quality = substr($quality, 0, 8);
			$dec_qual=eval hex("0x"."$quality");
			$dec_seq=eval hex("0x"."$sequence");

			print "truncated quality is $quality and dec equiv is $dec_qual\n";
			print "truncated sequence is $sequence and dec equiv is $dec_seq\n";
			$hex = dec2hex($dec_qual + $dec_seq);
			print "hex number is: $hex\n";


	}
}
close OUTFILE;
close INFILE;
exit;

sub dec2hex {
    # parameter passed to
    # the subfunction
    my $decnum = $_[0];
    # the final hex number
    my $hexnum="";
    my $tempval="";
    while ($decnum != 0) {
    # get the remainder (modulus function)
    # by dividing by 16
    $tempval = $decnum % 16;
    # convert to the appropriate letter
    # if the value is greater than 9
    if ($tempval > 9) {
    $tempval = chr($tempval + 55);
    }
    # 'concatenate' the number to 
    # what we have so far in what will
    # be the final variable
    $hexnum = $tempval . $hexnum ;
    # new actually divide by 16, and 
    # keep the integer value of the 
    # answer
    $decnum = int($decnum / 16); 
    # if we cant divide by 16, this is the
    # last step
    if ($decnum < 16) {
    # convert to letters again..
    if ($decnum > 9) {
    $decnum = chr($decnum + 55);
    }
    
    # add this onto the final answer.. 
    # reset decnum variable to zero so loop
    # will exit
    $hexnum = $decnum . $hexnum; 
    $decnum = 0 
    }
    }
    return $hexnum;
    } # end sub