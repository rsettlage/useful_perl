#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 9 April 2012
# fasta header reformatting to add length "_length"
#
#
#usage perl fasta_header_length.pl <fasta>
#example perl /groups/DAC/useful_perl/fasta_header_length.pl rh.ctg.fasta
#
#####################################################

$file_NAME = $ARGV[0];


print "working on file $file_NAME\n";

open (fastq_FILE, "$file_NAME")|| die "Couldn't open $file_NAME\n";

my $file_OUT = $ARGV[0].".length";

open (OUTFILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";

my $line=<fastq_FILE>;
chomp($line);
$header = $line;
$sequence="";
$first_line="T";
$length_SEQ=0;
while($line) {
		if (substr($line,0,1) eq '>' && $first_line eq 'F') {
			####get the header
			chomp($line);
			
			$header=$header."_".$length_SEQ;
			print OUTFILE "$header\n";
			print OUTFILE "$sequence";
			$header = $line;
			$sequence="";
			$length_SEQ=0;
		}else{
			if ($first_line eq 'T') {
				$line=<fastq_FILE>;
			}
			####get the sequence line
			$sequence .= $line;	
			chomp($line);
			$length_SEQ += length($line);
			$first_line='F';
		}
	$line=<fastq_FILE>;
}	
$header=$header."_".$length_SEQ;
print OUTFILE "$header\n";
print OUTFILE "$sequence";
close OUTFILE;
close fastq_FILE;
exit;