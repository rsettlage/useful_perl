#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 16 Oct 2013
#
#find and output all positions for sequence in fasta file
#
#usage perl find_coordinates_of_sequence.pl <file1> <sequence>
#example perl /groups/DAC/useful_perl/find_coordinates_of_sequence.pl file1.fasta ATC
#
#####################################################

$file_NAME = $ARGV[0];
$target_sequence = $ARGV[1];

print "working on file $file_NAME\n";
open (fasta_FILE, "$file_NAME")|| die "Couldn't open $file_NAME\n";

my $file_OUT = $ARGV[0].".map";
open (OUTFILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";

my $line=<fasta_FILE>;
chomp($line);
my $header="";
my $sequence="";
my $start="";
my $end="";
if ($line =~ /^>/) {
	####this is the header
	$header = $line;
	$sequence="";
	$start=0;
	$end=0;
	print "working on $header\n";
}
$line=<fasta_FILE>;
chomp($line);
while($line) {
		if ($line =~ /^>/) {
			####this is the header
			$header = $line;
			$sequence="";
			$start=0;
			$end=0;
			print "working on $header\n";
			$line=<fasta_FILE>;
			chomp($line);
		}
		while(substr($line,0,1) ne '>'){
			$sequence .= $line;	
			$line=<fasta_FILE>;
			if (!$line) {
				close OUTFILE;
				close fasta_FILE;
				exit;
			}
			chomp($line);
		}
		my $seqLen=length($sequence);
		my $offset = 0;
		my $result = index($sequence, $target_sequence, $offset);
		while ($result != -1) {
			$offset = $result + 1;
			if ($start == 0) {
				printf OUTFILE "-------------------------------$header\n";
				$start=$offset;
				$end=$offset;
			}elsif (($offset-$end)==1) {
				$end=$offset;
			}else{
				printf OUTFILE "$header $seqLen $start $end\n";
				$start=$offset;
				$end=$offset;
			}
			$result = index($sequence, $target_sequence, $offset);
		}
		if ($start>0) {
			printf OUTFILE "$header $start $end\n";
			$start=0;
			$end=0;
		}
}	


close OUTFILE;
close fasta_FILE;
exit;