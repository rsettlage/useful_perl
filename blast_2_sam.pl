#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 20 Nov 2012
# convert blast to sam, also create table of 'clipped' ends
#
#
#usage perl blast_2_sam.pl <blast table>
#example perl /groups/DAC/useful_perl/blast_2_sam.pl Para_SD1.it_3.sam
#
#####################################################

##assuming blast table format:
##qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore

$file_NAME = $ARGV[0];
$distance_END = $ARGV[1];

print "working on file $file_NAME\n";
open (blast_FILE, "$file_NAME")|| die "Couldn't open $file_NAME\n";

my $file_OUT = $file_NAME.".clipped";
open (OUTFILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";

my %ref_hash;

my $line=<blast_FILE>;
while ($line) {
	chomp($line);
	@hit_results = split(' ',$line);
	$query_name = $hit_results[0];
	$hit_name = $hit_results[1];
	$hit_pid = $hit_results[2];
	$qlen = $hit_results[3];
	$slen = $hit_results[4];
	$match_len = $hit_results[5];
	$mismatches = $hit_results[6];
	$gaps = $hit_results[7];
	$query_start = $hit_results[8];
	$query_end = $hit_results[9];
	$hit_start = $hit_results[10];
	$hit_end = $hit_results[11];
	$eval = $hit_results[12];

	if ($query_start>10) {
		printf OUTFILE "$hit_start\n";
		$ref_hash{$hit_start} ++;
	}
	if ($query_end<($qlen-10)){
		printf OUTFILE "$hit_end\n";
		$ref_hash{$hit_end} ++;
	}
	$line=<blast_FILE>;
}	

printf OUTFILE "summary time:\n";
foreach $element (keys %ref_hash){
	printf OUTFILE "$element\n";
}

close OUTFILE;
close blast_FILE;
exit;