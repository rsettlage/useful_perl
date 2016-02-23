#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 9 April 2012
# filter blast results for hits at the ends of contigs
#
#given a blast table, output results that hit within XXX of end of contigs, look on both strands, need to get lengths of blast targets if not present
#
#usage perl blast_filter_for_end_hits.pl <blast table> <length from end>
#example perl blast_filter_for_end_hits.pl rh.ctg.fasta.bln.table.out 200
#
# if you have a fasta file without reference length, use fasta_header_length.pl to put the length on the headers, then do blast
# if you already have a blast table with references that do not include the reference length, use fasta_header_length.pl, then grep '>' fasta.fa >references.txt 
# then gawk -F_ '{ print $1, $0 }' references.txt >references_wlength.txt then awk 'NR==FNR { A[$1]=$2; next }; A[$2] { $2=A[$2] } 1' blast.table references_wlength.txt > blast_wLength.txt
#####################################################

$file_NAME = $ARGV[0];
$distance_END = $ARGV[1];

print "working on file $file_NAME\n";
open (blast_FILE, "$file_NAME")|| die "Couldn't open $file_NAME\n";

my $file_OUT = $file_NAME.".end_hits";
open (OUTFILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";

my $line=<blast_FILE>;
chomp($line);
while ($line) {
	chomp($line);
	@hit_results = split(' ',$line);
	@hit_name_wlength = split('_',$hit_results[1]);
	$query_name = $hit_results[0];
	$hit_name = $hit_name_wlength[0];
	$hit_length = $hit_name_wlength[1];
	$distance_from_end = $hit_length - $distance_END;	
	if (($hit_results[8] < $distance_END) || ($hit_results[9] < $distance_END)) {
		printf OUTFILE "$query_name $hit_name $hit_results[6] $hit_results[7] $hit_results[8] $hit_results[9] $hit_length\n";
	}elsif(($hit_results[8] > $distance_from_end) || ($hit_results[9] > $distance_from_end)){
		printf OUTFILE "$query_name $hit_name $hit_results[6] $hit_results[7] $hit_results[8] $hit_results[9] $hit_length\n";
	}
	$line=<blast_FILE>;
}	

close OUTFILE;
close blast_FILE;
exit;