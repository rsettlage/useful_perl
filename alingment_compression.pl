#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 9 April 2012
# filter blast results for hits at the ends of contigs
#
#given a blast table, output results that hit within XXX of end of contigs, look on both strands, need to get lengths of blast targets if not present
#
#usage perl blast_filter_for_end_hits.pl <nucmer output> <max gap>
#example perl /groups/DAC/useful_perl/alingment_compression.pl turkey_2.1_turkey_4.0.coords_reformatted.166 10000
#
#####################################################

$file_NAME = $ARGV[0];
$gap_END = $ARGV[1];
$max_overlap = -1000;

print "working on file $file_NAME\n";
open (alignment_FILE, "$file_NAME")|| die "Couldn't open $file_NAME\n";

my $file_OUT = $file_NAME.".compressed_hits";
open (OUTFILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";

my $line=<alignment_FILE>;
chomp($line);
@current_hit_results = split(' ',$line);
$line=<alignment_FILE>;
while ($line) {
	$current_hit_orientation="F";  ###F means 3>4
	if ($current_hit_results[3]<$current_hit_results[4]) {
		$current_hit_orientation="R";
	}
	chomp($line);
	@next_hit_results = split(' ',$line);
	$next_hit_orientation="F";
	if ($next_hit_results[3]<$next_hit_results[4]) {
		#print "$next_hit_results[3] $next_hit_results[4]\n";
		#exit;
		$next_hit_orientation="R";
	}
	$gap_ok="F";
	if ($next_hit_orientation eq $current_hit_orientation) {
		if ($next_hit_orientation eq "F") {
			$gap_length=$current_hit_results[4]-$next_hit_results[3];
			if (($gap_length<$gap_END) && ($gap_length>$max_overlap)) {
				$gap_ok="T";
				}
		}else{
			$gap_length=$next_hit_results[3]-$current_hit_results[4];
			#print "$gap_length, $max_overlap\n";
			if (($gap_length<$gap_END) && ($gap_length>$max_overlap)) {
				$gap_ok="T";
				#print "$gap_ok\n";
			}
		}
	}
	if (($current_hit_results[17] ne $next_hit_results[17]) || ($current_hit_results[18] ne $next_hit_results[18]) || ($next_hit_orientation ne $current_hit_orientation) || ($gap_ok eq "F")) {
		printf OUTFILE "@current_hit_results\n";
		@current_hit_results=@next_hit_results;
	}else{
		$current_hit_results[1]=$next_hit_results[1];
		$current_hit_results[4]=$next_hit_results[4];
		$current_hit_results[6]=$next_hit_results[1]-$current_hit_results[0];
		$current_hit_results[7]=$current_hit_results[7]+$next_hit_results[7];
	}
	$line=<alignment_FILE>;
}	
printf OUTFILE "@next_hit_results\n";
close OUTFILE;
close alignment_FILE;
exit;
