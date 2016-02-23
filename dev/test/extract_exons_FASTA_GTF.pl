#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 9 April 2012
# extract exons from a fasta file and create another fasta file, do not assemble into transcripts
#
#given a GTF file and a fasta file, extract exons and create a new exon fasta file, keep gene_ID, transcript_ID and exon #
#
#usage perl extract_exons_FASTA_GTF.pl <FASTA> <GTF>
#example perl /groups/DAC/useful_perl/extract_exons_FASTA_GTF.pl rh.ctg.fasta merged.gtf
#
#####################################################

$fasta_NAME = $ARGV[0];
$gtf_NAME = $ARGV[1];

open (fasta_FILE, "$fasta_NAME")|| die "Couldn't open $fasta_NAME\n";
print "working on fasta_file $fasta_NAME\n";
open (gtf_FILE, "$gtf_NAME")|| die "Couldn't open $gtf_NAME\n";
print "working on gtf_file $gtf_NAME\n";

my $file_OUT = $fasta_NAME . ".exons.fasta";
open (out_FILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";
print out_FILE "creating $file_OUT for results\n";

while (<gtf_FILE>) {
	@exon_fields = split;
	if ($exon_fields[2] eq 'exon') {
		push @exon_lines, [ @exon_fields ];
	}  
}
$total_exons = @exon_lines;
$remaining_exons = $total_exons;
print "$total_exons exons to process\n";
$fasta_line = <fasta_FILE>;
while ($fasta_line) {
	chomp($fasta_line);
	if ($fasta_line =~ /^>/) {
		$fasta_line =~ s/^>//;
		$line = <fasta_FILE>;
		while (($line) && ($line !~ /^>/)) {
			chomp ($line);
			$sequence .= $line;
			$line = <fasta_FILE>;
		}
		$i = 0;
		while ($exon_lines[$i]) {   #need to go all the way through as several exons can come from a single fasta sequence
			#print "comparing $fasta_line with $exon_lines[$i][0]\n";
			if ($fasta_line eq $exon_lines[$i][0]) {  #ok, found a match, now pull out the exon and print it
				$exon_lines[$i][11] =~ s/\"//g;
				$exon_lines[$i][13] =~ s/\"//g;
				$exon_lines[$i][11] =~ s/;//g;
				$exon_lines[$i][13] =~ s/;//g;
				#print "$fasta_line\n";
				#print "\>$exon_lines[$i][11]_$exon_lines[$i][13] $exon_lines[$i][8] $exon_lines[$i][9] $exon_lines[$i][16] $exon_lines[$i][17] $exon_lines[$i][3] $exon_lines[$i][4] $exon_lines[$i][6] $exon_lines[$i][0]\n";
				printf out_FILE "\>$exon_lines[$i][11]_$exon_lines[$i][13] $exon_lines[$i][8] $exon_lines[$i][9] $exon_lines[$i][16] $exon_lines[$i][17] $exon_lines[$i][3] $exon_lines[$i][4] $exon_lines[$i][6] $exon_lines[$i][0]\n";
				$exon_start = $exon_lines[$i][3] -1;
				$exon_length = $exon_lines[$i][4] - $exon_lines[$i][3] +1;
				$exon_sequence = substr($sequence, $exon_start, $exon_length);
				printf out_FILE "$exon_sequence\n";
				#print "$sequence\n";
				splice @exon_lines, $i, 1;
			}
			$i++;
		}
		$remaining_exons = @exon_lines;
		$fasta_ref++;
		print "processed $fasta_ref references\n";
		print "$remaining_exons exons remaining out of $total_exons\n";
	}	
	$fasta_line = $line;
}	
#print exons you didnt find for some manual validation
my $file_OUT_singles = $fasta_NAME . ".notfound_exons.fasta";
open (out_FILE_singles, ">$file_OUT_singles") || die "Couldn't open $file_OUT_singles\n";
$i = 0;
	while ($exon_lines[$i]) {
			printf out_FILE_singles "$exon_lines[$i][0]\n";
			$i++;
	}

close out_FILE_singles;
close out_FILE;
close fasta_FILE;
close gtf_FILE;
exit;