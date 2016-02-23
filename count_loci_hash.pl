 #!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 9 April 2012
# count number of loci in a list of transcripts
#
#given a GTF file and a fasta file, extract exons and create a new exon fasta file, keep gene_ID, transcript_ID and exon #
#
#usage perl extract_exons_FASTA_GTF.pl <FASTA> <GTF>
#example perl /groups/DAC/useful_perl/count_loci_hash.pl transcripts.fa
#
#####################################################

$fasta_NAME = $ARGV[0];
$results_NAME=$fasta_NAME.".stats";

open (fasta_FILE, "$fasta_NAME")|| die "Couldn't open $fasta_NAME\n";
print "working on fasta_file $fasta_NAME\n";
open (out_FILE, ">$results_NAME")|| die "Couldn't open $results_NAME\n";
print "outputting data to $results_NAME\n";

print out_FILE "$results_NAME\n";
my @loci_totals;
my $total_loci=0;
my $total_transcripts;
my $current_loci=1;
my $current_count=1;
my $fasta_line = <fasta_FILE>;
while ($fasta_line) {
	chomp $fasta_line;
	if ($fasta_line =~ /^>/) {
		if ($fasta_line =~ /^>CUFF/) {
			$fasta_line =~ m/(\d+)/;
			$tmp_locus_fields[1] = $1;
		}else{
			@tmp_locus_fields = split('_', $fasta_line);
		}
		$total_transcripts += 1;
		if ($current_loci eq $tmp_locus_fields[1]) {
			$current_count+=1;
		}else{
			$loci_totals[$total_loci] = "$current_loci, $current_count";
			$current_count = 1;
			$current_loci = $tmp_locus_fields[1];
			$total_loci +=1;
		}
	} 
	$fasta_line = <fasta_FILE>;
}

$average_transcripts = $total_transcripts/$total_loci;

foreach(@loci_totals){
	 printf out_FILE "$_\n";
}

print out_FILE "$results_NAME\n";
print out_FILE "total loci = $total_loci\n"; 
print out_FILE "total transcripts = $total_transcripts\n"; 
print out_FILE "average transcripts per loci = $average_transcripts\n"; 

close out_FILE;
close fasta_FILE;
exit;