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

my %gtf_contigs;
my $temp = <gtf_FILE>;
while ($temp) {
	chomp $temp;
	@tmp_exon_fields = split(' ', $temp);
	if ($tmp_exon_fields[2] eq 'exon') {
		chomp $tmp_exon_fields[0];
		$temp .= "~~~";
		$gtf_contigs{$tmp_exon_fields[0]} .= $temp;
	} 
	$temp = <gtf_FILE>;
}
print "last exon parsed out is $tmp_exon_fields[0]\n";
print "this is combined as $gtf_contigs{$tmp_exon_fields[0]}\n";

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
			#print "sequence is $sequence\n";
			#exit;
		}
		#print "looking for $fasta_line\n";
		if (defined $gtf_contigs{$fasta_line}) {
			@exon_lines = split('~~~', $gtf_contigs{$fasta_line});
			$i=0;
			while ($exon_lines[$i]) { 
				#print "\n\nfirst line is \n$exon_lines[$i]\n\n";
				my @exon_fields = ();
				@exon_fields = split(" ",$exon_lines[$i]);
				$exon_fields[9] =~ s/\"//g;
				$exon_fields[11] =~ s/\"//g;
				$exon_fields[13] =~ s/\"//g;
				$exon_fields[9] =~ s/;//g;
				$exon_fields[11] =~ s/;//g;
				$exon_fields[13] =~ s/;//g;
				#print "exon 9 is $exon_fields[9] and 11 is $exon_fields[11]\n";
				#print "$fasta_line\n";
				#print "\>gID_$exon_fields[9]_$exon_fields[11]_$exon_fields[13] start $exon_fields[3] end $exon_fields[4] strand $exon_fields[6] $exon_fields[0]\n";
				printf out_FILE "\>gID_$exon_fields[9]_$exon_fields[11]_$exon_fields[13] start $exon_fields[3] end $exon_fields[4] strand $exon_fields[6] $exon_fields[0]\n";
				$exon_start = $exon_fields[3] -1;
				$exon_length = $exon_fields[4] - $exon_fields[3] +1;
				$exon_sequence = substr($sequence, $exon_start, $exon_length);
				printf out_FILE "$exon_sequence\n";
				$i++;
			}
			delete $gtf_contigs{$fasta_line};
		}	
		$fasta_line = $line;
		$sequence = "";
	}
}	
#print exons you didnt find for some manual validation
my $file_OUT_singles = $fasta_NAME . ".notfound_exons.fasta";
open (out_FILE_singles, ">$file_OUT_singles") || die "Couldn't open $file_OUT_singles\n";


foreach $element (keys %gtf_contigs){
	printf out_FILE_singles "$element\n";
}

close out_FILE_singles;
close out_FILE;
close fasta_FILE;
close gtf_FILE;
exit;