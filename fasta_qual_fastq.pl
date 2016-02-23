 #!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 9 April 2012
# combine fasta/qual to fastq                                       NOT FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#
#usage perl fasta_qual_fastq.pl <FASTA> <QUAL> <FASTQ>
#example perl /groups/DAC/useful_perl/fasta_qual_fastq.pl 454Trimmed.paired.2.fasta 454Trimmed.paired.2.qual 454Trimmed.paired.2.fastq
#
#####################################################

$fasta_NAME = $ARGV[0];
$qual_NAME = $ARGV[1];
$fastq_NAME= $ARGV[2];


open (fasta_FILE, "$fasta_NAME")|| die "Couldn't open $fasta_NAME\n";
print "working on fasta_file $fasta_NAME\n";
open (qual_FILE, "$qual_NAME")|| die "Couldn't open $qual_NAME\n";
print "working on fasta_file $qual_NAME\n";
open (out_FILE, ">$fastq_NAME")|| die "Couldn't open $fastq_NAME\n";
print "outputting data to $fastq_NAME\n";

print out_FILE "$fastq_NAME\n";

my @fasta_sequence;
my $fasta_line = <fasta_FILE>;
while ($fasta_line) {
	chomp $fasta_line;
	if ($fasta_line =~ /^clr/) {
		my $header=substr($fasta_line,0,15);
	}else{
		$fasta_sequence[$header].=$fasta_line;
	}
	$fasta_line = <fasta_FILE>;
}

my @qual_values;
my $qual_line = <qual_FILE>;
while ($qual_line) {
	chomp $qual_line;
	if ($qual_line =~ /^clr/) {
		my $header=substr($qual_line,0,15);
	}else{
		$qual_value[$header].=$qual_line;
	}
	$qual_line = <qual_FILE>;
}

foreach(@fasta_sequence){
	 printf out_FILE "$_\n";
}

print out_FILE "$results_NAME\n";
print out_FILE "total loci = $total_loci\n"; 
print out_FILE "total transcripts = $total_transcripts\n"; 
print out_FILE "average transcripts per loci = $average_transcripts\n"; 

close out_FILE;
close fasta_FILE;
exit;