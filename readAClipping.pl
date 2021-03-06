#!/usr/bin/perl -w

###################################################
#this script is readAClipping.pl
#script to remove trailing polyAs, all A's from 3' side
#

#this is the actual processing script
#ARGV[0]= file to be processed, use absolute path

# use perl readAClipping.pl read_file 
# perl /groups/DAC/useful_perl/readAClipping.pl Rhodococcus_3100_7_S2_L001_R1_001.fastq
###################################################

my $fastq_file = $ARGV[0];
my $fastqout_file_EXT = "_noA.fastq";
print "working on $fastq_file: ";

##make sure file is not zipped
my $file_EXT=$fastq_file;
$file_EXT=~ s/\.gz$// ;
if ($file_EXT ne $fastq_file) {
	`gzip -d $fastq_file` ;
	$fastq_file=$file_EXT;
	my $zipped = 1;
}

$temp=$fastq_file;
$temp =~ s/\.fastq$//;
$fastqout_file=$temp.$fastqout_file_EXT;

open (FASTQ, "$fastq_file") || die "Can't open $fastq_file $!\n";
open (OUTPUT, ">$fastqout_file") || die "Can't open $fastqout_file $!\n";

	while (<FASTQ>) {
		chomp;
		my $header = $_;
		chomp $header;
		my $sequence = <FASTQ>;
		chomp $sequence;
		my $optional = <FASTQ>;
		$optional ="+";
		my $score = <FASTQ>;
		chomp $score;

		$temp = $sequence;
		$temp =~ s/A+$// ;
		$tempLen = length($temp);
		$new_score = substr($score, 0, $tempLen);

		if ($tempLen>40) {
			print OUTPUT "$header\n";
			print OUTPUT "$temp\n";
			print OUTPUT "$optional\n";
			print OUTPUT "$new_score\n";
		}


}
close FASTQ;
close OUTPUT;

print " finished\n";
exit;