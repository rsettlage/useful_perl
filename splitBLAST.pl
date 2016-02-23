#!/usr/bin/perl -w

###################################################
#this script is splitBLAST.pl
#script to split a query file up into chunks to blast
#

#this is the actual processing script
#ARGV[0]= file to be processed, use absolute path
#ARGV[1]= % of Ns acceptable, 50% is default

# use perl splitBLAST.pl <partitions> <read_file> <blast options>
# perl /groups/DAC/useful_perl/splitBLAST.pl zz_test.fasta "blastn -task blastn -outfmt 6 -ungapped -num_threads 10 -evalue 1e-6 -num_descriptions 1"
###################################################

my $fasta_FILE = $ARGV[0];
my $prefix = "split_";
my $partitions = $ARGV[1]; ##should be 4-10x num of cores being used
my $minSequences = 25000;

if (!defined $ARGV[1]) {
	print "not enought args, use perl splitBLAST.pl <partitions> <read_file> <blast options>\n";
}

##get line count in file
my $line_count = `wc -l < $fasta_FILE` ;
my $temp = $line_count/$partitions/2;
if ($temp<$minSequences) {
	$partitions=$line_count/$minSequences;  ###need to make this even
}

open (FASTQ, "$fastq_file") || die "Can't open $fastq_file $!\n";
open (OUTPUT, ">$fastqout_file") || die "Can't open $fastqout_file $!\n";

	while (<FASTQ>) {
		chomp;
		my $header = $_;
		my $sequence = <FASTQ>;
		my $optional = <FASTQ>;
		my $score = <FASTQ>;
		chomp $sequence;
		$optional ="+";
		chomp $score;
		$temp = $sequence;
		$temp =~ s/N//g ;
		$tempLen = length($temp);
		$seqLenFraction = length($sequence) * $maxFraction;
		if ($temp ne $sequence) {
			if (($tempLen > $minLen) && ($tempLen > $seqLenFraction)) {
				$temp = $sequence;
				$temp =~ s/N+$// ;
				$tempLen = length($temp);
				$new_score = substr($score, 0, $tempLen);
				$temp =~ s/^N+// ;
				$tempLen = length($temp);
				$new_score =substr($score,-$tempLen); 
				print OUTPUT "$header\n";
				print OUTPUT "$temp\n";
				print OUTPUT "$optional\n";
				print OUTPUT "$new_score\n";
			}
		}else {
			print OUTPUT "$header\n";
			print OUTPUT "$sequence\n";
			print OUTPUT "$optional\n";
			print OUTPUT "$score\n";
		}
}
close FASTQ;
close OUTPUT;

print "$fastq_file finished\n";
exit;