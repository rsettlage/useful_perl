#!/usr/bin/perl -w

###################################################
#this script is AdaptorCounting
#perl AdapterCounting.pl filein fileout
###################################################

$fastq_file = $ARGV[0];
$fastqout_file = "test_count.txt"; #$ARGV[1];

my $start_reads;
my $start_length;
my $avg_start_length;
my $end_reads;
my $end_length;
my $avg_end_length;
#my $adaptor1 = "GATCGGA";
my $adaptor = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

if ($fastq_file =~ m/R2/) {
	$revcom= reverse $adaptor;
	$revcom=~ tr/ACGT/TGCA/;
	$adaptor=$revcom;
}

$adaptor_length = length($adaptor);

open (OUTPUT, ">>$fastqout_file") || die "Can't open $fastqout_file $!\n";
printf OUTPUT "$fastq_file\n";
printf OUTPUT "current_len current_adaptorL current_countL frequencyL current_adaptorR current_countR frequencyR\n";
print "$fastq_file\n";

for ($i=0;$i<$adaptor_length ;$i++) {
	my $current_adaptorL=substr $adaptor, 0, $adaptor_length-$i;
	my $current_countL=0;
	my $current_adaptorR=substr $adaptor, $i;
	my $current_countR=0;
	$current_len = length($current_adaptorL);
	$frequencyL=0;
	$frequencyR=0;
	open (FASTQ, "$fastq_file") || die "Can't open $fastq_file $!\n";
		while (<FASTQ>) {
			my $header = $_;
			my $sequence = <FASTQ>;
			my $optional = <FASTQ>;
			my $score = <FASTQ>;
			chomp $sequence;
			$read_len = length($sequence); # - $adaptor_length;
			$startL = index $sequence, $current_adaptorL;
			$startR = index $sequence, $current_adaptorR;
			if (($startL >= 0) && ($startL < $read_len-1)) {
				$current_countL++;
				$sequenceL=$sequence;
				$sequenceL =~ s/$current_adaptorL/1/g;
				$frequencyL+=($sequenceL =~ tr/1//);;
			}
			if (($startR >= 0) && ($startR < $read_len-1)) {
				$current_countR++;
				$sequenceR=$sequence;
				$sequenceR =~ s/$current_adaptorR/1/g;
				$frequencyR+=($sequenceR =~ tr/1//);;
			}
		}
	close FASTQ;
	printf OUTPUT "$current_len $current_adaptorL $current_countL $frequencyL $current_adaptorR $current_countR $frequencyR\n";
	print "$current_len \t $current_adaptorL \t $current_countL \t $frequencyL \t $current_adaptorR \t $current_countR $frequencyR\n";
}

close OUTPUT;

exit;

