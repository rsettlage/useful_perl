#!/usr/bin/perl -w

###################################################
#this script is readNClipping.pl
#script to remove reads that are all or mostly Ns, % acceptable Ns is preset to 30, can change using the pass in value
#

#this is the actual processing script
#ARGV[0]= file to be processed, use absolute path
#ARGV[1]= % of Ns acceptable, 50% is default

# use perl readNClipping.pl read_file <% Ns acceptable --optional>
# perl /groups/DAC/useful_perl/readNClipping.pl Rhodococcus_3100_7_S2_L001_R1_001.fastq
###################################################

my $fastq_file = $ARGV[0];
my $fastqout_file_EXT = "_mN.fastq";
my $minLen = 20;
my $maxFraction = 3;

if (defined $ARGV[1]) {
	$maxFraction = $ARGV[1];
}

##make sure file is not zipped
my $file_EXT=$fastq_file;
$file_EXT=~ s/\.gz$// ;
if ($file_EXT ne $fastq_file) {
	`gzip -d $fastq_file` ;
	$fastq_file=$file_EXT;
	my $zipped = 1;
}

$temp=$fastq_file;
$temp =~ s/\.[^.]*$//;
$fastqout_file=$temp.$fastqout_file_EXT;

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