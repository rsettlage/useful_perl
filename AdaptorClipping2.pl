#!/usr/bin/perl -w


###################################################
#this script is AdaptorClipping.pl
#script to automatically remove sequence from fastq files


#this is the actual processing script
#ARGV[0]= file
#ARGV[1]= sequence to trim
#ARGV[2]= minimum length of sequence to keep
#ARGV[2]= minimum length of sequence to trim

# use perl AdaptorClipping.pl 110526_HWUSI-EAS381R_00025_Stevens_Biedler_Holliday_SS_110526  GATCGGAAGAGCGGTTCAGCAGGAATGCCGA 40 10
###################################################

if(scalar(@ARGV)<3){
	die "use: perl AdaptorClipping.pl 110526_HWUSI-EAS381R_00025_Stevens_Biedler_Holliday_SS_110526 GATCGGAAGAG 40 10 \n";
}


my $fastqfile = $ARGV[0];

open(fastq_FILE1, $fastqfile)|| die "Couldn't open $fastqfile\n";

my $base=$fastqfile;
$base=~ s/\.fastq$//;
my $trimmed_out=$base . "_AT.fastq";  #####<----------------------change this
open(OUTPUT, ">$trimmed_out") || die "Couldn't open $trimmed_out\n";

print "working on $fastqfile\n";
print "creating $trimmed_out\n";

###############
##time to trim
##############

	my $adaptor1 = $ARGV[1];
	my $minLength = $ARGV[2];
	my $minClip = $ARGV[3];
	my $header = "";
	my $sequence = "";
	my $optional = "";
	my $score = "";
	my $new_sequence = "";
	my $new_score = "";

	print "min length is $minLength\n";
	print "min to clip is $minClip\n";
	print "looking for adaptor sequence $adaptor1\n";

	while (<fastq_FILE1>) {
		chomp;
		$header = $_;
		$sequence = <fastq_FILE1>;
		$optional = <fastq_FILE1>;
		$score = <fastq_FILE1>;
		chomp $sequence;
		chomp $optional;
		chomp $score;
		my $start = -1;
		$start = index $sequence, $adaptor1;
		# print "\nstart = $start\n";
		# print "$sequence\n";
		if ($start > 0) {
			$new_sequence = substr $sequence, 0, $start;
			$new_score = substr $score, 0, $start;
		}elsif ($start == -1) {
			$new_sequence = substr $sequence, 0, -$minClip;
			$new_score = substr $score, 0, -$minClip;
		}
		if (length($new_sequence) > $minLength) {
			print OUTPUT "$header\n";
			print OUTPUT "$new_sequence\n";
			print OUTPUT "$optional\n";
			print OUTPUT "$new_score\n";
		}
	}

close fastq_FILE1;
close OUTPUT;


exit;
