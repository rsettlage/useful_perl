#!/usr/bin/perl -w


###################################################
#this script is splitIlluminaMP_by_linker.pl
#script to automatically remove sequence from fastq files


#this is the actual processing script
#ARGV[0]= file
#ARGV[1]= minimum start position
#ARGV[2]= minimum end position

# use perl /groups/DAC/useful_perl/splitIlluminaMP_by_linker.pl DSilvestris_S1_L001_R2_001.fastq_full_linker.fastq 50 100
###################################################

if(scalar(@ARGV)<3){
	die "use: perl /groups/DAC/useful_perl/splitIlluminaMP_by_linker.pl DSilvestris_S1_L001_R2_001.fastq_full_linker.fastq 40 110 \n";
}


my $fastqfile = $ARGV[0];

open(fastq_FILE1, $fastqfile)|| die "Couldn't open $fastqfile\n";

##current MP linker, full length concatenated and a probable impurity, link1 should be in read1 while link2 should be in read2, might as well look for both and not worry about it
my $full_link="CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG";
my $impure_link1="CTGTCTCTTATACACATCAGATGTGTATAAGAGACAG";
my $impure_link2="CTGTCTCTTATACACATCTGATGTGTATAAGAGACAG";

my $sub_left="GTCTCTTATACACAT";
my $sub_right="GATGTGTATAAGAGAC";

my $base=$fastqfile;
$base=$fastqfile; ##~ s/\.fastq$//;
my $trimmed_out_left=$base . ".clip_link__R1.fastq";  #####<----------------------change this
my $trimmed_out_right=$base . ".clip_link__R2.fastq";  #####<----------------------change this


open(OUTPUT_LEFT, ">$trimmed_out_left") || die "Couldn't open $trimmed_out_left\n";
open(OUTPUT_RIGHT, ">$trimmed_out_right") || die "Couldn't open $trimmed_out_right\n";

print "working on $fastqfile\n";
print "creating $trimmed_out_left and $trimmed_out_right\n";

###############
##time to trim
##############

	my $minLeft = $ARGV[1];
	my $minRight = $ARGV[2] - length($impure_link2);
	my $header = "";
	my $sequence = "";
	my $optional = "";
	my $score = "";
	my $left_sequence = "";
	my $left_score = "";
	my $right_sequence = "";
	my $right_score = "";

	print "min left length is $minLeft\n";
	print "min right length is $minRight\n";
	print "looking for adaptor sequence $full_link $impure_link1 and $impure_link2\n";

	while (<fastq_FILE1>) {
		chomp;
		$header = $_;
		$sequence = <fastq_FILE1>;
		$optional = <fastq_FILE1>;
		$score = <fastq_FILE1>;
		chomp $sequence;
		chomp $optional;
		chomp $score;
		my $start_f = -1;
		$start_f = index $sequence, $full_link;
		if($start_f<0){
			$start_f = index $sequence, $impure_link1;
		}
		if($start_f<0){
			$start_2 = index $sequence, $impure_link2;
		}
		#print "\nstart_f = $start_f\n";
		# print "$sequence\n";
		
		if (($start_f > $minLeft) & ($start_f < $minRight)) {
			$left_sequence = substr $sequence, 0, $start_f;
			$left_score = substr $score, 0, $start_f;
			$right_sequence = substr $sequence, $start_f + length($full_link);
			$right_score = substr $score, $start_f + length($full_link);
			#print "got one\n";
		}
		if (($start_f > $minLeft) & ($start_f < $minRight)) {
			print OUTPUT_LEFT "$header\n";
			print OUTPUT_LEFT "$left_sequence\n";
			print OUTPUT_LEFT "$optional\n";
			print OUTPUT_LEFT "$left_score\n";

			print OUTPUT_RIGHT "$header\n";
			print OUTPUT_RIGHT "$right_sequence\n";
			print OUTPUT_RIGHT "$optional\n";
			print OUTPUT_RIGHT "$right_score\n";
		}
	}

close fastq_FILE1;
close OUTPUT_LEFT;
close OUTPUT_RIGHT;

exit;
