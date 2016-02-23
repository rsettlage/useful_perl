#!/usr/bin/perl -w


###################################################
#this script is HardClipLength.pl
#script to automatically remove sequence from fastq files


#this is the actual processing script
#ARGV[0]= file
#ARGV[2]= minimum length of sequence to keep
#ARGV[2]= minimum length of sequence to trim, use positive number to indicate final length, use negative number to trim XX bases from tail

# use perl HardClipLength.pl 110526_HWUSI-EAS381R_00025_Stevens_Biedler_Holliday_SS_110526 40 10
###################################################

if(scalar(@ARGV)<2){
	die "use: perl HardClipLength.pl 110526_HWUSI-EAS381R_00025_Stevens_Biedler_Holliday_SS_110526 40 10 \n";
}


my $fastqfile = $ARGV[0];
my $base=$fastqfile;
my $zippedFlag = 0;

if ($fastqfile =~ m/\.gz$/) {
	$zippedFlag = 1;
	open(fastq_FILE1, "gunzip -c $fastqfile |") || die "Couldn't open $fastqfile\n";
	$base=~ s/\.fastq.gz$//;
}else{
	open(fastq_FILE1, $fastqfile) || die "Couldn't open $fastqfile\n";
	$base=~ s/\.fastq$//;
}

my $trimmed_out=$base . "_CT.fastq";  #####<----------------------change this
open(OUTPUT, ">$trimmed_out") || die "Couldn't open $trimmed_out\n";

print "working on $fastqfile\n";
print "creating $trimmed_out\n";

	if ($ARGV[1] == 0) {
		print "looks like trimming/clipping is turned off\n";
		`cp -v $fastqfile $trimmed_out`;
		close fastq_FILE1;
		close OUTPUT;
		exit;
	}

###############
##time to trim
##############

	my $minLength = $ARGV[1];
	my $minClip = $ARGV[2];
	my $header = "";
	my $sequence = "";
	my $optional = "";
	my $score = "";
	my $new_sequence = "";
	my $new_score = "";

	print "min length is $minLength\n";
	print "min to clip is $minClip\n";

	while (<fastq_FILE1>) {
		chomp;
		$header = $_;
		$sequence = <fastq_FILE1>;
		$optional = <fastq_FILE1>;
		$score = <fastq_FILE1>;
		chomp $sequence;
		chomp $optional;
		chomp $score;
		if($minClip>0){
			if (length($sequence) > $minClip) {
				$new_sequence = substr $sequence, 0, $minClip;
				$new_score = substr $score, 0, $minClip;
			}else{
				$new_sequence=$sequence;
				$new_score=$score;
			}
		}else{
			if (length($sequence) + $minClip > 0) {
				$new_sequence = substr $sequence, 0, $minClip;
				$new_score = substr $score, 0, $minClip;
			}else{
				$new_sequence="A";
			}
		}
		if (length($new_sequence) > $minLength) {
			printf OUTPUT "$header\n";
			printf OUTPUT "$new_sequence\n";
			printf OUTPUT "$optional\n";
			printf OUTPUT "$new_score\n";
		}
	}

close fastq_FILE1;
close OUTPUT;

if ($zippedFlag eq 1) {
	`gzip -v $trimmed_out`;
}

exit;