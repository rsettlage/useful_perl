#!/usr/bin/perl -w
print "Content-Type: text/html\n\n";

###################################################
#this script is IllxTallyReads.pl
#script to count number of different reads.

#this is the actual processing script
#ARGV[0]= fastq file to parse


# use perl /data/pipeline_in/Runs/useful_perl/IllxTallyReads.pl test_R2.fastq
###################################################

use File::Find;
use File::Copy;
use Data::Dumper;
use warnings;

$fastqout_file=$ARGV[0].".out";
print "starting main script\n";
print "working on data in $ARGV[0]\n";
print "outputting to $fastqout_file\n";

open (FASTQ, "$ARGV[0]") || die "Can't open $ARGV[0]!\n";
open (OUTPUT, ">$fastqout_file") || die "Can't open $fastqout_file!\n";
#Do sort and tally reads in num of qseq files, grab random tiles from each lane
	my %hash;
	my %hashcount;
		while (<FASTQ>) {
			chomp;
			my $header = $_;
			my $sequence = <FASTQ>;
			my $optional = <FASTQ>;
			my $score = <FASTQ>;
			chomp $sequence;
			if ($sequence ne "") {
				#push (@seq, $temp[8]);
				$hash{$sequence}+=1;
			}
		}
	close FASTQ;
	# @sorted_lines =  sort @seq
	for my $key (keys(%hash)) { #this counts the frequence of occurences, ie we saw 10 sequences that appeared on the tile 50 times
		if ($hash{$key}>100) {
			print OUTPUT "$key $hash{$key}\n";
		}
	}
	close OUTPUT;
exit;