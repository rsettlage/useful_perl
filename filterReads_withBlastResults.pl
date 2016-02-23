#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2010
#
# filterReads_withBlastResults.pl
# This program filters sequencing reads with blast results
######################################

use strict;
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, @cutoffValues,  $blastFn, $seqFn, $outFn);
$cmd = "$0 @ARGV";


GetOptions(
	"h|?|help"		=> \$helpFlag,
	"cutoff=s" => \@cutoffValues,
	"blast=s"	=> \$blastFn,
	"sequence=s"	=> \$seqFn,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;

$blastFn = shift if(!defined $blastFn);
$seqFn = shift if(!defined $seqFn);
$outFn = shift if(!defined $outFn);

if($#cutoffValues == -1){ print STDERR "\nNeed cutoff values.!!\n\n"; help(1); }
if(!defined $blastFn){ print STDERR "\nNeed a blast table file!!\n\n"; help(1); }
if(!defined $seqFn){ print STDERR "\nNeed a sequence file!!\n\n"; help(1); }
if(!defined $outFn){ print STDERR "\nNeed a output file!!\n\n"; help(1); }


if(defined $outFn && ($blastFn eq $outFn || $seqFn eq $outFn))
{ print STDERR " Error) Input and output files are same \n"; exit; }


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  -b <blast table file>  -s <sequence file>  -c <(cutoff match length):(cutoff percentage)>  -c <(cutoff match length):(cutoff percentage)> ...  -o <out file>\n";
	print STDERR "  ex) $0 -b s_3_1.bln.table -s s_3_1.fastq -c 60:0 -c 50:90 -o s_3_1.filtered.fastq\n\n";
	exit($return);
}



my ($in, $out, @arr, $i, %acceptReads, @lenNum, $matchLen, @cutoffMatchLen, @cutoffPercentages);

for($i = 0; $i <= $#cutoffValues; $i++) {
	@arr = split /:/, $cutoffValues[$i];
	if(!defined $arr[1]){
		 print STDERR "\nWrong format of cutoff value!! : '-c $cutoffValues[$i]'\n\n"; help(1);
	}
	$cutoffMatchLen[$i] = $arr[0];
	$cutoffPercentages[$i] = $arr[1];
}


#### blast table output format.
###[0]query_id [1]db_id [2]percentage of identity [3]aligned_length(identity) [4]mimatches [5]gap_openings [6]q.start [7]q.end [8]s.start [9]s.end [10]e-value [11]hit
$in = openInput($blastFn);
while(<$in>){
	chomp;
	next if(/^#/ || /^\s*$/);
	
	@arr = split /\t/;

	$arr[0] =~ /\//;
	$arr[0] = $`;  ### it takes only 'HWUSI-EAS1737:1:1:1053:16879#0' from 'HWUSI-EAS1737:1:1:1053:16879#0/1'..
	##### to treat paired end reads as one read.
		
	next if(defined $acceptReads{$arr[0]});
	$matchLen = $arr[3]-$arr[4]-$arr[5];

	#### use multiple criteria ...
	for($i = 0; $i <= $#cutoffMatchLen; $i++) {
		if(	$matchLen >= $cutoffMatchLen[$i] && $arr[2] >= $cutoffPercentages[$i] ){
			$acceptReads{$arr[0]} = 1;
			last;
		}
	}
}

close($in);


my($name, $buf);
$in = openInput($seqFn);
$out = openOutput($outFn);

$buf = <$in>;
while($buf){
	if($buf =~ /^>([\w_+\-:\.\|#]+)/ || $buf =~ /^@([\w_+\-:\.\|#]+)/){ #### it does not read '/1' or '/2' from a read name. ex) HWUSI-EAS1737:1:1:1053:16879#0/1
		$name = $1;		

		##### ----------------------
		if(!defined $acceptReads{$name} || $acceptReads{$name} != 1){ $buf = <$in>; next; } ### if you need to remove reads satisfying the criteria, change this line..
		##### ----------------------

		print $out $buf;
		$buf = <$in>;
		while(defined $buf && ($buf !~ /^>/ && $buf !~ /^@/)){
			print $out $buf;
			$buf = <$in>;
		}
	}
	else{
		$buf = <$in>;
	}
}
close($in);
close($out) if(defined $outFn);


sub openInput
{
	my ($fn) = @_;

	return *STDIN unless defined $fn;

	my ($fd);
	open($fd, $fn =~ /\.gz/ ? "zcat $fn|" : ($fn =~ /\.bz2/ ? "bunzip2 -c $fn|" : $fn)) || die "Could not open '$fn' : $!\n";
	return $fd;
}

sub openOutput
{
	my ($fn) = @_;

	return *STDOUT unless defined $fn;

	my ($fd);
	open($fd, $fn =~ /.gz$/ ? "| gzip -c > $fn" : ($fn =~ /\.bz2/ ? "| bzip2 -c > $fn" : ">$fn")) || die "Could not write '$fn' : $!\n";
	return $fd;
}


sub rc
{
	my ($seq, $type) = @_;

	my $rc = reverse $seq;
	if(defined $type && ($type eq "rna" || $type eq "RNA")) # RNA
	{   $rc =~ y/acgtuACGTU/ugcaaUGCAA/;  }
	else ## DNA
	{   $rc =~ y/acgtuACGTU/tgcaaTGCAA/;  }

	return $rc;
}


