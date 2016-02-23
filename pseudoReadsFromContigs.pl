#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2011
#
# pseudoReadsFromContigs.pl
# It produces psedo reads from contig sequences.
######################################

use strict;
use warnings "all";
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $verboseFlag, $auto, $outFn);
my ($readLen, $stepSize, $copyNum, $contigFn);
$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..

$readLen = 1000;
$stepSize = 100;
$copyNum = 10;
$auto = "on";

GetOptions(
	"h|?|help"	=> \$helpFlag,
	"verbose"		=> \$verboseFlag,
	"auto"		=> \$auto,

	"len=i"	=> \$readLen,
	"step=i"	=> \$stepSize,
	"copyNum=s"	=> \$copyNum,
	"input=s"	=> \$contigFn,
	
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;


################################################

$contigFn = shift if(!defined $contigFn);

if(!defined $contigFn){ print STDERR "\nNeed a contig sequence!!\n\n"; help(1); }

if(defined $outFn && ($contigFn eq $outFn))
{ print STDERR " Error) Input and output files are same \n"; exit; }


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  -i <contig sequence file>  -o <out file>\n";

	print STDERR "  --------- other options (can be skipped) ---------\n";
	print STDERR "        [ -l <length of reads to be generated (default : $readLen)> ]\n";
	print STDERR "        [ -s <step size (default : $stepSize)> : if it is 100, the first cutting position is 100*k (k is k-th copy of a contig)]\n";
	print STDERR "        [ -c <copy number of a contig (default : $copyNum)> ]\n";
	print STDERR "        [ -auto <on | off (defalut : $auto)> : \n";
	print STDERR "              if a contig name has coverage information, the COPY NUMBER of the contig is estimated from the information.\n";
	print STDERR "              ex) 454 contig headers have 'numreads=NN'. coyp_number = NN*300/length of contig]\n";
	print STDERR "              ex) velvet contig headers have 'cov_NN']\n";

	print STDERR "  --------------------------------------------------\n\n";

	print STDERR "  ex) $0 -i contig.fa -o pseudoReads.fa -auto off -l 2000 -s 200 -c 20\n\n";
	exit($return);
}

################################################

my ($in, $out, @arr, $cid, $i, $k, $j, $pos, $seq, $sLen, $buf, $name, $rName);
$in = openInput($contigFn);
$out = openOutput($outFn);
#print $out "# $cmd\n";

$cid = 1;
$buf = <$in>;
while($buf){
	if($buf =~ /^>/){
		@arr = split /[\t\s\r\n]+/, $';
		$name = $arr[0];

		$seq = "";
		$buf = <$in>;
		while(defined $buf && $buf !~ /^>/){
			$buf =~ s/[\r\n]+//g;
			$seq .= $buf;
			$buf = <$in>;
		}
		
		$sLen = length($seq);

		if($auto eq "on"){
			if(defined $arr[2] && $arr[2] =~ /numreads=(\d+)/){
				$copyNum = $1*300/$sLen;
			}
			elsif($arr[0] =~ /cov_(\d+)/){
				$copyNum = $1+1;
			}
		}
		$copyNum = 3 if($copyNum < 3); ## the copy number should be bigger than 3		

		$name = $` if($name =~ /[:_]/);

		##### generate reads to forward direction
		for($i = 1, $j = 0; $i <= $copyNum/2; $i++, $j += $stepSize){
			if($j >= $sLen){
				$j = 0;
			}
			elsif($j > $readLen){
				$j = $j%$readLen;
			}


			$k = 1;
			$pos = $j;
			if($pos >= $stepSize){
				$rName = $name . "_$cid" . "_$i" . "_$k";				
				print $out ">$rName\n" . substr($seq, 0, $pos) . "\n";
				$k++;
			}
			for(; $pos < $sLen; $k++, $pos += $readLen){
				$rName = $name . "_$cid" . "_$i" . "_$k";					
				if($sLen - $pos >= $stepSize){
					print $out ">$rName\n" . substr($seq, $pos, $readLen) . "\n";
				}
			}
		}

		##### generate reads to reverse direction
		for($j = $sLen - $stepSize; $i <= $copyNum; $i++, $j -= $stepSize){
			if($j < 0){
				$j =  $sLen - $stepSize;
			}
			elsif($j < $sLen - $readLen){
				$j =  $sLen - $stepSize - $j%$readLen;;
			}

			$k = 1;
			$pos = $j;
			for(; ; $pos -= $readLen){
				$rName = $name . "_$cid" . "_$i" . "_$k";
				if($sLen - $pos >= $stepSize){
					print $out ">$rName\n" . substr($seq, $pos, $readLen) . "\n";
				}
				$k++;
				last if($pos < $readLen);
			}
			
			if($pos > $stepSize){
				$rName = $name . "_$cid" . "_$i" . "_$k";				
				print $out ">$rName\n" . substr($seq, 0, $pos) . "\n";
				$k++;
			}
		}
		$cid++;
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
