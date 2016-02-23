#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2010
#
# template.pl
# This file will be used as a template for other perl scripts.
######################################

use strict;
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $len, $geneNum, $gffFn, $seqFn, $outFn);
$cmd = "$0 @ARGV";

my $minCutoff = 0; ### minimum cutoff value


GetOptions(
	"h|?|help"	=> \$helpFlag,
	"min=i"			=> \$minCutoff,
	"s|sequence=s"	=> \$seqFn,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;

$seqFn = shift if(!defined $seqFn);
$outFn = shift if(!defined $outFn);

if(!defined $seqFn){ print STDERR "\nNeed a sequence file!!\n\n"; help(1); }


if(defined $outFn && ($seqFn eq $outFn))
{ print STDERR " Error) Input and output files are same \n"; exit 1; }


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  -s <sequence file> -o <out file>  [-m <minimum length of contigs>] \n";
	print STDERR "  ex) $0 -s sequence.fa -o out.txt\n\n";
	exit($return);
}

my($name, $buf, $seq, @arrSeq);
my $in = openInput($seqFn);
my $out = openOutput($outFn);
print $out "# $cmd\n";

$buf = <$in>;
while($buf){
	if($buf =~ /^>([\w_+\-:\.\|#\/]+)/){
		$name = $1;		
		$seq = "";
		$buf = <$in>;
		while(defined $buf && $buf !~ /^>/){
			$buf =~ s/[\r\n]+//g;
			$seq .= $buf;
			#@arr = split /\t/, $buf;
			$buf = <$in>;
		}
		my $len = length($seq);
		next if($len < $minCutoff);
		$arrSeq[$#arrSeq+1]{name} = $name;
		$arrSeq[$#arrSeq]{len} = $len;
		
	}
	else{
		$buf = <$in>;
	}
}

if($#arrSeq == -1){
	close($in);
	print STDERR "There is no sequence in the $seqFn file\n";
	exit 1;
}
my @sorted = sort {$b->{len} <=> $a->{len}} @arrSeq;



my ($totLen, $minLen, $maxLen, $midLen, $avgLen);
$totLen = 0;
$minLen = 100000000;
$maxLen = 0;
foreach my $elm (@sorted) {
	$totLen += $elm->{len};
	$minLen = $elm->{len} if($minLen > $elm->{len});
	$maxLen = $elm->{len} if($maxLen < $elm->{len});

	#print $out "$elm->{name}\t$elm->{len}\n";
}

$midLen = $totLen/2;

print $out "total Len : $totLen, num of seq : " . ($#arrSeq+1);
$avgLen = $totLen/($#arrSeq+1);

my $N50 = 0;
$totLen = 0;
foreach my $elm (@sorted) {
	$totLen += $elm->{len};	
	if($totLen > $midLen){
		$N50 = $elm->{len};
		last;
	}
	#print "totLen, $midLen\n";
}

print $out ", Min : $minLen, Max : $maxLen, Avg : $avgLen, N50 : $N50\n\n";

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


