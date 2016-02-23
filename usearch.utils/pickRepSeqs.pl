use strict;

my ($otuListFile, $rawSeqFile) = @ARGV or die "usage: $0 otuListFile rawSeqFile";

my %otu;
open F, $otuListFile or die "cannot open $otuListFile";
<F>; #skip first one ('NA')
while (<F>)
{
	my @r=split(/\t/, $_, 3);
	$otu{$r[1]}=$r[0];
}

my $otu='';
open OUT, ">rep_seq.fa";
open F, $rawSeqFile or die "cannot open $rawSeqFile";
while (<F>)
{
	if (/^>(\S+)/)
	{
		$otu=$otu{$1};
		print OUT ">$otu\n" if $otu;
	}
	else
	{
		print OUT if $otu;
	}
}
