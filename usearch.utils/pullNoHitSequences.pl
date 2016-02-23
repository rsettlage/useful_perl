use strict;
my $usage = "$0 uc_file fasta_file
pull sequences with no hits in uc file from fasta_file, write to stdout
";

my $ucfile= shift or die $usage;
my $fastaFile = shift or die $usage;

open F, $ucfile or die "cannot find uc file $ucfile";
#N       *       *       *       .       *       *       *       20__31970;size=1; 
my %nohitid;
while (<F>)
{
	if (/^N/)
	{
		my @r = split /\t/;
		$nohitid{$r[8]} = 1;
		#print STDERR "$r[8]\n";
	}
}

open F, $fastaFile or die "cannot open fasta file $fastaFile";
my $outputting = 0;
my $numSeqsInFasta=0;
my $numSeqsPrinted=0;
while (<F>)
{
	if (/>(\S+)/)
	{
		$numSeqsInFasta++;
		$outputting = 0;
		if ($nohitid{$1})
		{
			$outputting = 1; 
			$numSeqsPrinted++ if $outputting;
		}
	}
	print if $outputting;
}
print STDERR "Num seqs in $fastaFile: $numSeqsInFasta\nNum seqs printed: $numSeqsPrinted\nNum seqs without hits = ", scalar(keys %nohitid), "\n";
