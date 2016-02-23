use strict;

my $otuFastaFile = shift or die "specify otu fasta file, relabel seq ids as OTU_0001, OTU_0002, ...";
die "Cannot find $otuFastaFile" unless -f $otuFastaFile;

# first find out how many sequences to see what field width to use
open F, "wc -l $otuFastaFile | ";
$_=<F>;
my ($numOtus)=split;
$numOtus = int($numOtus/2);
print STDERR "Num reads = $numOtus";
my $fw = length("$numOtus");

print STDERR "Field width will be $fw\n";

my $otu = 1;
open F, $otuFastaFile or die "cannot open $otuFastaFile";
while (<F>)
{
	if (/^>/)
	{
		printf ">OTU_%0${fw}d\n", $otu;
		$otu++;
	}
	else
	{
		print;
	}
}
