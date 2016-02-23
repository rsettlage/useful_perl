use strict;

open F, shift or die "specify fasta file";
my %seq;
my $id;
while(<F>)
{
	chomp;
	if (/^>/)
	{
		$id=$_;
	}
	else
	{
		$seq{$id} .= $_;
	}
}
foreach my $id (sort {length($seq{$b})<=>length($seq{$a})} keys %seq)
{
	print "$id\n$seq{$id}\n";
}
