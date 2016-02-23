use strict;

my $usage = "$0 otu_by_sampleTable assigned_tax.txt seqID_to_otuID.list, output to STDOUT";
my ($otu_by_sampleTable, $assigned_tax, $seqID_to_otuID) = @ARGV or die $usage;

my %seqToOtu;
open F, $seqID_to_otuID or die "cannot open $seqID_to_otuID";
while (<F>)
{
	chomp;
	my ($seq, $otu) = split;
	$seqToOtu{$seq}=$otu;
}
printf STDERR "Number of otus read = %d\n", scalar keys %seqToOtu;

open F, $assigned_tax or die "cannot open taxonomy $assigned_tax";
my %tax;
while (<F>)
{
	chomp;
	my ($seq, $tax) = split;
	$tax{$seqToOtu{$seq}} = $tax;
}
printf STDERR "Number of taxa mapped to seqs to otus: %d\n", scalar keys %tax;

open F, $otu_by_sampleTable or die "cannot open otu_by_sample_table $otu_by_sampleTable";
$_=<F>;
my ($id, @r)= split /\t/;
print join("\t", $id, "Taxon", @r);
while (<F>)
{
	my ($id, @r)=split /\t/;
	my $taxon = $tax{$id};
	$taxon = "None" unless $taxon;
	print join("\t", $id, $taxon, @r);
}

