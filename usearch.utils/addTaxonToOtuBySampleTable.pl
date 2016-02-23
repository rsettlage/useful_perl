use strict;

my $usage = "$0 otu_by_sampleTable assigned_tax.txt, output to STDOUT";
my ($otu_by_sampleTable, $assigned_tax) = @ARGV or die $usage;

open F, $assigned_tax or die "cannot open taxonomy $assigned_tax";
my %tax;
while (<F>)
{
	chomp;
	my ($otu, $tax) = split /\t/;
	$tax =~ tr/ //d; # spaces added in qiime 1.8.0
	$tax =~ s/(;[a-z]__)+$//;
	$tax{$otu} = $tax;
}

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

