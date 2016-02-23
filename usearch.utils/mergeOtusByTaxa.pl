use strict;

my $otuBySampleTable = shift or die "specify otuBySampleTable";

my $resolution = shift;
if ($resolution)
{
	die "resolution not between 2 and 6" unless ($resolution >= 2 and $resolution <= 6);
}
my %tbs;
my @s;
my %ttot;

open F, $otuBySampleTable or die "cannot open $otuBySampleTable for reading";
$_=<F>;
chomp;
my ($otu, $taxon, @colHead) = split /\t/;
while (<F>)
{
	next if /^#/;
	my ($otu, $tax, @r) = split /\t/;
	if ($resolution)
	{
		my @t = split(/;/, $tax);
		$tax = join(";", splice(@t,0,$resolution));
	}
	for my $colHead (@colHead)
	{
		my $n = shift @r;
		$tbs{$tax}{$colHead} += $n;
		$ttot{$tax} += $n; 
	}
}

print join("\t", "Taxon", @colHead), "\n";
foreach my $t (sort {$ttot{$b} <=>$ttot{$a}} keys %ttot)
{
	print $t;
	foreach my $colHead (@colHead)
	{ print "\t$tbs{$t}{$colHead}";}
	print "\n";
}
