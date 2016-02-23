use strict;

my $usage = "specify tab-delimited table to be normalized by column [h [r [p]]]
	optional h = number of header rows to skip (default 1)
	optional r = number of row_id columns to skip (default 1)
	optional p: if specified, multiply proportions by 100 for percentages.";
my $file = shift or die  $usage;

my $headerRows = (@ARGV > 0) ? $ARGV[0] : 1;
my $rowIdColumns = (@ARGV > 1) ? $ARGV[1] : 1;
my $percentMultiplier = ($ARGV[2] eq 'p') ? 100 : 1;
print STDERR "Number of header rows = $headerRows\nNumber of rowId columns = $rowIdColumns\n";



my @d;
my @colSum;
my @hRow;
my @rowIds;
open F, $file or die "cannot open $file\n";
for my $i (1..$headerRows)
{
	$_ = <F>;
	push @hRow, $_;
}

while (<F>)
{
	chomp;
	my @r = split /\t/;
	my @rids = splice(@r,0,$rowIdColumns);
	push @rowIds, join("\t", @rids);
	my $idx = 0;
	foreach my $cellValue (@r)
	{
		$colSum[$idx] += $cellValue;
		$idx++;
	}
	push @d, \@r;
}

print @hRow;
for my $rowId (@rowIds)
{
	print $rowId;
	my $rowValues = shift @d;
	foreach my $idx (0..$#colSum)
	{
		printf "\t%.6f", $rowValues->[$idx]/$colSum[$idx]*$percentMultiplier;
	}
	print "\n";
}
