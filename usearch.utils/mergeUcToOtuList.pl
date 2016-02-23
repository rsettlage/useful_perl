use strict;

my ($uc, $otuList) = @ARGV or die "specify uc_file and otu_list_file to merge";

my $minLength = 240;
my $minPct = 97;

my %otu;
open F, $uc or die "cannot open uc_file $uc";
while (<F>)
{
	chomp;
	my @r=split;
	my $id=$r[8];
	my $len = $r[2];
	my $pct = $r[3];
	my $otu = ($len >= $minLength && $pct >= $minPct) ? $r[9] : 'NA';
	push @{$otu{$otu}}, $id;
}

#for my $otu (sort {scalar @{$otu{$b}} <=> scalar @{$otu{$a}} } keys %otu)
#{
#	print join("\t", $otu, @{$otu{$otu}}), "\n";
#}

open F, $otuList or die "cannot open otu_list_file $otuList";

open OUT, '>'. $uc . '_foldedInto_'. $otuList;
while (<F>)
{
	chomp;
	my ($otu, $rest) = split(/\t/, $_, 2);
	if ($otu{$otu})
	{
		print OUT join("\t", $otu, $rest, @{$otu{$otu}}), "\n";
	}
	else {print OUT;}
}
