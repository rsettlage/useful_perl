use strict;

my $otuMemberListFile = shift or die "usage: $0 otuMemberListFile [minCount default: 0]";
my $minOtuCount = shift;
my %obs;
my %s;
my %otu;
open F, $otuMemberListFile or die "cannot open $otuMemberListFile";

while (<F>)
{
	chomp;
	my ($otu, @list) = split;
	$otu{$otu} = scalar @list;
	foreach my $id (@list)
	{
		$id =~ /(\S+?)_/ or die "cannot parse id $id for otu $otu\n";
		my $sample = $1;
		$obs{$otu}{$sample}++;
		$s{$sample}++;
	}
}
my @sample = sort {$a<=>$b} keys %s;
#print join("\t", "OTU", @sample);
print ("OTU");
foreach my $sample (@sample)
{
	print "\tS$sample";
}
print "\n";
#print "Total";
#foreach my $sample (@sample)
#{
#	print "\t$s{$sample}";
#}
#print "\n";
foreach my $otu (sort {$otu{$b} <=> $otu{$a} or $a cmp $b} keys %otu)
{
	last if $otu{$otu} < $minOtuCount;
	print $otu;
	foreach my $sample (@sample)
	{
		print "\t", int($obs{$otu}{$sample});
	}
	print "\n";
}

