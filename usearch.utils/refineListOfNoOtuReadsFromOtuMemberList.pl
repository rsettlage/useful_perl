use strict;

#assume first line of otu member list table is 'NA', the ones we want, but that some of these have been actually re-assigned to later OTUs. Find just the ones that aren't re-assigned and spit them out.

open F, "nootu_v_repseq97.uc_foldedInto_otu_member_list.txt";
$_=<F>;
chomp;
my @r=split;
shift @r;
my %h;
foreach my $id (@r)
{
	$h{$id}=1;
}
while (<F>)
{
	chomp;
	@r=split;
	foreach my $id (@r)
	{
		delete $h{$id};
	}
}
print join("\n", sort keys %h);
