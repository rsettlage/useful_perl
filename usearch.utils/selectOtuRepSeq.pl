use strict;

my $usage = "$0 otu_member_list.txt derep.fa
	output to otu_rep_seq.fa
";

my ($otuMemberList, $derepFa) = @ARGV;
die $usage unless ($otuMemberList && $derepFa);

my %otuSeq;
my %seqOtu;
open F, $otuMemberList or die "Cannot open $otuMemberList";
open OUT, ">otu_rep_seq.fa";
while (<F>)
{
	my ($otu, $seqId) = split('\t', $_, 3);
	$seqOtu{$seqId} = $otu;
	#print OUT "$otu\t$seqId\n";
}
printf "Number of otu-seqs found = %d\n", scalar keys %seqOtu;
open F, $derepFa or die "cannot open $derepFa";
my $otu = undef;
while (<F>)
{
	if (/^>(\S+)/)
	{
		my $seqId = $1;
		$seqId =~ s/;size=\d+;//;
		$otu = $seqOtu{$seqId};
		print OUT ">$otu $seqId\n" if $otu;
	}
	elsif($otu)
	{
		print OUT $_;
	}
}
print "Representative otu seqs written to otu_rep_seq.fa\n";
