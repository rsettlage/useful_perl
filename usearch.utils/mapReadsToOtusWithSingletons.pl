use strict;

my $usage = "$0 identicalSeqIds derep_uc <singleton_uc optional>
   for each OTU, list each un-deprep seq_id (sample embeded)
   output to otu_member_list.txt
";
my ($identicalSeqIds, $derep_uc, $singleton_uc) = @ARGV;
die $usage unless ($identicalSeqIds && $derep_uc);

-f $derep_uc or die "cannot find uc file $derep_uc";
open F, $derep_uc or die "cannot open $derep_uc";
#H	2180	253	97.6	+	0	0	253M	4_24_2013__281240;size=503109;	5_4_2013__115003;size=2;
#H	2180	253	97.2	+	0	0	253M	4_24_2013__53;size=475403;	5_4_2013__115003;size=2;
#H	2180	253	98.0	+	0	0	253M	4_24_2013__468;size=616513;	5_4_2013__115003;size=2;
#or
#H	5	252	100.0	+	0	0	252M	4_24_2013__23;size=91600;	OTU_0006
#H	2	253	99.6	+	0	0	253M	4_24_2013__406;size=216833;	OTU_0003

my %derep2otu;
my %otuCount;
my $notu = 0;
my %notu;
while (<F>)
{
	#next if /^N/; # no hit
	chomp;
	my ($h, undef, $len, $pct, $orient, undef, undef, $mlen, $seqId, $otuId) =split /\t/;
	die "problem: no seqID:\n$_\n" unless $seqId;
	if ($otuId eq '*')
	{
		$notu++;
		$otuId = sprintf "NOTU_%04d", $notu;
	}
	$seqId =~ s/;size=\d.*//;
	$derep2otu{$seqId}=$otuId;
	$otuCount{$otuId}++;
}
printf STDERR "Number of distinct seqIDs read from $derep_uc = %d\n", scalar keys %derep2otu;
printf STDERR "Number of distinct OTUs = %d\n", scalar keys %otuCount;
printf STDERR "Highest NOTU otu is $notu\n";


my %otulist;
open F, $identicalSeqIds or die "cannot open $identicalSeqIds";
while (<F>)
{
	chomp;
	my @ids = split;
	my $otu = '';
	foreach my $id (@ids)
	{
		if (exists ($derep2otu{$id}))
		{
			$otu = $derep2otu{$id};
			last;
		} 
	}
	#die("Not otu: n". substr($_, 0, 70). "\n") unless $otu;
	unless ($otu)
	{
		$notu++;
		$otu = sprintf("NOTU_%04d", $notu);
	}
	push @{$otulist{$otu}}, @ids;
}
print STDERR "Number without otu = $notu\n";
printf STDERR "Number of OTUs with lists: = %d\n", scalar keys %otulist;
if ($singleton_uc)
{
	open F, $singleton_uc or die "cannot open $singleton_uc";
	while (<F>)
	{
		next if /^N/; # no hit
		chomp;
		my ($h, undef, $len, $pct, $orient, undef, undef, $mlen, $seqId, $otuId) =split /\t/;
		$seqId =~ s/;size=\d.*//;
		push @{$otulist{$otuId}}, $seqId;
	}
}

open F, ">otu_member_list_with_singletons.txt";
foreach my $otu (sort {scalar(@{$otulist{$b}}) <=> scalar(@{$otulist{$a}})} keys %otulist)
{
	print F join("\t", $otu, @{$otulist{$otu}}), "\n";
}

