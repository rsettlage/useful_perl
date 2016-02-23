use strict;

# read uc file, then read derep seq file, then read underep seq file
# for each un-reprep seq, print seq_id (sample embeded) and OTU

my ($ucfile, $identicalSeqIds) = @ARGV;
die "usage: $0 ucfile identicalSeqIds, output to otu_member_list.txt" unless ($ucfile && $identicalSeqIds);
-f $ucfile or die "cannot find uc file $ucfile";

open F, $ucfile or die "cannot open $ucfile";
#H	2180	253	97.6	+	0	0	253M	4_24_2013__281240;size=503109;	5_4_2013__115003;size=2;
#H	2180	253	97.2	+	0	0	253M	4_24_2013__53;size=475403;	5_4_2013__115003;size=2;
#H	2180	253	98.0	+	0	0	253M	4_24_2013__468;size=616513;	5_4_2013__115003;size=2;
#or
#H	5	252	100.0	+	0	0	252M	4_24_2013__23;size=91600;	OTU_0006
#H	2	253	99.6	+	0	0	253M	4_24_2013__406;size=216833;	OTU_0003

my $lineCount=0;
my %derep2otu;
my %otuCount;
my $notu = 0;
my %notu;
while (<F>)
{
	$lineCount++;
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
printf STDERR "Number of distinct seqIDs read from $ucfile = %d, numLines=$lineCount\n", scalar keys %derep2otu;
printf STDERR "Number of distinct OTUs = %d\n", scalar keys %otuCount;
printf STDERR "Highest NOTU otu is $notu\n";

open OUT, ">seqID_to_otuID.list";
foreach  my $id (keys %derep2otu)
{
	print OUT "$id\t$derep2otu{$id}\n";
}
close OUT;

my %otulist;
open F, $identicalSeqIds or die "cannot open $identicalSeqIds";
while (<F>)
{
	chomp;
	my @ids = split;
	my $otu = '';
	foreach my $id (@ids)
	{
		$otu = $derep2otu{$id} if exists $derep2otu{$id};
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

open F, ">otu_member_list.txt";
foreach my $otu (sort {scalar(@{$otulist{$b}}) <=> scalar(@{$otulist{$a}})} keys %otulist)
{
	print F join("\t", $otu, @{$otulist{$otu}}), "\n";
}

