use strict;

# read uc file, then read derep seq file, then read underep seq file
# for each un-reprep seq, print seq_id (sample embeded) and OTU

my ($ucfile, $derepSeqFile, $rawSeqFile) = @ARGV;
die "usage: $0 ucfile derepSeqFile rawSeqFile, output to otu_member_list_relabel.txt" unless ($ucfile && $derepSeqFile && $rawSeqFile);
-f $ucfile or die "cannot find uc file $ucfile";

open F, $ucfile or die "cannot open $ucfile";
my %derep2otu;
while (<F>)
{
	next if /^N/; # no hit
	chomp;
	my @r=split;
	$derep2otu{$r[8]}=$r[9];
}

open F, $derepSeqFile or die "cannot open $derepSeqFile";
my %derepSeq;
my $id;
while (<F>)
{
	if (/^>(\S+)/)
	{	
		$id = $1;
		warn("Hey, no otu for id $id in $derepSeqFile\nn") unless $derep2otu{$id};
	}
	else
	{
		chomp;
		$derepSeq{$id} .= $_;
	} 
}
open L, ">seqID_to_otuID.list";
my %seq2otu;
foreach  my $id (keys %derepSeq)
{
	$seq2otu{$derepSeq{$id}} = $derep2otu{$id};
	print L "$id\t$derep2otu{$id}\n";
}

my %nonderep2otu;
my %otumember;
open F, $rawSeqFile or die "cannot open $rawSeqFile";
my $seq;
while (<F>)
{
	        if (/^>(\S+)/)
        {
                $id = $1;
		$seq = '';
		$nonderep2otu{$id} = 'NA';
        }
        else
        {
                chomp;
		$seq .= $_;
                $nonderep2otu{$id} = $seq2otu{$seq} if $seq2otu{$seq};
        }
}


my %otulist;
foreach my $id (keys %nonderep2otu)
{
	push @{$otulist{$nonderep2otu{$id}}}, $id;
}

open F, ">otu_member_list.txt";
#my $otu = 1;
my $fieldWidth = int(log(keys %otulist)/log(10))+1;
foreach my $otu (sort {scalar(@{$otulist{$b}}) <=> scalar(@{$otulist{$a}})} keys %otulist)
{
	#my $otuLabel = $id;
	#if ($id ne "NA")
	#{
#		$otuLabel = sprintf "OTU_%0${fieldWidth}d", $otu;
#		$otu++;
#	}
	print F join("\t", $otu, @{$otulist{$otu}}), "\n";
#	print L "$id\t$otuLabel\n";
}

