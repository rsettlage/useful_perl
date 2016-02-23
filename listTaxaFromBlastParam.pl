use strict;
use DBI;
use Getopt::Std;
use vars qw($opt_p $opt_r $opt_s $opt_t $opt_h $opt_d $opt_l $opt_b $opt_c $opt_g $opt_d $opt_i);
$opt_p = 0.8;
$opt_t = 'nucl';
$opt_s = 0;
$opt_l = 0; # debugging output, zero for off
$opt_g = 0; # print lineage in higherTaxonCounts
$opt_d = "/groups/DAC/taxonomy"; 

my $usage = "Usage: $0 [-options] blast_hit_table 
   hit table format requires only gi numbers in second column, bit scores in last column
   Options: 
        -p num: proportion of low-level taxa required to be subtended by higher taxon [default $opt_p]
        -r rank: zero-out higher-taxa below this rank, when encountered [default none]
        -t prot or nucl: database to look in for gi identifiers [default $opt_t]
        -s slack: consider hits with bit-scores as low as max-slack [default $opt_s]
        -b: add this flag to turn on writing out the higher taxon assigned to each read [default off]
        -l: with -b writes out low-level taxa for each read as well (does nothing without -b)
        -c: flags printing of higher taxon summary at checkpoints
        -g: flags printing lineage in higher taxon summary
        -i: flags printing taxon_id in summary table
        -d path_to_taxonomy_sqlitedb [default $opt_d]
";

if (!@ARGV or !getopts('p:r:t:s:d:hblcg') or $opt_h)
{
    print $usage;
    exit(0);
}
my $scoreSlack = $opt_s;
my $dnaOrProt = $opt_t;
my $supportThreshold = $opt_p;
my $minRank = $opt_r;
my $printBestTaxonPerRead = $opt_b;
my $printLowLevelTaxaPerRead = $opt_l;
my $printAtCheckpoints = $opt_c;

my $taxonomyDir = $opt_d; 
my $giTaxidTable = "gi_taxid_$dnaOrProt";
my $gidbname = "$taxonomyDir/gi_taxid.db";
my $taxonomyDb = "$taxonomyDir/taxonomy.db";

my $blastOutput = shift or die $usage;
-f $blastOutput or die "cannot find file $blastOutput";

my $outfilebase = $blastOutput . "_p$supportThreshold";
$outfilebase .= "s$scoreSlack" if $scoreSlack;
$outfilebase .= $minRank if $minRank;


my $dbargs = {AutoCommit => 0,
	      PrintError => 1,
              RaiseError => 1};

my $dbh = DBI->connect(          
    "dbi:SQLite:dbname=$taxonomyDb", 
    "",                          
    "",                          
    $dbargs,         
) or die $DBI::errstr;
my $gidbh = DBI->connect(          
    "dbi:SQLite:dbname=$gidbname", 
    "",                          
    "",                          
    $dbargs,         
) or die $DBI::errstr;

$dbh->do("pragma synchronous = off");
$dbh->do("pragma journal_mode = off");
#$dbh->do("pragma locking_mode = exclusive");
$dbh->do("pragma temp_store = memory");
$dbh->do("PRAGMA PAGE_SIZE = 4096");
$dbh->do("PRAGMA default_cache_size=700000"); 
$dbh->do("PRAGMA cache_size=700000"); 

my %badtaxa;
my $getBadTaxa = $dbh->prepare("select tax_id from name where name_txt like 'unclassified%' or name_txt like 'environmental%'");
$getBadTaxa->execute();
while (my ($badtax) = $getBadTaxa->fetchrow_array())
{
    $badtaxa{$badtax} = 1;
}
$getBadTaxa->finish();

my $getRankParent = $dbh->prepare("select rank, parent_tax_id from node where tax_id = ?");

my $getTaxonName = $dbh->prepare("select name_txt from name where tax_id = ?");

my $getTaxon = $gidbh->prepare("select taxid from $giTaxidTable where gi = ?");

open F, $blastOutput or die "cannot open $blastOutput";

my %readTaxon;
my %giTaxon;
my %taxonName;
my %nodeParent;
my %nodeRank;
#my %nodeName;
my $limit = $opt_d+0;

my %higherTaxonCount;
my %bestTaxon;
my %taxonCount;
my $time = time();
my $curRead;
my $minScore;
my $lineCount = 0;
my $checkpoint = 320000;
my %taxa;
my %giLackingTaxa;

if ($printBestTaxonPerRead)
{
    print STDERR "Writing read higher taxon to ${outfilebase}_bestTaxon.txt\n";
    open BEST, ">${outfilebase}_bestTaxonPerRead.txt";    
}

print STDERR "now read $blastOutput\n";
open F, $blastOutput or die "cannot open $blastOutput";
while (<F>)
{
    /^(\S+)\s+(\S+).*\s([\d\.]+)$/ or die "cannot parse $_";
    my ($read, $gi, $score) = ($1, $2, $3);
    if ($read ne $curRead)
    {
	if (%taxa)
	{
	    my $bestTaxon = grabBestTaxon(\%taxa);
	    #$bestTaxon{$curRead} = $bestTaxon;
	    $higherTaxonCount{$bestTaxon}++;
	    if ($printBestTaxonPerRead)
	    {
		print BEST "$curRead\t$bestTaxon";
		print BEST '  ', join(' ', sort keys % taxa) if ($printLowLevelTaxaPerRead);
		print BEST "\n";
	    }
	    %taxa = ();
	}
	$curRead = $read;
	$minScore = $score-$scoreSlack;
	$limit-- if $limit;
    }
    $lineCount++;
    if ($lineCount == $checkpoint)
    {
	my $elapsed = time()-$time;
	my $rate = $elapsed ? $lineCount/$elapsed : 0;
	printf "lines=%10d  time=%6d  rate=%10d  giLackingTaxa=%10d\n", $lineCount, $elapsed, $rate, scalar keys %giLackingTaxa;
	$checkpoint *= 2;
	if ($printAtCheckpoints)
	{
	    writeHigherTaxaToFile();
	}
    }
    
    next unless ($score >= $minScore);
    $gi = $1 if $gi =~ /gi\|(\d+)/;
    my $taxon = $giTaxon{$gi};
    unless ($taxon)
    {
	$getTaxon->execute($gi);
	$taxon = $getTaxon->fetchrow_array();
	$giTaxon{$gi} = $taxon;
    }
    if ($taxon)
    {
	$taxa{$taxon} = 1 unless $badtaxa{$taxon};
    }
    else
    {
	$giLackingTaxa{$gi} = 1;
    }
}
if (%taxa)
{
    my $bestTaxon = grabBestTaxon(\%taxa);
    #$bestTaxon{$curRead} = $bestTaxon;
    $higherTaxonCount{$bestTaxon}++;
    %taxa = ();
}
my $elapsed = time()-$time;
my $rate = $elapsed ? $lineCount/$elapsed : 0;
printf "total_lines=%10d  time=%6d  rate=%10d  giLackingTaxa=%10d\n", $lineCount, $elapsed, $rate, scalar keys %giLackingTaxa;
if (0)
{
    open NOTAX, ">gi_lacking_taxon.txt";
    print NOTAX join("\n", keys %giLackingTaxa);
}
%giLackingTaxa = ();

printf STDERR "num low-level taxa = %d\n", scalar keys %taxonCount;
printf STDERR "num high-level taxa = %d\n", scalar keys %higherTaxonCount;

if (1)
{
    writeHigherTaxaToFile();
}

$getRankParent->finish();
$getTaxonName->finish();
$getTaxon->finish();

$dbh->disconnect();
$gidbh->disconnect();

sub writeHigherTaxaToFile()
{
    print STDERR "Writing higher taxon counts to ${blastOutput}_higherTaxonCounts.txt\n";
    open HTC, ">${outfilebase}_higherTaxonCounts.txt";
    
    foreach my $taxon (sort {$higherTaxonCount{$b} <=> $higherTaxonCount{$a}} keys %higherTaxonCount)
    {
	my $name = $taxonName{$taxon};
	unless ($name)
	{
	    $getTaxonName->execute($taxon);
	    $name = $getTaxonName->fetchrow_array();
	    $taxonName{$taxon} = $name;
	}
	print HTC $higherTaxonCount{$taxon};
        print HTC "\t$taxon" if $opt_i;
        print HTC "\t$name\t$nodeRank{$taxon}";
	if ($opt_g)
	{#print lineage
	    my @lineage;
	    my $anc = $taxon;
	    while ($anc > 1)
	    {
		$name = $taxonName{$anc};
		unless ($name)
		{
		    $getTaxonName->execute($anc);
		    $name = $getTaxonName->fetchrow_array();
		    $taxonName{$anc} = $name;
		}
	        last if $name =~ /cellular.organisms/;
		push @lineage, $name;
		$anc = $nodeParent{$anc}
	    }
	    print HTC "\t", join(";", reverse @lineage);
	}
	print HTC "\n";
    }
    close HTC;
}

sub grabBestTaxon($)
{
    my $taxa = shift; # reference to hash
    return 0 unless scalar keys %$taxa;
    print STDERR "grabBestTaxon($limit): ", join("  ", %$taxa), "\n" if $limit > 0;
    my $supportThresholdSum = $supportThreshold * scalar keys %$taxa;
    my %readSupport = ();
    foreach my $taxon (keys %$taxa)
    {
	my @trail = ();
	my $startTaxon = $taxon;
	my $startSupport = 1;

	while ($taxon > 1)
	{
	    my $anc = $nodeParent{$taxon};
	    unless ($anc)
	    {
		$getRankParent->execute($taxon);
		my $rank;
		($rank, $anc) = $getRankParent->fetchrow_array();
		$nodeParent{$taxon} = $anc;
		$nodeRank{$taxon} = $rank;
		if (0 && $limit > 0)
		{
		    printf STDERR "got rank, anc for %8d: %s  %d\n", $taxon, $rank, $anc; 
		    #exit(0);
		}
	    }
	    $readSupport{$taxon} += $startSupport;
	    if ($minRank and $nodeRank{$taxon} =~ /^$minRank/)
	    {# zero-out things below minimum rank
		while (my $t = pop @trail)
		{
		    $readSupport{$t} = 0;
		}
	    }
	    push @trail, $taxon;
	    $taxon = $anc;
	}
    }
    if ($limit > 0)
    {
	print STDERR "numtaxa=", scalar(keys %$taxa), ", readsupportSize=", scalar keys %readSupport, " need $supportThresholdSum\n";
	if (keys %$taxa < 10)
	{
	    print STDERR join("\t", sort keys %$taxa), "\n";
	}
    }
    foreach my $taxon (keys %readSupport)
    {
	if ($readSupport{$taxon} >= $supportThresholdSum)
	{# zero out parent line going up, leave only one node with support
	    $taxon = $nodeParent{$taxon};
	    while ($taxon > 1)
	    {
		$readSupport{$taxon} = 0;
		$taxon = $nodeParent{$taxon};
	    }
	}
    }
    my $retval = -9;
    foreach my $taxon (sort keys %readSupport)
    {
	if ($readSupport{$taxon} >= $supportThresholdSum)
	{
	    
	    if ($limit > 0)
	    {
		printf STDERR "acceptable($limit): $taxon: thresh=%.2f, got=%.2f, %s\n", $supportThresholdSum, $readSupport{$taxon}, $nodeRank{$taxon};
	    }
	    $retval = $taxon if $readSupport{$taxon} > $readSupport{$retval};
	}

    }
    print STDERR "returning $retval\n" if $limit > 0;
    return $retval;
}

