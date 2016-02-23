#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Description

## Usage:
## CHECK usage with option  -h

use File::Basename;
use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $execSff, $seqFn, $inSff, $outSff);

my $keyLen = 4;
$execSff = `which sfffile`;
$execSff =~ s/[\r\n]//;

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"exec=s"		=> \$execSff,	
	"fastaFile=s"		=> \$seqFn,
	"sff=s"		=> \$inSff,
	"output=s"		=> \$outSff,		## output dir
) || die "\n";


if ($helpFlag || !$seqFn || !$inSff)
{
	die("Usage : $0 -f <fasta file> -s <original sff file>  [-e 'fullpath of 454/bin/sfffile'] [-o out_sff_file]\n\n");
}

if(!defined $execSff  || $execSff eq ""){
	die("Need to a full path of 'sfffile' program\n\n");
}
(-e $execSff) || die "Could not find '$execSff' program.\n\n";

if(!$outSff) {
	$outSff = $seqFn;	
	$outSff =~ s/(seq|fa|fna|fasta)$/sff/;
	if($inSff eq $outSff){
		$outSff = $inSff;	
		$outSff =~ s/sff$/new.sff/g;
	}
}


my $binPath = dirname($0);


my %reads;
my $tmpFile = "$seqFn.tmp";


load($seqFn);
process($tmpFile);			

system("$execSff -i $tmpFile -t $tmpFile -o \"$outSff\" \"$inSff\"");
#unlink($tmpFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n";
	my $in = openInput($fileName);

	while (<$in>)
	{
		if (/^>\s*(\S+).*trimmed=(\d+)-(\d+)/)
		{
			my ($id, $beg, $end) = ($1, $2, $3);
			$id =~ s/:\d+$//;			 
			$reads{$id} = [$beg+$keyLen, $end+$keyLen] if !exists $reads{$id} || $end-$beg > $reads{$id}[1]-$reads{$id}[0];
		}
	}

	close($in) if defined $fileName;
}

sub process
{
	my ($fileName) = @_;

	my $out = openOutput($fileName);

	foreach my $id (keys %reads)
	{
		print $out join("\t", $id, @{$reads{$id}}), "\n";
	}

	close($out) if defined $fileName;
}

#-------------------------------------------------------------------------------

sub openInput
{
	my ($fileName) = @_;

	return *STDIN unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz(ip)?$/ ? "zcat $fileName |" : $fileName =~ /.bz(ip)?2$/ ? "bzcat $fileName |" : $fileName) || die("Open error: $fileName");
	return $fd;
}

sub openOutput
{
	my ($fileName) = @_;

	return *STDOUT unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

