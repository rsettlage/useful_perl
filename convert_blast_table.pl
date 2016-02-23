#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: 2006

## Convert a raw megablast output to tab delimeted output as -D 3

## Usage:
## CHECK usage with option  -h

## $Id: blast_table.pl,v 1.4 2007-09-09 03:25:30 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);
#use My::Bioinfo;

my ($helpFlag, $inFile, $outFile, $bl2seq, $verbose);
my (@types, $queFw, $subFw, $queBest, $subBest, $headerFlag, $quiet, $alignmentFlag, $allFlag, $lenFlag);
my $numReg = '[\d\+\-\.e]+';

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"inFile=s"		=> \$inFile,		## input file
	"outFile=s"		=> \$outFile,		## output file
	"bl2seq"			=> \$bl2seq,
	"v|verbose+"	=> \$verbose,		## verbose output
	"q|qforward"	=> \$queFw,
	"s|sforward"	=> \$subFw,
	"header!"		=> \$headerFlag,
	"type=s"			=> \@types,
	"qb|qbest"		=> \$queBest,
	"sb|sbest"		=> \$subBest,
	"a|alignment!"	=> \$alignmentFlag,
	"all!"			=> \$allFlag,
	"qu|quite"		=> \$quiet,
	"length"			=> \$lenFlag,
) || die "\n";

checkOptions();

my (%qBest, %sBest);

my ($query, @subjects, $x1, $x2, $y1, $y2, $score, $evalue, $length, $match, $gap, $pi, $qseq, $sseq, $frame, $qlength, $slength);
$x1 = 0; $x2 = 0;
$y1 = 0; $y2 = 0;

my $fd = openInput($inFile);
my $out = openOutput($outFile);

if (!defined $bl2seq)
{
	my $program = <$fd>;
	print $out "# $program" if $headerFlag;
}

print $out "# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n" if $headerFlag;

while (<$fd>)
{
	# Query= Dp21 P1-B112000FW52758 Description:  P1-B11(2000)FW(52>758);
	#          (707 letters)
	if (/^\s*Query=(.*)/)
	{
		onMatch() if @subjects;
		$query = _id($1, \@types);
		@subjects = ();
		while (<$fd>)
		{
			if (/([,\d]+)\s*letter/) { $qlength = $1; last; }
		}
		die "Unknown length in $_" if $lenFlag && !$qlength;
		$qlength =~ s/,//g;
	}
	#>Supercontig281 Name=Daphnia+pulex,9XSTL+Supercontig281;
	#                       Length = 151782
	elsif (/^>/)
	{
		onMatch() if @subjects;
		@subjects = ();
		push(@subjects, _id($', \@types));
		$score = -1;
		while (<$fd>)
		{
			if (/^\s+Length\s*=\s*(\d+)/i) { $slength = $1; last; }
			next if !/^( +)/ || length($1) > 3;
			push(@subjects, _id($', \@types));
		}
		die "Unknown length in $_" if $lenFlag && !$slength;
	}
	# Score = 5069 bits (2557), Expect = 0.0
	#elsif (/^\s*Score\s*=\s*($numReg) bits.*Expect(\(\d+\))*\s*=\s*($numReg)/)
	elsif (/^\s*Score\s*=\s*($numReg) bits \((\d+)\).*Expect(\(\d+\))*\s*=\s*($numReg)/)
	{
		onMatch() if $score != -1;
		$score = $2;
		$evalue = $4;
		$evalue =~ s/^e/1e/;
	}
	# Identities = 89/94 (94%), Gaps = 2/94 (2%)
	elsif (/^\s*Identities\s*=\s*(\d+)\/(\d+)\s*\(($numReg)\%\)/)
	{
		$length = $2;
		$match = $1;
		$pi = $3;
		if ($' =~ /.*Gaps\s*=\s*(\d+)\/(\d+)/)
		{
			die "Inconsistent matching length: $length $2\n" if $2 != $length;
			$gap = $1;
		}
		else
		{
			$gap = 0;
		}
	}
	elsif (/^\s*Frame\s*=\s*([\d\+\-]+)/)
	{
		$frame = $1;
	}
#	elsif (/^\s*Strand\s*=\s*(\w+)\s*\/\s*(\w+)/)
#	{
#		$x1 = 0; $x2 = 0;
#		$y1 = 0; $y2 = 0;
#	}
	#Query: 541707 agtcggttcggtcctccagttggttttacccaaccttcaacctgcccgtggctagatcac 541766
	elsif (/^\s*Query:\s*(\d+)\s+(\S+)\s+(\d+)/)
	{
		if ($x1 == 0) { $x1 = $1; }
		$x2 = $3;
		$qseq .= $2 if $alignmentFlag;
	}
	#Sbjct: 120830 agtcggttcggtcctccagttgatgttactcaaccttcaacctgctcatgactagctcaa 120771
	elsif (/^\s*Sbjct:\s*(\d+)\s+(\S+)\s+(\d+)/)
	{
		if ($y1 == 0) { $y1 = $1; }
		$y2 = $3;
		$sseq .= $2 if $alignmentFlag;
	}
	elsif (/^Database:/)
	{
		print $out '#', $_ if $headerFlag;
	}
}

onMatch() if defined $query;

close($fd);
close($out);


sub onMatch
{
	return if !@subjects;
	
	goto clean if ($queBest && exists $qBest{$query} && $qBest{$query} < $evalue)
	           || ($subBest && exists $sBest{$query}{$subjects[0]} && $sBest{$query}{$subjects[0]} < $evalue);

	$qBest{$query}               = $evalue if $queBest;
	$sBest{$query}{$subjects[0]} = $evalue if $subBest;

	foreach my $subject (@subjects)
	{
		my $qstrand = $x1 <= $x2 ? 1 : 0;
		my $sstrand = $y1 <= $y2 ? 1 : 0;
		if (($queFw && !$qstrand) || ($subFw && !$sstrand))
		{
			($x1, $x2, $y1, $y2) = ($x2, $x1, $y2, $y1);
			($qseq, $sseq) = (_rc($qseq), _rc($sseq)) if $alignmentFlag;
		}
		print $out join("\t", $query, $subject, $pi, $length, $length-$match-$gap, $gap, $x1, $x2, $y1, $y2, $evalue, $score);
		print $out "\t$frame" if defined $frame;
		print $out "\t$qlength\t$slength" if $lenFlag;
		print $out "\n";
		if ($alignmentFlag)
		{
			print $out join("\t", '#=Q', $query, $qseq), "\n";
			print $out join("\t", '#=S', $subject, $sseq), "\n";
		}

		last if !$allFlag;
	}

	die "Wrong format\n" if $score < 0;

clean:
	$x1 = 0; $x2 = 0;
	$y1 = 0; $y2 = 0;
	$qseq = $sseq = '' if $alignmentFlag;
	$frame = undef if defined $frame;
}

#-------------------------------------------------------------------------------

sub openInput
{
	my ($fileName) = @_;

	return STDIN unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /\.gz/ ? "zcat $fileName|" : $fileName) || die("Open error:$fileName");
	return $fd;
}

sub openOutput
{
	my ($fileName) = @_;

	return STDOUT unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag)
	{
		die("Arguments: <option> [[-i] in_file] [[-o] out_file] [-v]\n"
		  . "\t-b(l2seq) : input is the output of bl2seq\n"
		  . "\t-t(ype) string : how to extract ids from sequence header (gi, ref, 454, ...)\n"
		  . "\t    accept one or more types\n"
		  . "\t-qf(orward) | -sf(orward) : forward strand in either query or subject sequences\n"
		  . "\t-qb(est) : output only best matche for each query\n"
		  . "\t-sb(est) : output only best matche for each query and each subject\n"
		  . "\t-a(lignment) : output or not header\n"
		  . "\t-(no)he(ader) : output or not alignment\n"
		  );
	}

	$headerFlag = 0 if !defined $headerFlag;
}

sub _id
{
  my ($header, $types) = @_;

  my $id = '';

  foreach my $type (@$types)
  {
    if ($type eq '454')
    {
      last if ($id) = $header =~ /uaccno=(\S+)/;
    }
    elsif ($type =~ /^col_(\d+)/)
    {
      my $col = $1;
      $id = (split(/\s+/, $header))[$col-1];
    }
    elsif ($header =~ /$type\|([^|\s]+)/)
    {
      $id = $1;
      last;
    }
  }

  if (!$id)
  {
    $header =~ /^>*\s*(\S+)/;
    $id = $1;
  }

  return $id;
}

