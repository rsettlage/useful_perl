#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

use strict;
use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $dir, $seqName, $capital);
GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"capital"	=> \$capital,		
	"direction=s"	=> \$dir,		
	"name=s"	=> \$seqName
) || help(1);

help(0) if defined $helpFlag;


my $seqFile=shift;
my $pos = shift;
my $len = shift;

if(defined $pos && $pos =~ /-/){
	$pos = ($` < $' ? $` : $'); 
	$len = abs($' - $`) + 1 if(!defined $len);
}
elsif(defined $pos && $pos =~ /\.\./){
	$pos = ($` < $' ? $` : $'); 
	$len = abs($' - $`) + 1 if(!defined $len);
}
elsif(defined $pos && $pos =~ /:/){
	$pos = ($` < $' ? $` : $'); 
	$len = abs($' - $`) + 1 if(!defined $len);
}

$dir = "a" if(!defined $dir);

$pos = 1 if(!defined $pos);
$len = 0 if(!defined $len);

(defined $seqFile) || help(1);


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print " Usage: $0  sequence_file  position  length  [-d direction] [-n name of the sequence] [-c (print capital)]\n";
	print "  or : $0  sequence_file  start:end  [-d direction] [-n name of the sequence] [-c (print capital)]\n";
	print "  or : $0  sequence_file  start-end  [-d direction] [-n name of the sequence] [-c (print capital)]\n\n";
	exit($return);
}

(-e $seqFile) || die "Could not find the '$seqFile' file.\n";

#print "$pos, $len --  \n";
open(IN, $seqFile =~ /.gz$/ ? "zcat $seqFile |" : $seqFile);

my ($seq, $offset, $readLen, $subLen, $buf, $name, $i, $found);

$offset = 0;
$subLen = 0;
$buf = <IN>;

$buf =~ s/[\r\n]//g;
if($buf !~ /^>/){
	$readLen = length($buf);
	if($offset + $readLen > $pos){
		$seq = uc($len != 0 ? substr($buf, ($pos-$offset-1), $len) : substr($buf, ($pos-$offset-1)));;
		#$seq = substr($buf, ($pos-$offset-1), $len);
	}
	$offset += $readLen;
}

while(defined $buf){
	if($buf =~ /^>/){
		$offset = 0;
		$subLen = 0;
		
		last if(!defined $seqName && defined $name);
		my @arr = split /[\t\s\r\n]+/, $';
		$name = $arr[0];
#		$name =~ s/[\r\n]//g;
		
		if($name =~ /[ \t]/){ $name = $`;	}
	}
	
	$buf = <IN>;
	if(!defined $seqName || $name eq $seqName){
		while(defined $buf){
			$buf =~ s/[\r\n]//g;
			$readLen = length($buf);
			#print "$offset + $readLen > $pos && $offset < $pos \n";
			if($offset + $readLen >= $pos && $offset < $pos){				
				#$seq = substr($buf, ($pos-$offset-1), $len);
				$seq = uc($len != 0 ? substr($buf, ($pos-$offset-1), $len) : substr($buf, ($pos-$offset-1)));;
				$subLen = $len - length($seq);
			}
			elsif($offset >= $pos && $subLen != 0){				
				#$seq .= substr($buf, 0, $subLen);
				$seq .= uc($len != 0 ? substr($buf, 0, $subLen) : substr($buf, 0));
				$subLen = $len - length($seq);
			}
			elsif($offset > $pos && $subLen == 0) { $found= 1;  last; }
			$offset += $readLen; 
		
			$buf = <IN>;
			if(!defined $buf ||  $buf =~ /^>/){ 
				print STDOUT "The length of the sequence is '$offset'\n";
				if(defined $seq && $seq ne "") { $found = 1; }
			 	last; 
			}
		}
	}
	else { 
		while(defined($buf = <IN>)){ last if($buf =~ /^>/);}
	}
	last if(defined $found);
}

(defined $found) || die "Could not find '$seqName' \n";

#print ">$name\n" if(defined $name);
my $end = $pos -1 + length($seq);
if($dir eq "+" || $dir eq "a") 
{
	print ">$name:$pos:$end:+\n" . (defined $capital ? uc $seq : $seq) ."\n";
}
if($dir eq "-" || $dir eq "a") 
{
	$seq = rc($seq);
	print ">$name:$pos:$end:-\n" . (defined $capital ? uc $seq : $seq) . "\n";
}
 
sub rc
{
	my ($seq, $type) = @_;

	my $rc = reverse $seq;
	if (!$type) # DNA
	{
		$rc =~ y/acgtuACGTU/tgcaaTGCAA/;
	}
	else # RNA
	{
		$rc =~ y/acgtuACGTU/ugcaaUGCAA/;
	}

	return $rc;
}
