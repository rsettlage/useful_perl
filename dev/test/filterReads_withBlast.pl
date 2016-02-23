#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2011
#
# filterReads_withBlast.pl
# This program filters sequencing reads with blast results.
######################################

use strict;
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $evalFlag, @cutoff_E_values, @cutoff_L_values, @blastFn, @seqFn, $outFn);
$cmd = "$0 @ARGV";


GetOptions(
	"h|?|help"		=> \$helpFlag,
	"e=s" => \@cutoff_E_values,
	"l=s" => \@cutoff_L_values,
	"blast=s"	=> \@blastFn,
	"sequence=s"	=> \@seqFn
) || help(1);


help(0) if defined $helpFlag;

if($#blastFn < 0){ print STDERR "\nNeed one blast table file for single end  OR  two files for paired end!!\n\n"; help(1); }
if($#seqFn < 0){ print STDERR "\nNeed one fastq file for single end  OR two files for paired end!!\n\n"; help(1); }

if($#cutoff_E_values == -1 && $#cutoff_L_values == -1){ print STDERR "\nNeed cutoff (e-values) or (length:percentage) value!!\n\n"; help(1); }
if($#cutoff_E_values != -1 && $#cutoff_L_values != -1){ print STDERR "\nNeed only one of (cutoff e-values) and (length:percentage) values!!\n\n"; help(1); }


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  -b <blast table file 1>  -b <blast table file 2>  -s <sequence file 1>  -s <sequence file 2>\n";
	print STDERR "        ( -e <cutoff e-value[:reference]>  [-e <cutoff e-value[:reference]>  ....] )  or \n";
	print STDERR "        ( -l <(cutoff match length):(cutoff percentage)> [-l<(cutoff match length):(cutoff percentage)> ...] )\n\n";
	print STDERR "  ex 1) $0 -b s_1_1.bln.table -b s_1_2.bln.table -s s_1_1.fastq -s s_1_2.fastq  -e 1e-04:pTARBAC2.1 -e 1e-07:NC_010473\n";
	print STDERR "  ex 2) $0 -b s_1_1.bln.table -b s_1_2.bln.table -s s_1_1.fastq -s s_1_2.fastq  -e 1e-04\n\n";
	print STDERR "  ex 2) $0 -b s_1_1.bln.table -b s_1_2.bln.table -s s_1_1.fastq -s s_1_2.fastq  -l 60:10 -l 50:90 \n\n";
	exit($return);
}



my ($in, $out, @arr, $i, %hitReads, $val);

my $log10 = log(10);

#### different criteria for multiple references..
if($#cutoff_E_values > 0 || ($#cutoff_E_values == 0 && $cutoff_E_values[0] =~ /:/)){
	my %cutoff;
	for($i = 0; $i <= $#cutoff_E_values; $i++) {
		my ($eval, $ref) = $cutoff_E_values[$i] =~ /(.*):(.*)/;
		if(!defined $ref){
			print STDERR "[Error] The cutoff '$cutoff_E_values[$i]' format is wrong. If there two or more cutoff values, you have to use '(evalue):(reference)'\n"; help(1);
			exit;
		}		
		$cutoff{$ref} = $eval;
	}


	my ($query, $ref);
	my (%used);
	#### blast table output format.
	###[0]query_id [1]db_id [2]percentage of identity [3]aligned_length(identity) [4]mimatches [5]gap_openings [6]q.start [7]q.end [8]s.start [9]s.end [10]e-value [11]hit
	for(my $fi = 0; $fi <= $#blastFn; $fi++){
		$in = openInput($blastFn[$fi]);
		while(<$in>){
			s/[\r\n]+//g;
			next if(/^#/ || /^\s*$/);
			
			@arr = split /\t/;

			$arr[0] =~ /\//;
			$arr[0] = $`;  ### it takes only 'HWUSI-EAS1737:1:1:1053:16879#0' from 'HWUSI-EAS1737:1:1:1053:16879#0/1'..
			($query, $ref) = ($arr[0], $arr[1]);
			##### to treat paired end reads as one read.
				
			next if(defined $used{$ref}{$query});
			$used{$ref}{$query} = 1;

			#### use multiple criteria ...
			if(defined $cutoff{$ref} && $cutoff{$ref} >= $arr[10]){
				$hitReads{$query} = 1;
			}
		}
		close($in);
	}
}
####### for only one cutoff value..
elsif($#cutoff_E_values == 0)
{
	my $cutoff = $cutoff_E_values[0];	

	my ($query);
	my (%used);
	#### blast table output format.
	###[0]query_id [1]db_id [2]percentage of identity [3]aligned_length(identity) [4]mimatches [5]gap_openings [6]q.start [7]q.end [8]s.start [9]s.end [10]e-value [11]hit
	for(my $fi = 0; $fi <= $#blastFn; $fi++){
		$in = openInput($blastFn[$fi]);
		while(<$in>){
			s/[\r\n]+//g;
			next if(/^#/ || /^\s*$/);
			
			@arr = split /\t/;

			$arr[0] =~ /\//;
			$arr[0] = $`;  ### it takes only 'HWUSI-EAS1737:1:1:1053:16879#0' from 'HWUSI-EAS1737:1:1:1053:16879#0/1'..
			($query) = ($arr[0]);
			##### to treat paired end reads as one read.
				
			next if(defined $used{$query});
			$used{$query} = 1;

			if($cutoff >= $arr[10]){
				$hitReads{$query} = 1;
			}
		}
		close($in);
	}
}
###### for length and identity cutoff
else
{
	my (@cutoffMatchLen, @cutoffPercentages, $matchLen);

	for($i = 0; $i <= $#cutoff_L_values; $i++) {
		@arr = split /:/, $cutoff_L_values[$i];
		if(!defined $arr[1]){
			 print STDERR "\nWrong format of cutoff value!! : '-l $cutoff_L_values[$i]'\n\n"; help(1);
		}
		$cutoffMatchLen[$i] = $arr[0];
		$cutoffPercentages[$i] = $arr[1];
	}

	my ($query);
	my (%used);
	#### blast table output format.
	###[0]query_id [1]db_id [2]percentage of identity [3]aligned_length(identity) [4]mimatches [5]gap_openings [6]q.start [7]q.end [8]s.start [9]s.end [10]e-value [11]hit
	for(my $fi = 0; $fi <= $#blastFn; $fi++){
		$in = openInput($blastFn[$fi]);
		while(<$in>){
			s/[\r\n]+//g;
			next if(/^#/ || /^\s*$/);
			
			@arr = split /\t/;

			$arr[0] =~ /\//;
			$arr[0] = $`;  ### it takes only 'HWUSI-EAS1737:1:1:1053:16879#0' from 'HWUSI-EAS1737:1:1:1053:16879#0/1'..
			($query) = ($arr[0]);
			##### to treat paired end reads as one read.
				
			next if(defined $used{$query});
			$used{$query} = 1;
			$matchLen = $arr[3]-$arr[4]-$arr[5];

			#### use multiple criteria ...
			for($i = 0; $i <= $#cutoffMatchLen; $i++) {
				if(	$matchLen >= $cutoffMatchLen[$i] && $arr[2] >= $cutoffPercentages[$i] ){
					$hitReads{$query} = 1;
					last;
				}
			}
		}
		close($in);
	}
}


my($name, $buf, $type);
for(my $fi = 0; $fi <= $#seqFn; $fi++)
{
	$in = openInput($seqFn[$fi]);
	my $outNoHit = openOutput("$seqFn[$fi].nohit");
	my $outHit = openOutput("$seqFn[$fi].hit");

	print "Writing $seqFn[$fi].hit and $seqFn[$fi].nohit\n";

	$buf = <$in>;
	while($buf){
		if($buf =~ /^>([\w_+\-:\.\|#]+)/ || $buf =~ /^@([\w_+\-:\.\|#]+)/){ #### it does not read '/1' or '/2' from a read name. ex) HWUSI-EAS1737:1:1:1053:16879#0/1
			$name = $1;		

			if($buf =~ /^>/) { $type ="fa"; }
			else { $type ="fq"; }

			##### ----------------------
			#if(!defined $hitReads{$name} || $hitReads{$name} == 1){ $buf = <$in>; next; } ### if you need to remove reads satisfying the criteria, change this line..
			if(defined $hitReads{$name} && $hitReads{$name} == 1){ 
				$out = $outHit;
			} 
			else{
				$out = $outNoHit;
			}

			print $out $buf;
			$buf = <$in>;
			if($type eq "fq"){ 
				print $out $buf;
				$buf = <$in>;
				print $out $buf;
				$buf = <$in>;
				print $out $buf;
				$buf = <$in>;
			}
			else{
				while(defined $buf && ($buf !~ /^>/)){
					print $out $buf;
					$buf = <$in>;
				}
			}
		}
		else{
			$buf = <$in>;
		}
	}
	close($in);
	close($outNoHit);
	close($outHit);
}



sub openInput
{
	my ($fn) = @_;

	return *STDIN unless defined $fn;

	my ($fd);
	open($fd, $fn =~ /\.gz/ ? "zcat $fn|" : ($fn =~ /\.bz2/ ? "bunzip2 -c $fn|" : $fn)) || die "Could not open '$fn' : $!\n";
	return $fd;
}

sub openOutput
{
	my ($fn) = @_;

	return *STDOUT unless defined $fn;

	my ($fd);
	open($fd, $fn =~ /.gz$/ ? "| gzip -c > $fn" : ($fn =~ /\.bz2/ ? "| bzip2 -c > $fn" : ">$fn")) || die "Could not write '$fn' : $!\n";
	return $fd;
}


sub rc
{
	my ($seq, $type) = @_;

	my $rc = reverse $seq;
	if(defined $type && ($type eq "rna" || $type eq "RNA")) # RNA
	{   $rc =~ y/acgtuACGTU/ugcaaUGCAA/;  }
	else ## DNA
	{   $rc =~ y/acgtuACGTU/tgcaaTGCAA/;  }

	return $rc;
}


