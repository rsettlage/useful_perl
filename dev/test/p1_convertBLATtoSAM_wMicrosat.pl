#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2010
#
# convertBLATtoSAM_withMicrosat.pl
# convert pslx files to SAM files considering microsatellite positions.
######################################

use strict;
use Getopt::Long qw(:config no_ignore_case);
use IO::Handle;
STDOUT->autoflush(1);
STDERR->autoflush(1);

my ($verbose, $helpFlag, $pairedFlag, $pslFn, $microFlag, $qMicFn, $querySeqFn, $refFn, $micFn, $cutoffError, $outFn, $cmd, $allowdBestHitNum);
$cmd = "$0 @ARGV";

#$verbose = 1;
$cutoffError = 0.07;
$allowdBestHitNum = 1;


GetOptions(
	"verbose"		=> \$verbose,
	"h|?|help"		=> \$helpFlag,
	"paired"	=> \$pairedFlag,
	"Micro"		=> \$microFlag,
	"microsat=s"	=> \$micFn,
	"qMic=s"	=> \$qMicFn,
	"cutoff=i"	=> \$cutoffError,
	"p|psl=s"	=> \$pslFn,
	"seq=s"	=> \$querySeqFn,
	"reference=s"	=> \$refFn,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;


if(!defined $pslFn){ print STDERR "\nNeed a pslx file or directory containing splitted pslx files from BLAT!!\n\n"; help(1); }
if(!defined $querySeqFn){ print STDERR "\nNeed a (FASTQ|FASTA) query sequence file!!\n\n"; help(1); }
if(!defined $refFn){ print STDERR "\nNeed a reference file!!\n\n"; help(1); }
if(!defined $micFn){ print STDERR "\nNeed a microsatellite list file!!\n\n"; help(1); }
if(!defined $qMicFn){ print STDERR "\nNeed a query microsatellite list file!!\n\n"; help(1); }
if(!defined $outFn){ print STDERR "\nNeed a output file!!\n\n"; help(1); }


if(defined $outFn && ((defined $pslFn && $pslFn eq $outFn) || $querySeqFn eq $outFn || (defined $qMicFn && $qMicFn eq $outFn) || $refFn eq $outFn || $micFn eq $outFn))
{ print STDERR " Error) Input and output files are same \n"; exit(1); }

sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  -m <microsat list file>  -r <reference file>  -p <pslx file>  -s <FASTQ query sequence file>  [-qMic <query microsats list file>]  <[-c <rate of min. mismatch/gap>] [-M (show only hits related to microsats)]  -o <out file>\n";
	print STDERR "  ex) $0 -m ~/sequencing/targetEnrichment/design1_400.fa_mic.lst -r bwaidx/design1_400.fa -p tmp.MCF7.sureSelect.design1_400.unmapped.pslx -s MCF7.sureSelect.design1_400.unmapped.fa  -qMic  MCF7.sureSelect.unmappedRead.microsat.lst  -M  -o out.sam\n\n";
	print STDERR "  ex) $0 -m ../Misa_microsat/hg19_mic.lst -r ~/human/hg19.fa -p split_unmapped -s NC09.New.hg19.unmap.fa -qMic NC09.New.hg19.unmap.microsat.lst -o test.sam\n";
	print STDERR "  ex) $0 -m chr1_mic.lst -r ~/human/hg19/chr1.fa -p split_unmapped/chr1.pslx -s NC09.New.hg19.unmap.fa -qMic NC09.New.hg19.unmap.microsat.lst -o chr1.test.sam\n";
	print STDERR "  ex) $0 -m chr1_mic.lst -r chr1.test.fa -p chr1.test.pslx -s chr1.test.read.fa -qMic chr1.test.qMic_list -o test.sam \n";
	print STDERR "  ex) $0 -m ../../design3_around_90.fa.mic.lst -r ../../design3_around_90.fa -qMic GM13869_cot1.hg19_re.unmapped.mic.lst -p GM13869_cot1.hg19_re.unmapped.pslx -s GM13869_cot1.hg19_re.unmapped.fa -o test.sam\n";
	exit($return);
}


my ($in, $out, @arr, $i, $j, $ref, $name, %allMicInRef, %idxMicPos, @g_log_n);

for($i = 1; $i < 256; $i++) { $g_log_n[$i] = int(4.343 * log($i) + 0.5);	}

$in = openInput($micFn);
while(<$in>){
	s/[\r\n]+//g;
	next if(/^#/ || /^\s*$/);

	@arr = split /\t/;
	$ref = $arr[0];
	$i = $#{$allMicInRef{$ref}} + 1;
	$allMicInRef{$ref}[$i]{start} = $arr[1] - 1; ##### -1 : for blat.. blat is 0 based starting.
	$allMicInRef{$ref}[$i]{end} = $arr[1] + $arr[2] -1 -1;
	$allMicInRef{$ref}[$i]{len} = $arr[2];
	$allMicInRef{$ref}[$i]{str} = $_;

	$j = int($arr[1]/10000);
	if(!defined $idxMicPos{$ref}[$j]) { $idxMicPos{$ref}[$j] = $i; }
}

close($in);

my (%allMicInQuery, %skipQuery);
if($qMicFn){
	$in = openInput($qMicFn);
	while(<$in>){
		s/[\r\n]+//g;
		next if(/^#/ || /^\s*$/);
		
		if(/^>([\w_\-:\.\|#\/]+)\t([\-\d]+)/){
		
			if($2 != -1) { $name = $1; }
			else {
				$skipQuery{$1} = 1;
				$name = undef;
			}
		}
		elsif(defined $name){		
			if(/\t/){
				@arr = split /\t/;
				if(defined $allMicInQuery{$name}){ $allMicInQuery{$name} = $allMicInQuery{$name}.",$arr[0]-$arr[1]"; }
				else { $allMicInQuery{$name} = "$arr[0]-$arr[1]"; }
			}
			else  {	$allMicInQuery{$name} = $_; }
		}
	}

	close($in);
}

my(%refSeq);
$in = openInput($refFn);

while(<$in>){
	s/[\r\n]+//g;
	next if(/^#/ || /^\s*$/);
	
	if(/^>([\w_+\-:\.\|#\/]+)/){
		$ref = $1;
	}
	else{
		$refSeq{$ref} .= uc($_);
	}
}
close($in);



my ($inS, $buf,  $micI, $strand, $outMs,
	$arrMics, $gapNum, $gapLen, @mapList, $num, $matchLen, @pslList, $str, $rPos, $qPos, 
	$qName, $qStart, $qEnd, $qLen, $rName, $rStart, $rEnd, $isRelatedMic, $motifLen, $isMic,
	$len, $seq, $revSeq, $oriSeq, $qual, $qSeq, $templateSeq, $aln, $errNum, $minErrNum, $secMinErrNum, $isSkip);

if(defined $pslFn){
	if(-d $pslFn){
		for($i = 0; ; $i++){
			if(-e "$pslFn/$i.pslx"){
				$pslList[$i] = "$pslFn/$i.pslx";
			}
			elsif($i != 0){ last; }
		}
	}
	else { $pslList[0] = $pslFn;}
}


$minErrNum = 1000;
$secMinErrNum = 1000;
$qName = "";
$inS = openInput($querySeqFn);
$out = openOutput($outFn);
$outMs = openOutput($outFn . ".microsat") if(defined $microFlag);

print "# $cmd\n"  if(defined $verbose);

my $ii = 0; ### for test.

$buf = <$inS>;

my $bTestStop = 0;
my $maxUnaligned = 30;
my $cutoff;
for(my $p = 0; $p <= $#pslList; $p++){
	next if(!defined $pslList[$p]);

	$in = openInput($pslList[$p]);
	#<$in>;<$in>;<$in>;<$in>;<$in>;
	while(<$in>){
		s/[\r\n]+//g;
		next if(/^#/ || /^\s*$/ || /^[a-z\-\s]/);

		$str = $_;
		## [0]:match, [1]:mismatch, [2]rep. match [3]N, [4]query gap num., [5]query gap bases, 
		## [6]ref gap num., [7]ref gap bases, [8]strand, [9]q name, [10]q size, [11]q start, [12]q end, 
		## [13]r name, [14]r size, [15]r start, [16]r end, 
		## [17], block num. [18]block sizes, [19]q starts, [20]r starts 
		@arr = split /\t/;

		next if($#arr < 10 
			|| ($arr[11] > 10 && $arr[10]-$arr[12] > 10) #### both match ends are not near read ends.
			|| $skipQuery{$arr[9]});

		
		#if($arr[9] eq "6777:F:58028:58127:R:58138:58180:M:75_2"){ print "===-------------$arr[13] --- 6777:F:58028:58127:R:58138:58180:M:75_2 \n, ($arr[11] > $maxUnaligned)\n"; }
						
		if(!defined $qName || $qName ne $arr[9]){
			getBestMapping() if($qName ne "");
	
			$bTestStop = 0;

			$#mapList = -1;
			$minErrNum = 1000;
			$secMinErrNum = 1000;
					
			####### get a sequence for a query.
			while(defined $buf){
				if($buf =~ /^@([\w_+\-:\.\|#\/]+)/)
				{
					if($1 eq $arr[9]){				
						$oriSeq = <$inS>;
						$oriSeq =~ s/[\r\n]+//g;
						$revSeq = rc($oriSeq);
						<$inS>;  ### pass '+'
						$qual = <$inS>;
						$qual =~ s/[\r\n]+//g;
						$buf = <$inS>;
						last;
					}
				}
				elsif($buf =~ /^>([\w_+\-:\.\|#\/]+)/){
					if($1 eq $arr[9]){				
						$oriSeq = <$inS>;
						$oriSeq =~ s/[\r\n]+//g;
						$revSeq = rc($oriSeq);
						$buf = <$inS>;
						$qual = "*";
						last;
					}
				}
				$buf = <$inS>;
			}
			
			print "Testing $arr[9]\n" if(defined $verbose);
		}

		
		$strand = $arr[8];
		$qName = $arr[9];
		$qLen = $arr[10];
		### $arr[11] : starts from 0, $arr[12] : starts from 1
		($qStart, $qEnd) = ($arr[8] eq "+" ? ($arr[11], $arr[12]-1) : ($arr[10]-$arr[12], $arr[10]-$arr[11]-1));
		$ref = $rName = $arr[13];
		### $arr[15] : starts from 0, $arr[16] : starts from 1
		$rStart = $arr[15];
		$rEnd = $arr[16]-1;

		
		#if($qName =~ "R:6:69:2300:16514"){
		#	print "$_\n";exit(1);
		#}

		$errNum = $arr[1];
		$isRelatedMic = 0;
		$gapNum = 0;
		if($arr[8] eq "+") {$seq = $oriSeq;}
		else {$seq = $revSeq; } 

		if(!defined $arr[18]){
			exitMsg(__LINE__, "undefined \$arr[18]\n$_\n");
		}

		next if(!defined $allMicInRef{$rName});
		$arrMics = $allMicInRef{$rName};
		$j = int($rStart/10000);
		$micI = $idxMicPos{$rName}[$j];
		### note all microsatelite positions are startinng from 0....
		while(defined $micI && defined $arrMics->[$micI] && $arrMics->[$micI]{end} < $rStart){ #### before reference start..
			$micI++;
		}	
		next if(!defined $micI || !defined $arrMics->[$micI] || $arrMics->[$micI]{start} >= $rEnd);
		
		next if(
			($arr[11] > $maxUnaligned && 
						($strand eq "+" ? ($arr[15] < $arrMics->[$micI]{start} - 3)  ## head clipped..
						: ($arr[16] > $arrMics->[$micI]{end} + 3)) ) ## tail clipped
			|| 
			($arr[10]-$arr[12] > $maxUnaligned && 
						($strand eq "+" ? ($arr[16] > $arrMics->[$micI]{end} + 3) 
						:	($arr[15] < $arrMics->[$micI]{start} - 3)) ) ); 
		
		#if($arr[9] eq "6777:F:58028:58127:R:58138:58180:M:75_2"){ print "===------------- $ref --- 6777:F:58028:58127:R:58138:58180:M:75_2 \n, ".
		#	"($arr[11] > $maxUnaligned && $arr[15] < $arrMics->[$micI]{start} - 3) #|| ($arr[10]-$arr[12] > $maxUnaligned && $arr[16] > $arrMics->[$micI]{end} + 3)\n"; }

		my @blocks = split /,/, $arr[18];
		my @QStarts = split /,/, $arr[19];  #### their positions start from 0.
		my @RStarts = split /,/, $arr[20];  #### their positions start from 0.
		my @QSeq = split /,/, $arr[21];
		my @RSeq = split /,/, $arr[22];

		next if($#blocks > 2); #### if there are gaps more than 2... 
		
		my @micAddLengths;
		my @micList;

		$isSkip = 0;
		##########################################################################
		#### reads microsats in a query sequence.
		my (@qArrMic, $qm, $qmNum);
		
		if(defined $allMicInQuery{$qName}){	
			@arr = split /[\-,]/, $allMicInQuery{$qName}; 

			### note all microsatelite positions are startinng from 0....
			for($i = 0, $j = 0; $i <= $#arr; $i+= 2, $j++){
				($qArrMic[$j]{start}, $qArrMic[$j]{end}) = $strand eq "+" ? ($arr[$i]-1, $arr[$i+1]-1) : ($qLen-$arr[$i+1], $qLen-$arr[$i]);
				$qArrMic[$j]{len} = $qArrMic[$j]{end} - $qArrMic[$j]{start} + 1;
			}
			@qArrMic = reverse @qArrMic if($strand eq "+");
		}		

		### note all microsatelite positions are startinng from 0....
		#### only consider reads including at least two bases of flanking at the both ends.
		next if( ($#qArrMic == 0 && $qArrMic[0]{len} > 20 && ($qArrMic[0]{start} < 3 || $qLen - $qArrMic[0]{end} < 2)) || ### only one microsats and one of both ends does not end in the read.
			($#qArrMic == 1 && ($qArrMic[0]{start} < 3 && $qLen - $qArrMic[1]{end} < 2)) );  ### only two microsats and ends of both microsat do not end in the read.
#exit if($qName =~ "R:6:69:2300:16514");
		$qm = 0;
		$qmNum = $#qArrMic + 1;
		##########################################################################		
				
		$num = $#blocks + 1;	
		#if($qName =~ "R:6:69:2300:16514"){print "num : $num, rStart, $RStarts[0]\n"; }
		#### to prevent this kind of alignment -- aaaaaaaaaaaa, aaaaaaaaaaa, aaaaaaaacagtgcgaggctatg or caacagtgcgaggctatgtcaaaaaaaaaaaaaaa, aaaaaaaaaaa, aaaaaaaaaaa 		
		my $rPrevEnd;
		my $tmpI = $micI;
		for($i = 0; $i < $num; $i++){
			$qPos = $QStarts[$i] + $blocks[$i];
			$rPos = $RStarts[$i] + $blocks[$i];

			while($qm < $qmNum && $qArrMic[$qm]{end} <= $QStarts[$i]){ $qm++; }
			while(defined $arrMics->[$tmpI] && $i != $num-1 && $arrMics->[$tmpI]{end} < $RStarts[$i+1]){ $tmpI++;	} ##### skip microsats before $rStart(ref start)

			#if($qName eq "HWUSI-EAS381R:5:1:8369:972#0" && $rName eq "chr2_230252180_230253762")
			#{ print "qm : $qm, $qArrMic[$qm]{start}-1 <= $QStarts[$i] && $qPos-1 <= $qArrMic[$qm]{end}\n($i != 0 && $RStarts[$i-1]+$blocks[$i-1] < $RStarts[$i]-5) || ($i != $num-1 && $rPos < $RStarts[$i+1]-5)\n"; }
			
			#print "$rPos, $arrMics->[$tmpI]{start} <= $RStarts[$i+1]+1 && $RStarts[$i+1]+1 <= $arrMics->[$tmpI]{end} && $blocks[$i] <= 10)\n" if($qName =~ "R:6:69:2300:16514");
			#$qArrMic[$qm]{start} - position from 1, $QStarts[$i] - position from 0
			if(	
				(
					($num -1 == 0) 
					|| ($i != $num-1 && $rPos < $RStarts[$i+1]-5) ### this is for forward direction
				) &&
				(
					($qm < $qmNum && $qArrMic[$qm]{start} <= $QStarts[$i]+5 && $qPos-1 <= $qArrMic[$qm]{end}) ||  ###### $qArrMic[$qm]{start}-1 || end +1 is for 1 base extension of microsat.. ex) aaaaaaag vs aaaaaaag match
					(defined $arrMics->[$tmpI] && $i != $num-1 &&    #### if start pos. of next block is in a microsat.
						$arrMics->[$tmpI]{start} <= $RStarts[$i+1]+1 && $RStarts[$i+1]+1 <= $arrMics->[$tmpI]{end} && $blocks[$i] <= 10 
							&& ($rPos < $arrMics->[$tmpI]{start}-3 || $RStarts[$i] >= $arrMics->[$tmpI]{start}))
				)
			){  
				### if a block is included in a microsatellite at the query read, remove the block.. ### if distance between two blocks should be bigger than 5		
				for($j = $i; $j < $num; $j++){
					$blocks[$j] = $blocks[$j+1];
					$QStarts[$j] = $QStarts[$j+1];
					$RStarts[$j] = $RStarts[$j+1];
					$QSeq[$j] = $QSeq[$j+1];
					$RSeq[$j] = $RSeq[$j+1];
				}
				$#blocks--;	
				$num--;
			}
			else{
				last;
			}
		}
		#if($qName =~ "R:6:69:2300:16514"){print "num : $num, isSkip : $isSkip, rStart, $RStarts[0]\n"; }
		next if($isSkip == 1 || $num == 0);
		
		### for reverse direction.
		$qm = $qmNum-1;
		for($i = $num-1; $i >= 0; $i--){
			$qPos = $QStarts[$i] + $blocks[$i];
			$rPos = $RStarts[$i] + $blocks[$i];
			$rPrevEnd = $RStarts[$i-1] + $blocks[$i-1] if($i != 0);

			while($qm >= 0 && $qArrMic[$qm]{start} >= $qPos){ $qm--;  }
			while(defined $arrMics->[$tmpI] && $i != 0 && $arrMics->[$tmpI]{end} < $rPrevEnd){ $tmpI++;	} ##### skip microsats before $rStart(ref start)
			
			if(	
				(
					($num -1 == 0) 
					|| ($i > 0 && $rPrevEnd < $RStarts[$i]-5)  ### this is for reverse direction
				) &&
				(
					($qm >= 0 && $qm < $qmNum && $qArrMic[$qm]{start} <= $QStarts[$i] && $qPos-1-5 <= $qArrMic[$qm]{end}) ||  ###### $qArrMic[$qm]{start}-1 || end +1 is for 1 base extension of microsat.. ex) aaaaaaag vs aaaaaaag match
					(defined $arrMics->[$tmpI] && $i != 0 &&    #### if end pos. of previous block is in a microsat.
						$arrMics->[$tmpI]{start} <= $rPrevEnd && $rPrevEnd <= $arrMics->[$tmpI]{end} && $blocks[$i] <= 10 
							&& ($rPos <= $arrMics->[$tmpI]{end} || $RStarts[$i] > $arrMics->[$tmpI]{end}+3))
				)
			){
				### if a block is included in a microsatellite at the query read, remove the block.. ### if distance between two blocks is bigger than 5

				for($j = $i; $j < $num; $j++){
					$blocks[$j] = $blocks[$j+1];
					$QStarts[$j] = $QStarts[$j+1];
					$RStarts[$j] = $RStarts[$j+1];
					$QSeq[$j] = $QSeq[$j+1];
					$RSeq[$j] = $RSeq[$j+1];
				}
				$#blocks--;	
				$num--;
			}
			else{
				last;
			}
		}	
		#if($qName =~ "R:6:69:2300:16514"){print "num : $num, isSkip : $isSkip, rStart, $RStarts[0]\n"; }
		next if($isSkip == 1 || $num == 0);

		$qStart = $QStarts[0];
		$rStart = $RStarts[0];
		$qEnd = $QStarts[$#blocks] + $blocks[$#blocks] -1;
		$rEnd = $RStarts[$#blocks] + $blocks[$#blocks] -1;


		######################################################
		### need to add codes to figure out whether unaligned parts are microsats or not
		########################## for the left unaligned sequence.
		### note all microsatelite positions are startinng from 0....
		if($qStart > 2){ ### at least three...
			testClippedHead(\@blocks, \@QStarts, \@RStarts, \@QSeq, \@RSeq, \$num);	
		}
		elsif($qStart > 0){
			$qm = 0;
			if($qLen - $qStart < 12){ next; }# total match length should be at least 12bases.

			while($qm < $qmNum && $qArrMic[$qm]{end} < $qStart){ $qm++; } ##### skip query microsats before $qStart(query start)
			while(defined $arrMics->[$micI] && $arrMics->[$micI]{end} < $rStart){ $micI++;	} ##### skip microsats before $rStart(ref start)

			if(defined $arrMics->[$micI] && $arrMics->[$micI]{start} <= $rStart){ #### if 5' unaligned seq is releated to microsat...
				$matchLen = 0;
				
				if($qm < $qmNum && $qArrMic[$qm]{start} <= $qStart+1) {
					$matchLen = $qStart+1 - $qArrMic[$qm]{start};
				}
				
				$errNum += $qStart - $matchLen;
				$micAddLengths[$#micAddLengths+1] = $matchLen;
				$micList[$#micList+1] = $micI;
				$isRelatedMic = 1;
			}
			else { $errNum += $qStart; } 
			$gapNum++;
		}
		#if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2"){ print "-------------, $rName, $arrMics->[$micI]{start} <= $rStart ------err : $errNum, 6777:F:58028:58127:R:58138:58180:M:75_2 \n"; }
		
		my $qTail = $qLen - $qEnd - 1;
		if($qTail > 2){	### at least three...;
			testClippedTail(\@blocks, \@QStarts, \@RStarts, \@QSeq, \@RSeq, \$num);
		}
		elsif($qTail > 0){
			$qm = 0;
			if($qLen - $qStart < 12){ next; }# total match length should be at least 12bases.

			while($qm < $qmNum && $qArrMic[$qm]{end} < $qStart){ $qm++; } ##### skip query microsats before $qStart(query start)
			my $tmpI = $micI;
			while(defined $arrMics->[$tmpI] && $arrMics->[$tmpI]{end} < $rEnd){ $tmpI++;	} ##### skip microsats before $rStart(ref start)

			if(defined $arrMics->[$tmpI] && $arrMics->[$tmpI]{start} <= $rEnd){ #### if 5' unaligned seq is releated to microsat...
				$matchLen = 0;
				
				if($qm < $qmNum && $qArrMic[$qm]{start} <= $qEnd) {
					if(!defined $qArrMic[$qm]{end} || !defined $qEnd){
						exitMsg(__LINE__, "$qName, $qm, $qmNum, start : $qArrMic[$qm]{start}, $qArrMic[$qm]{end} - $qEnd\n");
					}
					$matchLen = $qArrMic[$qm]{end} - $qEnd;
				}

				$errNum += $qTail - $matchLen;
				$micAddLengths[$#micAddLengths+1] = $matchLen;
				$micList[$#micList+1] = $tmpI;
				$isRelatedMic = 1;
			}
			else { $errNum += $qTail; } 
			$gapNum++;
		}
		###################################
	
		#if($qName =~ "AS381R:6:33:1898:12992"){print "$qName, num : $num, $qStart, isSkip : $isSkip, errNum $errNum\n";}
		
		################################################## for middle gaps	
		$num = $#blocks + 1;
		my ($rGap, $qGap, $qGapEnd, $longerGap, $numLongGap);
		$numLongGap = 0; #### to allow only one long gap longer than 5.
		$qm = 0;
		for($i = 0; $i < $num; $i++){
			$isMic = 0;
			$micAddLengths[$#micAddLengths+1] = 0;
			$micList[$#micList + 1] = -1;
			$qPos = $QStarts[$i] + $blocks[$i];
			$rPos = $RStarts[$i] + $blocks[$i];
			($rGap, $qGap) = ($i < $num -1 ? ($RStarts[$i+1]-$rPos, $QStarts[$i+1]-$qPos) : (0, $qLen - $qPos) ) ;

			#if($qName =~ "R:6:69:2300:16514"){print "i $i, rStart, $RStarts[$i], rPos : $rPos\n";}
			
			$longerGap = max($rGap, $qGap);
			if($longerGap > 5){
				$numLongGap++;
				#### only 1 long gap is allowed...
				if($numLongGap > 1){
					$isSkip = 1; 
					last; ### exit For statement.
				}
			}

			while(defined $arrMics->[$micI] && $arrMics->[$micI]{end} < $RStarts[$i]){ $micI++;	} ##### skip microsats before $pos
			if(defined $arrMics->[$micI] && $arrMics->[$micI]{end} < $rPos){
				#if($i == 0 && $QStarts[$i] > 5){ print "$_\n"; exit(1);} #### TEST
				#if($qName =~ "AS381R:6:33:1898:12992"){print " test mic : $arrMics->[$micI]{start}-$arrMics->[$micI]{end}\n";}

				$micList[$#micList] = $micI;
				$isRelatedMic = 1;
				$micAddLengths[$#micAddLengths+1] = 0;
				$micList[$#micList + 1] = -1;
			}

			while($qm < $qmNum && $qArrMic[$qm]{end} < $qPos){ $qm++;	} ##### skip microsats before $pos
			while(defined $arrMics->[$micI] && $arrMics->[$micI]{end} < $rPos-1){ $micI++;	} ##### skip microsats before $pos
				
			if(defined $arrMics->[$micI] && $arrMics->[$micI]{start} <= $rPos){ ### it means that the block end position is in the middle of microsatellites.
					
				if($rGap < $qGap && $rGap <= 3){	 #### if this is insertion ...
					#### get the match length...
					#if($qName =~ "AS381R:6:33:1898:12992"){print " 1. errNum $errNum, rGap :$rGap\n";}
					$errNum += $rGap;

		#if($qName =~ "R:6:69:2300:16514"){print "i $i, 1. errNum : $errNum\n";}
					$matchLen = 0;
					$qGapEnd = ($i == $num -1 ? $qLen : $QStarts[$i+1]);
					
					### mic covers the gap completely..
					if($qm < $qmNum && $qArrMic[$qm]{start} <= $qPos && $qGapEnd <= $qArrMic[$qm]{end}){
						$matchLen = $qGap;
						$isMic = 1;
					}
					### overlapping... microsat - gap. matchLen is from the overlapping part
					elsif($qm < $qmNum && $qArrMic[$qm]{start} <= $qPos && $qPos < $qArrMic[$qm]{end}){
						$matchLen = $qArrMic[$qm]{end} - $qPos;
					}
					### overlapping... gap - microsat.
					elsif($qm < $qmNum && $qArrMic[$qm]{start} < $qGapEnd && $qGapEnd <= $qArrMic[$qm]{end}){
						$matchLen = $qGapEnd - $qArrMic[$qm]{start};
					}

					#$isMic = 1; ##### Need to test....

					$gapLen = $qGap - $matchLen;

					#if($qName =~ "AS381R:6:33:1898:12992"){print " 1. errNum $errNum, isMic : $isMic, min($gapNum, $gapLen), if($qm < $qmNum && $qArrMic[$qm]{start} <= $qPos && $qGapEnd <= $qArrMic[$qm]{end})\n";}
					if($isMic == 1){
						$errNum += min($gapNum, $gapLen);
					}
					else{
						$errNum += max($gapNum, $gapLen);
					}	

		#if($qName =~ "R:6:69:2300:16514" && $rName eq "chr5:93207559:93207762"){print "2. errNum : $errNum\n";}
					$micList[$#micList] = $micI; 
					$micAddLengths[$#micAddLengths] = ($i == $num -1 ? $matchLen : $qGap);					

				}  #### if($i < $num-1 && $rGap < $qGap)
				elsif($qGap < $rGap && $qGap <= 3){  ### deletion
					########################################		
					$errNum += $qGap;
					if($arrMics->[$micI]{start} <= $rPos && $RStarts[$i+1] <= $arrMics->[$micI]{end}){ ### if a microsat in a reference continues until the next block..
						$matchLen = $rGap;
						$isMic = 1;
					}
					elsif($arrMics->[$micI]{start} <= $rPos && $arrMics->[$micI]{end} < $RStarts[$i+1]){ ### if a microsat in a reference ends before the next block..
						$matchLen = $arrMics->[$micI]{end} - $rPos;
					}	

					#$isMic = 1; ##### Need to test....

					$gapLen = $rGap - $matchLen;
					
					$micAddLengths[$#micAddLengths] = -$rGap; 
					$micList[$#micList] = $micI; 
					if($isMic == 1){
						$errNum += min($gapNum, $gapLen);
					}
					else{
						$errNum += max($gapNum, $gapLen);
					}
					#print "---------------------m\n" if($rName eq "chr13_45941212");

				}
				else{ #### both gaps are bigger than 3
					#$errNum += $blocks[$i];
					$isSkip = 1; #print "both gaps are bigger than 3 $qGap, $rGap\n";
				}			
				$isRelatedMic = 1;
			} #### if(defined $arrMics->[$micI] && $arrMics->[$micI]{start} <= $rPos)

			#elsif(defined $arrMics->[$micI] && $rGap > 0 && $arrMics->[$micI]{start} <= $RStarts[$i+1]){   ### TEST
			elsif(defined $arrMics->[$micI] && $rGap > 0 && $qGap <= 3 && $arrMics->[$micI]{start} <= $RStarts[$i+1]){   ### if a microsat in a reference starts between the current block and the next block.
					
				$matchLen = 0;
				if($RStarts[$i+1] <= $arrMics->[$micI]{end}){	  ### if a microsat in a reference lasts until the next block..
					$matchLen = $RStarts[$i+1] - $arrMics->[$micI]{start};
				}
				else { ####  $rPos <  $arrMics->[$micI]{start}, $arrMics->[$micI]{end} < $RStarts[$i+1]  ##   ### if a microsat in a reference ends before the next block..
					$matchLen = $arrMics->[$micI]{end} - $arrMics->[$micI]{start} + 1;
				}

				$gapLen = $rGap - $matchLen;

				$micAddLengths[$#micAddLengths] = -$rGap; 
				$micList[$#micList] = $micI; 
				$errNum += $gapLen;
				$isRelatedMic = 1;
			}
			elsif($qGap > 3 && $rGap > 3){  				
				$isSkip = 1; 
			}
			else{ $errNum += max($qGap, $rGap);}

			$gapNum++;	
		}
		#if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2"){ print "---------------88-- isSkip $isSkip --err : $errNum, qStart $qStart, rStart $rStart, 6777:F:58028:58127:R:58138:58180:M:75_2 \n"; }
		
		#if($qName =~ "6777:F:58028:58127:R:58138:58180:M:75_2"){print "$qName, num : $num, $qStart, isSkip : $isSkip, errNum $errNum, \n$_\n"; exit(1); }
		next if($isSkip == 1);		
		
		$i = $#mapList + 1;
		$mapList[$i]{errNum} = $errNum;
		$mapList[$i]{str} = $str;
		$mapList[$i]{qStart} = $qStart+1;	###### change to 1 based starting..
		$mapList[$i]{qEnd} = $qEnd+1;	###### change to 1 based starting..
		$mapList[$i]{strand} = $strand;

		$mapList[$i]{rName} = $rName;
		$mapList[$i]{rStart} = $rStart+1; ###### change to 1 based starting..
		$mapList[$i]{rEnd} = $rEnd+1;	###### change to 1 based starting..
		$mapList[$i]{relatedMic} = $isRelatedMic;
		$mapList[$i]{blocks} = \@blocks;
		$mapList[$i]{QStarts} = \@QStarts;
		$mapList[$i]{RStarts} = \@RStarts;
		$mapList[$i]{QSeq} = \@QSeq;
		$mapList[$i]{RSeq} = \@RSeq;
		$mapList[$i]{micAddLengths} = \@micAddLengths;
		$mapList[$i]{micList} = \@micList;

		
		if($minErrNum > $errNum){
			$secMinErrNum = $minErrNum;
			$minErrNum = $errNum;
			#print "$qName, errNum : $errNum, minErrNum : $minErrNum, secMinErrNum : $secMinErrNum\n" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");			
		}

		

	}
	close($in);
}
getBestMapping();

close($inS);
close($out);
close($outMs) if(defined $microFlag);


sub testClippedHead
{
	my ($blocks, $QStarts, $RStarts, $QSeq, $RSeq, $num) = @_;
	my $sLen; #### length of a query part to be compared..
		

	my $bOver = 0;
	if($qStart > 5){				
		$sLen = $qStart;
	}
	else{
		$sLen = 6;
		$bOver = 1;
	}

	my $searchMax = 40 + $sLen;
	my $tmpI = $micI;
	while(defined $arrMics->[$tmpI] && $arrMics->[$tmpI]{end} < $rStart-$searchMax){ $tmpI++;	} ##### skip microsats before a search start pos.


	my $sSeq = substr($seq, 0, $sLen); ### query part
	#print "$tSeq\n$sSeq\n";
	my ($minI, $minMis, $mis) = (-1, 40, 0); ### number of mismatch.
	#my @sArr = split //, $sSeq;
	#my @tArr = split //, $tSeq;

	
	my ($i, $j, $k, $tmpSEnd);
	for($i = $rStart -1, $k = 0; $i >= 0 && $k < $searchMax; $i--, $k++){				
		#if($i+$sLen >= $rStart){ $mis = 0; }
		#elsif(defined $arrMics->[$tmpI]){
		$tmpSEnd = $sLen;

		##### if the gap bases are not overlapping with microsatellites, they will be counted as errors($mis).
		#### later, the gap bases(query bases) overlapping with microsatellites will be compared(from tmpSStart to tmpSEnd..)
		my $compEnd = $i+$sLen-1;
		if($bOver == 1 && $rStart <= $compEnd ){ $mis = 0; }
		elsif(defined $arrMics->[$tmpI]){			
			if($arrMics->[$tmpI]{start} <= $compEnd && $arrMics->[$tmpI]{start} - 3 <= $rStart 
				&& $rStart <= $arrMics->[$tmpI]{end}){ ### all bases are in a microsatellite  ## && $compEnd <= $arrMics->[$tmpI]{end}
				$mis = 0;
				$tmpSEnd = $arrMics->[$tmpI]{end} - $rStart if($arrMics->[$tmpI]{end} - $rStart < $tmpSEnd);
				#print "0, mis : $mis\n" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");
			}
			elsif($rStart <= $compEnd 
				|| ($arrMics->[$tmpI]{end} < $compEnd && $arrMics->[$tmpI]{end} < $rStart) 
				|| ($compEnd < $arrMics->[$tmpI]{start} && $rStart < $arrMics->[$tmpI]{start}) ){
				$mis = $rStart - $compEnd - 1;
				#print "1, mis : $mis\n" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");				
				if($mis < 0){
					$mis *= -1;
					$tmpSEnd = $sLen - $mis;
				}
			}
			elsif($compEnd <= $arrMics->[$tmpI]{start} && $arrMics->[$tmpI]{end} <= $rStart){
				$mis = $arrMics->[$tmpI]{start}- $compEnd - 1 + $rStart - $arrMics->[$tmpI]{end} - 1;
				$mis = 0 if($mis < 0);
				#print "2, mis : $mis\n" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");
			}
			elsif($compEnd <= $arrMics->[$tmpI]{start} && $arrMics->[$tmpI]{start} <= $rStart){
				$mis = $arrMics->[$tmpI]{start}- $compEnd - 1;
				$mis = 0 if($mis < 0);
				#print "3, mis : $mis\n" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");
			}
			elsif($compEnd <= $arrMics->[$tmpI]{end} && $arrMics->[$tmpI]{end} <= $rStart){
				$mis = $rStart - $arrMics->[$tmpI]{end} - 1;
				$mis = 0 if($mis < 0);
				#print "4, mis : $mis\n" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");
			}
			else{
				exitMsg(__LINE__, "$qName, $rName\n");
			}
		}
		else {			
			$mis = $rStart - $compEnd - 1;
			#print "5, mis : $mis\n";
			if($mis < 0){
				$mis *= -1;
				$tmpSEnd = $sLen - $mis;
			}
		}

		last if($mis > 10 && $compEnd < $rStart);

		
		#print "i :$i, rStart : $rStart, sLen : $sLen, tmpSEnd $tmpSEnd, mic($arrMics->[$tmpI]{start}-$arrMics->[$tmpI]{end}) " . 
		#	substr($refSeq{$ref}, $i, $tmpSEnd) . ", $sSeq, mis($mis->" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");

		for($j = 0; $j < $tmpSEnd; $j++){			
			#if($sArr[$j] ne $tArr[$i+$j]){
			if(substr($sSeq, $j, 1) ne substr($refSeq{$ref}, $i+$j, 1)){
				$mis++;
			}
			last if(($sLen > 20 && $mis > 2) || ($sLen <= 20 && $mis > 1));
		}
		if($mis == 0){ $minI = $i; $minMis = 0; last;}
		if( ($mis == 1 || ($sLen > 20 && $mis == 2)) && ($mis < $minMis || $minI == -1)){ $minI = $i;	$minMis = $mis;}
		#print "$mis)\n" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");
	}
	#print "\nmStart : " .($minI)."\n" if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2");

	#exit(1);
	#if($qName eq "6777:F:58028:58127:R:58138:58180:M:75_2"){ print "-------------, $minI ------err : $errNum, 6777:F:58028:58127:R:58138:58180:M:75_2 \n"; }

	if($minI != -1){
		#### found match...		
		### $tStart : 0 based
		my $mStart = $minI; ### match start at the reference..
		#print "if($qStart == $rStart - $mStart(".($rStart - $mStart).")\n";
		#### if extensions of query and reference are same..
		if($qStart == $rStart - $mStart){
			$blocks->[0] += $qStart;
			$QStarts->[0] = 0;
			$RStarts->[0] -= $qStart;
			$QSeq->[0] = substr($seq, 0, $qStart) . $QSeq->[0];
			$RSeq->[0] = substr($refSeq{$ref}, $mStart, $qStart) . $RSeq->[0];	
		}
		else{			
			#for($j = 0; $j < $$num; $j++){
			#	print "$j : q : $QStarts->[$j]-".($QStarts->[$j]+$blocks->[$j]).", r : $RStarts->[$j]-".($RStarts->[$j]+$blocks->[$j])."\n";
			#}

			for($j = $$num; $j > 0; $j--){
				$blocks->[$j] = $blocks->[$j-1];
				$QStarts->[$j] = $QStarts->[$j-1];
				$RStarts->[$j] = $RStarts->[$j-1];
				$QSeq->[$j] = $QSeq->[$j-1];
				$RSeq->[$j] = $RSeq->[$j-1];
			}
								
			my $rDiff = $rStart - $mStart;
			if($qStart < $rDiff){
				$blocks->[0] = $qStart;
				$QSeq->[0] = substr($seq, 0, $qStart);						
				$RSeq->[0] = substr($refSeq{$ref}, $mStart, $qStart);
			}
			else{
				$blocks->[0] = $rDiff;
				$QSeq->[0] = substr($seq, 0, $rDiff);
				$RSeq->[0] = substr($refSeq{$ref}, $mStart, $rDiff);
			}
			
			$$num++;
		}
		$qStart = $QStarts->[0] = 0;
		$rStart = $RStarts->[0] = $mStart;
		
		$errNum += $minMis;
		#for($j = 0; $j < $$num; $j++){
		#	print "$j : q : $QStarts->[$j]-".($QStarts->[$j]+$blocks->[$j]).", r : $RStarts->[$j]-".($RStarts->[$j]+$blocks->[$j])."\n";
		#}
	}
	### if it did not find match, then clipped head is counted as error
	else{
		$errNum += $qStart;
	}
}


sub testClippedTail
{
	my ($blocks, $QStarts, $RStarts, $QSeq, $RSeq, $num) = @_;
		
	my $sLen; #### length of a query part to be compared..
	
	my $bOver = 0;
	my $qTail = $qLen - $qEnd - 1;
	if($qTail > 5){				
		$sLen = $qTail;
	}
	else{
		$sLen = 6;
		$bOver = 1; ### if it is 1, allow to overlap..
	}
	
	my $refLen = length($refSeq{$ref});
	my $searchMax = 40 + $sLen;
	my $tmpI = $micI;
	while(defined $arrMics->[$tmpI] && $arrMics->[$tmpI]{end} < $rEnd - $sLen){ $tmpI++;	} ##### skip microsats before a search start pos.

	my $sSeq = substr($seq, $qLen-$sLen); ### query part

	my $tStart = $rEnd-$sLen+1;
	my $tLen = ($sLen+40) + $sLen;
	my $tSeq = substr($refSeq{$ref}, $tStart, $tLen); ### reference part
	$tLen = length($tSeq);

	#print "$sSeq\n$tSeq\n";
	my ($minI, $minMis, $mis) = (-1, 40, 0); ### number of mismatch.
	if(!defined $sSeq || !defined $tSeq){
		close($out);
		exitMsg(__LINE__, "No sequence .., $qName\n ref : $ref, length : " .(length($refSeq{$ref})). ", tStart : $tStart, tLen : $tLen\n undefined sSeq($sSeq) || tSeq($tSeq)\n");
	}
	
	#$arrMics->[$tmpI]{start} = 123;
	#$arrMics->[$tmpI]{end} = 132;##TEST
	#my @sArr = split //, $sSeq;
	#my @tArr = split //, $tSeq;
	my ($i, $j, $k, $tmpSStart, $tmpSEnd, $gStart, $gEnd);
	for($i = $rEnd-$sLen+2, $k = 0; $i < $refLen - $sLen && $k < $searchMax; $i++, $k++){
		$tmpSStart = 0;
		$tmpSEnd = $sLen;

		$mis = 0;

		##### if the gap bases are not overlapping with microsatellites, they will be counted as errors($mis).
		#### later, the gap bases(query bases) of overlapping with microsatellites will be compared(from tmpSStart to tmpSEnd..)
		if($bOver == 1 && $i <= $rEnd){
			$mis = 0;
		}
		elsif(defined $arrMics->[$tmpI]){
			if( $arrMics->[$tmpI]{start} <= $rEnd #&& $arrMics->[$tmpI]{start} <= $i 
				&& $rEnd - 3 <= $arrMics->[$tmpI]{end} && $i <= $arrMics->[$tmpI]{end} ) ### all bases are in a microsatellite				
			{ 
				$mis = 0;
				my $s = $sLen - ($rEnd - $arrMics->[$tmpI]{start});
				$tmpSStart = $sLen - ($rEnd - $arrMics->[$tmpI]{start}) if(($rEnd - $arrMics->[$tmpI]{start}) < $sLen);
				#print "0, mis : $mis\n";
			}
			elsif($i <= $rEnd 
				|| ($arrMics->[$tmpI]{end} < $rEnd && $arrMics->[$tmpI]{end} < $i) 
				|| ($rEnd < $arrMics->[$tmpI]{start} && $i < $arrMics->[$tmpI]{start}) ){
				$mis = $i - $rEnd - 1;
				#print "1, mis : $mis\n";
				if($mis < 0){
					$mis *= -1;
					$tmpSStart = $mis;
				}
			}
			elsif($rEnd <= $arrMics->[$tmpI]{start} && $arrMics->[$tmpI]{end} <= $i){
				$mis = $arrMics->[$tmpI]{start}- $rEnd - 1 + $i - $arrMics->[$tmpI]{end} - 1;
				$mis = 0 if($mis < 0);
				#print "2, mis : $mis\n";
			}
			elsif($rEnd <= $arrMics->[$tmpI]{start} && $arrMics->[$tmpI]{start} <= $i){
				$mis = $arrMics->[$tmpI]{start}- $rEnd - 1;
				$mis = 0 if($mis < 0);
				#print "3, mis : $mis\n";
			}
			elsif($rEnd <= $arrMics->[$tmpI]{end} && $arrMics->[$tmpI]{end} <= $i){
				$mis = $i - $arrMics->[$tmpI]{end} - 1;
				$mis = 0 if($mis < 0);
				#print "4, mis : $mis\n";
			}
			else{
				exitMsg(__LINE__, "$qName, $rName\n");
			}
		}
		else {			
			$mis = $i - $rEnd - 1;
			#print "5, mis : $mis\n";
			if($mis < 0){
				$mis *= -1;
				$tmpSStart = $mis;
			}
		}

		last if($mis > 10 && $rEnd < $i);


		#print "i :$i, rEnd : $rEnd, sLen : $sLen, tmpSStart $tmpSStart, " . (defined $arrMics->[$tmpI] ?"mic($arrMics->[$tmpI]{start}-$arrMics->[$tmpI]{end}), ":"") . 
		#	substr($refSeq{$ref}, $i+$tmpSStart, ($sLen-$tmpSStart)) . ", $sSeq, mis($mis->";
		for($j = $tmpSStart; $j < $sLen; $j++){
			
			#if($sArr[$j] ne $tArr[$i+$j]){
			if(substr($sSeq, $j, 1) ne substr($refSeq{$ref}, $i+$j, 1)){
				$mis++;
			}
			last if(($sLen > 20 && $mis > 2) || ($sLen <= 20 && $mis > 1));			
		}
		if($mis == 0){ $minI = $i; $minMis = 0; last;}
		if( ($mis == 1 || ($sLen > 20 && $mis == 2)) && ($mis < $minMis || $minI == -1)){ $minI = $i;	$minMis = $mis;}
		#print "$mis)\n";
	}
	#print "\nmStart : " .($minI).", minMis : $minMis\n";
	#if($minI != -1 && defined $arrMics->[$tmpI]){ print "$qName, $rName\n"; exit(1); }

	if($minI != -1){
		#### found match...

		### $tStart : 1 based .. note.. it is different from testClippedHead()
		my $mStart = $minI; ### match start at the reference..
		my $rDiff = $sLen - ($rEnd+1 - $mStart);
		#print "if($qTail == $sLen - ($rEnd - $mStart)($rDiff)\n";
		if($qTail == $rDiff){
			$blocks->[$$num-1] += $qTail;
			$QSeq->[$$num-1] .= substr($seq, $qEnd+1);
			$RSeq->[$$num-1] .= substr($refSeq{$ref}, $rEnd+1, $qTail);	
		}
		else{		
			##### deletion
			if($qTail < $rDiff){

				$QStarts->[$$num] = $qEnd + 1;
				$RStarts->[$$num] = $rEnd + ($rDiff-$qTail) + 1;
				$blocks->[$$num] = $qTail;
				$QSeq->[$$num] = substr($seq, $qEnd+1);						
				$RSeq->[$$num] = substr($refSeq{$ref}, $RStarts->[$$num], $qTail);
			}
			##### insertion
			else{
				$QStarts->[$$num] = $qEnd + ($qTail-$rDiff) + 1;
				$RStarts->[$$num] = $rEnd + 1;
				$blocks->[$$num] = $rDiff;
				$QSeq->[$$num] = substr($seq, $QStarts->[$$num]);
				$RSeq->[$$num] = substr($refSeq{$ref}, $rEnd+1, $rDiff);
			}
			
			$$num++;

			#print "$#$blocks : $QStarts->[$#$blocks], $RStarts->[$#$blocks], \n -- new, num : $$num\n";			
			#for($j = 0; $j < $$num; $j++){
			#	print "$j : q : $QStarts->[$j]-".($QStarts->[$j]+$blocks->[$j]).", r : $RStarts->[$j]-".($RStarts->[$j]+$blocks->[$j])."\n";
			#}
		}
		
		$qEnd = $QStarts->[$$num-1]+$blocks->[$$num-1] -1;
		$rEnd = $RStarts->[$$num-1]+$blocks->[$$num-1] -1;
		$errNum += $minMis;;
	}
	### if it did not find match, then clipped head is counted as error
	else{
		$errNum += $qStart;
	}
}


sub getBestMapping{
	my ($i, $num, $numHitMin, $numHitSecMin);
	my @relatedMic = ("N", "T");
	$num = $#mapList + 1;
	
	my $cutoff;
	if(defined $cutoffError && $cutoffError >= 1){ $cutoff = $cutoffError;}
	elsif(defined $cutoffError){ $cutoff =  $qLen * $cutoffError;}
	else { $cutoff = $qLen/20; }
	
	if($minErrNum == 100 || $num == 0){
		print "$qName : No satisfactory mapping..\n" if(defined $verbose); 
	}
	elsif($minErrNum > $cutoff){
		print "$qName : minErrNum is $minErrNum (bigger than $cutoff)\n" if(defined $verbose); 
		return;
	}	

	($numHitMin, $numHitSecMin) = (0,0);
	for($i = 0; $i < $num; $i++){
		if($mapList[$i]{errNum} == $minErrNum){ $numHitMin++;	}
		elsif($mapList[$i]{errNum} == $secMinErrNum){ $numHitSecMin++;	}
	}
	
	if($numHitMin > $allowdBestHitNum){
		print "$qName : numHitMin : $numHitMin > $allowdBestHitNum\n" if(defined $verbose); 
		return;
	}
	
	#### set mapping score
	my ($mapScore, $errDiff);
	$errDiff = $secMinErrNum - $minErrNum;
	if($numHitMin > 1) { $mapScore = 0; }
	elsif($errDiff > 3) { $mapScore = 37; }
	else{
		print "minErrNum : $minErrNum, secMinErrNum : $secMinErrNum, errDiff : $errDiff, numHitMin : $numHitMin, secMinErrNum : $secMinErrNum\n" if(defined $verbose);

		$secMinErrNum *= (3-$errDiff+1)*5;
		$secMinErrNum = 255 if($secMinErrNum > 255);
		
		if(!defined $g_log_n[$secMinErrNum]){
			exitMsg(__LINE__, "$qName, minErrNum : $minErrNum, secMinErrNum : $secMinErrNum, $g_log_n[$secMinErrNum]\n");
		}

		$mapScore = (23 < $g_log_n[$secMinErrNum])? 0 : 23 - $g_log_n[$secMinErrNum];
		print "mapScore : $mapScore\n" if(defined $verbose); 
	}


	my $hitType = ($numHitMin == 1 ? "U" : "M");
	for($i = 0; $i < $num; $i++){
		if($mapList[$i]{errNum} == $minErrNum){

			writeSAM($qName, $mapList[$i]{rName}, $mapList[$i]{strand}, $oriSeq, $qual, $mapList[$i]{QStarts}, $mapList[$i]{RStarts}, $mapList[$i]{blocks},
				$mapList[$i]{QSeq}, $mapList[$i]{RSeq}, $hitType, $mapScore);

			if($microFlag && $mapList[$i]{relatedMic} == 1){
				print $outMs ">$qName\t$hitType\t$mapList[$i]{errNum}\t$qLen\t$mapList[$i]{qStart}\t$mapList[$i]{qEnd}\t$mapList[$i]{rName}\t$mapList[$i]{rStart}\t$mapList[$i]{rEnd}\n";

				my $j = 0;
				foreach my $mic ( @{$mapList[$i]{micList}} ) {
					print $outMs "$j-th, mic : $mic, $mapList[$i]{micAddLengths}[$j], $allMicInRef{$mapList[$i]{rName}}[$mic]{str}\n" if($mic != -1); 
					$j++;
				}
			}
		}
	}
}

sub writeSAM
{
	## <QNAME> <FLAG> <RNAME> <POS in ref> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [<TAG>..]
	my ($qName, $rName, $st, $seq, $qual, $qStarts, $rStarts, $blocks, $qSequences, $rSequences, $hitType, $mapScore) = @_;
	
	my ($flag, @arrTemp);
		
	$flag = 0;

	if($qName =~ /\/1/){
		$flag |= 0x1 if($pairedFlag);
		#$flag |= 0x40;
		$qName = $`;
	}
	elsif($qName =~ /\/2/ || $qName =~ /\/3/){
		$flag |= 0x1;
		$flag |= 0x80;
		$qName = $`;
	}

	if($hitType eq "M"){
		$flag |= 0x100;
	}

=head
	if(defined $in2) #### if pairs exist.
	{
		$flag |= 0x1;	## paired.
		if($arr[17] eq "" && $arr[13] ne "" && $arr[20] ne "" && $arr[13] ne $arr[20]){
			$flag |= 0x2;	## proper pair.
		}
		if($arr[10] eq "" || $arr[13] eq ""){
			$flag |= 0x4;	## itself is unmapped.
		}
		if($arr[20] eq ""){
			$flag |= 0x8;	## mate is unmapped.
		}
		if($pST ne "" && $pST eq "-"){
			$flag |= 0x20;
		}

		if($pair == 1) {
			$flag |= 0x40;
		}
		elsif($pair == 2){
			$flag |= 0x80;
		}
	}
=cut

	if($st eq "-"){
		$flag |= 0x10;
		$seq = rc($seq);
		$qual = reverse($qual);
		$arr[14] = revMatchInfo($arr[14]);
	}
	my $qLen = length($seq);

=head	### to chage quality score from illumina scale to phred scale
	@arrTemp = split //, $qual;
	$qual = "";
	foreach my $qu (@arrTemp) {
		$qual .= pack("c", ord($qu) - 31);  ### to convert Illumina quality(ASCII+64) to standard quality(ASCII+33). 
	}
=cut

	my ($rStart, $isINDEL, $mat, $str, $k, @cigar, $md, $mn, $gapLen, $mLen, $bNum, $i, $j, $delSeq, @qArr, @rArr, $qGap, $rGap, $blockLen);
	$isINDEL = 0;
	$mat = 0;
	$md = "";
	$mn = 0;
	$rStart = $rStarts->[0];

	$bNum = $#$blocks + 1;
	$mLen = 0;
	if($qStarts->[0] > 0) { 
		$cigar[0] = "$qStarts->[0]S"; 
		#$bTestStop = 1 if($qStarts->[0] > 2);
	}

	$blockLen = 0;
	for($i = 0; $i < $bNum; $i++){
		compareTwoSubStr($blocks->[$i], \$qSequences->[$i], \$rSequences->[$i], \$mLen, \$md, \$mn);
=head
		@qArr = split //, $qSequences->[$i];
		@rArr = split //, $rSequences->[$i];
		
		for($j = 0; $j < $blocks->[$i]; $j++){
			
			if($qArr[$j] eq $rArr[$j]){ $mLen++; }
			else{ 
				$md .= $mLen . uc($rArr[$j]);
				$mLen = 0; 
				$mn++;
			}
		}
=cut
		$blockLen += $blocks->[$i];
		if($i == $bNum -1){
			
			push(@cigar, $blockLen."M");
			$md .= $mLen if($mLen != 0);
			last; ######
		}

		($qGap, $rGap) = ($qStarts->[$i+1] - ($qStarts->[$i] + $blocks->[$i]), $rStarts->[$i+1] - ($rStarts->[$i] + $blocks->[$i]));
		if($qGap < $rGap ){	### deletion
			
			if($qGap != 0){
				$str = substr($refSeq{$rName}, $rStarts->[$i] + $blocks->[$i], $qGap); 				
				compareTwoSubStr(length($str), \substr($seq, $qStarts->[$i]+$blocks->[$i], $qGap), \$str, \$mLen, \$md, \$mn);
				#print "----- substr($refSeq{$rName}, $rStarts->[$i] + $blocks->[$i], $qGap) :  " . substr($refSeq{$rName}, $rStarts->[$i] + $blocks->[$i], $qGap) . "\n";
				$md .= $mLen;
				push(@cigar, ($blockLen+$qGap)."M");
				
			}
			else{
				push(@cigar, $blockLen."M");
				$md .= $mLen;
			}
			$gapLen = $rGap - $qGap;
			$delSeq = substr($refSeq{$rName}, $rStarts->[$i] + $blocks->[$i], $gapLen);
			$md .= "^$delSeq";
			$mLen = 0;
			
			push(@cigar, $gapLen ."D");
			$mn += $rGap;
			$blockLen = 0;
		}
		elsif($qGap > $rGap){ #### insertion
			if($rGap != 0){
				$str = substr($refSeq{$rName}, $rStarts->[$i] + $blocks->[$i], $rGap); 				
				compareTwoSubStr(length($str), \substr($seq, $qStarts->[$i]+$blocks->[$i], $rGap), \$str, \$mLen, \$md, \$mn);

				push(@cigar, ($blockLen+$rGap)."M");
				$md .= $mLen;
			}
			else{
				push(@cigar, $blockLen."M");				
			}

			$gapLen = $qGap - $rGap;
			push(@cigar, $gapLen."I");
			
			$mn += $qGap;
			#$mLen = 0; ####### insertions are not included in the MD count..
			$blockLen = 0;
		}
		else{  
			#my $qstr = ;
			#my $rstr = ;
			compareTwoSubStr($qGap, \substr($seq, $qStarts->[$i]+$blocks->[$i], $qGap), \substr($refSeq{$rName}, $rStarts->[$i]+$blocks->[$i], $rGap), \$mLen, \$md, \$mn);
			#print STDERR "Error in writing cigar!!!, $qGap, $rGap, $qName($qStarts->[$i]+$blocks->[$i] -> $qStarts->[$i+1]), $rName($rStarts->[$i]+$blocks->[$i] -> $rStarts->[$i+1])\n"; exit -1;
			$blockLen += $qGap;
		}
	}
	if($qStarts->[$i] + $blocks->[$i] < $qLen) { 
		push(@cigar, ($qLen - $qStarts->[$i]- $blocks->[$i]) . "S"); 		
	}


	## <QNAME> <FLAG> <RNAME> <POS in ref> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [<TAG>..]
	
	$rStart++; ### $rStart is 0 base. to make it 1 base.
	print $out "$qName\t$flag\t$rName\t$rStart\t$mapScore\t" . join("", @cigar) . "\t" .
			"*". #(defined $in2  && $arr[20] ne "N" ? ($arr[17] eq "" ? "=":$arr[17]):"*").															## <MRNM> reference for pair.
			"\t". "0" . #(defined $in2 ? ($arr[17] eq ""  && $arr[20] ne "N" ? $arr[12]+$arr[19] : $arr[19]) : "0") .		## <MPOS> reference for pair.
			"\t". "0" . #(defined $in2 && $arr[17] eq "" ? $arr[19] : "0") .													## <ISIZE> reference for pair.
			"\t$seq\t$qual" .																															## sequence , quality
			#"\tAS:i:$arr[15]" . 
			"\tNM:i:$mn" . "\tMD:Z:$md\tYB:Z:BL\n";	

	if($bTestStop == 1){
		print "Exit for test..\n";
		close($out); exit(1);
	}
}

sub compareTwoSubStr{
	my ($len, $qstr, $rstr, $mLen, $md, $mn) = @_;
	
	my @qArr = split //, $$qstr;
	my @rArr = split //, $$rstr;
		
	for(my $j = 0; $j < $len; $j++){
		
		if($qArr[$j] eq $rArr[$j]){ $$mLen++; }
		else{ 
			$$md .= $$mLen . uc($rArr[$j]);
			$$mLen = 0; 
			$$mn++;
		}			
	}
}


sub revMatchInfo
{
	my ($info) = @_;
	return undef if(!defined $info);
	
	my ($alt, $str, $mat);
	$alt = "";
	while($info =~ /(\d+)|([A-Z]+)|\^|\$/g){
		$str = $&;
		if($str eq "^"){ ## indel start
			$alt = "\$" . $alt;
		}
		elsif($str eq "\$"){ ## indel end
			$alt = "^" . $alt;
		}
		elsif($str =~ /[A-Z]/){
			$alt = rc($str) . $alt;
		}
		else{  ### numbers.
			$alt = $str . $alt;
		}
	}

	return $alt;
}

sub exitMsg
{
	my ($line, $msg) = @_;
	$|=1;
	print STDERR "[Error] at the line $line : $msg\n" if($msg);
	exit(1);
}



sub min
{
	my ($i, $min);
	$min = $_[0];
	for($i = 1; $i <= $#_; $i++){
		$min = $_[$i]  if($min > $_[$i]);
	}
	return $min;
}


sub max
{
	my ($i, $max);
	$max = $_[0];
	for($i = 1; $i <= $#_; $i++){
		$max = $_[$i]  if($max < $_[$i]);
	}
	return $max;
}


sub openInput
{
	my ($fn) = @_;

	#return *STDIN unless defined $fn;
	
	print "reading '$fn'\n";
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



