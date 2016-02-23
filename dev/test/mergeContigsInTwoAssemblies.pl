#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2010
#
# mergeContigsByHomology.pl
# merges overlapping contigs
=head
1. If headers of two contig sets are same, you have to change one of them.
ex) If names of contigs from newbler start with 'contig' like 'contig00001' and names of contigs from CLCbio start w
ith 'contig' also, then you need to change names of contigs from CLCbio to 'CLC_contig'.

2. merge two fasta files(from two different assemblies) into one fasta file. Let's assume the merged fasta file as '
merged.fa'.

3. blast for itself. (it needs -m 8 -F F -e 1e-10 ). Let's assume output as 'merged.bln.raw'

4. ./mergeContigsInTwoAssemblies.pl -b merged.bln.table  -s merged.fa  -o merged.new.fa
-- If you need to add scaffold information, then you can use '-S' option.
-- If some contigs are not involved to the merge process, you can decide to keep contigs from one assembly. If you w
ant to keep contigs from CLCbio(their names start with 'CLC_contig'), then use '-keep CLC_contig'. But if you set sc
affold '-S' option, it will keep newbler contigs ('contig00~').
=cut
######################################


use strict;
use warnings "all";
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $verboseFlag, $lenMinOverlap, $identityMin, $endClippingAllowed, $blnFn, $seqFn, $scaffoldFn, $outFn,
	$indelPenalty, $mismatchPenalty);
$cmd = "$0 @ARGV";

my $keepHead;# = "CLC"; #### if there are two different heads (ex. "contig", "Contig_"), contigs with this head will be kept in a new sequence file.


#$verboseFlag = 1;
$endClippingAllowed = 4;
$lenMinOverlap = 40;
$identityMin = 90;
$indelPenalty = 1;   #### penalty for INDEL  ... considering hybrid assembly between illumina and 454
$mismatchPenalty = 2;   #### penalty for mismatch..

GetOptions(
	"h|?|help"		=> \$helpFlag,
	"verbose"		=> \$verboseFlag,
	"Length=i"	=> \$lenMinOverlap,
	"Identity=i"	=> \$identityMin,
	"blast=s"	=> \$blnFn,
	"seq=s"	=> \$seqFn,
	"Scaffold=s" => \$scaffoldFn,
	"keep=s"	=>	\$keepHead,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;

if(!defined $blnFn){ print STDERR "\n Error ..... Need a blast table file!!\n\n"; help(1); }
if(!defined $seqFn){ print STDERR "\n Error ..... Need a merged sequence file!!\n\n"; help(1); }
#if(!defined $scaffoldFn){ print STDERR "\nNeed a '454Scaffolds.txt' file!!\n\n"; help(1); }
if(!defined $outFn){ print STDERR "\n Error ..... Need a output file!!\n\n"; help(1); }


if(defined $outFn && ($blnFn eq $outFn || $seqFn eq $outFn || (defined $scaffoldFn && $scaffoldFn eq $outFn)))
{ print STDERR " Error) Input and output files are same \n"; exit; }

sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  -b <blast table file>  -s <contig FASTA file>   -o <output FASTA file> \n\n";
	print STDERR "  --------- other options (can be skipped) ---------\n";
	print STDERR "        [ -keep <header of contigs to be kept (default : contig) (ex. velvet contigs : 'NODE_')> ]\n";
	print STDERR "        [ -S <'454Scaffolds.txt' file> ]\n";
	print STDERR "        [ -L <mininum overlapping length (default : $lenMinOverlap) ]\n";
	print STDERR "        [ -I <minimum identiity (default : $identityMin) ]\n";
	print STDERR "  --------------------------------------------------\n\n";
	print STDERR "  ex) $0 -b megedContigs.NoMT.NoInc.bln.table.gz -s megedContigs.NoMT.NoInc.fa -S Good.Step6Len12/assembly/454Scaffolds.txt -o megedContigs.NoMT.NoInc.new.fa\n\n";
	exit($return);
}

my @strandTxt = ("+", "-");
my ($F, $R, $left, $right) = (0, 1, 5, 3);
my ($in, $out, @arr, $i);


$keepHead = "contig" if(defined $scaffoldFn || !defined $keepHead);


my ($seq, $name, %contigLen);
$in = openInput($seqFn);

$seq = "";
while(<$in>){	
	s/[\r\n]+//g;
	if(/^>/){
		if(defined $name){ $contigLen{$name} = length($seq);		}
		@arr = split /[\t\s]+/, $';
		$name = $arr[0];
		$seq = "";
	}
	else{
		s/[\s\t]+//g;
		$seq .= $_;
	}
}
if(defined $name){ $contigLen{$name} = length($seq); }
close($in);


my (@scaffold, %scaffoldContigs, $sid, $sNum, $j, $str);
if(defined $scaffoldFn){
	$in = openInput($scaffoldFn);

	while(<$in>){	
		s/[\r\n]+//g;
		@arr = split /\t/;

		$arr[0] =~ /(\d+)/;
		if(!defined $sid || $sid != $1){
			#### if scaffold has only one contig, delete it..
			if(defined $sid && scalar @{$scaffold[$sid]{contig}} == 1){
				delete $scaffoldContigs{$scaffold[$sid]{contig}[0]{name}};
				$scaffold[$sid] = undef;
			}
			$j = 1; #### order of a contig in a scaffold
		}
		$sid = $1;
		if($arr[5] =~ /^contig/){			
			
			$scaffoldContigs{$arr[5]}{id} = $sid;
			$scaffoldContigs{$arr[5]}{order} = $j;
			$j++;

			$i = $#{$scaffold[$sid]{contig}} + 1;
			$scaffold[$sid]{name} = $arr[0];
			$scaffold[$sid]{contig}[$i]{name} = $arr[5];
			$scaffold[$sid]{contig}[$i]{start} = $arr[1];
			$scaffold[$sid]{contig}[$i]{end} = $arr[2];
			$scaffold[$sid]{contig}[$i]{strand} = $arr[8];
		}
		else{
			$scaffold[$sid]{gap}[$#{$scaffold[$sid]{gap}} + 1]{start} = $arr[1];
			$scaffold[$sid]{gap}[$#{$scaffold[$sid]{gap}}]{end} = $arr[2];
			$scaffold[$sid]{gap}[$#{$scaffold[$sid]{gap}}]{len} = $arr[2] - $arr[1] + 1;
		}
	}
	$sNum = $sid+1;
	close($in);
}

my($qName, $rName, $identity, $matchLen, $mismatch, $indel, $qStart, $qEnd, $rStart, $rEnd, $qLen, $rLen, $rStrand,
	%endMatches, %includeLowMark, %includeLow, %includeHigh, %deleteHash, $numInclude, $numEndMatch);

($numInclude, $numEndMatch) = (0, 0);
my (%leftRepeatHit, %rightRepeatHit, $leftHitMaxLen, $rightHitMaxLen); ### to see an end is repeat...

$in = openInput($blnFn);

while(<$in>){
	s/[\r\n]+//g;
	next if(/^#/ || /^\s*$/);
	
	@arr = split /\t/;
	if($#arr < 10){
		print STDERR "The blast file does not have enough information.\n";
		close($in);
		exit;
	}
	next if($arr[0] eq $arr[1]);
	$name = $arr[0];  #### keep original name..
	
	if(defined $qName && $qName ne $arr[0]){
		
		if(defined $endMatches{$qName}){
			if(defined $leftRepeatHit{num}){ 
				print "$qName 5' end is repeat..\n"  if(defined $verboseFlag);
				foreach $rName (keys %{$endMatches{$qName}{5}}) {
					if($endMatches{$qName}{5}{$rName}{rStrand} == $F){
						if(defined $endMatches{$rName} && defined $endMatches{$rName}{3}{$qName}){
							delete $endMatches{$rName}{3}{$qName};
							if(keys %{$endMatches{$rName}{3}} == 0){
								delete $endMatches{$rName}{3};
								delete $endMatches{$rName} if(!defined $endMatches{$rName}{5});
							}
						}
					}
					else{ ## if($endMatches{$qName}{5}{$rName}{rStrand} == $R)
						if(defined $endMatches{$rName} && defined $endMatches{$rName}{5}{$qName}){
							delete $endMatches{$rName}{5}{$qName};
							if(keys %{$endMatches{$rName}{5}} == 0){
								delete $endMatches{$rName}{5};
								delete $endMatches{$rName} if(!defined $endMatches{$rName}{3});
							}
						}
					}					
				}
				delete $endMatches{$qName}{5}; 
			}			
			if(defined $rightRepeatHit{num}){
				print "$qName 3' end is repeat..\n"  if(defined $verboseFlag);
				foreach $rName (keys %{$endMatches{$qName}{3}}) {
					if($endMatches{$qName}{3}{$rName}{rStrand} == $F){
						if(defined $endMatches{$rName} && defined $endMatches{$rName}{5}{$qName}){
							delete $endMatches{$rName}{5}{$qName};
							if(keys %{$endMatches{$rName}{5}} == 0){
								delete $endMatches{$rName}{5};
								delete $endMatches{$rName} if(!defined $endMatches{$rName}{3});
							}
						}
					}
					else{ ## if($endMatches{$qName}{3}{$rName}{rStrand} == $R)
						if(defined $endMatches{$rName} && defined $endMatches{$rName}{3}{$qName}){
							delete $endMatches{$rName}{3}{$qName};
							if(keys %{$endMatches{$rName}{3}} == 0){
								delete $endMatches{$rName}{3};
								delete $endMatches{$rName} if(!defined $endMatches{$rName}{5});
							}
						}						
					}					
				}

				delete $endMatches{$qName}{3}; 
			}

			if(!defined $endMatches{$qName}{5} && !defined $endMatches{$qName}{3}){
				delete $endMatches{$qName};
				print "$qName does not have end pairs anymore.\n"  if(defined $verboseFlag);
			}

			#if($qName eq "Contig_833"){
			#	print "Contig_833, $endMatches{$qName}, $leftRepeatHit{num}, $rightRepeatHit{num}\n"; exit;
			#}
		}
		(%leftRepeatHit, %rightRepeatHit, $leftHitMaxLen, $rightHitMaxLen) = ((), (), 0, 0); ## , $leftHitMaxLen, $rightHitMaxLen not using
	}

	($qName, $rName, $identity, $matchLen, $mismatch, $indel, $qStart, $qEnd, $rStart, $rEnd, $qLen, $rLen) 
		= ($arr[0], $arr[1], $arr[2], $arr[3], $arr[4], $arr[5], $arr[6], $arr[7], $arr[8], $arr[9], $contigLen{$arr[0]}, $contigLen{$arr[1]});
		
	$rStrand = $F;
	my ($rOriStart, $rOriEnd) = ($rStart, $rEnd);
	if($rStart > $rEnd){
		($rStart, $rEnd) = ($rLen - $rStart + 1, $rLen - $rEnd + 1);
		$rStrand = $R;
	}

	$identity = ($matchLen - ($mismatch*$mismatchPenalty) - ($indel*$indelPenalty)) / $matchLen * 100;
	
	#next if($identity < $identityMin);
	
	#if($qName eq "Contig_185663" && $rName eq "Contig_92550"){
	#	print " $qName, $rName --) , $_\n ".
	#		"if($rStart <= $endClippingAllowed+1 && $qEnd >= $qLen -$endClippingAllowed){\nelsif($qStart <= $endClippingAllowed+1 && $rEnd >= $rLen -$endClippingAllowed){		\n";
	#}
	
	if($identity >= $identityMin && $rStart <= $endClippingAllowed+1 && $rEnd >= $rLen-$endClippingAllowed && 
		(!defined $includeLow{$rName} || !defined $includeLow{$rName}{$qName}) && (!defined $includeLow{$qName} || !defined $includeLow{$qName}{$rName}) ){		
		($includeLow{$rName}{$qName}{lStart}, $includeLow{$rName}{$qName}{lEnd}) = ($rStrand == $F ? ($rOriStart, $rOriEnd) : ($rOriEnd, $rOriStart));
		($includeLow{$rName}{$qName}{hStart}, $includeLow{$rName}{$qName}{hEnd}) = ($rStrand == $F ? ($qStart, $qEnd) : ($qEnd, $qStart));
		$includeLow{$rName}{$qName}{strand} = $rStrand;
		
		#if($qName =~ /^NODE_163_/ && $rName =~ /^NODE_271_/){ print "2. $qName, $rName, \n"; } ##TEST

		if($rName =~ /^$keepHead/){
			if(!defined $includeLowMark{$rName}){
				$includeLowMark{$rName} = 1; ### the 'rName' contig is for keeping 
			}
			if($qName =~ /^$keepHead/){
				$includeLowMark{$rName} = 2; ### high contig including 'rName' contig is also for keeping 
			}
		}
		#print "\$includeLow{$rName}{$qName}, incMatch->{strand} : $includeLow{$rName}{$qName}{strand}\n";
		
		$includeHigh{$qName}{$rName} = 1; #####
		$numInclude++;
	}
	elsif($identity >= $identityMin && $qStart <= $endClippingAllowed+1 && $qEnd >= $qLen -$endClippingAllowed && 
		(!defined $includeLow{$rName} || !defined $includeLow{$rName}{$qName}) && (!defined $includeLow{$qName} || !defined $includeLow{$qName}{$rName}) ){
		$includeLow{$qName}{$rName}{lStart} = $qStart;
		$includeLow{$qName}{$rName}{lEnd} = $qEnd;
		$includeLow{$qName}{$rName}{hStart} = $rOriStart;
		$includeLow{$qName}{$rName}{hEnd} = $rOriEnd;
		$includeLow{$qName}{$rName}{strand} = $rStrand;
		
		#if($qName =~ /^NODE_163_/ && $rName =~ /^NODE_271_/){ print "3. $qName, $rName, \n"; } ##TEST

		if($qName =~ /^$keepHead/){
			if(!defined $includeLowMark{$qName}){
				$includeLowMark{$qName} = 1; ### the 'rName' contig is for keeping 
			}
			if($rName =~ /^$keepHead/){
				$includeLowMark{$qName} = 2; ### high contig including 'rName' contig is also for keeping 
			}
		}
		

		#print "\$includeLow{$qName}{$rName}, incMatch->{strand} : $includeLow{$qName}{$rName}{strand}\n";
		
		$includeHigh{$rName}{$qName} = 1; #####
		$numInclude++;
	}
	elsif(!defined $includeLow{$qName} && !defined $includeLow{$rName}){
		if($identity >= $identityMin && $matchLen >= $lenMinOverlap){
			
			#( ($rStrand == $F && $rEnd >= $rLen -$endClippingAllowed) || ($rStrand == $R && $rEnd <= $endClippingAllowed+1) )
			if($qStart <= $endClippingAllowed+1 && $rEnd >= $rLen -$endClippingAllowed){		
				### it is for 5' end match of $qName contig.
				#if($qName =~ /^NODE_163_/ && $rName =~ /^NODE_271_/){ print "4. $qName, $rName, \n$_\n"; } ##TEST
			
				# if the repeat length is longer than match-20, skip...
				next if(defined $leftRepeatHit{num} && $leftRepeatHit{len} > $qEnd-$qStart-20);

				### if the both ends match with same contig.. treat it as a repeat..
				if(defined $endMatches{$qName} && defined $endMatches{$qName}{3}{$rName}){  
					$leftRepeatHit{num}++;
					my $len = $qEnd-$qStart+1;
					$leftRepeatHit{len} = $len if(!defined $leftRepeatHit{len} || $leftRepeatHit{len} < $len);
										
					$rightRepeatHit{num}++;
					$len = $endMatches{$qName}{3}{$rName}{qEnd}-$endMatches{$qName}{3}{$rName}{qStart}+1;
					$endMatches{$qName}{3}{$rName}{qLen} = $len;
					$endMatches{$qName}{3}{$rName}{identity} = $identity;
					$rightRepeatHit{len} = $len if(!defined $rightRepeatHit{len} || $rightRepeatHit{len} < $len);
					#if($qName =~ /^NODE_179/) { print "$qName right repeat 1, rName : $rName, len : $len\n" }; ##TEST

					next;
				}
				 ### if the both ends match with same contig.. treat it as a repeat..
				if($rStrand == $F && defined $endMatches{$rName}{5}{$qName}) { ### if F, it should be only 3'end of contig{rName}
					print "strand is $strandTxt[$rStrand]. qName : $qName.. $qName and $rName pair should be $rName(3)-$qName(5). But $qName(3|5)-$rName(5) exists... Removing pair\n" if(defined $verboseFlag);	
					delete $endMatches{$rName}{5}{$qName};
					next;
				}
				if($rStrand == $R && defined $endMatches{$rName}{3}{$qName}) { ### if R, it should be only 5'end of contig{rName}
					print "strand is $strandTxt[$rStrand]. qName : $qName.. $qName and $rName pair should be $qName(5)-$rName(5). But $rName(3)-$qName(5|3) exists... Removing pair\n" if(defined $verboseFlag);	
					delete $endMatches{$rName}{3}{$qName};
					next;
				}
				delete $leftRepeatHit{num};
				

				$leftRepeatHit{noRepeatEnd} = $qEnd if(!defined $leftRepeatHit{noRepeatEnd} || $leftRepeatHit{noRepeatEnd} < $qEnd);;
				#$leftHitMaxLen, $rightHitMaxLen

				putAtMatchHash(\%endMatches, $qName, $rName, $left, $rStrand, $qStart, $qEnd, $rOriStart, $rOriEnd, $identity);

			}
			elsif($rStart <= $endClippingAllowed+1 && $qEnd >= $qLen -$endClippingAllowed){
				### it is for 3' end match of $qName contig.
				#if($qName =~ /^NODE_163_/ && $rName =~ /^NODE_271_/){ print "5. $qName, $rName, \n"; } ##TEST

				#if($qName eq "Contig_185663" && $rName eq "Contig_92550"){ print "entered.. at right hit\n";}
				next if(defined $rightRepeatHit{num} && $rightRepeatHit{len} > $qEnd-$qStart-20);
				
				#if($qName =~ /^NODE_179/) { print "$qName hits rName : $rName, len : " .($qEnd-$qStart-20). "\n" }; ##TEST

				if(defined $endMatches{$qName} && defined $endMatches{$qName}{5}{$rName}){
					$leftRepeatHit{num}++;
					my $len = $endMatches{$qName}{5}{$rName}{qEnd}-$endMatches{$qName}{5}{$rName}{qStart}+1;
					$endMatches{$qName}{5}{$rName}{qLen} = $len;
					$endMatches{$qName}{5}{$rName}{identity} = $identity;
					$leftRepeatHit{len} = $len if(!defined $leftRepeatHit{len} || $leftRepeatHit{len} < $len);

					
					$rightRepeatHit{num}++;
					$len = $qEnd-$qStart+1;
					$rightRepeatHit{len} = $len if(!defined $rightRepeatHit{len} || $rightRepeatHit{len} < $len);
					#if($qName =~ /^NODE_179/) { print "$qName right repeat 2, rName : $rName, len : $len\n" }; ##TEST

					next;
				}
				### if the both ends match with same contig.. treat it as a repeat..
				if($rStrand == $F && defined $endMatches{$rName}{3}{$qName}) { ### it should be only 5'end of contig{rName}
					print "strand is $strandTxt[$rStrand]. qName : $qName.. $qName and $rName pair should be $qName(3)-$rName(5). But $rName(3)-$qName(5|3) exists... Removing pair\n" if(defined $verboseFlag);	
					delete $endMatches{$rName}{3}{$qName};
					next;
				}
				if($rStrand == $R && defined $endMatches{$rName}{5}{$qName}) { ### it should be only 3'end of contig{rName}
					print "strand is $strandTxt[$rStrand]. qName : $qName.. $qName and $rName pair should be $qName(3)-$rName(3). But $qName(3)-$rName(5|3) exists... Removing pair\n" if(defined $verboseFlag);	
					delete $endMatches{$rName}{5}{$qName};
					next;
				}
				delete $rightRepeatHit{num};

				$rightRepeatHit{noRepeatEnd} = $qStart if(!defined $rightRepeatHit{noRepeatEnd} || $rightRepeatHit{noRepeatEnd} > $qStart);
				putAtMatchHash(\%endMatches, $qName, $rName, $right, $rStrand, $qStart, $qEnd, $rOriStart, $rOriEnd, $identity);

			}
			elsif($qStart <= $endClippingAllowed+1 && #### to check the left end is repeat...
				(!defined $leftRepeatHit{noRepeatEnd} || $leftRepeatHit{noRepeatEnd}-20 < $qEnd))
			{ 
				#if($qName =~ /^NODE_163_/ && $rName =~ /^NODE_271_/){ print "6. $qName, $rName, \n"; } ##TEST
				$leftRepeatHit{num}++;
				my $len = $qEnd-$qStart+1;
				$leftRepeatHit{len} = $len if(!defined $leftRepeatHit{len} || $leftRepeatHit{len} < $len);
			}
			elsif($qEnd >= $qLen -$endClippingAllowed && #### to check the right end is repeat...
				(!defined $rightRepeatHit{noRepeatEnd} || $rightRepeatHit{noRepeatEnd}+20 > $qStart))
			{
				#if($qName =~ /^NODE_163_/ && $rName =~ /^NODE_271_/){ print "7. $qName, $rName, \n"; } ##TEST
				my $len = $qEnd-$qStart+1;				
				$rightRepeatHit{num}++;
				$rightRepeatHit{len} = $len if(!defined $rightRepeatHit{len} || $rightRepeatHit{len} < $len);
				#if($qName =~ /^NODE_179/) { print "$qName right repeat 3, rName : $rName, len : $len, \$rightRepeatHit{noRepeatEnd} : $rightRepeatHit{noRepeatEnd}, $qStart\n" }; ##TEST

			}
		}
		elsif(!defined $includeLow{$qName}){  ##### for weak pairs..
		}
	}
	
	#exit if($qName eq "contig137159");
	#if(defined $endMatches{contig137157}{3}{contig20550}){
	#	print STDERR "$qName, $rName -- defined $endMatches{contig137157}{3}{contig20550}\n";
	#}
	#else{
	#	print STDERR "$qName, $rName -- undefined \$endMatches{contig137157}{3}{contig20550}\n";
	#}

}

#### for the last contig...
if(defined $qName && defined $endMatches{$qName}){
	if(defined $leftRepeatHit{num}){ 
		print "$qName 5' end is repeat..\n"  if(defined $verboseFlag);
		delete $endMatches{$qName}{5}; 
	}			
	if(defined $rightRepeatHit{num}){
		print "$qName 3' end is repeat..\n"  if(defined $verboseFlag);
		delete $endMatches{$qName}{3}; 
	}
	delete $endMatches{$qName} if(!defined $endMatches{$qName}{5} && !defined $endMatches{$qName}{3}); 
}

close($in);


sub putAtMatchHash
{
	my ($matchHash, $qName, $rName, $direction, $strand, $qStart, $qEnd, $rStart, $rEnd, $identity) = @_;
	$matchHash->{$qName}{qName} = $qName;
	$matchHash->{$qName}{$direction}{$rName}{rName} = $rName; 
	$matchHash->{$qName}{$direction}{$rName}{qStart} = $qStart; # = $_;	
	$matchHash->{$qName}{$direction}{$rName}{qEnd} = $qEnd;
	$matchHash->{$qName}{$direction}{$rName}{qLen} = abs($qEnd-$qStart)+1;
	$matchHash->{$qName}{$direction}{$rName}{identity} = $identity;
	$matchHash->{$qName}{$direction}{$rName}{rStart} = $rStart;
	$matchHash->{$qName}{$direction}{$rName}{rEnd} = $rEnd;
	$matchHash->{$qName}{$direction}{$rName}{rStrand} = $strand;
	
	$matchHash->{$rName}{qName} = $rName;
	if($rStrand == $F){
		$direction = ($direction == 5 ? 3 : 5);
		$matchHash->{$rName}{$direction}{$qName}{rName} = $qName; 
		$matchHash->{$rName}{$direction}{$qName}{qStart} = $rStart; # = $_;	
		$matchHash->{$rName}{$direction}{$qName}{qEnd} = $rEnd;		
		$matchHash->{$rName}{$direction}{$qName}{qLen} = abs($rEnd-$rStart)+1;
		$matchHash->{$rName}{$direction}{$qName}{identity} = $identity;
		$matchHash->{$rName}{$direction}{$qName}{rStart} = $qStart;
		$matchHash->{$rName}{$direction}{$qName}{rEnd} = $qEnd;
		$matchHash->{$rName}{$direction}{$qName}{rStrand} = $strand;
	}
	else{
		$matchHash->{$rName}{$direction}{$qName}{rName} = $qName; 
		$matchHash->{$rName}{$direction}{$qName}{qStart} = $rEnd; # = $_;	
		$matchHash->{$rName}{$direction}{$qName}{qEnd} = $rStart;		
		$matchHash->{$rName}{$direction}{$qName}{qLen} = abs($rEnd-$rStart)+1;		
		$matchHash->{$rName}{$direction}{$qName}{identity} = $identity;
		$matchHash->{$rName}{$direction}{$qName}{rStart} = $qEnd;
		$matchHash->{$rName}{$direction}{$qName}{rEnd} = $qStart;
		$matchHash->{$rName}{$direction}{$qName}{rStrand} = $strand;
	}
}



#foreach my $aa (keys %{$endMatches{Contig_159061}{5}}) {
#	print "5 : $aa\n";
#}
#foreach my $aa (keys %{$endMatches{Contig_159061}{3}}) {
#	print "3 : $aa\n";
#}
#exit;

my %scaffoldHighContigs; #### to store a list of contigs including(covering) another contigs in a scaffold.

foreach $rName (keys %includeLow) {
	if(defined $scaffoldContigs{$rName}){
		
		#### if(high level contig(qName, contig inclucing rName) is also member of scaffold... break their matching...
		foreach $qName (keys %{$includeLow{$rName}}) {	
			if(defined $scaffoldContigs{$qName}){
				delete $includeLow{$rName}{$qName};
				delete $includeHigh{$qName}{$rName};
			}
			else{
				print "$rName(length:$contigLen{$rName}) is included in $qName(length:$contigLen{$qName}), but is a member of a scaffold" .
					"($includeLow{$rName}{$qName}{lStart}-$includeLow{$rName}{$qName}{lEnd} vs " . 
					"$includeLow{$rName}{$qName}{hStart}-$includeLow{$rName}{$qName}{hEnd} : $strandTxt[$includeLow{$rName}{$qName}{strand}]). " .
					"It will not be removed..\n"  if(defined $verboseFlag);
				
				
				$scaffoldHighContigs{$qName}[$#{$scaffoldHighContigs{$qName}}+1] = $rName;
			}
		}
		my $size = keys %{$includeLow{$rName}}; 
		delete $includeLow{$rName} if($size == 0); #### if both rName and qName contigs are member of scaffold... break pair...
	}
	else{
		foreach $qName (keys %{$includeLow{$rName}}) {				
			$deleteHash{$rName} = 1; #[$#{$deleteHash{$rName}}+1] = $includeLow{$qName}{$rName};
			print "$rName is included in $qName .. deleted\n"  if(defined $verboseFlag);
			#delete $endMatches{$rName} if($endMatches{$rName}); ###### delete node here....
			$i++;
		}
	}
}


my ($recodeType, %contigSeq);
$in = openInput($seqFn);
$out = openOutput($outFn);
$recodeType = 0; ### 0:not recorde, 1:endmatch, 2:file out.
while(<$in>){	
	s/[\r\n]+//g;
	if(/^>/){
		@arr = split /[\t\s]+/, $';
		$name = $arr[0];
				
		if(defined $deleteHash{$name}){
			$recodeType = 0;
		}
		elsif(defined $endMatches{$name} || defined $scaffoldContigs{$name} || defined $scaffoldHighContigs{$name} || defined $includeHigh{$name}){
			$recodeType = 1;
			#print $out "$_\n";
		}
		elsif($name =~ /^$keepHead/){
			$recodeType = 2;
			print $out "$_\n";
		}
		else{ ### ignore..
			$recodeType = 0; 
		}
	}
	else{
		if($recodeType == 1){
			$contigSeq{$name} .= $_;
			#print $out "$_\n";
		}
		elsif($recodeType == 2){
			print $out "$_\n";
		}
	}
}
close($in);




############################### 
#### build graphs for paired contigs...
my @names = sort {length($a) <=> length($b) || $a cmp $b} keys %endMatches;

$str = $outFn;
$str =~ s/\.(fa|fna|seq)$/.merged_graphs.txt/;
$str .= ".merged_graphs.txt" if($str eq $outFn);
my $graphOut = openOutput($str);


$numEndMatch = int($numEndMatch / 2);
#print "# end : $#names, # include : $numInclude\n"; 

my (@mergedGroup, $gi, $FNum, $RNum, $leftNum, $rightNum, $nodeNum, $node, $match, $start,
 %groupContigs);
$gi = 1;

foreach $qName (@names) {
	next if(!defined $endMatches{$qName} || defined $includeLow{$qName} || defined $deleteHash{$qName}); ###### if it is deleted....
	my (@nodeList, @strandList, @leftNodeList, @leftStrandList, @rightNodeList, @rightStrandList);

	($FNum, $RNum) = (1, 0); ### ($FNum = 1) includes the current contig...
	print "\n------- Testing pairs from $qName\n" if(defined $verboseFlag);
	#searchLeftNode($endMatches{$qName}, $F);
	print "(left)\n" if(defined $verboseFlag);

	searchNode($endMatches{$qName}, $F, $left, \@leftNodeList, \@leftStrandList);		

	print "(right)\n" if(defined $verboseFlag);
	searchNode($endMatches{$qName}, $F, $right, \@rightNodeList, \@rightStrandList);

	$leftNum = scalar @leftNodeList;
	$rightNum = scalar @rightNodeList;
	next if($leftNum == 0 && $rightNum == 0);
		
	@nodeList = reverse @leftNodeList;	
	@strandList = reverse @leftStrandList;

	push(@nodeList, $qName);
	push(@strandList, $F);
	push(@nodeList, @rightNodeList);
	push(@strandList, @rightStrandList);
	$nodeNum = $leftNum + $rightNum + 1;	### +1 : current contig..

	if($RNum > $FNum || ($RNum == $FNum && $strandList[0] == $R && $strandList[$#strandList] == $R) ){
		@nodeList = reverse @nodeList;	
		@strandList = reverse @strandList;
		for($i = 0; $i < $nodeNum; $i++){
			$strandList[$i] = ($strandList[$i] == $F ? $R : $F);
		}
	}

	
	my ($minLen, $len, $weak_i, $dir, %prevNodes);
	$minLen = 10000000;
	for($i = 1; $i < $nodeNum; $i++){ 
		if($strandList[$i] == $F){ $dir = 5; }
		else { $dir = 3; }
		$len = $endMatches{$nodeList[$i]}{$dir}{$nodeList[$i-1]}{qLen}*$endMatches{$nodeList[$i]}{$dir}{$nodeList[$i-1]}{identity};
		if($minLen > $len){
			$minLen = $len;
			$weak_i = $i;
		}
		
		if($nodeList[$i] eq $nodeList[0] && $strandList[$i] == $strandList[0]){
			if($strandList[$i] == $F){ $dir = 5; }
			else { $dir = 3; }

			if($#nodeList < $i*2){
				for($j = 1; $j <= $i; $j++){
					$nodeList[$i+$j] = $nodeList[$j];
					$strandList[$i+$j] = $strandList[$j];
				}
			}

			if(defined $verboseFlag){
				print " ** A recursive link was detected.\n ** " ;
				for($j = 0; $j < $nodeNum; $j++){
					print " - $nodeList[$j]($strandTxt[$strandList[$j]], $j)";
				}
				print "\n";
				print " ** Cutting $nodeList[$weak_i-1]-$nodeList[$weak_i] link. \n" ;
			}

			$weak_i = 0 if($weak_i == $i);
			my $end = ($i+$weak_i-1);
			@nodeList = @nodeList[$weak_i..$end];
			@strandList = @strandList[$weak_i..$end];
			$nodeNum = $i;

			last;
		}

		if(defined $prevNodes{$nodeList[$i]}){
			if(defined $verboseFlag){
				print " ** A recursive link was detected.\n ** " ;
				for($j = 0; $j < $nodeNum; $j++){
					print " - $nodeList[$j]($strandTxt[$strandList[$j]], $j)";
				}
				print "\n";
				print " ** Cutting $nodeList[$weak_i-1]-$nodeList[$weak_i] link. \n" ;
			}

			$j = $prevNodes{$nodeList[$i]}{i} + 1;
			if($prevNodes{$nodeList[$j]}{len_withPrev} > $len){
				print " ** Cutting $nodeList[$i-1]-$nodeList[$i] link. \n" if(defined $verboseFlag);
				@nodeList = @nodeList[0..$i-1];
				@strandList = @strandList[0..$i-1];
				$nodeNum = $i;
			}
			else{
				print " ** Cutting $nodeList[$i-1]-$nodeList[$i] link. \n" if(defined $verboseFlag);
				@nodeList = @nodeList[$j..$#nodeList];
				@strandList = @strandList[$j..$#nodeList];
				$nodeNum = $#nodeList-$j+1;
			}
			
			last;
			
		}
		$prevNodes{$nodeList[$i]}{i} = $i;
		$prevNodes{$nodeList[$i]}{len_withPrev} = $len;

	}

	$mergedGroup[$gi]{name} = sprintf("merged%06d", $gi);
	print "Group $mergedGroup[$gi]{name} [graph] : " if(defined $verboseFlag);	

	
	if(defined $verboseFlag){
		for($i = 0; $i < $nodeNum; $i++){
			print " - $nodeList[$i]($strandTxt[$strandList[$i]])";
		}
		print "\n";
	}
	my ($rvSeq, @seqPos, @relativePos, @mergedName, $rStart, $rEnd);

	$name = $nodeList[0];
	$groupContigs{$nodeList[0]}{gi} = $gi;
	$groupContigs{$nodeList[0]}{order} = 0;
	
	$rStart = 1; #### to show position in the merged sequence.
	$start = ($strandList[0] == $F ? 1 : $contigLen{$nodeList[0]});
	$seqPos[0]{start} = $start;
	for($i = 0; $i < $nodeNum - 1; $i++){
		$name .= ".$nodeList[$i+1]";
		$groupContigs{$nodeList[$i+1]}{gi} = $gi;
		$groupContigs{$nodeList[$i+1]}{order} = $i+1;


		$node = $endMatches{$nodeList[$i]};

		#$match = ($strandList[$i] == $F ? $node->{3}{$nodeList[$i+1]} : $node->{5}{$nodeList[$i+1]});
		
		if($strandList[$i] == $F) { 
			$match = $node->{3}{$nodeList[$i+1]};
			
			if(!defined $match){
				print STDERR "ERROR at line ".__LINE__.": No \$match object. $nodeList[$i]($strandTxt[$strandList[$i]]:$contigLen{$nodeList[$i]}) vs $nodeList[$i+1]($strandTxt[$strandList[$i+1]]:$contigLen{$nodeList[$i+1]})\n";
				foreach my $aa (values %{$endMatches{$nodeList[$i]}{3}}) {
					print STDERR "($aa->{rName}), $aa->{qEnd}, $aa->{rEnd}\n";
				}
				exit;
			}

			print " match : $nodeList[$i]($strandTxt[$strandList[$i]]:$contigLen{$nodeList[$i]}) vs $nodeList[$i+1]($strandTxt[$strandList[$i+1]]:$contigLen{$nodeList[$i+1]}), " .
				"$match->{qStart}..$match->{qEnd}, $match->{rStart}..$match->{rEnd} \n" if(defined $verboseFlag); 
			
			if(!defined $nodeList[$i] || !defined $start || !defined $match->{qStart}){
				print STDERR "ERROR at line ".__LINE__.": \n";
				exit;
			}

			if($start <= $match->{qStart}){
				$seqPos[$i]{end} = $match->{qStart};
				$rEnd = $rStart+abs($start-$match->{qStart});
				print "+ $nodeList[$i]($strandTxt[$strandList[$i]]), $start..$match->{qStart} ($rStart-$rEnd)\n" if(defined $verboseFlag);	
				$rStart = $rEnd+1;

				if( $match->{qEnd}-$match->{qStart} <= abs($match->{rEnd}-$match->{rStart}) ){ #### choose shorter contig.
					$seqPos[$i]{end} = $match->{qEnd};
					$rEnd = $rStart+abs($match->{qStart}+1-$match->{qEnd});
					print "+ $nodeList[$i]($strandTxt[$strandList[$i]]), " . 
						($match->{qStart}+1). "..$match->{qEnd} ($rStart-$rEnd)\n" if(defined $verboseFlag); 
					$rStart = $rEnd+1;
				}
				else{
					if($strandList[$i+1] == $F){
						$seqPos[$i+1]{start} = $match->{rStart}+1;
						$seqPos[$i+1]{end} = $match->{rEnd};
					}
					else{
						$seqPos[$i+1]{start} = $match->{rStart}-1;
						$seqPos[$i+1]{end} = $match->{rEnd};
					}
					$rEnd = $rStart+abs($seqPos[$i+1]{start}-$seqPos[$i+1]{end});
					print "+ $nodeList[$i+1]($strandTxt[$strandList[$i+1]]), $seqPos[$i+1]{start} ..$seqPos[$i+1]{end} ($rStart-$rEnd)\n" if(defined $verboseFlag); 
					$rStart = $rEnd+1;
				}
			}
			else{  ### if($start <= $match->{qStart}){
				$seqPos[$i]{end} = $match->{qEnd};
				$rEnd = $rStart+abs($start-$match->{qEnd});
				print "+ $nodeList[$i]($strandTxt[$strandList[$i]]), $start..$match->{qEnd} ($rStart-$rEnd)\n" if(defined $verboseFlag); 
				$rStart = $rEnd+1;
			}
			$start = ($strandList[$i+1] == $F ? $match->{rEnd}+1 : $match->{rEnd}-1);
		}
		else { #### for $R
			$match = $node->{5}{$nodeList[$i+1]};
			if(!defined $match){
				print STDERR "ERROR at line ".__LINE__.": No \$match object. $nodeList[$i]($strandTxt[$strandList[$i]]:$contigLen{$nodeList[$i]}) vs $nodeList[$i+1]($strandTxt[$strandList[$i+1]]:$contigLen{$nodeList[$i+1]})\n";
				foreach my $aa (values %{$endMatches{$nodeList[$i]}{5}}) {
					print STDERR "($aa->{rName}), $aa->{qEnd}, $aa->{rEnd}\n";
				}
				exit;
			}

			print " match : $nodeList[$i]($strandTxt[$strandList[$i]]:$contigLen{$nodeList[$i]}) vs $nodeList[$i+1]($strandTxt[$strandList[$i+1]]:$contigLen{$nodeList[$i+1]}), " .
				"$match->{qStart}..$match->{qEnd}, $match->{rStart}..$match->{rEnd} \n" if(defined $verboseFlag); 
	
			if(!defined $match->{qEnd}){
				print STDERR "ERROR at line ".__LINE__.": \n";
				exit; 
			}
			if($start >= $match->{qEnd}){				
				$seqPos[$i]{end} = $match->{qEnd};
				$rEnd = $rStart+abs($start-$match->{qEnd});
				print "+ $nodeList[$i]($strandTxt[$strandList[$i]]), $start..$match->{qEnd} ($rStart-$rEnd)\n" if(defined $verboseFlag);				
				$rStart = $rEnd+1;

				if( $match->{qEnd}-$match->{qStart} <= abs($match->{rEnd}-$match->{rStart}) ){ #### choose shorter contig.
					$seqPos[$i]{end} = $match->{qStart};
					$rEnd = $rStart+abs($match->{qEnd}-1-$match->{qStart});
					print "+ $nodeList[$i]($strandTxt[$strandList[$i]]), ".($match->{qEnd}-1)."..$match->{qStart} ($rStart-$rEnd)\n" if(defined $verboseFlag); 
					$rStart = $rEnd+1;
				}
				else{
					if($strandList[$i+1] == $F){
						$seqPos[$i+1]{start} = $match->{rEnd}+1;
						$seqPos[$i+1]{end} = $match->{rStart};
					}
					else{
						$seqPos[$i+1]{start} = $match->{rEnd}-1;
						$seqPos[$i+1]{end} = $match->{rStart};
					}
					$rEnd = $rStart+abs($seqPos[$i+1]{start}-$seqPos[$i+1]{end});
					print "+ $nodeList[$i+1]($strandTxt[$strandList[$i+1]]), $seqPos[$i+1]{start}..$seqPos[$i+1]{end} ($rStart-$rEnd)\n" if(defined $verboseFlag); 					
					$rStart = $rEnd+1;
				}
			}
			else{
				$seqPos[$i]{end} = $match->{qStart};
				$rEnd = $rStart+abs($start-$match->{qStart});
				print "+ $nodeList[$i]($strandTxt[$strandList[$i]]), $start..$match->{qStart} ($rStart-$rEnd)\n" if(defined $verboseFlag); 
				$rStart = $rEnd+1;
			}

			$start = ($strandList[$i+1] == $F ? $match->{rStart}+1 : $match->{rStart}-1);
		}		
		$seqPos[$i+1]{start} = $start if(!defined $seqPos[$i+1]{start});

		delete $endMatches{$nodeList[$i]};
		$deleteHash{$nodeList[$i]} = 1;
	}
	if($strandList[$i] == $F) { 
		$seqPos[$i]{end} = $contigLen{$nodeList[$i]};
		$rEnd = $rStart+abs($start-$contigLen{$nodeList[$i]});
		print "+ $nodeList[$i]($strandTxt[$strandList[$i]]), $start..$contigLen{$nodeList[$i]} ($rStart-$rEnd)\n" if(defined $verboseFlag); 
		$rStart = $rEnd+1;
	}
	else{ 
		$seqPos[$i]{end} = 1;
		$rEnd = $rStart+abs($start-1);
		print "+ $nodeList[$i]($strandTxt[$strandList[$i]]), $start..1 ($rStart-$rEnd)\n" if(defined $verboseFlag); 
		$rStart = $rEnd+1;
	}
	delete $endMatches{$nodeList[$i]};
	$deleteHash{$nodeList[$i]} = 1;
		
	my $prevLen = 0;

	for($i = 0; $i <= $#seqPos; $i++){
		$seqPos[$i]{rStart} = $prevLen + 1;
		$seqPos[$i]{rEnd} = $prevLen + 1 + abs($seqPos[$i]{end} - $seqPos[$i]{start});
		$prevLen += 1 + abs($seqPos[$i]{end} - $seqPos[$i]{start});		
		print $graphOut "$mergedGroup[$gi]{name}\t$seqPos[$i]{rStart}\t$seqPos[$i]{rEnd}\t$nodeList[$i]\t$seqPos[$i]{start}\t$seqPos[$i]{end}\t$strandTxt[$strandList[$i]]\n";
	}

	
	$mergedGroup[$gi]{nameAdd} = $name;	
	$mergedGroup[$gi]{nodeList} = \@nodeList;
	$mergedGroup[$gi]{strandList} = \@strandList;
	$mergedGroup[$gi]{seqPos} = \@seqPos;
	$mergedGroup[$gi]{num} = scalar @nodeList;
	$mergedGroup[$gi]{len} = $seqPos[$#seqPos]{rEnd};

	print "\nGroup $mergedGroup[$gi]{name} contig list\n" if(defined $verboseFlag);	
	for($i = 0; $i <= $#nodeList; $i++){
		print "$nodeList[$i]\n" if(defined $verboseFlag);	
	}
	print "\n\n" if(defined $verboseFlag);

	$gi++;	
}


#print "From here, single contig Group...\n" if(defined $verboseFlag);
#### to consider contigs covering other contigs in scaffolds.
#### If a contig includes another contig in a scaffold, they will be merged.
foreach $qName (keys %scaffoldHighContigs) {
	next if(defined $includeLow{$qName} || defined $deleteHash{$qName});
	#print "Group $gi\n" if(defined $verboseFlag);	
	#print "+ $qName($strandTxt[$F]), 1..$contigLen{$qName}\n" if(defined $verboseFlag); 
	$mergedGroup[$gi]{name} = $qName;
	$mergedGroup[$gi]{nodeList}[0] = $qName;
	$mergedGroup[$gi]{strandList}[0] = $F;
	$mergedGroup[$gi]{seqPos}[0]{start} = 1;
	$mergedGroup[$gi]{seqPos}[0]{end} = $contigLen{$qName};
	$mergedGroup[$gi]{seqPos}[0]{rStart} = 1;
	$mergedGroup[$gi]{seqPos}[0]{rEnd} = $contigLen{$qName};
	$mergedGroup[$gi]{num} = 1;
	$mergedGroup[$gi]{len} = $contigLen{$qName};

	#print $graphOut "$mergedGroup[$gi]{name}\t1\t$contigLen{$qName}\t$qName\t$seqPos[$i]{start}\t$seqPos[$i]{end}\t$strandTxt[$strandList[$i]]\n";

	$deleteHash{$qName} = 1;
	$gi++;
}
my $gNum = $gi;
#print "\n";

close($graphOut);

if(!defined $scaffoldFn){
	printMergedSeq();
	exit(0);
}

################################################
### build graphs for merged contig groups and scaffolds
my ($conflictStr, $conflictOut);

$str = $outFn;
$str =~ s/\.(fa|fna|seq)$/.scaffolds.conflict.txt/;
$str .= ".scaffolds.conflict.txt" if($str eq $outFn);
$conflictOut = openOutput($str);


my $ci;
for($sid = 1; $sid < $sNum; $sid++){
	next if(!defined $scaffold[$sid]);
	
	$conflictStr = "";
	my $cNum = $#{$scaffold[$sid]{contig}} + 1;
	print $conflictOut "--> $scaffold[$sid]{name}, num : $cNum\n";
	my @contigGroupInfo;
	for($ci = 0; $ci < $cNum; $ci++){
		$qName = $scaffold[$sid]{contig}[$ci]{name}; #### $ci th contig in scaffold_$sid..
		if(defined $includeLow{$qName}){			
			foreach $rName (keys %{$includeLow{$qName}} ) {
				if(defined $groupContigs{$rName}){  ### if a contig(rName) covering the qName contig and it is a part of a merged contig.. $groupContigs{$rName}
					$gi = $groupContigs{$rName}{gi};
					$contigGroupInfo[$ci]{$gi}{order} = $groupContigs{$rName}{order};
					$contigGroupInfo[$ci]{$gi}{rName} = $rName;
					print $conflictOut " - $qName in $scaffold[$sid]{name} is covered by $rName in $mergedGroup[$gi]{name}" 
						. "($includeLow{$qName}{$rName}{lStart}-$includeLow{$qName}{$rName}{lEnd} vs " 
						. "$includeLow{$qName}{$rName}{hStart}-$includeLow{$qName}{$rName}{hEnd}:$strandTxt[$includeLow{$qName}{$rName}{strand}])\n";
				}
			}			
		}
		elsif(defined $groupContigs{$qName}){
			$gi = $groupContigs{$qName}{gi};
			$contigGroupInfo[$ci]{$gi}{order} = $groupContigs{$qName}{order};
		}
	}

	
	### $contigGroupInfo[$ci] has merged contig ID(gi)s which the $ci-th contig is included in..
	### if(defined $groupContigs{$qName}), there is only one 'gi', but if(defined $includeLow{$qName}) there could be several 'gi'..
	my (%continuousRecord, %giAtPrev, @arrGI, @overGroup); 
	for($ci = 0; $ci < $cNum; $ci++){		
		addMergedGroupToScaffoldOverlapList($sid, $ci, $cNum, \@arrGI, \%giAtPrev, \%continuousRecord, \@contigGroupInfo, \@overGroup);

		%giAtPrev = ();
		@arrGI = keys %{$contigGroupInfo[$ci]};
		foreach $gi (@arrGI) {
			#print " adding $gi at contig $ci\n";
			$continuousRecord{$gi}++;
			$giAtPrev{$gi} = $continuousRecord{$gi};
		}
	}
	addMergedGroupToScaffoldOverlapList($sid, $ci, $cNum, \@arrGI, \%giAtPrev, \%continuousRecord, \@contigGroupInfo, \@overGroup);
	$scaffold[$sid]{overGroup} = \@overGroup;

=head
	$overGroup[$j]{gi}
	$overGroup[$j]{startCid}	
	$overGroup[$j]{endCid}
	$overGroup[$j]{strand}
=cut
	#### if the first merged group overlappig this scaffold is including the first contig of this scaffold.. 
	if(scalar @overGroup > 0){
		my $isIncluded = 0;
		if($overGroup[0]{startCid} == 0) { 		
			### 						
			if($overGroup[$#overGroup]{endCid} == $cNum-1){
				$isIncluded = 1;
				for($i = 1; $i <= $#overGroup; $i++){
					if($overGroup[0]{gi} != $overGroup[$i]{gi}) {
						$isIncluded = 0;
						last;
					}
				}
			}
		
			if($isIncluded == 1){
				my $gi = $overGroup[0]{gi};
				

				$scaffold[$sid]{leftGroup} = $gi;
				$scaffold[$sid]{rightGroup} = $gi;
				$scaffold[$sid]{leftStrand} = $overGroup[0]{strand};
				$scaffold[$sid]{rightStrand} = $overGroup[0]{strand};
				
				print STDERR "ERROR at line ".__LINE__.": ($overGroup[0]{strand}), sid : $sid, !defined \$scaffold[12]{leftStrand}, gi : $scaffold[12]{leftGroup}\n" if(!defined $scaffold[12]{leftStrand});

				$i = $#{$mergedGroup[ $gi ]{includeScaffold}} + 1;

				$mergedGroup[ $gi ]{includeScaffold}[$i]{sid} = $sid;
				#### get contig index in group from the matching contigs(startCid, endCid) in a scaffold..
				#### check start and stop.....
				my ($start, $end) = ($contigGroupInfo[ $overGroup[0]{startCid} ]{$gi}{order}, $contigGroupInfo[ $overGroup[0]{endCid} ]{$gi}{order});
				($start, $end) = ($end, $start) if($start > $end);

				#print "group : $mergedGroup[$gi]{name}. start in Group($start), end in Group($end), scaffold : $scaffold[$sid]{name}, start in scaffold($overGroup[0]{startCid}), end in scaffold($overGroup[0]{endCid})\n";
				#exit;
				$mergedGroup[ $gi ]{includeScaffold}[$i]{inStartCid} = $start;
				$mergedGroup[ $gi ]{includeScaffold}[$i]{inEndCid} = $end;
				$mergedGroup[ $gi ]{includeScaffold}[$i]{strand} = $overGroup[0]{strand};
			}
			else {  #### for the left side of the scaffold...
				#### this is a number of contigs in a scaffold overlapping with a merged group..
				my $overNum = $overGroup[0]{endCid} - $overGroup[0]{startCid} + 1;
				my $gi = $overGroup[0]{gi};
				my $strand = $overGroup[0]{strand};

				my ($tmpSid, $tmpOverNum);

				if($strand == $F){
					if(defined $mergedGroup[ $overGroup[0]{gi} ]{rightScaffold}){
						$tmpSid = $mergedGroup[ $gi ]{rightScaffold};
						$conflictStr .=	 " *** The right scaffold for '$mergedGroup[$gi]{name}' " .
							"has been already assigned to $scaffold[$tmpSid]{name}.\n";
						
						#### get the overlapping number of contigs at the previously assigned scaffold..
						if($mergedGroup[ $gi ]{rightStrand} == $F){
							$tmpOverNum = $scaffold[$tmpSid]{overGroup}[0]{endCid} - $scaffold[$tmpSid]{overGroup}[0]{startCid} + 1;
						}
						else{
							$tmpOverNum = $scaffold[$tmpSid]{overGroup}[$#{$scaffold[$tmpSid]{overGroup}}]{endCid} - 
								$scaffold[$tmpSid]{overGroup}[$#{$scaffold[$tmpSid]{overGroup}}]{startCid} + 1;
						}
						if($overNum > $tmpOverNum){
							$conflictStr .=	 "  -> $overNum contigs in $scaffold[$sid]{name} vs $tmpOverNum contigs in $scaffold[$tmpSid]{name} -> replaced to $scaffold[$sid]{name}\n";
							assignMergedGroupToScaffold($left, $strand, $gi, $sid);
						}
					}
					else{
						assignMergedGroupToScaffold($left, $strand, $gi, $sid);
					}
				}
				else{ ### if($strand == $R){
					if(defined $mergedGroup[ $gi ]{leftScaffold}){
						$tmpSid = $mergedGroup[ $gi ]{leftScaffold};
						$conflictStr .=	 " *** The left scaffold for '$mergedGroup[$gi]{name}' " .
							"has been already assigned to $scaffold[$mergedGroup[ $gi ]{leftScaffold}]{name}.\n";
						
						#### get the overlapping number of contigs at the previously assigned scaffold..
						if($mergedGroup[ $gi ]{leftStrand} == $F){
							$tmpOverNum = $scaffold[$tmpSid]{overGroup}[$#{$scaffold[$tmpSid]{overGroup}}]{endCid} - 
								$scaffold[$tmpSid]{overGroup}[$#{$scaffold[$tmpSid]{overGroup}}]{startCid} + 1;
						}
						else{
							$tmpOverNum = $scaffold[$tmpSid]{overGroup}[0]{endCid} - $scaffold[$tmpSid]{overGroup}[0]{startCid} + 1;
						}
						if($overNum > $tmpOverNum){
							$conflictStr .=	 "  -> $overNum contigs in $scaffold[$sid]{name} vs $tmpOverNum contigs in $scaffold[$tmpSid]{name} -> replaced to $scaffold[$sid]{name}\n";
							assignMergedGroupToScaffold($left, $strand, $gi, $sid);
						}
					}					
					else{
						assignMergedGroupToScaffold($left, $strand, $gi, $sid);
					}
				}
			}
		}

		if($isIncluded == 0 && $overGroup[$#overGroup]{endCid} == $cNum-1) { 			
			my $overNum = $overGroup[$#overGroup]{endCid} - $overGroup[$#overGroup]{startCid} + 1;
			my $gi = $overGroup[$#overGroup]{gi};
			my $strand = $overGroup[$#overGroup]{strand};
				
			my ($tmpSid, $tmpOverNum);
			if($strand == $F){
				if(defined $mergedGroup[ $gi ]{leftScaffold}){
					$tmpSid = $mergedGroup[ $gi ]{leftScaffold};
					$conflictStr .=	 " *** The left scaffold for '$mergedGroup[$gi]{name}'" .
						"has been already assigned for $scaffold[$mergedGroup[ $gi ]{leftScaffold}]{name}.\n";
					
					#### get the overlapping number of contigs at the previously assigned scaffold..
					if($mergedGroup[ $gi ]{leftStrand} == $F){
						$tmpOverNum = $scaffold[$tmpSid]{overGroup}[$#{$scaffold[$tmpSid]{overGroup}}]{endCid} - 
							$scaffold[$tmpSid]{overGroup}[$#{$scaffold[$tmpSid]{overGroup}}]{startCid} + 1;
					}
					else{
						$tmpOverNum = $scaffold[$tmpSid]{overGroup}[0]{endCid} - $scaffold[$tmpSid]{overGroup}[0]{startCid} + 1;
					}
					if($overNum > $tmpOverNum){
						$conflictStr .=	 "  -> $overNum contigs in $scaffold[$sid]{name} vs $tmpOverNum contigs in $scaffold[$tmpSid]{name} -> replaced to $scaffold[$sid]{name}\n";
						assignMergedGroupToScaffold($right, $strand, $gi, $sid);
					}
				}
				else{
					assignMergedGroupToScaffold($right, $strand, $gi, $sid);
				}
			}			
			else{
				if(defined $mergedGroup[ $gi ]{rightScaffold}){
					$tmpSid = $mergedGroup[ $gi ]{rightScaffold};
					$conflictStr .=	 " *** The right scaffold for '$mergedGroup[$gi]{name}'" .
						"has been already assigned for $scaffold[$mergedGroup[ $gi ]{rightScaffold}]{name}.\n";
				
					#### get the overlapping number of contigs at the previously assigned scaffold..
					if($mergedGroup[ $gi ]{rightStrand} == $F){
						$tmpOverNum = $scaffold[$tmpSid]{overGroup}[0]{endCid} - $scaffold[$tmpSid]{overGroup}[0]{startCid} + 1;
					}
					else{
						$tmpOverNum = $scaffold[$tmpSid]{overGroup}[$#{$scaffold[$tmpSid]{overGroup}}]{endCid} - 
							$scaffold[$tmpSid]{overGroup}[$#{$scaffold[$tmpSid]{overGroup}}]{startCid} + 1;
					}
					if($overNum > $tmpOverNum){
						
						$conflictStr .=	 "  -> $overNum contigs in $scaffold[$sid]{name} vs $tmpOverNum contigs in $scaffold[$tmpSid]{name} -> replaced to $scaffold[$sid]{name}\n";
						assignMergedGroupToScaffold($right, $strand, $gi, $sid);
					}
				}
				else{
					assignMergedGroupToScaffold($right, $strand, $gi, $sid);
				}
			}
		}
	}
	print $conflictOut $conflictStr . "\n";#  if(defined $verboseFlag);
}



sub assignMergedGroupToScaffold
{
	my ($direction, $strand, $gi, $sid) = @_;
	if($direction == $left){ #### merged contig group is connected to left side of scaffold..
		$scaffold[$sid]{leftGroup} = $gi;
		$scaffold[$sid]{leftStrand} = $strand;
		if($strand == $F){
			$mergedGroup[$gi]{rightScaffold} = $sid;
			$mergedGroup[$gi]{rightStrand} = $strand;
		}
		else{
			$mergedGroup[$gi]{leftScaffold} = $sid;
			$mergedGroup[$gi]{leftStrand} = $strand;
		}
	}
	else{
		$scaffold[$sid]{rightGroup} = $gi;
		$scaffold[$sid]{rightStrand} = $strand;
		if($strand == $F){
			$mergedGroup[$gi]{leftScaffold} = $sid;
			$mergedGroup[$gi]{leftStrand} = $strand;
		}
		else{
			$mergedGroup[$gi]{rightScaffold} = $sid;
			$mergedGroup[$gi]{rightStrand} = $strand;
		}
	}
}



#print "\nsNum : $sNum\n";
$str = $outFn;
$str =~ s/\.fa$/.scaffolds.fa/;
$str =~ s/\.fna$/.scaffolds.fna/;
$str =~ s/\.seq$/.scaffolds.seq/;
$str .= ".scaffolds.fna" if($str eq $outFn);
my $scaOut = openOutput($str);

$str =~ s/\.(fa|fna|seq)$/.txt/;
my $scaLogOut = openOutput($str);
print $scaLogOut "## M : a merged contig in a new scaffold\n" . 
			"## C : a contig in the original scaffold\n" . 
			"## N : a gap in a new scaffold\n" . 
			"#### E : a merged contig at the end of a scaffold\n" . 
			"#### MIS : a merged contig included in a scaffold completely\n" . 
			"#### SIM : a scaffold included in a merged contig completely\n" . 
			"#### J : a merged contig connecting two scaffolds\n";

$str =~ s/\.txt$/.status/;
my $scaStatusOut = openOutput($str);

my $newSid = 1;
my ($gStrand, $newName, $scaLog);
for($sid = 1; $sid < $sNum; $sid++){
	#print "... sid : $sid\n";
	if(!defined $scaffold[$sid]){
		print $scaStatusOut sprintf("scaffold%05d\tSingle contig scaffold\n", $sid);
		next;
	}
	next if(defined $scaffold[$sid]{used});
	
	
	my $prevGi;
	if(defined $scaffold[$sid]{overGroup} && $#{$scaffold[$sid]{overGroup}} > -1){
		#print "--- sid : $sid, overGroup : ".($#{$scaffold[$sid]{overGroup}}+1)."\n" if(defined $verboseFlag);;
		$newName = sprintf("newScaffold%04d", $newSid);
		print $scaOut ">$newName\n";
		
			
		my $headSid = $sid;
		$gi = $scaffold[$sid]{leftGroup};
		my $strand = $F;
		### add code to prevent recusive...

		while(defined $gi){
			#### $strand : strand of scaffold...$scaffold[$headSid]{leftStrand/rightStrand} : strand of merged contig group.
			### ex. of $strand == $R && $scaffold[$headSid]{rightStrand} == $R : ?--s(next)--? --g--> <--s(current)--
			#print " headsid : $headSid, strand : $strandTxt[$strand], group : $gi, left : $strandTxt[$scaffold[$headSid]{leftStrand}], right : $strandTxt[$scaffold[$headSid]{rightStrand}]\n";
			if(($strand == $F && defined $scaffold[$headSid]{leftStrand} && $scaffold[$headSid]{leftStrand} == $F) || 
				($strand == $R && defined $scaffold[$headSid]{rightStrand} && $scaffold[$headSid]{rightStrand} == $R)){
				if(defined $mergedGroup[$gi]{leftScaffold}){
					#print " next sid : defined $mergedGroup[$gi]{leftScaffold}\n";
					$strand = $mergedGroup[$gi]{leftStrand}; #### strand of the next scaffold..
					$headSid = $mergedGroup[$gi]{leftScaffold};		
					$gi = ($strand ==$F ? $scaffold[$headSid]{leftGroup} : $scaffold[$headSid]{rightGroup});
				}
				else{ last; }
			}
			else{
				if(defined $mergedGroup[$gi]{rightScaffold}){					
					#print " next sid : defined $mergedGroup[$gi]{rightScaffold}\n";
					$strand = ($mergedGroup[$gi]{rightStrand} == $F ? $R : $F);
					$headSid = $mergedGroup[$gi]{rightScaffold};				
					#$gi = $scaffold[$headSid]{rightGroup};					
					$gi = ($strand == $F ? $scaffold[$headSid]{leftGroup} : $scaffold[$headSid]{rightGroup});
				}
				else{ last; }
			}
			if($headSid == $sid){
				print " \n ** the current headSid is same with seed sid $sid. Searching was stopped.\n" if(defined $verboseFlag);
				last;
			}
		}
		my $tmpSid = $headSid;
		#print "--head scaffold : $scaffold[$tmpSid]{name}, $strandTxt[$strand]\n" if(defined $verboseFlag);
		
		$scaLog = "";
		print "$newName graph : " if(defined $verboseFlag);
		print $scaLogOut "#$newName graph : ";

		my $startInNew = 1;
		if($strand == $F && defined $scaffold[$tmpSid]{leftGroup}){
			#print "first GI :\n";
			print " $mergedGroup[$scaffold[$tmpSid]{leftGroup}]{name}($strandTxt[$scaffold[$tmpSid]{leftStrand}]) -" if(defined $verboseFlag);
			print $scaLogOut " $mergedGroup[$scaffold[$tmpSid]{leftGroup}]{name}($strandTxt[$scaffold[$tmpSid]{leftStrand}]) -";

			if(!defined $scaffold[$tmpSid]{leftStrand}){
				print STDERR "ERROR at line ".__LINE__.": !defined \$scaffold[$tmpSid]{leftStrand}\n";
				exit;
			}

			#print $scaLogOut sprintf("mergedScaffold%d\tGrouping\n", $newSid);
			$gStrand = ($strand == $scaffold[$tmpSid]{leftStrand} ? $F : $R);
			$startInNew = printMergedGroup($scaffold[$tmpSid]{leftGroup}, $gStrand, $startInNew, $newName);
		}
		elsif($strand == $R && defined $scaffold[$tmpSid]{rightGroup}){			
			print " $mergedGroup[$scaffold[$tmpSid]{rightGroup}]{name}($strandTxt[$scaffold[$tmpSid]{rightStrand}]) -" if(defined $verboseFlag);
			print $scaLogOut " $mergedGroup[$scaffold[$tmpSid]{rightGroup}]{name}($strandTxt[$scaffold[$tmpSid]{rightStrand}]) -";
			#print $scaLogOut sprintf("mergedScaffold%d\tGrouping\n", $newSid);
			$gStrand = ($strand == $scaffold[$tmpSid]{rightStrand} ? $F : $R);
			$startInNew = printMergedGroup($scaffold[$tmpSid]{rightGroup}, $gStrand, $startInNew, $newName);
		}
		print " $scaffold[$tmpSid]{name}($strandTxt[$strand]) -" if(defined $verboseFlag);
		print $scaLogOut " $scaffold[$tmpSid]{name}($strandTxt[$strand]) -";
		$startInNew = printScaffold($tmpSid, $strand, $startInNew, $newName);

		$prevGi = $gi;
		$gi = ($strand == $F ? $scaffold[$tmpSid]{rightGroup} : $scaffold[$tmpSid]{leftGroup});		
				
		while(defined $gi){	
			#print STDERR "$strand, right : $scaffold[$tmpSid]{rightStrand}, left : $scaffold[$tmpSid]{leftStrand}\n";
			$gStrand = ($strand == $F ? ($strand == $scaffold[$tmpSid]{rightStrand} ? $F : $R) : 
				($strand == $scaffold[$tmpSid]{leftStrand} ?	$F : $R) ); ## for strand == $R

			print " $mergedGroup[$gi]{name}($strandTxt[$gStrand]) -" if(defined $verboseFlag);
			print $scaLogOut " $mergedGroup[$gi]{name}($strandTxt[$gStrand]) -";
			if(!defined $prevGi || $prevGi != $gi)
			{
				if($strand == $F) { 				
					$startInNew = printMergedGroup($scaffold[$tmpSid]{rightGroup}, $gStrand, $startInNew, $newName);
				}
				else{
					$startInNew = printMergedGroup($scaffold[$tmpSid]{leftGroup}, $gStrand, $startInNew, $newName);
				}
			}
			if(($strand == $F && defined $scaffold[$tmpSid]{rightStrand} && $scaffold[$tmpSid]{rightStrand} == $F) || 
				($strand == $R && defined $scaffold[$tmpSid]{leftStrand} && $scaffold[$tmpSid]{leftStrand} == $R)){
				#print "\n--- 1, $strandTxt[$strand]\n";
				if(defined $mergedGroup[$gi]{rightScaffold}){					
					$strand = $mergedGroup[$gi]{rightStrand}; #### strand of the next scaffold..
					$tmpSid = $mergedGroup[$gi]{rightScaffold};
					$prevGi = $gi;
					$gi = ($strand ==$F ? $scaffold[$tmpSid]{rightGroup} : $scaffold[$tmpSid]{leftGroup});
					$gStrand = $F;
					$scaffold[$tmpSid]{used} = 1;
				}
				else{ $tmpSid = undef; last; }
			}
			else{			
				#print "\n--- 2, $strandTxt[$strand]\n";
				if(defined $mergedGroup[$gi]{leftScaffold}){
					$strand = ($mergedGroup[$gi]{leftStrand} == $F ? $R : $F);
					$tmpSid = $mergedGroup[$gi]{leftScaffold};
					$prevGi = $gi;
					$gi = ($strand ==$F ? $scaffold[$tmpSid]{rightGroup} : $scaffold[$tmpSid]{leftGroup});
					$gStrand = $R;
					$scaffold[$tmpSid]{used} = 1;
					#print "selected gi : $gi\n";
				}
				else{ $tmpSid = undef; last; }
				
			}	
			if(defined $tmpSid){
				print " $scaffold[$tmpSid]{name}($strandTxt[$strand]) -" if(defined $verboseFlag);
				print $scaLogOut " $scaffold[$tmpSid]{name}($strandTxt[$strand]) -";
				$startInNew = printScaffold($tmpSid, $strand, $startInNew, $newName);
			}
			if($headSid == $tmpSid){
				print " \n**** the current node sid is same with head sid $headSid. Searching was stopped.\n" if(defined $verboseFlag);
				last;
			}
		}
		
		print "\n\n" if(defined $verboseFlag);
		print $scaOut "\n";
		print $scaLogOut "\n$scaLog";
		print $scaStatusOut "$scaffold[$sid]{name}\tExtended by $newName\n";


		$newSid++;
	}
	else{
		print $scaStatusOut "$scaffold[$sid]{name}\tNot extended\n";
	}

	my $cNum = $#{$scaffold[$sid]{contig}} + 1;
}

close($scaOut);
close($scaLogOut);
close($scaStatusOut);

printMergedSeq();


sub printMergedSeq{

	print "Total merged groups : ".($gNum-1)."\n";

	my $rvSeq;
	for($gi = 1; $gi < $gNum; $gi++){
		print $out ">$mergedGroup[$gi]{name}\n";
		$nodeNum = $#{$mergedGroup[$gi]{nodeList}} + 1;
		for($i = 0; $i < $nodeNum; $i++){	
			if(!defined $contigSeq{$mergedGroup[$gi]{nodeList}[$i]}){
				print STDERR "No seq : $mergedGroup[$gi]{nodeList}[$i], !defined $contigSeq{$mergedGroup[$gi]{nodeList}[$i]} || !defined $mergedGroup[$gi]{seqPos}[$i]{start} || !defined $mergedGroup[$gi]{seqPos}[$i]{end}\n"; 
				exit;
			}
			next if(!defined $mergedGroup[$gi]{seqPos}[$i]{start} || !defined $mergedGroup[$gi]{seqPos}[$i]{end});

			if($mergedGroup[$gi]{strandList}[$i] == $F) {			
				print $out substr($contigSeq{$mergedGroup[$gi]{nodeList}[$i]}, 
					$mergedGroup[$gi]{seqPos}[$i]{start}-1, $mergedGroup[$gi]{seqPos}[$i]{end}-$mergedGroup[$gi]{seqPos}[$i]{start}+1);
			}
			else {
				$rvSeq = rc($contigSeq{$mergedGroup[$gi]{nodeList}[$i]});
				print $out substr($rvSeq, 
					$contigLen{$mergedGroup[$gi]{nodeList}[$i]} - $mergedGroup[$gi]{seqPos}[$i]{start}, 
					$mergedGroup[$gi]{seqPos}[$i]{start}-$mergedGroup[$gi]{seqPos}[$i]{end}+1);
			}
		}
		print $out "\n";
	}


	foreach $name (keys %endMatches) {
		print $out ">$name\n$contigSeq{$name}\n" if(!defined $deleteHash{$name} && defined $contigSeq{$name} && $name =~ /$keepHead/);
	}


	foreach $name (keys %includeHigh) {
		if(!defined $deleteHash{$name} && !defined $includeLow{$name} && $name !~ /^$keepHead/ ){ ## !defined $endMatches{$name} && 
			my $inName;
			foreach $qName (keys %{$includeHigh{$name}}) {

				### $includeLowMark{$qName} == 2 : a contig including qName contig is also for keeping..
				### so.. this $name contig does not need to be printed.
				if($qName =~ /^$keepHead/ && $includeLowMark{$qName} != 2){ 
					$inName = $qName;
					last;
				}
			}
			if(defined $inName && !defined $deleteHash{$name} && defined $contigSeq{$name})
			{ 
				print " - $name is including $inName, print sequence..\n" if(defined $verboseFlag);
				print $out ">$name\n$contigSeq{$name}\n";	
			}
		}
		elsif(!defined $endMatches{$name} && !defined $deleteHash{$name} && !defined $includeLow{$name} && $name =~ /^$keepHead/ && defined $contigSeq{$name}){
			print " - $name, print sequence..\n" if(defined $verboseFlag);
			print $out ">$name\n$contigSeq{$name}\n";	
		}
	}


	close($out) if(defined $outFn);

	print "'$cmd' has been finished.\n" if(defined $verboseFlag);
}


sub printScaffold{
	my ($sid, $strand, $startInNew, $newName) = @_;

	my ($i, $j, $cNum, $gNum, $loopEnd, $totlen, $len, $seq);
	$totlen = 0;
		
	$cNum = $#{$scaffold[$sid]{contig}} + 1;
	$gNum = $#{$scaffold[$sid]{overGroup}} + 1;	
	if($strand == $F){
		for($ci = 0, $j = 0; $ci < $cNum; ){
			#print STDERR " out : $ci < $cNum, j : $j\n";:$/

			if($j < $gNum){ $loopEnd = $scaffold[$sid]{overGroup}[$j]{startCid};	}
			else{ $loopEnd = $cNum-1;	}
			
			for(;$ci <= $loopEnd; $ci++){

				#print STDERR " in : $ci <= $loopEnd, j : $j\n";
				if($ci != 0)	{
					$len = $scaffold[$sid]{gap}[$ci-1]{len};
					#print STDERR "$ci gap len $len\n";

					if($j < $gNum && defined $scaffold[$sid]{overGroup}[$j] && $scaffold[$sid]{overGroup}[$j]{startCid} == $ci){
						$len -= $scaffold[$sid]{overGroup}[$j]{leftAdditionalLen};
						#print STDERR "left : $scaffold[$sid]{overGroup}[$j]{leftAdditionalLen}, $ci gap len $len\n";
						if($len < 0){
							print STDERR "left : $scaffold[$sid]{overGroup}[$j]{leftAdditionalLen}, $ci gap len $len\n";
							print STDERR " ------ length of gap is 0, $scaffold[$sid]{name}\n"; exit;
						}
					}
					if($j > 0 && defined $scaffold[$sid]{overGroup}[$j-1] && $scaffold[$sid]{overGroup}[$j-1]{endCid} == $ci-1){
						$len -= $scaffold[$sid]{overGroup}[$j-1]{rightAdditionalLen};
						
						#if($len < 0){
						#	print STDERR "right : $scaffold[$sid]{overGroup}[$j-1]{rightAdditionalLen}, $ci gap len $len\n";
						#	print STDERR " ------ length of gap is 0, $scaffold[$sid]{name}\n"; exit;
						#}
					}
					
					#print STDERR "$ci gap len $len\n";
					$len = 20 if($len < 20);
					print $scaOut "N"x$len;
					$scaLog .= "$newName\tN\t$scaffold[$sid]{name}\t$startInNew\t" . ($startInNew+$len-1) . "\t$strandTxt[$strand]\n";
					$startInNew += $len;
					$totlen += $len;
				} #### for if($ci != 0)..

				#if(($j == $gNum && ($j == 0 || $ci != $scaffold[$sid]{overGroup}[$j-1]{endCid} ) ) || 
				#	($j < $gNum && $ci < $scaffold[$sid]{overGroup}[$j]{startCid})){
				if($j == $gNum || $ci < $scaffold[$sid]{overGroup}[$j]{startCid}){
					
					if(!defined $contigSeq{$scaffold[$sid]{contig}[$ci]{name}}){
						print STDERR "ERROR at line ".__LINE__.": No sequence for $scaffold[$sid]{contig}[$ci]{name}, newName : $newName\n";
						exit;
					}

					$seq = \$contigSeq{$scaffold[$sid]{contig}[$ci]{name}}; ### for speeed, assigned a point instead of copying a sequence.
					$len = length($$seq);					
					print $scaOut ($strand == $F ? $$seq : rc($$seq));
					$scaLog .= "$newName\tC\t$scaffold[$sid]{name}\t$scaffold[$sid]{contig}[$ci]{name}\t$startInNew\t" . ($startInNew+$len-1) . "\t$strandTxt[$strand]\n";
					$startInNew += $len;
					$totlen += $len;
				}
			}
			
			#### if the merged group is not overlapping with the first or the last contig in a scaffold
			if($j < $gNum && $scaffold[$sid]{overGroup}[$j]{startCid} != 0 && $scaffold[$sid]{overGroup}[$j]{endCid} != $cNum -1){
				my $gi = $scaffold[$sid]{overGroup}[$j]{gi};
				my $gStrand = ($strand == $scaffold[$sid]{overGroup}[$j]{strand} ? $F : $R);
				$seq = getMergedGroupSequence($gi, $gStrand);
				$len = length($seq);
				print $scaOut $seq;
				$scaLog .= "$newName\tM\t$mergedGroup[$gi]{name}\t$startInNew\t" . ($startInNew+$len-1) . "\t$strandTxt[$gStrand]\n";
				
				#my ($startStr, $endStr) = ($gStrand == $F ? ("startCid", "endCid") : ("endCid", "startCid"));
				$scaLog .= "\tMIS\t$mergedGroup[$gi]{name}\t$scaffold[$sid]{name}\t". 
					"$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{startCid}]{name}(". 
					($scaffold[$sid]{overGroup}[$j]{startCid}+1) . 
					"-th)\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{endCid}]{name}(". 
					($scaffold[$sid]{overGroup}[$j]{endCid}+1) .
					"-th)\t$strandTxt[$gStrand]\n";				
				$startInNew += $len;
				$totlen += $len;
			}	
			elsif($j < $gNum && $scaffold[$sid]{overGroup}[$j]{endCid} == $cNum -1){ last; }

			$j++;

			#print STDERR "$ci, -----\n";
		}
	}
	else{ ### if $strand == $R
		for($ci = $cNum-1, $j = $gNum -1; $ci >= 0; ){					
			if($j >= 0){ $loopEnd = $scaffold[$sid]{overGroup}[$j]{endCid};	}
			else{ $loopEnd = 0;	}

			if(!defined $loopEnd){
				print STDERR "ERROR at line ".__LINE__.": $j-th overGroup of $scaffold[$sid]{name} is undefined. gNum : $gNum, $scaffold[$sid]{overGroup}[$j], $scaffold[$sid]{overGroup}[$j]{startCid};\n";
				exit;
			}

			for(;$ci >= $loopEnd; $ci--){
				if($ci != $cNum-1)	{
					$len = $scaffold[$sid]{gap}[$ci]{len};
					
					if($j < $gNum-1 && defined $scaffold[$sid]{overGroup}[$j+1] && $scaffold[$sid]{overGroup}[$j+1]{startCid} == $ci+1){
						$len -= $scaffold[$sid]{overGroup}[$j]{leftAdditionalLen};
						if($len < 0){
							#print STDERR " ------ length of gap is 0, $scaffold[$sid]{name}\n"; exit;
						}
					}
					if($j >= 0 && defined $scaffold[$sid]{overGroup}[$j] && $scaffold[$sid]{overGroup}[$j]{endCid} == $ci){
						$len -= $scaffold[$sid]{overGroup}[$j]{rightAdditionalLen};
						if($len < 0){
							#print STDERR " ------ length of gap is 0, $scaffold[$sid]{name}\n"; exit;
						}
					}
					$len = 20 if($len < 20);
					print $scaOut "N"x$len;
					$scaLog .= "$newName\tN\t$scaffold[$sid]{name}\t$startInNew\t" . ($startInNew+$len-1) . "\t$strandTxt[$strand]\n";
					$startInNew += $len;
					$totlen += $len;
				} #### for if($ci != 0)..

				if($j < 0 || $ci > $scaffold[$sid]{overGroup}[$j]{endCid}){
					$seq = \$contigSeq{$scaffold[$sid]{contig}[$ci]{name}}; ### for speeed, assigned a point instead of copying a sequence.
					$len = length($$seq);
					print $scaOut ($strand == $F ? $$seq : rc($$seq));
					$scaLog .= "$newName\tC\t$scaffold[$sid]{name}\t$scaffold[$sid]{contig}[$ci]{name}\t$startInNew\t" . ($startInNew+$len-1) . "\t$strandTxt[$strand]\n";				
					$startInNew += $len;
					$totlen += $len;
				}
			}

			#### if the merged group is not overlapping with the first or the last contig in a scaffold
			if($j >= 0 && $scaffold[$sid]{overGroup}[$j]{startCid} != 0 && $scaffold[$sid]{overGroup}[$j]{endCid} != $cNum -1){ 
				my $gi = $scaffold[$sid]{overGroup}[$j]{gi};
				my $gStrand = ($strand == $scaffold[$sid]{overGroup}[$j]{strand} ? $F : $R);
				$seq = getMergedGroupSequence($gi, $gStrand);
				$len = length($seq);
				print $scaOut $seq;
				$scaLog .= "$newName\tM\t$mergedGroup[$gi]{name}\t$strandTxt[$gStrand]\t$startInNew\t" . ($startInNew+$len-1) . "\n";
				$scaLog .= "\tMIS\t$mergedGroup[$gi]{name}\t$scaffold[$sid]{name}\t". 
					"$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{endCid}]{name}(". 
					($scaffold[$sid]{overGroup}[$j]{endCid}+1) . 
					"-th)\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{startCid}]{name}(". 
					($scaffold[$sid]{overGroup}[$j]{startCid}+1) .
					"-th)\t$strandTxt[$gStrand]\n";				
				$startInNew += $len;
				$totlen += $len;
			}	
			$j--;
		}
	}
	return $startInNew;
}

sub printMergedGroup{
	my ($gi, $strand, $startInNew, $newName) = @_;
	my ($seq, $len);
	$seq = getMergedGroupSequence($gi, $strand);
	$len = length($seq);
	print $scaOut $seq;
	my ($sid, $cNum, $sNum, $i, $j);
	$scaLog .= "$newName\tM\t$mergedGroup[$gi]{name}\t$startInNew\t" . ($startInNew+$len-1) . "\t$strandTxt[$strand]\n";
	if($strand == $F){ #### if strand of this group is forward...
		if(defined $mergedGroup[$gi]{leftScaffold}){
			$sid = $mergedGroup[$gi]{leftScaffold};
			$j = ($mergedGroup[$gi]{leftStrand} == $F ? $#{$scaffold[$sid]{overGroup}} : 0);
			###mergedGroup[$gi]{nodeList}[$scaffold[$sid]{overGroup}[$j]{startCid}]
			$scaLog .= "\t" . (defined $mergedGroup[$gi]{rightScaffold} ? "J" : "E") . 
				"\t$mergedGroup[$gi]{name}\t$scaffold[$sid]{name}\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{startCid}]{name}(". 
				($scaffold[$sid]{overGroup}[$j]{startCid}+1) . "-th)\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{endCid}]{name}(". 
				($scaffold[$sid]{overGroup}[$j]{endCid}+1) .
				"-th)\t".($strand == $scaffold[$sid]{overGroup}[$j]{strand} ? $strandTxt[$F] : $strandTxt[$R])."\n";
		}
			
		if(defined $mergedGroup[$gi]{includeScaffold}){
			$sNum = $#{$mergedGroup[$gi]{includeScaffold}} + 1;		
			my @includeScaffold = sort {$a->{inStartCid} <=> $b->{inStartCid}} @{$mergedGroup[$gi]{includeScaffold}};
			$j = 0;
			for($i = 0; $i < $sNum; $i++){
				$sid = $includeScaffold[$i]{sid};				
				$scaLog .= "\tSIM\t$mergedGroup[$gi]{name}\t$scaffold[$sid]{name}\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{startCid}]{name}(". 
					($scaffold[$sid]{overGroup}[$j]{startCid}+1) . "-th)\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{endCid}]{name}(". 
					($scaffold[$sid]{overGroup}[$j]{endCid}+1) .
					"-th)\t".($strand == $scaffold[$sid]{overGroup}[$j]{strand} ? $strandTxt[$F] : $strandTxt[$R])."\n";				
			}
		}
		if(defined $mergedGroup[$gi]{rightScaffold}){
			$sid = $mergedGroup[$gi]{rightScaffold};
			$j = ($mergedGroup[$gi]{rightStrand} == $F ? 0 : $#{$scaffold[$sid]{overGroup}});
			$scaLog .= "\t" . (defined $mergedGroup[$gi]{leftScaffold} ? "J" : "E") . 
				"\t$mergedGroup[$gi]{name}\t$scaffold[$sid]{name}\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{startCid}]{name}(". 
				($scaffold[$sid]{overGroup}[$j]{startCid}+1) . "-th)\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{endCid}]{name}(". 
				($scaffold[$sid]{overGroup}[$j]{endCid}+1) .
				"-th)\t".($strand == $scaffold[$sid]{overGroup}[$j]{strand} ? $strandTxt[$F] : $strandTxt[$R])."\n";				
			#$scaLog .= "\t$scaffold[$sid]{name}\t$scaffold[$sid]{overGroup}[$j]{startCid}\t$scaffold[$sid]{overGroup}[$j]{endCid}\t$strandTxt[$scaffold[$sid]{overGroup}[$j]{strand}]\n";
		}
	}
	else{  #### if strand of this group is reverse...
		if(defined $mergedGroup[$gi]{rightScaffold}){
			$sid = $mergedGroup[$gi]{rightScaffold};
			$j = ($mergedGroup[$gi]{rightStrand} == $F ? 0 : $#{$scaffold[$sid]{overGroup}});
			$scaLog .= "\t" . (defined $mergedGroup[$gi]{leftScaffold} ? "J" : "E") . 
				"\t$mergedGroup[$gi]{name}\t$scaffold[$sid]{name}\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{startCid}]{name}(". 
				($scaffold[$sid]{overGroup}[$j]{startCid}+1) . "-th)\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{endCid}]{name}(". 
				($scaffold[$sid]{overGroup}[$j]{endCid}+1) .
				"-th)\t".($strand == $scaffold[$sid]{overGroup}[$j]{strand} ? $strandTxt[$F] : $strandTxt[$R])."\n";
			#$scaLog .= "\t$scaffold[$sid]{name}\t$scaffold[$sid]{overGroup}[$j]{startCid}\t$scaffold[$sid]{overGroup}[$j]{endCid}\t$strandTxt[$scaffold[$sid]{overGroup}[$j]{strand}]\n";
		}
		if(defined $mergedGroup[$gi]{includeScaffold}){
			$sNum = $#{$mergedGroup[$gi]{includeScaffold}} + 1;		
			my @includeScaffold = sort {$b->{inStartCid} <=> $a->{inStartCid}} @{$mergedGroup[ $gi ]{includeScaffold}};
			$j = 0;
			for($i = 0; $i < $sNum; $i++){
				$sid = $includeScaffold[$i]{sid};
				$scaLog .= "\tSIM\t$mergedGroup[$gi]{name}\t$scaffold[$sid]{name}\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{startCid}]{name}(". 
					($scaffold[$sid]{overGroup}[$j]{startCid}+1) . "-th)\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{endCid}]{name}(". 
					($scaffold[$sid]{overGroup}[$j]{endCid}+1) .
					"-th)\t".($strand == $scaffold[$sid]{overGroup}[$j]{strand} ? $strandTxt[$F] : $strandTxt[$R])."\n";
				#$scaLog .= "\t$scaffold[$sid]{name}\t$scaffold[$sid]{overGroup}[0]{endCid}\t$scaffold[$sid]{overGroup}[0]{startCid}\t$strandTxt[$scaffold[$sid]{overGroup}[0]{strand}]\n";
			}
		}
		if(defined $mergedGroup[$gi]{leftScaffold}){
			$sid = $mergedGroup[$gi]{leftScaffold};
			$j = ($mergedGroup[$gi]{leftStrand} == $F ? $#{$scaffold[$sid]{overGroup}} : 0);
			$scaLog .= "\t" . (defined $mergedGroup[$gi]{rightScaffold} ? "J" : "E") . 
				"\t$mergedGroup[$gi]{name}\t$scaffold[$sid]{name}\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{startCid}]{name}(". 
				($scaffold[$sid]{overGroup}[$j]{startCid}+1) . "-th)\t$scaffold[$sid]{contig}[$scaffold[$sid]{overGroup}[$j]{endCid}]{name}(". 
				($scaffold[$sid]{overGroup}[$j]{endCid}+1) .
				"-th)\t".($strand == $scaffold[$sid]{overGroup}[$j]{strand} ? $strandTxt[$F] : $strandTxt[$R])."\n";
			#$scaLog .= "\t$scaffold[$sid]{name}\t$scaffold[$sid]{overGroup}[$j]{startCid}\t$scaffold[$sid]{overGroup}[$j]{endCid}\t$strandTxt[$scaffold[$sid]{overGroup}[$j]{strand}]\n";
		}		
	}

	return $startInNew+$len;
}

sub getMergedGroupSequence
{
	my ($gi, $strand) = @_;
	my ($i, $nodeNum, $seq);
	$seq = "";
	$nodeNum = $#{$mergedGroup[$gi]{nodeList}} + 1;

	for($i = 0; $i < $nodeNum; $i++){	
		if(!defined $contigSeq{$mergedGroup[$gi]{nodeList}[$i]} || !defined $mergedGroup[$gi]{seqPos}[$i]{start} || !defined $mergedGroup[$gi]{seqPos}[$i]{end}){
			print STDERR "No seq : $mergedGroup[$gi]{nodeList}[$i], !defined $contigSeq{$mergedGroup[$gi]{nodeList}[$i]} || !defined $mergedGroup[$gi]{seqPos}[$i]{start} || !defined $mergedGroup[$gi]{seqPos}[$i]{end}\n"; 
			exit;
		}
		if($mergedGroup[$gi]{strandList}[$i] == $F) {			
			$seq .= substr($contigSeq{$mergedGroup[$gi]{nodeList}[$i]}, 
				$mergedGroup[$gi]{seqPos}[$i]{start}-1, $mergedGroup[$gi]{seqPos}[$i]{end}-$mergedGroup[$gi]{seqPos}[$i]{start}+1);
		}
		else {
			my $rvSeq = rc($contigSeq{$mergedGroup[$gi]{nodeList}[$i]});
			$seq .= substr($rvSeq, 
				$contigLen{$mergedGroup[$gi]{nodeList}[$i]} - $mergedGroup[$gi]{seqPos}[$i]{start}, 
				$mergedGroup[$gi]{seqPos}[$i]{start}-$mergedGroup[$gi]{seqPos}[$i]{end}+1);
		}
	}

	return ($strand == $F ? $seq : rc($seq));
}


sub addMergedGroupToScaffoldOverlapList
{
	my ($sid, $ci, $cNum, $arrGI, $giAtPrev, $continuousRecord, $contigGroupInfo, $overGroup) = @_;
	my ($gi, $startCid, $endCid, $gStrand, $isStop, $rName, $i);

	my $logForCurrent_str = ""; #### for log...
	if($ci < $cNum){
		@$arrGI = keys %{$contigGroupInfo->[$ci]};
		if(scalar @$arrGI > 1){
			$conflictStr .= " ** $scaffold[$sid]{name}, (".($ci+1)."-th contig is included in ".(scalar @$arrGI)." groups. : @$arrGI \n";
		}
		$logForCurrent_str .= " merged groups at (".($ci+1)."-th contig($scaffold[$sid]{contig}[$ci]{name}) : ";
		foreach $gi (@$arrGI) {
			$logForCurrent_str .= "$mergedGroup[$gi]{name}, ";
			if(defined $giAtPrev->{$gi}){
				delete $giAtPrev->{$gi};  #### as a result, %giAtPrev will have only 'gi' ending at the previous contig..
			}
		}
		$logForCurrent_str .= "\n";
	}
	
	### @overGroup : for overlapping merged contigs..
	#### if there is any merged contig(contig group) ends at the previous($ci-1) contig...
	foreach $gi (keys %{$giAtPrev}) {					
		$rName = "";			
		$gStrand = $mergedGroup[$gi]{strandList}[$contigGroupInfo->[$ci-1]{$gi}{order}]; ### strand of the merged contig group
		$startCid = $ci-$continuousRecord->{$gi}; #### contig index at sid..			
		$endCid = $ci-1; #### contig index at sid..

		#print " -- Group $gi end($contigGroupInfo->[$ci-1]{$gi}{order}) at scaffold_$sid\[" .($ci-1)."] contig and (" . ($giAtPrev->{$gi}) . ") contigs are included \n" if(defined $verboseFlag);

		##### check.. order is correct...			
		$isStop = 0;
		$conflictStr .= " - $mergedGroup[$gi]{name} covers $scaffold[$sid]{name} from $scaffold[$sid]{contig}[$startCid]{name}(".($startCid+1)."-th) to $scaffold[$sid]{contig}[$endCid]{name}(".($endCid+1)."-th)\n";
		if($gStrand == $F){
			for($i = $startCid; $i < $endCid; $i++){
				if(!defined $contigGroupInfo->[$i]{$gi}{order} || !defined $contigGroupInfo->[$i+1]{$gi}{order}){
					print STDERR  "ERROR at line ".__LINE__.": i:$i, gi:$gi, $contigGroupInfo->[$i]{$gi}, defined $contigGroupInfo->[$i]{$gi}{order} || !defined $contigGroupInfo->[$i+1]{$gi}{order}\n"; 
					exit;
				}
				if($contigGroupInfo->[$i]{$gi}{order} > $contigGroupInfo->[$i+1]{$gi}{order}){
					$conflictStr .= " ** Order of contigs in $scaffold[$sid]{name} is not consistent with order of contigs in '$mergedGroup[$gi]{name}'. '$mergedGroup[$gi]{name}' was not merged.\n" if(defined $verboseFlag);
					last;
				}
			}
		}
		else{
			for($i = $startCid; $i < $endCid; $i++){
				if($contigGroupInfo->[$i]{$gi}{order} < $contigGroupInfo->[$i+1]{$gi}{order}){
					$conflictStr .= " ** Order of contigs in $scaffold[$sid]{name} is not consistent with order of contigs '$mergedGroup[$gi]{name}. '$mergedGroup[$gi]{name}' was not merged.'\n" if(defined $verboseFlag);
					last;
				}
			}
		}
		
		if($isStop == 1){ 
			delete $continuousRecord->{$gi};
			next;
		}

		##### check the length is about right...
		my ($leftAdditionalLen, $rightAdditionalLen) = (0, 0);
		my $gStrand;
		
		my $incMatch;
		
		if($ci == $cNum){
			$j = $contigGroupInfo->[$ci-1]{$gi}{order};
			if(defined $contigGroupInfo->[$ci-1]{$gi}{rName}){ ### if ($ci-1)th contig in the scaffold is included in another contig in a merged group..
				$qName = $scaffold[$sid]{contig}[$ci-1]{name};
				$rName = $contigGroupInfo->[$ci-1]{$gi}{rName};
				$incMatch = $includeLow{$scaffold[$sid]{contig}[$ci-1]{name}}{$rName};
				$gStrand = ($mergedGroup[$gi]{strandList}[$j] == $incMatch->{strand} ? $F : $R);
			}
			else{
				$rName = $scaffold[$sid]{contig}[$ci-1]{name};
				$gStrand = $mergedGroup[$gi]{strandList}[$j];
			}
		}
		if($ci < $cNum){  #### if the current contig is not the last contig in the scaffold...						
			$j = $contigGroupInfo->[$ci-1]{$gi}{order}; #### index of the contig in the merged contig group.
			if(defined $contigGroupInfo->[$ci-1]{$gi}{rName}){ ### if ($ci-1)th contig in the scaffold is included in another contig in a merged group..
				$qName = $scaffold[$sid]{contig}[$ci-1]{name};
				$rName = $contigGroupInfo->[$ci-1]{$gi}{rName};
=head
				$includeLow{$qName}{$rName}{lStart}
				$includeLow{$qName}{$rName}{lEnd}
				$includeLow{$qName}{$rName}{hStart}
				$includeLow{$qName}{$rName}{hEnd}
				$includeLow{$qName}{$rName}{strand}
=cut
				$incMatch = $includeLow{$scaffold[$sid]{contig}[$ci-1]{name}}{$rName};
				$gStrand = ($mergedGroup[$gi]{strandList}[$j] == $incMatch->{strand} ? $F : $R);
				if($includeLow{$qName}{$rName}{strand} == $F){
					$rightAdditionalLen = $contigLen{$rName} - $includeLow{$qName}{$rName}{hEnd};
				}
				else{ #### in the blast match, if strand == R, hEnd < hStart
					$rightAdditionalLen = $includeLow{$qName}{$rName}{hEnd} -1;
				}
				#print STDERR "rightAdditionalLen : $rightAdditionalLen, sid $scaffold[$sid]{name}, $scaffold[$sid]{contig}[$ci-1]{name}, gi $mergedGroup[$gi]{name}\n";				
			}
			else{
				$rName = $scaffold[$sid]{contig}[$ci-1]{name};
				$gStrand = $mergedGroup[$gi]{strandList}[$j];
			}
			#print STDERR " ----  rightAdditionalLen : $rightAdditionalLen, sid $scaffold[$sid]{name}, $scaffold[$sid]{contig}[$ci-1]{name}, gi $mergedGroup[$gi]{name}\n";
			if($gStrand == $F && $j != $mergedGroup[$gi]{num}-1){
				$rightAdditionalLen += $mergedGroup[$gi]{len} - 
					($mergedGroup[$gi]{seqPos}[$j]{rEnd} + 
					(!defined $incMatch || $incMatch->{strand} == $F ? 
							($contigLen{$rName}-$mergedGroup[$gi]{seqPos}[$j]{end}) : $mergedGroup[$gi]{seqPos}[$j]{end} ) ); 
					### considering --group--> --hContig--> --lContig--> or --group--> <--hContig-- --lContig-->
					### hContig : contig covering the lContig, lContig : contig in a scaffold
					#### if --group--> <--hContig--, {start} > {end}..

					##### $mergedGroup[$gi]{seqPos}[$j]{end} is an end position of the alignment.
					
					#if(!defined $incMatch || $incMatch->{strand} == $F ){
					#	print STDERR "rightAdditionalLen : $rightAdditionalLen, rName : $rName, $mergedGroup[$gi]{len} - ($mergedGroup[$gi]{seqPos}[$j]{rEnd} + ($contigLen{$rName}-$mergedGroup[$gi]{seqPos}[$j]{end}))\n";
					#}
					#else{
					#	print STDERR "rightAdditionalLen : $rightAdditionalLen, $mergedGroup[$gi]{len} - ($mergedGroup[$gi]{seqPos}[$j]{rEnd} + $mergedGroup[$gi]{seqPos}[$j]{end} );\n";
					#}
					#exit;
			}
			elsif($gStrand == $R && $j != 0){
				#####
				$rightAdditionalLen += $mergedGroup[$gi]{seqPos}[$j]{rStart} - 
					(!defined $incMatch || $incMatch->{strand} == $F ? 
						($contigLen{$rName}-$mergedGroup[$gi]{seqPos}[$j]{start}) : $mergedGroup[$gi]{seqPos}[$j]{start});
					### considering <--group-- --hContig--> --lContig--> or  <--group-- <--hContig-- --lContig-->
					### hContig : contig covering the lContig, lContig : contig in a scaffold
			}
			
			#if($rightAdditionalLen < 0){
			#	print STDERR "rightAdditionalLen : $rightAdditionalLen, $mergedGroup[$gi]{name}, $scaffold[$sid]{name}, $scaffold[$sid]{contig}[$ci-1]{name}\n";
			#	exit;
			#}
				
			##### if extention of contig is longer than gap...
			if($rightAdditionalLen > $scaffold[$sid]{gap}[$ci-1]{len}){
				$conflictStr .= " ** right side extention($rightAdditionalLen bases) from $scaffold[$sid]{contig}[$ci-1]{name} is longer than ($ci-th gap size:$scaffold[$sid]{gap}[$ci-1]{len}) of the $scaffold[$sid]{name}. $mergedGroup[$gi]{name} was not merged.\n" if(defined $verboseFlag);
				$isStop = 1;
			}
		}

		$incMatch = undef;
		
		if($isStop != 1 && $giAtPrev->{$gi} != $ci) { #### it means the merged contig is not overlapping from the first contig of the current scaffold..		
			$j = $contigGroupInfo->[$ci-$giAtPrev->{$gi}]{$gi}{order}; #### index of the contig in the merged contig group.			
					
			if(!defined $j){
				print STDERR "ERROR at line ".__LINE__.": ci:$ci, cNum:$cNum, gi:$gi, \$giAtPrev->{\$gi}:$giAtPrev->{$gi}, $contigGroupInfo->[$ci-$giAtPrev->{$gi}]{$gi}\n"; 
				exit;
			}
			if(defined $contigGroupInfo->[$ci-$giAtPrev->{$gi}]{$gi}{rName}){ ### if ($ci-1)th contig in the scaffold is included in another contig in a merged group..
				$qName = $scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}]{name};
				$rName = $contigGroupInfo->[$ci-$giAtPrev->{$gi}]{$gi}{rName};
=head
				$includeLow{$qName}{$rName}{lStart}
				$includeLow{$qName}{$rName}{lEnd}
				$includeLow{$qName}{$rName}{hStart}
				$includeLow{$qName}{$rName}{hEnd}
				$includeLow{$qName}{$rName}{strand}
=cut
				$incMatch = $includeLow{$scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}]{name}}{$rName};
				#print "$incMatch, \$includeLow{$scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}]{name}}{$rName}, incMatch->{strand} : $incMatch->{strand}\n";
				if(!defined $mergedGroup[$gi]{strandList}[$j] || !defined $incMatch->{strand}){
					print STDERR "ERROR at line ".__LINE__.": !defined (\$mergedGroup[$gi]{strandList}[$j]) || !defined (\$incMatch->{strand})\n"; 
					exit;
				}
				$gStrand = ($mergedGroup[$gi]{strandList}[$j] == $incMatch->{strand} ? $F : $R);
				if($includeLow{$qName}{$rName}{strand} == $F){
					$leftAdditionalLen = $includeLow{$qName}{$rName}{hStart} -1;
				}
				else{ #### in the blast match, if strand == R, hEnd < hStart
					$leftAdditionalLen = $contigLen{$rName} - $includeLow{$qName}{$rName}{hStart};
				}
			}
			else{
				$rName = $scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}]{name};
				$gStrand = $mergedGroup[$gi]{strandList}[$j];
			}
			if($gStrand == $F && $j != 0){
				$leftAdditionalLen += $mergedGroup[$gi]{seqPos}[$j]{rStart} - 
					(!defined $incMatch || $incMatch->{strand} == $F ? 
						$mergedGroup[$gi]{seqPos}[$j]{start} : ($contigLen{$rName}-$mergedGroup[$gi]{seqPos}[$j]{start}) );
					### considering --group--> --hContig--> --lContig--> or --group--> <--hContig-- --lContig-->
					### hContig : contig covering the lContig, lContig : contig in a scaffold

				#if(!defined $incMatch || $incMatch->{strand} == $F){
				#	print STDERR "leftAdditionalLen : $leftAdditionalLen, $mergedGroup[$gi]{seqPos}[$j]{rStart}-$mergedGroup[$gi]{seqPos}[$j]{start}\n";
				#}
				#else{
				#	print STDERR "leftAdditionalLen : $leftAdditionalLen, $mergedGroup[$gi]{seqPos}[$j]{rStart}-($contigLen{$rName}-$mergedGroup[$gi]{seqPos}[$j]{start})\n";
				#}
				
			}
			elsif($gStrand == $R && $j != $mergedGroup[$gi]{num}-1){
				$leftAdditionalLen += $mergedGroup[$gi]{len} - 
					($mergedGroup[$gi]{seqPos}[$j]{rEnd} + 
					(!defined $incMatch || $incMatch->{strand} == $F ? 
						$mergedGroup[$gi]{seqPos}[$j]{end} : ($contigLen{$rName}-$mergedGroup[$gi]{seqPos}[$j]{end}) ) ); 
					### considering <--group-- --hContig--> --lContig--> or <--group-- <--hContig-- --lContig-->
					### hContig : contig covering the lContig, lContig : contig in a scaffold
					#### if --group--> <--hContig--, {start} > {end}..
				##### $mergedGroup[$gi]{seqPos}[$j]{end} is an end position of the alignment.

				#if(!defined $incMatch || $incMatch->{strand} == $F){
				#	print STDERR "leftAdditionalLen : $leftAdditionalLen, $mergedGroup[$gi]{len}-($mergedGroup[$gi]{seqPos}[$j]{rEnd} + $mergedGroup[$gi]{seqPos}[$j]{end}) )\n";
				#}
				#else{
				#	print STDERR "leftAdditionalLen : $leftAdditionalLen, $mergedGroup[$gi]{len}-($mergedGroup[$gi]{seqPos}[$j]{rEnd} + ($contigLen{$rName}-$mergedGroup[$gi]{seqPos}[$j]{end})) )\n";
				#}
			}
			
			#if($leftAdditionalLen < 0){
			#	print STDERR "leftAdditionalLen : $leftAdditionalLen, $mergedGroup[$gi]{name}, $scaffold[$sid]{name}, $scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}]{name}\n";
			#	exit;
			#}

			#print " ** left side extention($leftAdditionalLen) from $scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}]{name} is not longer than (".($ci-$giAtPrev->{$gi})."-th gap size:$scaffold[$sid]{gap}[$ci--$giAtPrev->{$gi}-1]{len}) of the $scaffold[$sid]{name}\n";
				
			if($leftAdditionalLen > $scaffold[$sid]{gap}[$ci-$giAtPrev->{$gi}-1]{len}){
				$conflictStr .= " ** left side extention($leftAdditionalLen bases) from $scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}]{name} is longer than (".($ci-$giAtPrev->{$gi})."-th gap size:$scaffold[$sid]{gap}[$ci-$giAtPrev->{$gi}-1]{len}) of the $scaffold[$sid]{name}. '$mergedGroup[$gi]{name}' was not merged.\n" if(defined $verboseFlag);
				$isStop = 1;
			}
			elsif($#$overGroup > -1 && $overGroup->[$#$overGroup]{endCid} == $ci-$giAtPrev->{$gi}-1 &&
				$leftAdditionalLen+$overGroup->[$#$overGroup]{rightAdditionalLen} > $scaffold[$sid]{gap}[$ci-$giAtPrev->{$gi}-1]{len}){
				$conflictStr .= " ** sum of right side extention($overGroup->[$#$overGroup]{rightAdditionalLen} bases) from $scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}-1]{name} " .
					"and left side extention($leftAdditionalLen bases) from $scaffold[$sid]{contig}[$ci-$giAtPrev->{$gi}]{name} " .
					"is longer than (".($ci-$giAtPrev->{$gi})."-th gap size:$scaffold[$sid]{gap}[$ci-$giAtPrev->{$gi}-1]{len}) of the $scaffold[$sid]{name}. '$mergedGroup[$gi]{name}' was not merged.\n" if(defined $verboseFlag);
				$isStop = 1;
			}
		}
		if($isStop == 1){ 
			delete $continuousRecord->{$gi};
			next;
		}
		#print " -- rName : $rName, Group $gi passed the test.\n" if(defined $verboseFlag);

	#######
		$j = $#$overGroup+1;
		
		#### if two merged contig groups are overlapping.. choose longer one.
		if($j > 0 && ($overGroup->[$j-1]{endCid} >= $startCid) ){
			if($overGroup->[$j-1]{endCid} - $overGroup->[$j-1]{startCid} < $endCid - $startCid)
			{	
				$conflictStr .= "  ++ $mergedGroup[$overGroup->[$j-1]{gi}]{name} was replaced by $mergedGroup[$gi]{name}.\n";
				$j--;	
			}
			else	
			{	
				$conflictStr .= "  ++ $mergedGroup[$gi] was deleted.\n";
				next;	
			}
		}
		
		$overGroup->[$j]{gi} = $gi;
		$overGroup->[$j]{startCid} = $startCid; #### contig index at sid..			
		$overGroup->[$j]{endCid} = $endCid; #### contig index at sid..
		$overGroup->[$j]{strand} = $gStrand;
		$overGroup->[$j]{leftAdditionalLen} = $leftAdditionalLen;
		$overGroup->[$j]{rightAdditionalLen} = $rightAdditionalLen;


=head
for($i = 0; $i <= $#seqPos; $i++){
	$seqPos[$i]{rStart} = $prevLen + 1;
	$seqPos[$i]{rEnd} = $prevLen + 1 + abs($seqPos[$i]{end} - $seqPos[$i]{start});
}

$mergedGroup[$gi]{name} = $name;
$mergedGroup[$gi]{nodeList} = \@nodeList;
$mergedGroup[$gi]{strandList} = \@strandList;
$mergedGroup[$gi]{seqPos} = \@seqPos;
=cut
		#if($scaffold[$sid]{contig}[$i-1]{name}
		#print "delete \$continuousRecord->{$gi};\n";
		delete $continuousRecord->{$gi};
	}		

	print $conflictOut $logForCurrent_str;#  if(defined $verboseFlag);
}


sub searchNode{
	my ($node, $strand, $direction, $nodeList, $strandList) = @_;  #### $direction : left, right
	
	print "  qName : $node->{qName}, strand : $strandTxt[$strand]\n" if(defined $verboseFlag && $#$nodeList != -1);
	my (@matchList, $qName);
	$qName = $node->{qName};
	
	return if(scalar @$nodeList > 1 && $nodeList->[0] eq $qName); ### to prevent recursive paring. A-B-C-A-B...

	my $dir;
	if( ($strand == $F && $direction == $left) || ($strand == $R && $direction == $right)){
		$dir = "5";
	}
	else{
		$dir = "3";
	}


	my $stop; ### if data is not consistent... stop searching..
	my $rDir; ### end direction(5' or 3') of match.

	$stop = 0;
	###### look at the ($dir)' pairs...
	my @rNames = sort {$node->{$dir}{$b}{qLen} <=> $node->{$dir}{$a}{qLen}} keys %{$node->{$dir}};
	

#	if($qName =~ /NODE_318976_length_5613_cov_8.772849/){ # && $dir == 3
#		print " **** qName : $node->{qName}, dir : $dir, strand : $strandTxt[$strand]\n";
#		foreach $rName (@rNames) {
#			print "  - **** rName : $rName, $node->{$dir}{$rName}{qLen}\n";
#		}		
#	}

	foreach $rName (@rNames) {		

		if(!defined $includeLow{$rName} && !defined $deleteHash{$rName} && defined $endMatches{$rName}){
			#print "end direction $dir, ($qName), $rName, strand :$node->{$dir}{$rName}{rStrand}, $endMatches{$rName}{5}{$qName}, $endMatches{$rName}{3}{$qName}\n"; #exit;
			return if(scalar @$nodeList > 1 && $nodeList->[0] eq $rName); ### to prevent recursive paring. A-B-C-A-B...

			##### skip when qname has pair information with rName but rName does not have it.
			if(	($dir eq "5" && $node->{$dir}{$rName}{rStrand} == $F)||
					($dir eq "3" && $node->{$dir}{$rName}{rStrand} == $R)
			)	{ $rDir = "3"; }
			elsif( ($dir eq "5" && $node->{$dir}{$rName}{rStrand} == $R) ||
						($dir eq "3" && $node->{$dir}{$rName}{rStrand} == $F)
			)	{ $rDir = "5"; }
			next if(!defined $endMatches{$rName}{$rDir}{$qName});

			#### to see both contigs have common pairs...
			foreach my $tmpName (keys %{$endMatches{$rName}{$rDir}}) {
				next if($tmpName eq $qName);
				if(!defined $node->{5}{$tmpName} && !defined $node->{3}{$tmpName}){
					if($node->{$dir}{$rName}{qLen}*$node->{$dir}{$rName}{identity} > 
						$endMatches{$rName}{$rDir}{$tmpName}{qLen}*$endMatches{$rName}{$rDir}{$tmpName}{identity} + 40){
						unlink_endMatch($rName, $rDir, $tmpName);
					}
					else{
						print " *** Breaking pair : $qName($dir') paired with $rName($rDir'). $qName does not have $tmpName, which $rName paired with at the ($rDir)' end\n" if(defined $verboseFlag);
						$stop = 1;
					}					
					#return; ##TEST
				}
			}
			next if($stop == 1); #### if qName node does not have a pair rName has...
			
			print "   $rName, $strandTxt[$node->{$dir}{$rName}{rStrand}]\n" if(defined $verboseFlag);
			$matchList[$#matchList+1] = $node->{$dir}{$rName};
		}
	}

	if(scalar @matchList == 1){ #### if there is only one pair at the end.. search the paired node.		
		push(@{$nodeList}, $matchList[0]->{rName}); 
		push(@{$strandList}, ($matchList[0]->{rStrand} == $strand ? $F : $R));
		if($matchList[0]->{rStrand} == $strand){$FNum++;} else {$RNum++;}
		searchNode($endMatches{$matchList[0]->{rName}}, ($matchList[0]->{rStrand} == $strand ? $F : $R), 
			$direction, $nodeList, $strandList); 
	}
	elsif(scalar @matchList > 1){
		my @sortedList;
		
		if($dir eq "5") { @sortedList= sort{$a->{qEnd} <=> $b->{qEnd}} @matchList; }
		else { @sortedList= sort{$a->{qStart} <=> $b->{qStart}} @matchList; }
		### compare the nearest node with other nodes at 5'end..
		my ($i, $j, $num, $middleNode);
		$num = scalar @sortedList;
		

		my ($middleNode_i, $nodeDir_i, $middleNode_j, $nodeDir_j, $isInconsistant);
		$isInconsistant = 0;
		if(($strand == $F && $direction == $left) || ($strand == $R && $direction == $right)){ ##### this strand is decided by the first contig.. If it is $F, the current contig is same direction with the first contig.

			for($i = $num-1; $i > 0; $i--){
				$middleNode_i = $endMatches{$sortedList[$i]->{rName}};
				$nodeDir_i = ($sortedList[$i]->{rStrand} == $F ? "5" : "3");
				#print "   nodeDir : $nodeDir\n" if(defined $verboseFlag);				
				$j = $i-1;
				$middleNode_j = $endMatches{$sortedList[$j]->{rName}};
				$nodeDir_j = ($sortedList[$j]->{rStrand} == $F ? "3" : "5");
								
				#print " F $node->{qName} : $sortedList[$i]->{rName} , $nodeDir,  $sortedList[$j]->{rName}\n";
				if(!defined $middleNode_i->{$nodeDir_i}{$sortedList[$j]->{rName}} || 
						!defined $middleNode_j->{$nodeDir_j}{$sortedList[$i]->{rName}}){
					print " $node->{qName} : $sortedList[$i]->{rName} is not paired with $sortedList[$j]->{rName}. It will not be extended.\n" if(defined $verboseFlag);
					$isInconsistant = 1;
					last;
				}
			}
			##### all paired contigs should be paired with the $sortedList[$#sortedList] contig.
			if($isInconsistant == 0){				
				for($i = $num - 1; $i >= 0; $i--){
					push(@{$nodeList}, $sortedList[$i]->{rName}); 
					push(@{$strandList}, ($sortedList[$i]->{rStrand} == $strand ? $F : $R)); 
					if($matchList[$i]->{rStrand} == $strand){$FNum++;} else {$RNum++;}
				}
				
				searchNode($endMatches{$sortedList[0]->{rName}}, ($sortedList[0]->{rStrand} == $strand ? $F : $R), 
					$direction, $nodeList, $strandList); 
			}
		}
		else { ### if($strand == $R) 			
			for($i = $num-1; $i > 0; $i--){
				$middleNode_i = $endMatches{$sortedList[$i]->{rName}};
				$nodeDir_i = ($sortedList[$i]->{rStrand} == $F ? "5" : "3");
				$j = $i-1;
				$middleNode_j = $endMatches{$sortedList[$j]->{rName}};
				$nodeDir_j = ($sortedList[$j]->{rStrand} == $F ? "3" : "5");				

				#print " R $node->{qName} : $sortedList[$i]->{rName} , $nodeDir,  $sortedList[$j]->{rName}\n";
				if(!defined $middleNode_i->{$nodeDir_i}{$sortedList[$j]->{rName}} ||
					!defined $middleNode_j->{$nodeDir_j}{$sortedList[$i]->{rName}}){
					print " $node->{qName} : $sortedList[$i]->{rName} is not paired with $sortedList[$j]->{rName}. It will not be extended.\n" if(defined $verboseFlag);
					$isInconsistant = 1;
					last;
				}
			}

			if($isInconsistant == 0){
				for($i = 0; $i < $num; $i++){
					push(@{$nodeList}, $sortedList[$i]->{rName}); 
					push(@{$strandList}, ($sortedList[$i]->{rStrand} == $strand ? $F : $R)); 
					if($matchList[$i]->{rStrand} == $strand){$FNum++;} else {$RNum++;}
				}
				searchNode($endMatches{$sortedList[$#sortedList]->{rName}}, ($sortedList[$#sortedList]->{rStrand} == $strand ? $F : $R), 
					$direction, $nodeList, $strandList); 
			}
		}
	}
}

sub unlink_endMatch
{
	my ($qName, $qDir, $rName) = @_;

	my $rDir = (
		($qDir == 5 && $endMatches{$qName}{$qDir}{$rName}{rStrand} == $F) ||
		($qDir == 3 && $endMatches{$qName}{$qDir}{$rName}{rStrand} == $R)
		? 3 : 5);

	if(defined $endMatches{$rName}{$rDir}{$qName}){
		delete $endMatches{$rName}{$rDir}{$qName};
		if(keys %{$endMatches{$rName}{$rDir}} == 0){
			delete $endMatches{$rName}{$rDir};
			my $revDir = ($rDir == 5 ? 3 : 5); 
			#delete $endMatches{$rName} if(!defined $endMatches{$rName}{$revDir});
		}
	}

	if(defined $endMatches{$qName}{$qDir}{$rName}){
		delete $endMatches{$qName}{$qDir}{$rName};
		if(keys %{$endMatches{$qName}{$qDir}} == 0){
			delete $endMatches{$qName}{$qDir};
			my $revDir = ($qDir == 5 ? 3 : 5); 
			#delete $endMatches{$qName} if(!defined $endMatches{$qName}{$revDir});
		}
	}
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


