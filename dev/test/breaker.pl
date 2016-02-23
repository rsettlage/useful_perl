#!/usr/bin/perl

#AUTHOR
# Nagarjun Vijay (c) 2011
# mate maker.pl - program to break scaffolds into mate pairs
#

#DATE
# 3rd March 2011
# 

use strict;

my $scafile=$ARGV[0];

#reading the scaffold sequences
open SCAF, "<",$scafile or die $!;

#open contigs file
open CONTIG, ">>",$scafile . ".contigs.fa" or die $!;


#open contigs file
open LIB, ">>",$scafile . ".lib" or die $!;

my ($scaffoldcount,$seqst_temp,$contigcount,$notscaffoldcount,$z,$head,$i,$j,$libcount);
$scaffoldcount=0,$seqst_temp="",$contigcount=0,$notscaffoldcount=0,$libcount=0;
my $header = ;#reading the first header

while($z = ){
if($z=~ />/)
{
#split scaffold into contigs
if($seqst_temp=~ m/[N]+/){#is really a scaffold
$scaffoldcount++;
$contigcount +=SplitScaffold($seqst_temp);
}#end of scaffold check if
else{#not a scaffold - so not writing into breaker file
$head=">unbreak";
$notscaffoldcount++;
print CONTIG $head . "" . $notscaffoldcount ."\n";
print CONTIG "$seqst_temp\n";
}#end of scaffold check else

#next sequence
$seqst_temp="";
}
else
{
$z =~ s/[\n\t\f\r_0-9\s]//g;
$seqst_temp .= $z;
}
}#end of File while loop

if($seqst_temp=~ m/[N]+/){#is really a scaffold
$scaffoldcount++;
$contigcount +=SplitScaffold($seqst_temp);
}#end of scaffold check if
else{#not a scaffold - so not writing into breaker file
$head=">unbreak";
$notscaffoldcount++;
print CONTIG $head . "" . $notscaffoldcount ."\n";
print CONTIG "$seqst_temp\n";
}#end of scaffold check else
close SCAF;
close CONTIG;

#Split scaffolds, write contigs and breaker file
sub SplitScaffold{
my $seq = shift;
my ($head,$num,$contigs,$n);
$head=">breaker$scaffoldcount";
$num=1,$contigs=1;
my @values=split(/([N]+)/,$seq);
my $vallen= scalar (@values);
for($i=0;$i<$vallen;$i++) {
my $val=$values[$i];
if(($num % 2) == 0){
$n = length $val;
}
else{
print CONTIG $head . "" . $contigs ."\n";
print CONTIG "$val\n";
$contigs++;
my $matedist=(length $values[$i])+(length $values[$i+1])+(length $values[$i+2]);
for($j=$i+2;$j<$vallen;$j=$j+2){
print LIB "lib".$libcount." l".$libcount."s".$scaffoldcount."b".$i.".fasta l".$libcount."s".$scaffoldcount."b".$j.".fasta ".$matedist." 0.75 0\n";
$matedist=$matedist+(length $values[$j+1])+(length $values[$j+2]);
open READ1, ">>","l".$libcount."s".$scaffoldcount."b".$i.".fasta" or die $!;
open READ2, ">>","l".$libcount."s".$scaffoldcount."b".$j.".fasta" or die $!;
print READ1 ">l".$libcount."s".$scaffoldcount."b".$i."\n";
print READ1 $values[$i]."\n";
print READ2 ">l".$libcount."s".$scaffoldcount."b".$j."\n";
print READ2 reverseComplement($values[$j])."\n";
close READ1;
close READ2;
$libcount++;
}

}
$num++;
}

return $contigs;
}
#reverse complement a sequence
sub reverseComplement{
$_ = shift;
tr/ATGC/TACG/;
return (reverse());
}
#end of Breaker 
