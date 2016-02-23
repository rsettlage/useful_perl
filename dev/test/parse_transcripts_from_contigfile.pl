#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $gfffile=$ARGV[0];
my $contigfile=$ARGV[1];
my $cdnafile=$ARGV[2];

open (GFF, $gfffile)|| die "Couldn't open $gfffile\n";
open (CONTIG, $contigfile) || die "Couldn't open $contigfile\n";
open (OUT, ">$cdnafile") || die "couldn't open $cdnafile\n"; 

my %contighash;
my $cline=<CONTIG>;
my $linenum=0;

#First, let's get the coordinates for all the ocntigs
while ($cline){
    if($cline=~/>/){
	my $name=$cline;
	chomp $name;
	$name=~s/>//g;
	$name=~s/\s+//g;
	if(!exists($contighash{$name})){
	    $contighash{$name}=$linenum;
	    $cline=<CONTIG>;
	    $linenum++;
	}else{
	    print "I've seen $name before\n";
	    $cline=<CONTIG>;
	    $linenum++;
	}
    }else{
	$cline=<CONTIG>;
	$linenum++;
    }
}

#Now we can go into the GFF file, and combine the cDNA sequences into one long cDNA sequence.

my %cdnahash;

my ($contigname, $strand, $name, $other, @arr, @start, @end);


while(<GFF>){
    if($_!~/gff/){
	@arr=split(/\t/);
	my ($tempname3, @other)=split(/\s+/, $arr[8]);
	my ($other, $newname)=split(/\;/, $tempname3);
	$newname .= "_" . $arr[0];
	#If we're at a new hit

	if($name ne $newname){
	    #Get the contig for the current hit
	    my $contig="";
	    if(exists($contighash{$contigname})){
		my $linenumber=$contighash{$contigname};
		my $contigline=0;
		open(CON, $contigfile)|| die "Couldn't open $contigfile\n";
		my $line=<CON>;
		while($line){
		    if($contigline eq $linenumber+1){
			$contig=$line;
			last;
		    }else{
			$line=<CON>;
			$contigline++;
		    }
		}
		
	    }

	    my @exonseqs;
	    
	    for(my $j=0; $j<scalar(@start); $j++){
		my $exonseq=substr($contig, $start[$j]-1, $end[$j]-$start[$j]+1);
		chomp $exonseq;
		if($strand=~/\-/){
		    $exonseq=revcomp($exonseq);
		    unshift(@exonseqs, $exonseq);
		}else{
		    push(@exonseqs, $exonseq);
		}
	    }

	    my $cdna=join("", @exonseqs);
	    print OUT ">" . "$name\n$cdna\n";
	    $cdnahash{$name}=$cdna;
	    @start=();
	    @end=();
	    $name=$newname;
	    $contigname=$arr[0];
	    $strand=$arr[6];
	}


	push(@start, $arr[3]);
	push(@end, $arr[4]);
	
	
    }
}

	    



sub revcomp {

#Got this from http://www.perlmonks.org/?node_id=197793
    my $dna = shift; 
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}
