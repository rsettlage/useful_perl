#!/usr/bin/perl -w

###Bob Settlage DAC Jan 2012
### perl script to create genome file for BEDtools
### need to take fasta file as input and output
### tab-deliminited file containing sequence name -tab- sequence length
### perl fasta_to_oneline_bed_summary.pl rh_scf.fa

my $file=$ARGV[0];
my $out=$file . ".summary";

open(IN, $file);
open(OUT, ">$out")|| die "Couldn't open $out\n";

my $line=<IN>;

while($line){

    if($line=~/>/){

		my $name=$line;
		chomp $name;
		$name=~s/>//g;
		my $seq="";
		my $seq_len=0;
		$line=<IN>;

		while(($line) && ($line!~/>/)){
			chomp $line;
			$line=~s/\r//g;
			$seq.=$line;
			$line=<IN>;
			}
		$seq_len=length($seq);
		$delim="\t";
		$outlist=$name.$delim.$seq_len;
		print OUT "$outlist\n";
		#print "$outlist\n";
	#	$line=<IN>;
    }else{
	$line=<IN>;
    }
}

exit;