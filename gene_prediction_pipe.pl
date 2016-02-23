#!/usr/bin/perl -w
use strict;

my $contigfile=$ARGV[0];
my $similarfile=$ARGV[1];
my $cdnafile=$ARGV[2];
chomp ($contigfile, $similarfile);

print "Begining blat\n";
system("blat $contigfile $similarfile blatresults2.psl -minIdentity=80");

print "Beginning conversion to gff\n";
system("perl blat2gff.pl blatresults2.psl > blatresults2.gff");

print "Beginning parse of cDNA\n";
system("perl parse_trancript_from_contigfile.pl blatresults2.gff $contigfile $cdnafile");
