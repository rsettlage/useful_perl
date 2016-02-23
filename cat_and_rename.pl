#!/usr/bin/perl -w
use strict;

#open sample sheet

open (SS, $ARGV[0]) || die "couldn't open $ARGV[0]\n";

my $firstline=<SS>;
my $line=<SS>;
while($line){

    chomp $line;
    $line=~s/\r\n//g;
    my @data=split(',', $line);
    my $project=$data[9];
    my $sampleid=$data[2];
    my $userid=$data[5];
    chop $project;
    my $lastid=".\/Project_" . $project . "\/Sample_$sampleid\/*.fastq.gz";
    system("cat $lastid > .\/samples\/$project" . "_" . "$userid.fastq.gz");
#    print "lastid is $lastid\n";
    $line=<SS>;
}
close SS;
    
#put together the address for each sample using sample sheet information
#make samples/project directories
#cat all into each samples/project directory



