#!/usr/bin/perl
use strict;
use warnings;

###############################################################################
###      @author: Robert E. Settlage, Ph.D.                                                                                                  ###
###        Data Analysis Core @ the Virginia Bioinformatics Institute                                                  ###   
###        November 2011                                                                                                                            ###                                          
###                                                                                                                                                               ###                        ###                                                                                                                                                               ### 
###        use: perl /groups/DAC/useful_perl/blast_count_reformat.pl <count_file>                                                                                                                              ###
###                                                                                                                                                                  ###
###expecting 3 columns of data in csv file: sample, reference, count
###reformatting to reference with columns being sample and count
###############################################################################


$file_NAME=ARGV[0];
$new_NAME=$file_NAME."reformatted.txt";

open (original_FILE, "$file_NAME")|| die "Couldn't open $file_NAME\n";
print "working on $file_NAME\n";
open (new_FILE, ">$new_NAME")|| die "Couldn't open $new_NAME\n";
print "writting to $new_NAME\n";


my %samples;
my %references;
my %expression;
my $rows=0;
my $temp = <original_FILE>;
while ($temp) {
	chomp $temp;
	@tmp_fields = split(',', $temp); 
	$samples{$tmp_fields[0]}=$tmp_fields[0];
	$references{$tmp_fields[1]}=$tmp_fields[1];
	$expression[$rows]=$temp;
	$temp = <original_FILE>;
}


