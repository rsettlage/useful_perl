#!/usr/bin/perl

#####################################################
# debugging stuff
#
#####################################################

$fasta_NAME = $ARGV[0];
$gtf_NAME = $ARGV[1];

open (fasta_FILE, "$fasta_NAME")|| die "Couldn't open $fasta_NAME\n";
print "working on fasta_file $fasta_NAME\n";
open (gtf_FILE, "$gtf_NAME")|| die "Couldn't open $gtf_NAME\n";
print "working on gtf_file $gtf_NAME\n";

my $file_OUT = $fasta_NAME . ".exons.fasta";
open (out_FILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";
print out_FILE "creating $file_OUT for results\n";
print "creating $file_OUT for results\n";

exit;