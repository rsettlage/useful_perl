#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 16 Oct 2013
#
#remove lines from file 1 found in file 2
#
#usage perl remove_lines.pl <file1> <file2>
#example perl /groups/DAC/useful_perl/remove_lines.pl rh.ctg.fasta.bln.txt test.bln.txt
#
#####################################################

$file1_NAME = $ARGV[0];
$file2_NAME = $ARGV[1];

open (target_FILE, "$file1_NAME")|| die "Couldn't open $file1_NAME\n";
print "working on file1 $file1_NAME\n";
open (exclusion_FILE, "$file2_NAME")|| die "Couldn't open $file2_NAME\n";
print "getting include/exclude list from $file2_NAME\n";

my $file_OUT = $file1_NAME . "_reduced";
open (out_FILE, ">$file_OUT") || die "Couldn't open $file_OUT\n";
print out_FILE "creating $file_OUT for results\n";

print "starting file scan\n";
my %target_hash;
my $temp = <exclusion_FILE>;
while ($temp) {
	chomp $temp;
	@tmp_fields = split('\t', $temp);
	$target_hash{$tmp_fields[0]}=1;
	$temp = <exclusion_FILE>;
}

$target_line = <target_FILE>;
while ($target_line) {
	chomp($target_line);
	@tmp_target_fields = split('\t', $target_line);
	printf out_FILE "$target_line\n" if exists $target_hash{$tmp_target_fields[0]};
	$target_line = <target_FILE>;
}	

close target_FILE;
close exclusion_FILE;
close out_FILE;
exit;