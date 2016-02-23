#!/usr/bin/perl -w

#####################################################
# written by Bob Settlage, DAC, 23 March 2012
# fasta_fastq_read_compare
#
#given a gtf, fasta or fastq file, program creates a .common and .file1diff and .file2diff file,
#makes no assumptions on which (fasta or fastq) the input file is, but need to be the same
#
#usage perl fasta_fastq_read_compare.pl <gtf|fasta|fastq> <gtf|fasta|fastq>
#example perl fasta_fastq_read_compare.pl ./merged_asm_1/merged_test_1.gtf ./merged_asm_1/merged_test_2.gtf
#
#####################################################

$file_ext1 = ($ARGV[0] =~ m/([^.]+)$/)[0];
$file_ext2 = ($ARGV[1] =~ m/([^.]+)$/)[0];

if ($file_ext1 ne $file_ext2) {
	print "files are not the same type, please change the file extensions to match or choose another file pair\n";
	exit;
}

open (FILE_1, $ARGV[0])|| die "Couldn't open $ARGV[0]\n";
open (FILE_2, $ARGV[1])|| die "Couldn't open $ARGV[0]\n";

my $file_common = $ARGV[0].".common";
my $file_1_diff = $ARGV[0].".diff";
my $file_2_diff = $ARGV[1].".diff";
open (FILE_common, ">$file_common") || die "Couldn't open $file_common\n";
open (FILE_diff_1, ">$file_1_diff") || die "Couldn't open $file_1_diff\n";
open (FILE_diff_2, ">$file_2_diff") || die "Couldn't open $file_2_diff\n";

my $line_1=<FILE_1>;
my $line_2=<FILE_2>;
$read_FLAG = 1;

$line_num = 1;

if ($file_ext1 eq "gtf") {
	while($line_1 && $line_2) {
		if($line_1 eq $line_2){
			printf FILE_common $line_num." ".$line_1;
		}else{
			printf FILE_diff_1 $line_num." ".$line_1;
			printf FILE_diff_2 $line_num." ".$line_2;
		}
		$line_1=<FILE_1>;
		$line_2=<FILE_2>;
		$line_num ++;
	}
}

close FILE_1;
close FILE_2;
close FILE_common;
close FILE_diff_1;
close FILE_diff_2;
exit;