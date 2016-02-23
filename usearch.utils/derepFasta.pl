use strict;

my $usage = "$0 [file_of_fasta_files default: $ENV{'file_of_fasta_files'}]
output is set of repeated sequences in derep.fa
";

my $minLength = 10;
my $minCount = 2;
my %seq;
my $seq;
my $id;
my $fileOfFastaFiles;# = shift;
unless ($fileOfFastaFiles) 
{
	$fileOfFastaFiles = $ENV{'file_of_fasta_files'};
}
print "ENV{'file_of_fasta_files'} = $ENV{'file_of_fasta_files'}\n";
die "Cannot find file of fasta files" unless $fileOfFastaFiles;
print "File of file names = $fileOfFastaFiles\n";

my $fasta_directory = $ENV{'fasta_directory'};
print "fasta_directory = $fasta_directory\n";

open FOFF, $fileOfFastaFiles or die "cannot open $fileOfFastaFiles";
my %readsPerFile;
while (<FOFF>)
{	
	my $file = $_;
	chomp $file;
	print STDERR "open $fasta_directory/$file\n";
	open F, "$fasta_directory/$file" or die "cannot open $file";
	while (<F>)
	{
		if (/^>(\S*)/)
		{
			if ($seq)
			{
				$readsPerFile{$file}++;
				warn "short seq: $seq" unless (length($seq) >= $minLength);
				push @{$seq{$seq}}, $id;
			}	
			$id = $1;
			$seq = '';
		}
		else
		{
			chomp;
			tr/ //d;
			$seq .= $_;
		}
	}#while <F>
	#process last sequence of file
	push @{$seq{$seq}}, $id;
}
open IDENT_SEQ_IDS, ">identical_seq_ids.txt";
open DEREP, ">derep.fa";
#open SINGLE, ">singleton_seqs.fa";
foreach my $seq (sort {scalar(@{$seq{$b}}) <=> scalar(@{$seq{$a}})} keys %seq)
{#sort by num reads
	my $id = $seq{$seq}->[0];
	my $cnt = scalar(@{$seq{$seq}});
	if ($cnt >= $minCount)
	{
		print DEREP ">$id;size=$cnt;\n$seq\n";
		print IDENT_SEQ_IDS join(" ", @{$seq{$seq}}), "\n";
	}
	else
	{
	#		print SINGLE ">$id;size=$cnt;\n$seq\n";
	}
}
open F, ">raw_reads_per_file.txt";
foreach my $file (sort keys %readsPerFile)
{
	print F $file, "\t", $readsPerFile{$file}, "\n";
}
