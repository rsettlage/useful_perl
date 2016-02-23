use strict;

open F, $ENV{'file_of_fasta_files'} or die "cannot find $ENV{'file_of_fasta_files'}";
my @files;
while (<F>)
{
	chomp;
	push @files, $_;
}
my $identity = $ENV{'otu_id'};
$identity = 0.97 unless $identity; #default

my $usearch=$ENV{'usearch'};

my $number_of_jobs = $ENV{'number_of_usearch_jobs'};
my $filesPerJob = int(@files/$number_of_jobs);
my $job = 1;
my $fasta_directory = $ENV{fasta_directory};
while (@files)
{
	open F, ">usearch_fasta_vs_OTUs_job$job.sh";
	for my $line (1..$filesPerJob)
	{
		my $file = shift @files;
		print F "$usearch -usearch_global $fasta_directory/$file -db otu_and_notu.fa -id $identity -strand plus -uc $file.vs.otu_and_notu.uc\n"; 
	}	
	close F;
	$job++;
}
