use LWP::Simple;
use strict;

my $startUrl = "http://brcdownloads.vbi.vt.edu/patric2/genomes.May2012/";
my $list = get($startUrl);

my $usage = "$0 target_genus (or any prefix for directories on $startUrl\n";
my $target = shift or die $usage;
my $minPrefixLength = 5;
die "target genus (or prefix) must be at least $minPrefixLength long\n" unless length($target) >= $minPrefixLength;

$target;

while ($list =~ /href="([^">]*)">/g)
{
    my $dir = $1;
    if ($dir =~ /^$target\_/)
    {
	$dir =~ tr/\///d;
	print "$dir\n";
	my $subdir = get($startUrl . $dir);
	my $found = 0;
	my @files;
	while ($subdir =~ /$dir\.fna/gs)
	{
	   print $&, "\n";
	   push @files, $&;
	   $found++;
	}
	my %downloaded;
	for my $file (@files)
	{
	    next if $downloaded{$file};
	    next if -e $file; # do not re-download or clobber existing genomes of same name
	    $downloaded{$file}=1;
	    print STDERR "Getting $file\n";

	    my $status = getstore($startUrl . $dir . '/' . $file, $file);
	    warn "Error $status on $file" unless is_success($status);
	    #open F, ">$file";
	    #print F get($startUrl . $dir . '/' . $file);
	    #close F;
	}
    }
}
