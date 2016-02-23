use LWP::Simple;
use strict;

my $startUrl = "http://brcdownloads.vbi.vt.edu/patric2/genomes.May2012/";
my $usage = "$0 target_genus (or any prefix for directories on $startUrl\n";
my $target = shift or die $usage;
my $minPrefixLength = 5;
die "target genus (or prefix) must be at least $minPrefixLength long\n" unless length($target) >= $minPrefixLength;

my $list = get($startUrl);
while ($list =~ /href="([^">]*)">/g)
{
    my $dir = $1;
    if ($dir =~ /$target/)
    {
	$dir =~ tr/\///d;
	print "$dir\n";
	my $subdir = get($startUrl . $dir);
	my $found = 0;
	my @files;
	while ($subdir =~ /$dir\.PATRIC.faa/gs)
	{
	   print $&, "\n";
	   push @files, $&;
	   $found++;
	}
	unless ($found)
	{
	    while ($subdir =~ /$dir\.RefSeq.faa/gs)
	    {
		print "  $&\n";
		push @files, $&;
		$found++;
	    }
	}
	unless ($found)
	{
	    while ($subdir =~ /$dir\.\S+faa/gs)
	    {
		print "      $&\n";
		push @files, $&;
		$found++;
	    }
	}

	my %downloaded;
	for my $file (@files)
	{
	    next if $downloaded{$file};
	    $downloaded{$file}=1;
	    #open F, ">$file";
	    print F getstore($startUrl . $dir . '/' . $file, $file);
	    #close F;
	}
    }
}
