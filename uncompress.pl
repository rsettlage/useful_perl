#!/usr/bin/perl
use strict;
use warnings;

###############################################################################
###      @author: Robert E. Settlage, Ph.D.                                                                                                  
###
###        Data Analysis Core @ the Virginia Bioinformatics Institute          
###   
###        August 29, 2012                                                                                                                            
###                                          
###        uncompress.pl   
###                        
###        this script takes a compressed file and uncompresses it, it returns the basefile name
### 
###        use:    perl  uncompress.pl file.gz                                                                                                                        
###
###############################################################################


$file=$ARGV[0];

my $file_EXT=$file;

####check for gzip
$file_EXT=~ s/\.gz$// ;
if ($file_EXT ne $file) {
	`gzip -d $file` ;
	$exit_value = $file_EXT;
}

####check for bzip2
$file_EXT=~ s/\.bz2$// ;
if ($file_EXT ne $file) {
	`bunzip2 -v $file` ;
	$exit_value = $file_EXT;
}



exit $exit_value;