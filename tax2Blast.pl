#!/usr/bin/env perl

#Author: Bela Tiwari, 2005,and 2008, with additions provided by Frank Tanoh, 2006
#Note the VERY USEFUL facility for building and testing URLs at 
#http://www.ncbi.nlm.nih.gov/Class/wheeler/eutils/eu.cgi

##Note the NASTY hack at line 151

#Warning: any file with the same filename as the one about to be created is removed...this is 
#because of the append function being used to generate the file. 


##Addtions by Frank Tanoh July 2006 
# - ability to download sequences according to a taxon id instead of a taxon name.
# - ability to download sequences using "genus species" rather than just "genus"
# - made command line arguments more intuitive
# - added testing facility through command line flag so testing files downloaded into "testing" directories
# - added some error checking 

#Added to script, 2008 by B. Tiwari - ability to read in a file of taxids.
#Script assumes that file is a list of taxids - one per line in the file. The taxid MUST be in the first column of the file.
#It doesn't matter what else is in the file as long as the taxid is in the first line.


=head1 Name

tax2Blast.pl

=head1 Synopsis

This script is used to download sequences from the NCBI for particular genus or species. By default, the script will download all nucleotide sequences to a file and all peptide sequences to another file. If the -makeblast flag is provided on the command line, blast databases of these sequences will be made.

The names of the resulting fasta files and blast files are reported to screen when the script runs. 

=head2 Warning!

This script is really meant for organisms without vast quantities of sequence data in the NCBI databases. For example, please do NOT run this script to fetch all the mouse sequences. There are much better ways to fetch large volumes of data!


=head1 Usage

There are three modes this script can be used in. All the example commands below include the -makeblast flag
that will create a blast database as well as downloading sequences in fasta format. Omit this flag if you do not wish to create blast databases. 

1) download sequences from a given taxon (genus or species) and make a blast database of them by providing the genus or "genus species" name on the command line. 

tax2Blast.pl -taxon genusName  -makeblast

tax2Blast.pl -taxon "genus species"  -makeblast

2) download sequences and make a blast database of them using a particular taxonomic id from the NCBI. This is a numeric code which you can find in a number of ways. One is by visiting the NCBI Taxonomy databasee website and looking for the organsim or genus you are interested in:  http://www.ncbi.nlm.nih.gov/sites/entrez?db=Taxonomy

tax2Blast.pl  -taxid  taxidNumber  -makeblast

3) download sequences using a list of taxonomic ids from the NCBI. The script expects a tab or space delimited text file as input, where a list of taxonomic ids is given in the first column. 

tax2Blast.pl  -taxfile  inputFileName -makeblast


=head1 Assumptions

This script was written specifically for users of NEBC Bio-Linux (http://nebc.nox.ac.uk/biolinux.html), and thus there are assumptions made about default file locations and what programs area already available on the system. 

Note that the script assumes you have blastall and formatdb installed and on your path. 
The default location to download fasta files to is /home/db/downloads. 
The default location to build blast databases is /home/db/blastdb, which is also the default location for the BLASTDB environmental variable on Bio-Linux systems. 

These settings can be easily changed by editing $dbloc and $blastdbloc near the beginning of the script.

=head1 Additional usage information

You can edit the file to download only nucleotide or only peptide databases by editing the line 
my @dbtypes = ("nucleotide", "protein"); in the script.

There is an additional flag  -testing  that you can set on the command line. Within the script you can set a different location to download files to and build blast databases when this flag is used. 


=head1 Authors

Bela Tiwari and Frank Tanoh, NEBC

=cut



#Things that could be done - e.g. get uids from last weeks file, compare to a list of uids from esearch
#for this week - just get those new ids and append to db file. Rebuild blast databases and alert sysadmin that 
#databases have been updated and which ids have been added

#The list of ids can be used to look for links to other databases and a similar plan for updating and alerting the
#sysadmin could be done.
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=Taxonomy&id=6192&report=brief&mode=text

use warnings;
use strict;
use LWP::Simple;
use Cwd;
use Getopt::Long; 
use File::Basename;
use XML::Simple;

if(scalar @ARGV < 2)
{
        die "\nUsage: $0 -taxon <taxonName> -makeblast\nOR\n$0 -taxid <taxid> -makeblast\nOR\n$0 -taxfile <filename> -makeblast
	\nCommand line options are:
	-taxon  <taxonName>
	-taxid <taxID>  give a taxon id OR a taxon name, not both
	-taxfile <filename> give the name of a file containing a list of the taxids (not names) that you wish to collect sequence for
	-makeblast   (optional) if this flag is present, blast databases are made
	-testing   used for testing only - sets up alternative database locations
	\n\n\n";
}

###############Editable bits#################
#You might want to changes these to be command line options
my $dbloc = "/home/db/downloads";
my $testloc = "/home/manager/wagstaff/sept2008/taxtest";
my $blastdbloc = "/home/db/blastdb";
my $testblastdb = "/home/manager/wagstaff/sept2008/taxtest";

my @dbtypes = ("nucleotide", "protein");    #note that these are picked up as exact string matches in the mkblast function
my $format = "fasta";
my $file_ext = "_seq.fasta";


#############################################


my $base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
my $uid_count = "Null"; 
my $Webenv;
my $querkey;
#The max number of records you can get with a single URL is 500
my $retmax = 500;
my $tdatabase = "Taxonomy";
my $command = "&term=";
my $blast_dbname;

#Step 1 - Get the taxon id from the taxonomy database using the taxon name to search with
# declare the variable to use with GetOption 
my $taxon;
my $makeblast;
my $taxid;
my $taxfile;
my $testing;
my $taxon_change;  
my $gettaxid=0; 
my $gettaxon=0; 
my $gettaxfile=0;
my $makeblastid=0; 
my $blastdatadir;
my $continue = 0;
my $efetch_url;

#argument with the GetOption 
GetOptions ( "taxon:s" =>\$taxon, 
            "taxid:i" =>\$taxid, 
            "taxfile:s" =>\$taxfile,
	    "makeblast!" =>\$makeblast, 
	    "testing!" => \$testing
	   );
	   				#taxon, taxid      #taxonid,taxid_id

if ($taxid && $taxon) 
   {
     $gettaxid=1; 
     print "\n\nYou gave both a taxon id and a taxon name. We will search for sequeces associated with the taxon name $taxon only\n\n"; 
   } 
   elsif($taxon)
   { 
     $gettaxid=1;
   }
   
   elsif($taxfile)
   {
     $gettaxfile=1;
   }

elsif($taxid)
   { 
      $gettaxon=1; 
   }
else
   {
    print "\nThis script requires a valid taxon id or a taxon name to be given on the command line. Please try again\n\n";
    exit;
   } 
 
 
   
####################### get the taxid when TAXON was supplied by user
if ($gettaxid==1)
   {
     ($taxid, $continue) = getTaxidFromTaxon($taxon);
     

      if ($continue == 0)   { print "\n\nNo taxid found for taxon $taxon.\n\n"; exit; }     #exit if you found no taxon for a single taxid

  }
     
############### get the taxon when TAXID was supplied by user
if ($gettaxon ==1)
{
      $taxon = getTaxonFromTaxid($taxid);
      my @taxoncheckarray = split(/ /, $taxon);   #this array should have only two things
       #(species, genus) associated with it if it is a valid taxon.

      if (scalar(@taxoncheckarray) < 3)
     {
     	#print "Your taxon name is \t $taxon ...\n";
	      $continue = 1;
     }
     else
     {
     	print "\nThe taxid $taxid does not appear to be associated with a taxon in the NCBI taxonomy\n. Please try again.\n\n";
	    $continue = 0;
	    exit;
     }
}
############### get a list of the taxons when a taxidfile was supplied by user  - report how many taxa could be aligned
#with tax id.

my @taxidCollection = ();
my @taxonCollection = ();   #ideally we would do this with a hash, but it's a bit of an addin, so array it is.
if ($gettaxfile == 1)
{
   #open file, do a while, push results on array. Maybe we should think of writing a summary to a log file
   
   open(TAXFILE, $taxfile) or die "Couldn't open file $taxfile. Exiting. $!\n";
   
   while ( <TAXFILE> )
      {
      #take first column. Split on any white space in case they didn't use tabs
      my @taxidLineBits = split(/\s/);
      $taxid = $taxidLineBits[0];
      unless ($taxid =~ m/\d*/ ) {warn "The first column in $taxid does not look like a taxid....ignoring this line in $taxfile\n"; next;}
      #$taxid = $_;
      #chomp($taxid);
      $taxon = getTaxonFromTaxid($taxid);  #how is this working? $_ seems to be fine, but not $taxid....
      my @taxoncheckarray = split(/ /, $taxon);   #this array should have only two things
       #(species, genus) associated with it if it is a valid taxon.

      if (scalar(@taxoncheckarray) < 3)
     {
	      $continue = 1;
	      push (@taxidCollection, $taxid);
	      push (@taxonCollection, $taxon);
	      print "Will collect sequences from \t $taxon - taxid $taxid.\n";

     }
     else
     {
     	print "\nThe taxid $taxid does not appear to be associated with a taxon in the NCBI taxonomy. Ignoring this line of $taxfile.\n";
     }
   }

  close TAXFILE;
  my $taxCollectionSize = @taxidCollection;
  if ( $taxCollectionSize == 0 ) 
   {
      print "No valid taxids in $taxfile. Quitting program now.\n"; 
      exit;
    }

  }
   
 ####### use different directory if testing is yes .....

if ($continue) 
{  
 if($testing )
   {
     $dbloc = $testloc;
     $blastdatadir = $testblastdb;
   }
   
  else
   {
      $dbloc = $dbloc;
      $blastdatadir = $blastdbloc;
   }
   
   createdirs($dbloc, $blastdatadir);
     
   #######here we want to make a reasonable sequence and database filenames
	#e.g. if no spaces, take name as is. If spaces, then remove - take first letter of first word
	#all letters of last word 
	############################
########################################### i hav to check this part tomorrow .....	
##############################################################################################

if ($gettaxid == 1  or $gettaxon == 1)
{
      if ($taxon =~/\w+\s+\w+/)
         {
	   $blast_dbname = $taxon ;
	   $blast_dbname =~s/(^\w)(\w+)\s+(\w+)/$1_$3/ ;  #remove the space in the taxon name 
	 }
     else 
         {
	   $blast_dbname = $taxon;
	 }
}

else            #presumably gettaxfile == 1 is the only option left!
{
  #take base name of file as the database name - files created all have extensions so original file will not be overwritten.
  my ($file1, $longDir1, $ext1) = fileparse($taxfile, qr/\..*/);
  $blast_dbname = $file1;

}

###################################################################	


#Now get the nuc and prot seqs for organism and create fasta files

foreach my $database (@dbtypes)
{

  my $seq_filename = $blast_dbname . "_" . $database . $file_ext;
 if ( -e "$dbloc/$seq_filename") { unlink("$dbloc/$seq_filename"); }    #remove old files of same name

  if ($gettaxid == 1  or $gettaxon == 1)
    {
      $uid_count = getSeqsFromNCBI($taxid, $base_url, $database, $dbloc, $seq_filename);
     	#print "Using this URL:\n$efetch_url\n to fetch some of your sequences\n\n";
       if ( $uid_count == 0 ) { print "\nThere are no $database sequences for taxid:$taxid\ttaxon:$taxon\n"; }
	else { print "\nThere are $uid_count $database sequences for taxid:$taxid\ttaxon:$taxon\n"; }
       
    }

  else            #presumably gettaxfile == 1 is the only option left!
  {
    #run a loop so that additional sequences are added to the sequence file for each taxa
    my $count = 0;
    my $sum_uid_count=0;
    my @zeroCountTaxa;
    foreach my $taxa (@taxidCollection)
    {
          $uid_count = getSeqsFromNCBI($taxa, $base_url, $database,$dbloc, $seq_filename);
	  print "There are $uid_count $database sequences for $taxonCollection[$count] ($taxa) at the NCBI.\n";
          $count++;
	  if ($uid_count == 0) { push(@zeroCountTaxa, $taxa); }
	  $sum_uid_count = $sum_uid_count + $uid_count;  
	  
    }
    $uid_count = $sum_uid_count;   #need this because otherwise if the last taxa in a list returns 0 sequences, no blast
				   #database would be made. 
    print "$uid_count $database sequences have been fetched for $count taxa from $taxfile\n";

    my $zeroCounts = @zeroCountTaxa;
    unless ($zeroCounts == 0) 
    {
	print "Taxids with entries in NCBI, but for which no $database sequences were retrieved:\n";
	#if we'd done the taxon and taxid link via hashes, we could print out the name of the organism as well as the id
        #maybe a fix for another time.
	foreach my $noinf (@zeroCountTaxa) 
	{
		print "$noinf\n";
	}


    }

}



#I'm getting an error I shouldn't be getting in the nuc file:  Error: Error: failed to retrieve docsum
#I can't seem to remove it by working on the $efetch_out string, even when converted to an array.
	#The following lines are a nasty hack to remove these lines from the file using the perl command line
	#This should be fixed someday....perhaps something wrong with my URL, or something wrong with the db records themselves?
#UGLY HACK!!!
  {
	my $uncorrectedfilename = "$dbloc/$seq_filename";
	my $correctedfilename = "$uncorrectedfilename" . "_corrected";
	if ( -e "$dbloc/$seq_filename" ) #if your uid_count was 0, this file shouldn't exist
	{
		system ("perl -nle 'print unless /Error/' $uncorrectedfilename > $correctedfilename");
#delete uncorrected file
		unlink ($uncorrectedfilename);
		rename($correctedfilename, $uncorrectedfilename);
		print "\nYour $database sequences have been written to $uncorrectedfilename\n\n";
	}
  }

	if (($makeblast) && ($uid_count > 0)) 
	{
	makeblastdb($database, $seq_filename, $dbloc, $blast_dbname);

		 if ($gettaxid == 1  or $gettaxon == 1)
		 {
			print "A blast database called $blast_dbname containing $database sequences from $taxon has been made in directory $dbloc.\n\n";
		 }
		else 
		{
			print "A blast database called $blast_dbname containing $database sequences from taxa in $taxfile has been made in directory $dbloc.\n\n";
		}

	}

 }
}
#}

#################
sub getTaxidFromTaxon
{
 my $taxon = $_[0];

   	$taxon_change = $taxon;
      if ($taxon_change =~/(\w+)\s+(\w+)/)
         {
	     $taxon_change = $1 . "%20" . $2;
         }

     my $tesearch_url = $base_url . "esearch.fcgi?db=$tdatabase&term=$taxon_change";

      #print "Using this URL:\n $tesearch_url\n to get your taxon id...";

      my $tinit_output = get($tesearch_url);

      #print "$init_output\n";

      my @tinit_lines = split(/^/, $tinit_output);

          foreach my $line (@tinit_lines)
        {
          if ($line=~/<Id>(\d+)<\/Id>/)
             {
                $taxid = $1;
		$continue = 1;
              }
         }

       return ($taxid, $continue);
} 


#####################################################
sub getTaxonFromTaxid
{
  my $taxid = $_[0];
     $efetch_url = $base_url."efetch.fcgi?db=$tdatabase&id=$taxid&report=brief&mode=text";

    # print "Using this URL:\n $efetch_url\n to get your taxon ...\n";

     $taxon = get($efetch_url);
     chomp($taxon);

    return $taxon;
	
}


##############################################

sub getSeqsFromNCBI
{

my $taxid =  $_[0];
my $base_url = $_[1];
my $database = $_[2];
my $dbloc  = $_[3];
my $seq_filename = $_[4];
my $uid_count = 0;

#sample URL that works and returns a file
my $esearch_url = $base_url . "esearch.fcgi?db=" . $database . "&dbfrom=taxonomy&usehistory=y&term=txid" . $taxid . "[Organism:exp]";


my $init_output = get($esearch_url);

#it would be better if I'd done this with proper DOM handling, but I didn't
#I used XML::Simple...because it was easier!
my $xmlref = XMLin($init_output);
$Webenv = $xmlref->{'WebEnv'};
$uid_count = $xmlref->{'Count'};
$querkey = $xmlref->{'QueryKey'};

if ($uid_count>0)
{
#inefficient to be opening and closing this file within the loop.
  open (SEQFILE, ">>$dbloc/$seq_filename") || die "I could not open the file $dbloc/$seq_filename to write sequences to: $!\n";
   my $efetch_url;
  for (my $retstart = 0; $retstart < $uid_count; $retstart += $retmax)
  {
	$efetch_url = $base_url  ."efetch.fcgi?db=$database&WebEnv=$Webenv&query_key=$querkey&retmode=text";

	$efetch_url .= "&retstart=$retstart&retmax=$retmax&rettype=$format";

	my $efetch_out = get($efetch_url);
	print SEQFILE "$efetch_out";

  }

 close SEQFILE;
 }
return ($uid_count);

}

#############################################

sub makeblastdb 
{
	my $database = $_[0];
	my $seq_filename = $_[1];
        my $dbloc = $_[2];
        my $blast_dbname = $_[3];
	
	#check which settings to use for formatdb
	my $pep = "T";   #default setting for -p is T (peptide). Set to F for nucleotide:
	if ($database =~ /[Nn][Uu][Cc]/)  { $pep = "F"; }
	
	my $dbtype = substr($database,0,1);
	
	my $bldb_name = $blast_dbname. "_" . $dbtype;
	
	my $working_dir = getcwd();
	
	if ($blastdatadir)   { chdir $blastdatadir; }

	system ("formatdb -i $dbloc/$seq_filename -p $pep -n $bldb_name");

	chdir $working_dir;
	
	return 0;

}

#####################################
sub createdirs 
{
  foreach my $i (@_) 
  {
   unless( -d $i ) { mkdir ($i, 0755) or die "Cannot create directory $i. Quitting program because: $!\n"; } 
  }
}  
  

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	############################Notes about the script########################
	
# Need to get a list of accessions that match a species for
# 	a) nucs
# 	b) peps

# For each of the nuc and pep lists:	
#	Find out the number of accessions involved
################################
#To do the above - run Esearch in the Web Environment mode
#This will retrieve the total number of UIDs that match the query
#Is returned in the <Count> tag in the ESearch output

#Save the value of <Count> to $count
#The values of WebEnv into $Webenv
#The query_key into $key

######################################

#	set up a loop so that a URL is generated to retrieve entries in batches of 500
#Retrieve the entries

###################################

#Do this using EFetch to retrieve each batch of a size $retmax, e.g. $retmax=500

#Format as needed

