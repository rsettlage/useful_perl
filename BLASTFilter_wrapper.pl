#! /usr/bin/perl -w

###Bob Settlage VBI 2011

###perl blastFiltering on single or paired reads

###using Hongseok's script
###filterReads_withBlast.pl  -b <blast table file 1>  -b <blast table file 2>  -s <sequence file 1>  -s <sequence file 2>
###        ( -e <cutoff e-value[:reference]>  [-e <cutoff e-value[:reference]>  ....] )  or 
###        ( -l <(cutoff match length):(cutoff percentage)> [-l<(cutoff match length):(cutoff percentage)> ...] )

###  ex 1) filterReads_withBlast.pl -b s_1_1.bln.table -b s_1_2.bln.table -s s_1_1.fastq -s s_1_2.fastq  -e 1e-04:pTARBAC2.1 -e 1e-07:NC_010473
###  ex 2) filterReads_withBlast.pl -b s_1_1.bln.table -b s_1_2.bln.table -s s_1_1.fastq -s s_1_2.fastq  -e 1e-04
###  ex 3) filterReads_withBlast.pl -b s_1_1.bln.table -b s_1_2.bln.table -s s_1_1.fastq -s s_1_2.fastq  -l 60:10 -l 50:90 


### for this script, we will use a 
###
### find ./ -type f \( -iname 's_?_1*.bln.raw.PBS' \) -exec perl BLASTFilter_wrapper '{}' s -e 1e-04 \;
###
### to do the launching, so we will take in a single file argument and launch Hongseok's scrips

### perl filterBlastWrapper <s|p> <blast table file> <e|l> <cutoffs>


###launch appropriate 
	
	$numArgs = @ARGV;
	if ($ARGV[0] eq 's'){
		$i = 3;
		if ($ARGV[2] eq 'e') {
			while ($i<$numArgs){
				$limitList .= ' -e '.$ARGV[$i];
				$i++;
			}
		}else{
			while ($i<$numArgs+1){
				$limitList .= ' -l '.$ARGV[$i];
				$i++;
			}
		}
	
		my @files = glob($ARGV[1]);
		foreach my $file (@files) {
			$command_string = "/common/groups/dac/useful_perl/filterReads_withBlast.pl";
			my @tempArray = split(/\./,$files[$j]);
			$base_FileName = $tempArray[0];
			print "base file name is $base_FileName with a total of $numArgs arguments\n";
			$blastTable = $files[$j];
			$sequenceFile = $base_FileName.'.fastq';
			$command_string .= " -b ". $blastTable . " -s " . $sequenceFile . " " . $limitList;
			print "launching perl $command_string\n";
			system("perl $command_string");
			$j++;
		}
		
	}else{
		$i = 4;
		if ($ARGV[3] eq 'e') {
			while ($i<$numArgs){
				$limitList .= ' -e '.$ARGV[$i];
				$i++;
			}
		}else{
			while ($i<$numArgs+1){
				$limitList .= ' -l '.$ARGV[$i];
				$i++;
			}
		}
		
		my @files = glob($ARGV[1]);
		foreach my $file (@files) {
			$command_string = "/common/groups/dac/useful_perl/filterReads_withBlast.pl";
			my @tempArray1 = split(/\./,$files[$j]);
			$base_FileName1 = $tempArray1[0];
			print "base file name is $base_FileName with a total of $numArgs arguments\n";
			$blastTable1 = $files[$j];
			$sequenceFile1 = $base_FileName1.'.fastq';
			$j++;
			my @tempArray2 = split(/\./,$files[$j]);
			$base_FileName2 = $tempArray2[0];
			print "base file name is $base_FileName with a total of $numArgs arguments\n";
			$blastTable2 = $files[$j];
			$sequenceFile2 = $base_FileName2.'.fastq';
			
			$command_string .= " -b ". $blastTable1 . " -b " . $blastTable2 . " -s " . $sequenceFile1 . " -s " . $sequenceFile2 . " " . $limitList;
			print "launching perl $command_string\n";
			system("perl $command_string");
			$j++;
		}
		
	}
	




