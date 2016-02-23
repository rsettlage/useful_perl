#!/usr/bin/perl -w
###############################################################################
###      @author: Christopher Caldwell                                      ###
###                                                                         ###
###                      align_fastq.pl:                                    ###
###                              removes differences between fastq files    ###
###                                                                         ###
###                                                                         ###
###############################################################################

my $file1 = "";
my $file2 = "";
my $outfile1 = "";
my $outfile2 = "";
my $errfile1 = "";
my $errfile2 = "";

###############################################################################
###                                                                         ###
###                                                                         ###
###                               checkopt()                                ###
###                                                                         ###
###                                                                         ###
###############################################################################
sub checkopt(){
    $argc = @ARGV;
    if($argc < 2){
        print("Usage: perl align_fastq.pl <file1> <file2>\n");
        exit(1);
    }
    else{
        $file1 = $ARGV[0];
        $file2 = $ARGV[1];
        
        ($file1 =~ /(.*\/|)(.*)\.(.*)/);    
        $outfile1 = $2 . "_aligned." . $3;
        $errfile1 = $2 . "_errored." . $3;
        ($file2 =~ /(.*\/|)(.*)\.(.*)/);    
        $outfile2 = $2 . "_aligned." . $3;
        $errfile2 = $2 . "_errored." . $3;
    }
    
}

###############################################################################
###                                                                         ###
###                                                                         ###
###                  run_grep_check()                                       ###
###                                                                         ###
###                                                                         ###
###############################################################################
sub run_grep_check(){
    $file = shift;
    $line = shift;
    @result = `grep \"$line\" $file`;
    $lines = @result;
    return ($lines > 0);
}

###############################################################################
###                                                                         ###
###                                                                         ###
###                           align_files()                                 ###
###                                                                         ###
###                                                                         ###
###############################################################################
sub align_files(){
     
    $running_1 = 1;
    $running_2 = 1;
    $read1 = 1;
    $read2 = 1;
    $written_1 = 0;
    $written_2 = 0;
    $write_f1 = "";
    $write_f2 = "";
    $line_f1 = "";
    $line_f2 = "";
 
    ### OPEN FILES
    open(FILE1, "< $file1");
    open(OUTFILE1, "> $outfile1");
    open(ERRFILE1, "> $errfile1");
    open(FILE2, "< $file2");
    open(OUTFILE2, "> $outfile2");
    open(ERRFILE2, "> $errfile2");

    while(($running_1 + $running_2) != 0){
        if($read1 && $running_1){
            if(!($line_f1 = <FILE1>)){
                $running_1 = 0;
                next;
            }
            chomp($line_f1);
            $write_f1 = $line_f1 . "\n" . <FILE1> . <FILE1> . <FILE1>; 
        }
        if($read2 && $running_2){
            if(!($line_f2 = <FILE2>)){
                $running_2 = 0;
                next;
            }
            chomp($line_f2);
            $write_f2 = $line_f2 . "\n" . <FILE2> . <FILE2> . <FILE2>;
        }
        $f1_end = chop $line_f1;
        $f2_end = chop $line_f2;
        
        if($line_f1 eq $line_f2){
            $read1 = 1;
            $read2 = 1;
            if(!$written_1){
                print OUTFILE1 $write_f1;
            }
            if(!$written_2){
                print OUTFILE2 $write_f2;
            }
            $written_1 = 0;
            $written_2 = 0;
           # print "$line_f1 and $line_f2 are aligned...\n";
        }
        else{
            $grep1 = &run_grep_check($file2, $line_f1);
            $grep2 = &run_grep_check($file1, $line_f2);
                
            if(!$grep1 && !$grep2){ # both lines are not in files
                $read1 = 1;
                $read2 = 1;
                $written_1 = 0;
                $written_1 = 0;
                print ERRFILE1 $write_f1;
                print ERRFILE2 $write_f2;
              #  print "$line_f1 and $line_f2 are unique...\n";
            }
            elsif(!$grep1){ # line_f1 is not in file_2
                $read1 = 1;
                $read2 = 0;
                print ERRFILE1 $write_f1;
                if(!$written_2){
                    print OUTFILE2 $write_f2;
                    $written_2 = 1;
                }
                $written_1 = 0;
              #  print "$line_f1 is unique...but $line_f2 is shared...\n";
            }
            elsif(!$grep2){ # line_f2 is not in file_1
                $read2 = 1;
                $read1 = 0;   
                print ERRFILE2 $write_f2;
                if(!$written_1){
                    print OUTFILE1 $write_f1;
                    $written_1 = 1;
                }
                $written_2 = 0;
              #  print "$line_f2 is unique...but $line_f1 is shared...\n";
            }
            else{ # lines are in files just not in the right order.
                $read1 = 1;
                $read2 = 1;
                if(!$written_1){
                    print OUTFILE1 $write_f1;
                }
                if(!$written_2){
                    print OUTFILE2 $write_f2;
                }
                $written_1 = 0;
                $written_2 = 0;
              #  print "$line_f1 and $line_f2  are shared...\n";
            }
             
            $line_f1 = $line_f1 . $f1_end;
            $line_f2 = $line_f2 . $f2_end;
        }
    }
 
    ### CLOSE FILES
    close(FILE1);
    close(OUTFILE1);
    close(ERRFILE1);
    close(FILE2);
    close(OUTFILE2);
    close(ERRFILE2);
}

###############################################################################
###                                                                         ###
###                                                                         ###
###                               main()                                    ###
###                                                                         ###
###                                                                         ###
###############################################################################
sub main(){
    checkopt();
    align_files();
}

main();
