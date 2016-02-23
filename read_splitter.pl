#!/bin/bash

#####################################################
# written by Bob Settlage, DAC, 17 Sep 2012
# converting fastq to fasta and splitting in N references per subfile
#
#
#usage perl fasta_header_length.pl <fasta>
#example perl /groups/DAC/useful_perl/read_splitter.sh <fastq> <references per subfile>
#
#####################################################

$temp_NAME = $ARGV[0];
$target_NUM = $ARGV[1];

file_NAME=${temp_NAME##*/}
file_DIR=${temp_NAME%/*}
cd $file_DIR
pwd

echo splitting on $file_NAME into $target_NUM references in directory $file_DIR

mkdir $file_NAME
cp $file_NAME ./$file_NAME
cd ./$file_NAME

file_EXT=${file_NAME##*.}
if [ ${file_EXT} == 'gz' ];then
	fasta_suffix=.fasta
	new_file_NAME=$file_NAME$fasta_suffix
	cat $file_NAME | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){s/\./N/g;print;$i=-3}$i++;}' > $new_file_NAME
	file_NAME=$new_query_FILE
fi
file_EXT=${file_NAME##*.}
if [ ${file_EXT} == 'fastq' ] || [ ${file_EXT} == 'fq' ];then
	fasta_suffix=.fasta
	new_file_NAME=$file_NAME$fasta_suffix
	cat $file_NAME | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){s/\./N/g;print;$i=-3}$i++;}' > $new_file_NAME
	file_NAME=$new_query_FILE
fi

lines=`wc -l $f | cut -f1 -d' '` 







exit;