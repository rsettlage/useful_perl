#! /bin/blast

FILEIN=$1
FILEOUT=$FILEIN".fasta"

echo parsing $FILEIN to $FILEOUT

`cat $FILEIN | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){s/\./N/g;print;$i=-3}$i++;}' > $FILEOUT`


