#!/usr/bin/sh

#simple script for formatting databases

#####################################
## ARGV[0] = file to format
#####################################

database = ARGV[0];

#a = F means fasta format, use T for ASN.1
#b = F means text, T for binary
#p = F means nucleotide, use T for protein


exec /common/shared/blast/bin/formatdb -i database -a F -b F -p F -V



