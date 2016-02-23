#!/bin/bash
# sample Sun Grid Engine script to execute blast

# Specify the shell for this job
#$ -S /bin/bash

# Merge the standard output and standard error
#$ -j y

# Specify the location of the output
#$ -o $HOME/

#
# Explanation of blast options:
#
# Arguments:
#
# -p blastn                     Use the blastn program
# -d nt                         Use the nt database
# -i ${seqpath}/nt.123          Use the nt.123 query file
# -e 0.1                        Expectation value
# -o ${outpath}/out.123.blastn  Output file (final report)
# -a (number)                   number of cores to use

# Location of the BLAST formatted databases
export BLASTDB=/common/groups/dac/blastdbs

# Location of the BLAST matrices
export BLASTMAT=/common/shared/blast/data

# Location of BLAST executables
progpath=/common/shared/blast/bin

# Location of BLAST query files
seqpath=/common/groups/dac/blastdbs

# Location of BLAST output (final report) file
outpath=/common/groups/dac/

echo -n "---> begin: " ; date
${progpath}/blastall -p blastn \
    -F F \
    -e 1e-10 \
    -a 5 \
    -d /common/groups/dac/CLC_Velvet_SIPES_PE_PF_LIPES3k/combined_CLC_Velvet_SIPES_PE_PF_LIPES3k.fa \
    -i /common/groups/dac/CLC_Velvet_SIPES_PE_PF_LIPES3k/combined_CLC_Velvet_SIPES_PE_PF_LIPES3k.fa \
    -o /common/groups/dac/CLC_Velvet_SIPES_PE_PF_LIPES3k/combined_blastn_out.txt
echo -n "---> end: " ; date

