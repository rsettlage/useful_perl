#Author: Allan Dickerman
#Takes list of fastq files and finds all unique read sequences
# Outputs unique sequences with count of how many times each is found
# Output format suitable for downstream analysis by usearch tools (Robert Edgar's software)
import argparse
import sys
import math
import fileinput

parser = argparse.ArgumentParser()
parser.add_argument('-m', default=2, action='store', dest='minCount', type=int, help='minimum count of output sequences, def=2')
parser.add_argument('files', metavar='FILE', nargs='*', help='fastq files to process (can be fastq.gz)')
args = parser.parse_args()

seqCount={}
for line in fileinput.input(args.files):
  if fileinput.isfirstline():
    sys.stderr.write('reading '+ fileinput.filename()+" numUniq=%d\n"%len(seqCount))
    index=0
  if fileinput.filelineno() % 4 == 2:
    seq = line.rstrip()
    if seq not in seqCount:
      seqCount[seq]=0
    seqCount[seq]+=1
sys.stderr.write("Done reading input files, numUniq=%d\n"%len(seqCount))

idDigits = int(math.log(len(seqCount))/math.log(10)+1.1)
OUT=open("derep.fa", 'w')
i = 1
for seq in sorted(seqCount, key=seqCount.get, reverse=True):
  if seqCount[seq] < args.minCount:
    break
  OUT.write(">DerepSeq_%0*d;size=%d;\n%s\n"%(idDigits, i, seqCount[seq], seq))
  i+=1
