import sys
import fileinput
maxReadsPerOTU=20
otuReads={}
readOtu={}
curFile=''
readId=''
for line in fileinput.input():
  if (fileinput.filename() != curFile):
    curFile = fileinput.filename()
    print "Now reading from "+curFile+" lines.read=%d"%fileinput.lineno()
    print "readId = "+readId
  f=line.split('\t')
  readId = f[0]
  otu = f[1]
  if not otuReads.has_key(otu):
    otuReads[otu] = {}
  if len(otuReads[otu]) < maxReadsPerOTU:
    otuReads[otu][readId] = 1 
    readOtu[readId]=''


#3_S3_L001_R1_001_AT_QT.paired_matched.fastq.overlapped.fasta
#../orig_data/1_S1_L001_R1_001_AT_QT.fastq
#fileBase="%d_S%d_L001_"%(sample, sample)
readId=""
numFound=0
for sample in range(1,17):
  FASTQ = open ("../orig_data/%d_S%d_L001_R1_001_AT_QT.fastq"%(sample, sample))
  lineNo = 0
  for line in FASTQ:
    if (lineNo % 4 == 0):
      readId = line.split()[0][1:]
    elif (lineNo % 4 == 1 and readOtu.has_key(readId)):
      readOtu[readId]=line.strip()
      numFound+=1
    lineNo+=1
print "Numfound = %d"%numFound
    
numFound = 0
OUT = open("SeqsPerOTU_max%d.txt"%maxReadsPerOTU, 'w')
for otu in sorted(otuReads):
  for read in otuReads[otu]:
    if readOtu[read]:
      OUT.write(">"+otu+"\n"+readOtu[read]+"\n")
      numFound+=1
OUT.close()
print "Numseqs written = %d"%numFound


"""
S01.vs.otusmin10.usearch
M00720:166:000000000-AF64C:1:1101:15742:1340	OTU_6;size=211773;	100.0	0
M00720:166:000000000-AF64C:1:1101:16813:1403	OTU_8;size=488536;	98.9	0
M00720:166:000000000-AF64C:1:1101:15678:1417	OTU_55;size=7575;	100.0	0
M00720:166:000000000-AF64C:1:1101:15758:1424	OTU_2;size=1959063;	100.0	0
M00720:166:000000000-AF64C:1:1101:14790:1435	OTU_5;size=260315;	100.0	0
M00720:166:000000000-AF64C:1:1101:14060:1435	OTU_108;size=2865;	99.2	0
M00720:166:000000000-AF64C:1:1101:14840:1439	OTU_8;size=488536;	100.0	0
M00720:166:000000000-AF64C:1:1101:16531:1442	OTU_2;size=1959063;	99.4	0
M00720:166:000000000-AF64C:1:1101:16322:1445	OTU_3420;size=76;	97.0	0
M00720:166:000000000-AF64C:1:1101:15421:1447	OTU_2;size=1959063;	99.4	0
"""
