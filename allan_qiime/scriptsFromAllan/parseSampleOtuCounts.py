import sys
import fileinput
import re

getSamplePattern=re.compile("S(\d+)")
getOTUPattern=re.compile("(OTU_\d+)")
sampleOtuCount={}
otuCount={}
sampleCount={}
sampleList=[]
curFile=''

sample=''
for line in fileinput.input():
  if (fileinput.filename() != curFile):
    #print "changing file from %s\t%d" % (sample, sampleCount[sample])
    curFile = fileinput.filename()
    m = getSamplePattern.search(curFile)
    sample = m.group(1)
    sample = "S%02d"%int(sample)
    sampleOtuCount[sample]={}
    sampleCount[sample]=0
    sampleList.append(sample)
    print "Now reading from "+curFile+", sample "+sample+" lines.read=%d"%fileinput.lineno()
  f=line.split('\t')
  m=getOTUPattern.search(line)
  if (m):
    otu = m.group(1)
    if otu not in otuCount:
      otuCount[otu]=0
    if otu not in sampleOtuCount[sample]:
      sampleOtuCount[sample][otu] = 0
    sampleOtuCount[sample][otu] += 1
    otuCount[otu]+=1
    sampleCount[sample]+=1

for curSample in sorted(sampleCount, key=sampleCount.__getitem__, reverse=True):
  print "%s\t%d" % (curSample, sampleCount[curSample])
sampleList=sorted(sampleList)

OUT = open("CountsBySampleAndOTU.txt", 'w')
OUT.write("OTU\t"+"\t".join(sampleList)+"\n")
for otu in sorted(otuCount, key=otuCount.__getitem__, reverse=True):
  OUT.write(otu)
  for sample in sampleList:
    outNum=0
    if otu in sampleOtuCount[sample]:
      outNum=sampleOtuCount[sample][otu]
    OUT.write("\t%d"%outNum)
  OUT.write("\n")
OUT.close()
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
