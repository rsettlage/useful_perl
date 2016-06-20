import sys
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-m', default=0.8, action='store', dest='propThresh', type=float, help='proportion of reads must map to taxon (default: 0.8)')
parser.add_argument('--taxonomy', default="/home/allan/gg/97_otu_taxonomy.txt", action='store', dest='taxonomy', help='GreenGenes taxonomy (or similar)')
parser.add_argument('filename', action='store', help='fastq file')

args = parser.parse_args()

taxFile = open(args.taxonomy)
idToTax={}
for line in taxFile:
  f=line.split("\t")
  id = f[0]
  taxArray=f[1].rstrip().split('; ')
  #while taxArray[len(taxArray)-1].endswith('__'):
  while len(taxArray[-1]) < 4:
    #print "trimming: "+taxArray[-1]
    taxArray=taxArray[0:-1]
  idToTax[id] = taxArray

getOtuPattern = re.compile('OTU_(\d+)')
otuTaxonCount={}
otuCount={}
otuCountTotal={} # include non-matching lines
otuPctIdent={}
usearchFile = open(args.filename)
#print "opening "+sys.argv[1]
for line in usearchFile:
  f=line.rstrip().split()
  otu=f[8]
  m=getOtuPattern.match(otu)
  if not m:
    raise Exception("Cannot parse OTU: "+otu)
  otu="OTU_%04d"%int(m.group(1))
  taxId=f[9]
  pct=f[3]
  if not otuCountTotal.has_key(otu):
    otuCountTotal[otu]=0
  otuCountTotal[otu] += 1
  if line.startswith('N'):
    continue # non-match
  if not otuCount.has_key(otu):
    otuCount[otu]=0
    otuPctIdent[otu]=0
  otuCount[otu] += 1
  otuPctIdent[otu]+=float(pct)
  if not otuTaxonCount.has_key(otu):
    otuTaxonCount[otu]={}
  taxString = ''
  for i in range(1,len(idToTax[taxId])):
    taxString = "; ".join(idToTax[taxId][0:i])
    if not otuTaxonCount[otu].has_key(taxString):
      otuTaxonCount[otu][taxString] = 0
    otuTaxonCount[otu][taxString] += 1
  
otuTaxon={}
for otu in sorted(otuTaxonCount.keys()):
  taxonCatLen=0
  countThreshold = int(otuCount[otu]*args.propThresh)
  for taxonString in otuTaxonCount[otu]:
    if (otuTaxonCount[otu][taxonString] > countThreshold):
      if (taxonString.count(';') > taxonCatLen):
        otuTaxon[otu]=taxonString
        taxonCatLen=taxonString.count(';')

for otu in sorted(otuTaxon):
  taxonString = otuTaxon[otu]
  pctIdent = otuPctIdent[otu]/otuCount[otu]
  print otu+"\t"+taxonString+"\t%d/%d(%d)\t%.1f"%(otuTaxonCount[otu][taxonString], otuCount[otu], otuCountTotal[otu], pctIdent)
    


"""
/home/allan/gg/97_otu_taxonomy.txt
367523	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__
187144	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
836974	k__Bacteria; p__Cyanobacteria; c__Chloroplast; o__Cercozoa; f__; g__; s__
310669	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__
823916	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Enhydrobacter; s__
878161	k__Bacteria; p__Acidobacteria; c__Acidobacteriia; o__Acidobacteriales; f__Acidobacteriaceae; g__Terriglobus; s__
3064251	k__Bacteria; p__Verrucomicrobia; c__Opitutae; o__Puniceicoccales; f__Puniceicoccaceae; g__Puniceicoccus; s__
1138555	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Caldicoprobacteraceae; g__Caldicoprobacter; s__
3918	k__Bacteria; p__Spirochaetes; c__Spirochaetes; o__Spirochaetales; f__Spirochaetaceae; g__Treponema; s__
339472	k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Rhodospirillaceae; g__; s__

seqs20Per_vs_gg.usearch
OTU_1000;size=23;	301498	98.4	0
OTU_1000;size=23;	258202	99.4	0
OTU_1000;size=23;	301498	97.4	0
OTU_1004;size=43;	261084	97.1	0
OTU_1004;size=43;	173883	98.2	0
OTU_1004;size=43;	276478	97.1	0
OTU_1004;size=43;	173883	99.0	0
OTU_1004;size=43;	322840	98.2	0
OTU_1004;size=43;	173883	97.1	0
OTU_1004;size=43;	276478	97.4	0
H	32083	159	91.8	+	0	0	556I159M654I	OTU_1000;size=2;	265106
H	32083	185	92.4	+	0	0	545I185M639I	OTU_1000;size=2;	265106
H	32083	185	90.8	+	0	0	545I185M639I	OTU_1000;size=2;	265106
H	32083	185	92.4	+	0	0	545I185M639I	OTU_1000;size=2;	265106
"""
