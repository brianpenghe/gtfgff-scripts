##################################
#                                #
# Last modified 03/02/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import math
import random
import string
from cistematic.core.geneinfo import geneinfoDB
from cistematic.genomes import Genome
from sets import Set

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for i in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-i-1]]
    return sequence

try:
	import psyco
	psyco.full()
except:
	pass

def main(argv):

    if len(argv) < 4:
        print 'usage: python %s genome gtf outfilename [-polyA length]' % argv[0]
        sys.exit(1)

    genome = argv[1]
    gtf=argv[2]
    outputfilename = argv[3]
    doPolyA=False
    if '-polyA' in argv:
        doPolyA=True
        tailsize=int(argv[argv.index('-polyA')+1])
        tail=''
        for i in range(tailsize):
            tail=tail+'A'
        print 'will add a polyA tail of ', tailsize, 'nt'

    outfile = open(outputfilename, 'w')

    hg = Genome(genome)

    j=0
    lineslist = open(gtf)
    TranscriptDict={}
    for line in lineslist:
        j+=1
        if j % 100000 == 0:
            print j, 'lines processed'
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        if 'transcript_name "' in fields[8]:
            TranscriptID=fields[8].split('transcript_name "')[1].split('";')[0]
        else:
            TranscriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if TranscriptDict.has_key(TranscriptID):
            pass
        else:
            TranscriptDict[TranscriptID]=[]
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        orientation=fields[6]
        TranscriptDict[TranscriptID].append((chr,left,right,orientation))

    g=0 
    print 'Found', len(TranscriptDict.keys()), 'transcripts'
    for transcript in TranscriptDict.keys():
        g+=1
        if g % 1000 == 0:
            print g, 'transcripts sequences processed'
        sequence=''
        leftEnds=[]
        rightEnds=[] 
        TranscriptDict[transcript].sort()
        orientation = TranscriptDict[transcript][0][3]
        if orientation=='+' or orientation=='F':
            for (chr,left,right,orientation) in TranscriptDict[transcript]:
                leftEnds.append(left)
                rightEnds.append(right)
                try:
                    sequence=sequence+hg.sequence(chr[3:len(chr)],left,right-left)
                    print "can't retrieve sequence"
                except:
                    for p in range(left,right-left):
                        try:
                            sequence=sequence+hg.sequence(chr[3:len(chr)],p,1)
                        except:
                            sequence=sequence+'N'
                            missed+=1
            sense='plus_strand'
        if orientation=='-' or orientation=='R':
            for (chr,left,right,orientation) in reversed(TranscriptDict[transcript]):
                leftEnds.append(left)
                rightEnds.append(right)
                try:
                    exonsequence=hg.sequence(chr[3:len(chr)],left-1,right-left+1)
                    sequence=sequence+getReverseComplement(exonsequence)
                except:
                    for p in range(left-1,right-left+1):
                        try:
                            sequence=sequence+getReverseComplement(hg.sequence(chr[3:len(chr)],p,1))
                        except:
                            sequence=sequence+'N'
                            missed+=1
            sense='minus_strand'
        LeftEnd=min(leftEnds)
        RightEnd=max(rightEnds)
        outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense
        outfile.write(outline+'\n')
        if doPolyA:
            outfile.write(sequence+tail+'\n')
        else:
            outfile.write(sequence+'\n')

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

