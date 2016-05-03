##################################
#                                #
# Last modified 09/02/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import math
import random
import string

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','X':'X','a':'t','t':'a','g':'c','c':'g','n':'n','x':'x','R':'R','r':'r','M':'M','m':'m','Y':'Y','y':'y','S':'S','s':'s','K':'K','k':'k','W':'W','w':'w'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s fasta gtf outfilename [-polyA length] [-fastaChrFieldID ID] [-addChrToGTFchrID] [-noPositionalInformation] [-KeepUundeterminedStrand]' % argv[0]
        print '\t use the -KeepUundeterminedStrand if the file contains strands specified with a dot; those will be considered to be on the plus strand if the option is on, otherwise they will be skipped'
        sys.exit(1)

    fasta = argv[1]
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

    doKeepUundeterminedStrand = False
    if '-KeepUundeterminedStrand' in argv:
        doKeepUundeterminedStrand = True

    doAddChr=False
    if '-addChrToGTFchrID' in argv:
        doAddChr=True

    doFastaChrFieldID = False
    if '-fastaChrFieldID' in argv:
        doFastaChrFieldID = True
        FastaChrFieldID = int(argv[argv.index('-fastaChrFieldID')+1])

    doNoPosInfo=False
    if '-noPositionalInformation' in argv:
        doNoPosInfo=True

    outfile = open(outputfilename, 'w')

    GenomeDict={}
    sequence=''
    inputdatafile = open(fasta)
    for line in inputdatafile:
        if line[0]=='>':
            if sequence != '':
                GenomeDict[chr] = ''.join(sequence)
            chr = line.strip().split('>')[1]
            if doFastaChrFieldID:
                chr = chr.split('\t')[FastaChrFieldID]
            print chr
            sequence=[]
            Keep=False
            continue
        else:
            sequence.append(line.strip())
    GenomeDict[chr] = ''.join(sequence)

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
        chr=fields[0]
        if GenomeDict.has_key(chr):
            pass
        else:
            continue
        if 'gene_name "' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName=fields[8].split('gene_id "')[1].split('";')[0]
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        if 'transcript_name "' in fields[8]:
            transcriptName=fields[8].split('transcript_name "')[1].split('";')[0]
        else:
            transcriptName=fields[8].split('transcript_id "')[1].split('";')[0]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        transcript = (geneID, geneName, transcriptName, transcriptID)
        if TranscriptDict.has_key(transcript):
            pass
        else:
            TranscriptDict[transcript]=[]
        if doAddChr:
            if chr == 'MtDNA':
                chr = 'chrM'
            else:
                chr = 'chr' + chr
        left=int(fields[3])
        right=int(fields[4])
        orientation=fields[6]
        TranscriptDict[transcript].append((geneName,chr,left,right,orientation))

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
        orientation = TranscriptDict[transcript][0][4]
        if orientation == '.':
            if doKeepUundeterminedStrand:
                orientation = '+'
            else:
                continue
        if orientation=='+' or orientation=='F':
            for (geneName,chr,left,right,orientation) in TranscriptDict[transcript]:
                leftEnds.append(left)
                rightEnds.append(right)
                sequence=sequence+ GenomeDict[chr][left-1:right]
            sense='plus_strand'
        if orientation=='-' or orientation=='R':
            for (geneName,chr,left,right,orientation) in TranscriptDict[transcript]:
                leftEnds.append(left)
                rightEnds.append(right)
                sequence=sequence+GenomeDict[chr][left-1:right]
            sense='minus_strand'
            sequence = getReverseComplement(sequence)
        LeftEnd=min(leftEnds)
        RightEnd=max(rightEnds)
        (geneID, geneName, transcriptName, transcriptID) = transcript
        if doNoPosInfo:
            outline='>'+geneName + ':' + geneID + ':' + transcriptName + ':' + transcriptID
        else:
            outline='>'+geneName + ':' + geneID + ':' + transcriptName + ':' + transcriptID + ':' + chr + ':' + str(LeftEnd) + '-' + str(RightEnd) + '-' + sense
        outfile.write(outline+'\n')
        if doPolyA:
            sequence = sequence+tail
        for b in range(0,len(sequence ),50):
            outfile.write(sequence[b:min(b+50, len(sequence))] + '\n')

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

