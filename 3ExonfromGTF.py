##################################
#                                #
# Last modified 03/08/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
from sets import Set

def overlap(left1,right1,left2,right2):

    overlapOrNot = False
 
    if (right1 >= left2 and right1 <= right2) or (left1 >= left2 and left1 <= right2):
        overlapOrNot = True

    return overlapOrNot

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s gtf [single | merge] outfilename' % sys.argv[0]
        print '       Note: monoexonic genes will not be outputted'
        print '       merge option not done yet'
        sys.exit(1)

    gtf = sys.argv[1]
    MergeOrNot = sys.argv[2]
    outfilename = sys.argv[3]

    GeneDict={}

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        chr=fields[0]
        strand=fields[6]
        left=int(fields[3])
        right=int(fields[4])
        transcriptID = fields[8].split('transcript_id "')[1].split('"')[0]
        geneID = fields[8].split('gene_id "')[1].split('"')[0]
        if 'gene_name "' in fields[8]:
            geneName = fields[8].split('gene_name "')[1].split('"')[0]
        else:
            geneName = geneID
        if 'gene_name "' in fields[8]:
            transcriptName = fields[8].split('transcript_name "')[1].split('"')[0]
        else:
            transcriptName = transcriptID
        if GeneDict.has_key((geneID,geneName)):
            pass
        else:
            GeneDict[(geneID,geneName)]={}
        if GeneDict[(geneID,geneName)].has_key((transcriptID,transcriptName)):
            pass
        else:
            GeneDict[(geneID,geneName)][(transcriptID,transcriptName)] = []
        GeneDict[(geneID,geneName)][(transcriptID,transcriptName)].append((chr,left,right,strand))

    print 'finished inputting gtf'

    UTR5DictCoverage={}
    UTRDict={}

    outfile = open(outfilename, 'w')

    outfile.write('#chr\tleft\tright\tstrand\tGeneID\tGeneName\tTranscriptID\tTranscriptName\tOverlap_with_5UTR\n')

    for (geneID,geneName) in GeneDict:
        for (transcriptID,transcriptName) in GeneDict[(geneID,geneName)].keys():
            GeneDict[(geneID,geneName)][(transcriptID,transcriptName)].sort()
            chr = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][0][0]
            strand = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][0][3]
            if strand == '+':
                UTR5Left = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][0][1]
                UTR5Right = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][0][2]
            if strand == '-':
                UTR5Left = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][-1][1]
                UTR5Right = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][-1][2]
            if UTR5DictCoverage.has_key(chr):
                pass
            else:
                UTR5DictCoverage[chr]={}
            for i in range(UTR5Left,UTR5Right):
                UTR5DictCoverage[chr][i]=0
            if UTRDict.has_key(chr):
                pass
            else:
                UTRDict[chr]={}
        if len(GeneDict[(geneID,geneName)].keys()) == 1 and len(GeneDict[(geneID,geneName)][GeneDict[(geneID,geneName)].keys()[0]]) == 1:
            continue
        UTRs={}
        for (transcriptID,transcriptName) in GeneDict[(geneID,geneName)].keys():
            if strand == '+':
                UTR3Left = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][-1][1]
                UTR3Right = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][-1][2]
            if strand == '-':
                UTR3Left = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][0][1]
                UTR3Right = GeneDict[(geneID,geneName)][(transcriptID,transcriptName)][0][2]
            UTR = (chr,UTR3Left,UTR3Right,strand)
            if UTRs.has_key(UTR):
                pass
            else:
                UTRs[UTR]=[]
            UTRs[UTR].append((transcriptID,transcriptName))
        if MergeOrNot == 'single':
            pass
        if MergeOrNot == 'merge':
            UTRkeys = UTRs.keys()
            UTRkeys.sort()
#            if len(UTRkeys) == 1:
#                continue
#            for i in range(len(UTRkeys)-1):
#                if overlap(UTRkeys[i][1],UTRkeys[i][2],UTRkeys[i+1][1],UTRkeys[i+1][2]):
        for (chr,UTR3Left,UTR3Right,strand) in UTRs.keys():
            outline = chr + '\t' +str(UTR3Left) + '\t' + str(UTR3Right) + '\t' + strand + '\t' + geneID + '\t' + geneName + '\t'
            for (transcriptID,transcriptName) in UTRs[(chr,UTR3Left,UTR3Right,strand)]:
                outline = outline + transcriptID+ ','
            outline = outline[0:-1] + '\t'
            for (transcriptID,transcriptName) in UTRs[(chr,UTR3Left,UTR3Right,strand)]:
                outline = outline + transcriptName + ','
            outline = outline[0:-1]
            UTR5Overlap=False
            for i in range(UTR3Left,UTR3Right):
                if UTR5DictCoverage[chr].has_key(i):
                    UTR5Overlap=True
                    break
            if UTR5Overlap:
                outline = outline + '\tYes'
            else:
                outline = outline + '\t-'
            outfile.write(outline + '\n')

    outfile.close()
        
run()

