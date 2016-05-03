##################################
#                                #
# Last modified 03/02/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
from sets import Set

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s gtf repeatmasker outfilename' % argv[0]
        print '\tNote: the following repeatMasker format is assumed:'
        print '\t585     2643    40      30      0       chr2LHet        0       328     -368544 -       Copia2_I        LTR     Copia   -2478   1684    1347    1'
        print '\tNote: it is also assumed that the left coordinates of repetitive elements are unique'
        sys.exit(1)

    gtf = argv[1]
    RM = argv[2]
    outfilename = argv[3]

    RepeatDict = {}
    lineslist  = open(RM)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        chr = fields[5]
        left = int(fields[6])
        right = int(fields[7])
        if RepeatDict.has_key(chr):
            pass
        else:
            RepeatDict[chr]={}
        RepeatDict[chr][left]=fields

    GeneDict={}

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        chr= fields[0]
        left = int(fields[3])
        right = int(fields[4])
        strand = fields[6]
        geneID = fields[8].split('gene_id "')[1].split('";')[0]
        if 'gene_name' in fields[8]:
            geneName = fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName = fields[8].split('gene_id "')[1].split('";')[0]
        if 'transcript_name' in fields[8]:
            transcriptName = fields[8].split('transcript_name "')[1].split('";')[0]
        else:
            transcriptName = fields[8].split('transcript_id "')[1].split('";')[0]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if GeneDict.has_key((geneName,geneID)):
            pass
        else:
            GeneDict[(geneName,geneID)]={}
        if GeneDict[(geneName,geneID)].has_key((transcriptID,transcriptName)):
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)].append((chr,left,right,strand))
        else:
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]=[]
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)].append((chr,left,right,strand))

    outfile = open(outfilename, 'w')

    outline='#GeneID\tGeneName\tTranscriptID\tTranscriptName\tIntronLeft\tIntronRight\tIntronNumber\tNumberIntrons\tRepeatFields\n'
    outfile.write(outline)

    print len(GeneDict.keys())

    for (geneName,geneID) in GeneDict.keys():
        longestLength=0
        for (transcriptID,transcriptName) in GeneDict[(geneName,geneID)].keys():
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]=list(Set(GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]))
            NumberIntrons = len(GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]) - 1
            strand = GeneDict[(geneName,geneID)][(transcriptID,transcriptName)][0][3]
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)].sort()
            if strand == '-':
                GeneDict[(geneName,geneID)][(transcriptID,transcriptName)].reverse()
            IN = 0
            for E in range(NumberIntrons):
                (chr,left,right,strand) = GeneDict[(geneName,geneID)][(transcriptID,transcriptName)][E]
                (nextExonchr,nextExonleft,nextExonright,nextstrand) = GeneDict[(geneName,geneID)][(transcriptID,transcriptName)][E+1]
                IN += 1
                for i in range(right,nextExonleft):
                    if RepeatDict.has_key(chr):
                        if RepeatDict[chr].has_key(i):
                            outline = geneID + '\t' + geneName + '\t' + transcriptID + '\t' + transcriptName + '\t' + str(right) + '\t' + str(nextExonleft) + '\t' + str(IN) + '\t' + str(NumberIntrons)
                            for fields in RepeatDict[chr][i]:
                                outline = outline + '\t' + fields
                            outfile.write(outline + '\n')

    outfile.close()
        
if __name__ == '__main__':
    main(sys.argv)

