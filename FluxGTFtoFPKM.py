##################################
#                                #
# Last modified 09/15/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import math
import random
import string

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s gtf outfilename' % argv[0]
        sys.exit(1)

    gtf=argv[1]
    outputfilename = argv[2]

    TranscriptDict={}

    lineslist = open(gtf)
    for line in lineslist:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'transcript':
            continue
        if '; RPKM ' not in fields[8]:
            continue
        geneID = fields[8].split('gene_id "')[1].split('";')[0]
        transcriptID = fields[8].split('transcript_id "')[1].split('";')[0]
        RPKM = fields[8].split('; RPKM ')[1]
        TranscriptDict[(geneID,transcriptID)] = RPKM

    outfile = open(outputfilename, 'w')
    outline = '#GeneID\tTranscriptName\tRPKM'   
    outfile.write(outline + '\n')

    keys = TranscriptDict.keys()
    keys.sort()

    for (geneID,transcriptID) in keys:
        outline = geneID + '\t' + transcriptID + '\t' + str(TranscriptDict[(geneID,transcriptID)])
        outfile.write(outline + '\n')

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

