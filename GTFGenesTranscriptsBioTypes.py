##################################
#                                #
# Last modified 09/07/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
from sets import Set

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s gtf outfile [-field ID]' % sys.argv[0]
        print '\use the -field option if the bioType is not specified in the optional attributes as gene_type and transcript_type'
        sys.exit(1)

    gtf = sys.argv[1]
    outfilename = sys.argv[2]

    doField = False
    if '-field' in sys.argv:
        doField = True
        BTfieldID = int(sys.argv[sys.argv.index('-field') + 1])

    GeneDict={}

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        transcriptID = fields[8].split('transcript_id "')[1].split('"')[0]
        geneID = fields[8].split('gene_id "')[1].split('"')[0]
        if 'transcript_name "' in fields[8]:
            transcriptName = fields[8].split('transcript_name "')[1].split('"')[0]
        else:
            transcriptName = transcriptID
        if 'gene_name "' in fields[8]:
            geneName = fields[8].split('gene_name "')[1].split('"')[0]
        else:
            geneName = geneID
        if doField:
            geneType = fields[BTfieldID]
            transcriptType = 'n/a'
        else:
            geneType = fields[8].split('gene_type "')[1].split('"')[0]
            transcriptType = fields[8].split('transcript_type "')[1].split('"')[0]
        if GeneDict.has_key((geneID,geneName,geneType)):
            pass
        else:
            GeneDict[(geneID,geneName,geneType)] = {}
        GeneDict[(geneID,geneName,geneType)][(transcriptID,transcriptName,transcriptType)] = 1

    outfile = open(outfilename,'w')
    outline = '#GeneID\tGeneName\tGeneType\tTranscriptID\tTranscriptName\tTranscriptType'
    outfile.write(outline + '\n')

    keys = GeneDict.keys()
    keys.sort()

    for (geneID,geneName,geneType) in keys:
        transcripts = GeneDict[(geneID,geneName,geneType)].keys()
        transcripts.sort()
        for (transcriptID,transcriptName,transcriptType) in transcripts:
            outline = geneID + '\t' + geneName + '\t' + geneType + '\t' + transcriptID + '\t' + transcriptName + '\t' + transcriptType
            outfile.write(outline + '\n')

    outfile.close()
        
run()

