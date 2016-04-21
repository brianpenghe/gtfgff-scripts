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

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s gtf outfilename [-longestOnly]' % sys.argv[0]
        sys.exit(1)

    gtf = sys.argv[1]
    outfilename = sys.argv[2]

    doLongest=False
    if '-longestOnly' in sys.argv:
        doLongest=True

    GeneDict={}

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        if 'gene_name' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName=fields[8].split('gene_id "')[1].split('";')[0]
        if 'transcript_name' in fields[8]:
            transcriptName=fields[8].split('transcript_name "')[1].split('";')[0]
        else:
            transcriptName=fields[8].split('transcript_id "')[1].split('";')[0]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if GeneDict.has_key((geneName,geneID)):
            pass
        else:
            GeneDict[(geneName,geneID)]={}
        if GeneDict[(geneName,geneID)].has_key((transcriptID,transcriptName)):
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)].append((chr,left,right))
        else:
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]=[]
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)].append((chr,left,right))

    outfile = open(outfilename, 'w')

    outline='#GeneID\tGeneName\tLongest_Isoform_Length\tTranscriptID\tTranscriptName\tTranscript_Length\n'
    outfile.write(outline)

    for (geneName,geneID) in GeneDict.keys():
        longestLength=0
        for (transcriptID,transcriptName) in GeneDict[(geneName,geneID)].keys():
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]=list(Set(GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]))
            length=0
            for (chr,left,right) in GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]:
                length+=(right-left)
            if length > longestLength:
                longestLength = length
        for (transcriptID,transcriptName) in GeneDict[(geneName,geneID)].keys():
            outline=geneID + '\t' + geneName + '\t' + str(longestLength)
            GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]=list(Set(GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]))
            length=0
            for (chr,left,right) in GeneDict[(geneName,geneID)][(transcriptID,transcriptName)]:
                length+=(right-left)
            if doLongest:
                if length == longestLength:
                    outline= outline + '\t' + transcriptID + '\t' + transcriptName + '\t' + str(length)
                    outfile.write(outline+'\n')
                    break
            else:
                outline= outline + '\t' + transcriptID + '\t' + transcriptName + '\t' + str(length)
                outfile.write(outline+'\n')

    outfile.close()
        
run()

