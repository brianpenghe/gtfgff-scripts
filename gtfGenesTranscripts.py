##################################
#                                #
# Last modified 03/11/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
from sets import Set

try:
	import psyco
	psyco.full()
except:
	pass

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s inputfilename outfilename' % sys.argv[0]
        sys.exit(1)

    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]

    TranscriptDict={}
    listoflines = open(inputfilename)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        geneID=fields[8].split('gene_id "')[1].split('"')[0]
        transcriptID=fields[8].split('transcript_id "')[1].split('"')[0]
        if 'gene_name "' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('"')[0]
        else:
            geneName=geneID
        if 'transcript_name "' in fields[8]:
            transcriptName=fields[8].split('transcript_name "')[1].split('"')[0]
        else:
            transcriptName = transcriptID
        if TranscriptDict.has_key((geneID,geneName,transcriptID,transcriptName)):
            continue
        else:
            TranscriptDict[(geneID,geneName,transcriptID,transcriptName)] = 0

    outfile = open(outputfilename, 'w')
    outfile.write('#GeneID\tGeneName\tTranscriptID\tTranscriptName\n')

    for (geneID,geneName,transcriptID,transcriptName) in TranscriptDict.keys():
        outline= geneID + '\t' + geneName + '\t' + transcriptID + '\t' + transcriptName
        outfile.write(outline+'\n')
            
    outfile.close()

run()

