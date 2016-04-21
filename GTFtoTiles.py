##################################
#                                #
# Last modified 09/13/2014       # 
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
        print 'usage: python %s gtf outprefix' % sys.argv[0]
        sys.exit(1)

    inputfilename = sys.argv[1]
    outprefix = sys.argv[2]

    GeneDict = {}

    listoflines = open(inputfilename)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        geneType = fields[8].split('gene_type "')[1].split('"')[0]
        geneID = fields[8].split('gene_id "')[1].split('"')[0]
        if 'gene_name "' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('"')[0]
        else:
            geneName = geneID
        transcriptID = fields[8].split('transcript_id "')[1].split('"')[0]
        if 'transcript_name "' in fields[8]:
            transcriptName=fields[8].split('transcript_name "')[1].split('"')[0]
        else:
            transcriptName = transcriptID
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        strand=fields[6]
        if GeneDict.has_key((geneType,strand)):
            pass
        else:
            GeneDict[(geneType,strand)] = {}
        if GeneDict[(geneType,strand)].has_key((geneID,geneName)):
            pass
        else:
            GeneDict[(geneType,strand)][(geneID,geneName)]={}
        if GeneDict[(geneType,strand)][(geneID,geneName)].has_key((transcriptID,transcriptName)):
            pass
        else:
            GeneDict[(geneType,strand)][(geneID,geneName)][(transcriptID,transcriptName)] = []
        GeneDict[(geneType,strand)][(geneID,geneName)][(transcriptID,transcriptName)].append((chr,left,right))

    print GeneDict

    for (geneType,strand) in GeneDict.keys():
        print geneType,strand
        outfileGenes = open(outprefix + '.' + geneType + '.genes.' + strand.replace('+','plus').replace('-','minus') + '.tile','w')
        outfileGenesList = []
        outfileTranscripts = open(outprefix + '.' + geneType + '.transcripts.' + strand.replace('+','plus').replace('-','minus') + '.tile','w')
        outfileTranscriptsList = []
        outfileExons = open(outprefix + '.' + geneType + '.exons.' + strand.replace('+','plus').replace('-','minus') + '.tile','w')
        outfileExonsList = []
        for (geneID,geneName) in GeneDict[(geneType,strand)].keys():
            Gpositions = []
            for (transcriptID,transcriptName) in GeneDict[(geneType,strand)][(geneID,geneName)].keys():
                Tpositions = []
                for (chr,left,right) in GeneDict[(geneType,strand)][(geneID,geneName)][(transcriptID,transcriptName)]:
                    outfileExonsList.append((chr,left,right))
                    Tpositions.append(left)
                    Tpositions.append(right)
                    Gpositions.append(left)
                    Gpositions.append(right)
                outfileTranscriptsList.append((chr,min(Tpositions),max(Tpositions),'id=' + transcriptName))
            outfileGenesList.append((chr,min(Gpositions),max(Gpositions),'id=' + geneName))
        outfileExonsList.sort()
        for (chr,left,right) in outfileExonsList:
            outline = chr + '\t' + str(left) + '\t' + str(right)
            outfileExons.write(outline + '\n')
        outfileTranscriptsList.sort()
        for (chr,left,right,ID) in outfileTranscriptsList:
            outline = chr + '\t' + str(left) + '\t' + str(right) + '\t' + ID
            outfileTranscripts.write(outline + '\n')
        outfileGenesList.sort()
        for (chr,left,right,ID) in outfileGenesList:
            outline = chr + '\t' + str(left) + '\t' + str(right) + '\t' + ID
            outfileGenes.write(outline + '\n')
        outfileExons.close()        
        outfileGenes.close()
        outfileTranscripts.close()

run()

