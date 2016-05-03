##################################
#                                #
# Last modified 11/22/2011       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s gtf outputfilename ' % argv[0]
        sys.exit(1)
    
    inputfilename = argv[1]
    outfilename = argv[2]

    outfile = open(outfilename, 'w')
    outline = '#chr\tleft\tright\tstrand\tGeneID(s)\tGeneName(s)\tTranscriptID(s)\tTranscriptName(s)'
    outfile.write(outline + '\n')

    ExonDict = {}

    linelist = open(inputfilename)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        chr = fields[0]
        start = int(fields[3])
        stop = int(fields[4])
        strand = fields[6]
        exon = (chr,start,stop,strand)
        if ExonDict.has_key(exon):
            pass
        else:
            ExonDict[exon] = {}
            ExonDict[exon]['geneIDs'] = []
            ExonDict[exon]['transcriptIDs'] = []
            ExonDict[exon]['geneNames'] = []
            ExonDict[exon]['transcriptNames'] = []
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
        ExonDict[exon]['geneIDs'].append(geneID)
        ExonDict[exon]['transcriptIDs'].append(transcriptID)
        ExonDict[exon]['geneNames'].append(geneName)
        ExonDict[exon]['transcriptNames'].append(transcriptName)


    exons = ExonDict.keys()
    exons.sort()

    for exon in exons:
        (chr,start,stop,strand) = exon
        outline = chr + '\t' + str(start)  + '\t' + str(stop) + '\t' + strand + '\t'
        ExonDict[exon]['geneIDs'] = list(Set(ExonDict[exon]['geneIDs']))
        ExonDict[exon]['geneIDs'].sort()
        for geneID in ExonDict[exon]['geneIDs']:
            outline = outline + geneID + ','
        outline = outline[0:-1] + '\t'
        ExonDict[exon]['geneNames'] = list(Set(ExonDict[exon]['geneNames']))
        ExonDict[exon]['geneNames'].sort()
        for geneName in ExonDict[exon]['geneNames']:
            outline = outline + geneName + ','
        outline = outline[0:-1] + '\t'
        ExonDict[exon]['transcriptIDs'] = list(Set(ExonDict[exon]['transcriptIDs']))
        ExonDict[exon]['transcriptIDs'].sort()
        for transcriptID in ExonDict[exon]['transcriptIDs']:
            outline = outline + transcriptID + ','
        outline = outline[0:-1] + '\t'
        ExonDict[exon]['transcriptNames'] = list(Set(ExonDict[exon]['transcriptNames']))
        ExonDict[exon]['transcriptNames'].sort()
        for transcriptName in ExonDict[exon]['transcriptNames']:
            outline = outline + transcriptName + ','
        outline = outline[0:-1] + '\t'
        outfile.write(outline+'\n')
   
    outfile.close()
   
if __name__ == '__main__':
    main(sys.argv)
