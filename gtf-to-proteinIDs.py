##################################
#                                #
# Last modified 05/19/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
from sets import Set

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s gtf outfilename' % argv[0]
        sys.exit(1)

    gtf = argv[1]
    outputfilename = argv[2]

    outfile = open(outputfilename, 'w')

    GeneTranscriptProteinDict={}

    listoflines = open(gtf)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        geneID = fields[8].split('gene_id "')[1].split('"')[0]
        geneName = fields[8].split('gene_name "')[1].split('"')[0]
        transcriptID = fields[8].split('transcript_id "')[1].split('"')[0]
        transcriptName = fields[8].split('transcript_name "')[1].split('"')[0]
        if 'protein_id "' in line:
            proteinID = fields[8].split('protein_id "')[1].split('"')[0]
        else:
            proteinID = 'N\A'
        GeneTranscriptProteinDict[(geneID,geneName,transcriptID,transcriptName,proteinID)]=1

    keys = GeneTranscriptProteinDict.keys()
    keys.sort()

    outline='#geneID\tgeneName\ttranscriptID\ttranscriptName\tproteinID'
    outfile.write(outline+'\n')

    for (geneID,geneName,transcriptID,transcriptName,proteinID) in keys:
        outline = geneID + '\t' + geneName + '\t' + transcriptID + '\t' + transcriptName + '\t' + proteinID
        outfile.write(outline+'\n')
            
    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

