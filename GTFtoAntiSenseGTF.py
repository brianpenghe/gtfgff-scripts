##################################
#                                #
# Last modified 12/10/2014       # 
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
        print 'usage: python %s gtf gene_sufffix outfile' % argv[0]
        sys.exit(1)

    gtf = argv[1]
    suff = argv[2]
    outfilename = argv[3]

    outfile = open(outfilename, 'w')

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
        strand = fields[6]
        transcriptID = fields[8].split('transcript_id "')[1].split('"')[0]
        geneID = fields[8].split('gene_id "')[1].split('"')[0]
        if 'gene_name "' in fields[8]:
            geneName = fields[8].split('gene_name "')[1].split('"')[0]
        else:
            geneName = geneID
        if 'transcript_name "' in fields[8]:
            transcriptName = fields[8].split('transcript_name "')[1].split('"')[0]
        else:
            transcriptName = geneID
        outline = fields[0] + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5] + '\t'
        if strand == '+':
            outline = outline + '-\t'
        if strand == '-':
            outline = outline + '+\t'
        outline = outline + fields[7] + '\t' + fields[8] + '\n'
        outline = outline.replace('"' + geneID + '";','"' + geneID + '-' + suff + '";')
        outline = outline.replace('"' + transcriptID + '";','"' + transcriptID + '-' + suff + '";')
        outline = outline.replace('"' + transcriptName + '";','"' + transcriptName + '-' + suff + '";')
        outline = outline.replace('"' + geneName + '";','"' + geneName + '-' + suff + '";')
        outfile.write(outline)
   

    outfile.close()
        
if __name__ == '__main__':
    main(sys.argv)

