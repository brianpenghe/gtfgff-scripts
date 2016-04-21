##################################
#                                #
# Last modified 06/11/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import math
import random
import string

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s annotaiton_gtf gtf outfilename' % sys.argv[0]
        print '\t use - to read the gtf file from stdin'
        sys.exit(1)

    annotation_gtf = sys.argv[1]
    gtf=sys.argv[2]
    outputfilename = sys.argv[3]

    lineslist = open(annotation_gtf)
    GeneDict={}
    for line in lineslist:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        if 'gene_name "' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName=fields[8].split('gene_id "')[1].split('";')[0]
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        GeneDict[geneID]=geneName

    outfile = open(outputfilename, 'w')
    outline = '#GeneID\tGeneName\tRPKM1\tRPKM2'   
    outfile.write(outline + '\n')

    skipped = 0
    lineslist = open(gtf)
    for line in lineslist:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        geneID = fields[8].split('gene_id "')[1].split('";')[0]
        if 'RPKM1 "' in fields[8] and 'RPKM2 "':
            RPKM1 = fields[8].split('RPKM1 "')[1].split('";')[0]
            RPKM2 = fields[8].split('RPKM2 "')[1].split('";')[0]
        else:
            RPKM1 = fields[8].split('RPKM "')[1].split('";')[0]
            RPKM2 = 'n/a'
        if GeneDict.has_key(geneID):
            pass
        else:
            print 'skipping geneID', geneID, ', not found in annotation'
            skipped+=1
            continue
        outline = geneID + '\t'+ GeneDict[geneID] + '\t'+ RPKM1 + '\t'+ RPKM2
        outfile.write(outline + '\n')

    print 'skipped', skipped, 'entries'

    outfile.close()

run()

