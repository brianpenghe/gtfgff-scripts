##################################
#                                #
# Last modified 09/16/2011       # 
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
        print 'usage: python %s gtf old-to-new-names outfilename [-genePredictionToGeneName]' % argv[0]
        print '       old-to-new-names format: old <tab> new'
        print '       [-genePredictionToGeneName] will convert attritbutes from "GenePrediction NV18831-RA" to "gene_name "NV18831-RA"; transcript_name "NV18831-RA"; " and keep notes as formatted'
        sys.exit(1)

    gtf = argv[1]
    old_to_new = argv[2]
    outfilename = argv[3]

    doGPtoGN=False
    if '-genePredictionToGeneName' in argv:
        doGPtoGN=True

    outfile = open(outfilename, 'w')

    chrDict={}
    lineslist  = open(old_to_new)
    for line in lineslist:
        fields=line.strip().split('\t')
        chrDict[fields[0]] = fields[1]

    print chrDict

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            outfile.write(line)
            continue
        fields=line.strip().split('\t')
        if doGPtoGN:
            GP = fields[8].split('GenePrediction ')[1].split(' ')[0]
            if 'Note "' in fields[8]:
                note = fields[8].split('Note "')[1].split('"')[0]
                outline = chrDict[fields[0]] + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + '0' + '\t' + fields[6] + '\t' + fields[7] + '\t' + 'gene_id "' + GP + '-gene"; transcript_id "' + GP + '"; GenePrediction "' + GP + '"; Note "' + note + '";' + '\n'
            else:
                outline = chrDict[fields[0]] + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + '0' + '\t' + fields[6] + '\t' + fields[7] + '\t' + 'gene_id "' + GP + '-gene"; transcript_id "' + GP + '"; GenePrediction "' + GP + '";' + '\n'
        else:
            outline = chrDict[fields[0]] + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5] + '\t' + fields[6] + '\t' + fields[7] + '\t' + fields[8] + '\n'
        outfile.write(outline)

    outfile.close()
        
if __name__ == '__main__':
    main(sys.argv)

