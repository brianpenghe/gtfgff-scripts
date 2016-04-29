##################################
#                                #
# Last modified 04/03/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import numpy
from sets import Set

def main(argv):

    if len(argv) < 5:
        print 'usage: python %s input GeneID/NameFieldID name|geneID gtf outfilename' % argv[0]
        sys.exit(1)

    input = argv[1]
    fieldID = int(argv[2])
    nameOrID = argv[3]
    gtf = argv[4]
    outfilename = argv[5]

    WantedDict={}

    listoflines = open(input)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        ID = fields[fieldID]
        WantedDict[ID] = {}

    outfile = open(outfilename, 'w')

    listoflines = open(gtf)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if nameOrID == 'name':
            geneID = fields[8].split('gene_name "')[1].split('";')[0]
        if nameOrID == 'geneID':
            geneID = fields[8].split('gene_id "')[1].split('";')[0]
        if WantedDict.has_key(geneID):
            outfile.write(line)
    outfile.close()

if __name__ == '__main__':
    main(sys.argv)
