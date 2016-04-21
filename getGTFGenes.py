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

def run():

    if len(sys.argv) < 5:
        print 'usage: python %s input GeneID/NameFieldID name|geneID gtf outfilename' % sys.argv[0]
        sys.exit(1)

    input = sys.argv[1]
    fieldID = int(sys.argv[2])
    nameOrID = sys.argv[3]
    gtf = sys.argv[4]
    outfilename = sys.argv[5]

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

run()
