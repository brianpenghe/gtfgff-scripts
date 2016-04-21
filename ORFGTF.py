##################################
#                                #
# Last modified 09/28/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s inputfilename outfilename'
        print '\tThis script will take a GTF file and keep the CDS entries only, converting the "CDS" label to "exon"'
        sys.exit(1)

    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]

    outfile = open(outputfilename, 'w')

    listoflines = open(inputfilename)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'CDS':
            continue
        else:
            outline = fields[0] + '\t' + fields[1] + '\t' + 'exon' + '\t' + fields[3] + '\t' + fields[4] + '\t' + fields[5] + '\t' + fields[6] + '\t' + fields[7] + '\t' + fields[8]
            outfile.write(outline + '\n')
 
    outfile.close()

run()

