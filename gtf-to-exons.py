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

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s gtf outputfilename ' % sys.argv[0]
        sys.exit(1)
    
    inputfilename = sys.argv[1]
    outfilename = sys.argv[2]

    outfile = open(outfilename, 'w')
    outline = '#chr\tleft\tright\tstrand\tattributes\t'
    outfile.write(outline+'\n')

    linelist = open(inputfilename)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        chr=fields[0]
        start=fields[3]
        stop=fields[4]
        strand=fields[6]
        outline = chr + '\t' + start  + '\t' + stop + '\t' + strand + '\t' + fields[8]
        outfile.write(outline+'\n')
   
    outfile.close()
   
run()
