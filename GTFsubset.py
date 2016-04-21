##################################
#                                #
# Last modified 09/02/2011       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
from sets import Set

def run():

    if len(sys.argv) < 5:
        print 'usage: python %s gtf chr left right outfilename' % sys.argv[0]
        print '     Note: only elements entirely in the region will be outputed'
        sys.exit(1)

    gtf = sys.argv[1]
    TheChr = sys.argv[2]
    outfilename = sys.argv[5]
    start = int(sys.argv[3])
    stop = int(sys.argv[4])

    outfile = open(outfilename, 'w')

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        if chr == TheChr:
            if (left > start and left < stop) and (right > start and right < stop):
                outfile.write(line)

    outfile.close()
        
run()

