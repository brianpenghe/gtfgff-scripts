##################################
#                                #
# Last modified 03/11/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
from sets import Set

def main(argv):

    if len(argv) < 5:
        print 'usage: python %s gtf wanted_list IDfield <genes | transcripts> outfilename' % argv[0]
        print '     Note: only elements entirely in the region will be outputed'
        sys.exit(1)

    gtf = argv[1]
    wanted = argv[2]
    fieldID = int(argv[3])
    GenesOrTranscripts = argv[4]

    outfile = open(argv[5], 'w')

    WantedDict={}

    lineslist  = open(wanted)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        WantedDict[fields[fieldID]] = 0

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            outfile.write(line)
            continue
        fields=line.strip().split('\t')
        if GenesOrTranscripts == 'genes':
             ID = fields[8].split('gene_id "')[1].split('"')[0]
        if GenesOrTranscripts == 'transcripts':
             ID = fields[8].split('transcript_id "')[1].split('"')[0]
        if WantedDict.has_key(ID):
             outfile.write(line)

    outfile.close()
        
if __name__ == '__main__':
    main(sys.argv)

