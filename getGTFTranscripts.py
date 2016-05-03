##################################
#                                #
# Last modified 10/15/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import numpy
from sets import Set

def main(argv):

    if len(argv) < 4:
        print 'usage: python %s input TranscriptIDFieldID gtf outfilename' % argv[0]
        sys.exit(1)

    input = argv[1]
    fieldID = int(argv[2])
    gtf = argv[3]
    outfilename = argv[4]

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
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if WantedDict.has_key(transcriptID):
            outfile.write(line)
    outfile.close()

if __name__ == '__main__':
    main(sys.argv)
