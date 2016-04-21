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

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s input TranscriptIDFieldID gtf outfilename' % sys.argv[0]
        sys.exit(1)

    input = sys.argv[1]
    fieldID = int(sys.argv[2])
    gtf = sys.argv[3]
    outfilename = sys.argv[4]

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

run()
