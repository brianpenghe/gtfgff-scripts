##################################
#                                #
# Last modified 02/05/2011       # 
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
        print 'usage: python %s gtf outfilename' % argv[0]
        sys.exit(1)

    gtf = argv[1]
    outfilename = argv[2]

    outfile = open(outfilename, 'w')

    outfile.write('#chr\tleft\tright\tstrand\n')


    TranscritpDict={}
    JunctionList=[]

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields = line.strip().split('\t')
        if fields[2]!='exon':
            continue
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        strand = fields[6]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if TranscritpDict.has_key(transcriptID):
            pass
        else:
            TranscritpDict[transcriptID]=[]
        TranscritpDict[transcriptID].append((chr,left,right,strand))

    for transcriptID in TranscritpDict.keys():
        if len(TranscritpDict[transcriptID])==1:
            continue
        TranscritpDict[transcriptID]=list(Set(TranscritpDict[transcriptID]))
        TranscritpDict[transcriptID].sort()
        strand = TranscritpDict[transcriptID][0][3]
        chr = TranscritpDict[transcriptID][0][0]
        if strand == '-':
            TranscritpDict[transcriptID].reverse()
        for i in range(len(TranscritpDict[transcriptID])-1):
            left=TranscritpDict[transcriptID][i][2]
            right=TranscritpDict[transcriptID][i+1][1]
            JunctionList.append((chr,left,right,strand))
           

    JunctionList=list(Set(JunctionList))
    JunctionList.sort()

    for (chr,left,right,strand) in JunctionList:
        outline=chr+'\t'+str(left)+'\t'+str(right)+'\t'+strand
        outfile.write(outline+'\n')

    outfile.close()
        
if __name__ == '__main__':
    main(sys.argv)

