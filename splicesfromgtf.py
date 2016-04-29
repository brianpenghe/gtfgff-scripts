##################################
#                                #
# Last modified 07/19/2010        # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s  gtf_filename outputfilename ' % argv[0]
        sys.exit(1)
    
    inputfilename = argv[1]
    outfilename = argv[2]

    linelist  = open(inputfilename)
    transcriptDict={}
    for line in linelist:
        if line.startswith('#'):
            continue
        fields=line.split('\n')[0].split('\t')
        if fields[2]!='exon':
            continue
        transcriptID=fields[8].split('transcript_id "')[1].split('"')[0]
        chr=fields[0]
        start=int(fields[3])
        stop=int(fields[4])
        orientation=fields[6]
        if transcriptDict.has_key(transcriptID):
            transcriptDict[transcriptID].append((chr,start,stop,orientation))
        else:
            transcriptDict[transcriptID]=[]
            transcriptDict[transcriptID].append((chr,start,stop,orientation))
    splicesList=[]
    outfile = open(outfilename, 'w')
    for transcriptID in transcriptDict.keys():
        transcriptDict[transcriptID].sort()
        elements = len(transcriptDict[transcriptID])
        if elements == 1:
            continue
        (currentChr,currentStart,currentStop,orientation)=transcriptDict[transcriptID][0]
        for i in range(1,elements):
            (chr,start,stop,orientation)=transcriptDict[transcriptID][i]
            if currentStop+1<start:
                splice=chr+':'+str(currentStop)+'-'+str(start)+orientation
                splicesList.append(splice)
            (currentChr,currentStart,currentStop,orientation)=transcriptDict[transcriptID][i]
 
    FinalSplicesList = list(Set(splicesList))
    FinalSplicesList.sort()
    for splice in FinalSplicesList:
        outfile.write(splice+'\n')
   
    outfile.close()

if __name__ == '__main__':
    main(sys.argv)
