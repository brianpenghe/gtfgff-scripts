##################################
#                                #
# Last modified 11/06/2011       # 
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
        print 'usage: python %s gtf length outfilename' % argv[0]
        sys.exit(1)

    gtf = argv[1]
    NewLength = int(argv[2])
    outfilename = argv[3]

    TranscriptDict={}
    ExonDict={}

    lineslist  = open(gtf)
    for line in lineslist:
        if line[0]=='#':
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        chr=fields[0]
        strand=fields[6]
        left=int(fields[3])
        right=int(fields[4])
        transcriptID = fields[8].split('transcript_id "')[1].split('"')[0]
        geneID = fields[8].split('gene_id "')[1].split('"')[0]
        geneName = fields[8].split('gene_name "')[1].split('"')[0]
        transcriptName = fields[8].split('transcript_name "')[1].split('"')[0]
        if TranscriptDict.has_key((transcriptID,transcriptName,geneID,geneName)):
            pass
        else:
            TranscriptDict[(transcriptID,transcriptName,geneID,geneName)]=[]
        TranscriptDict[(transcriptID,transcriptName,geneID,geneName)].append((chr,left,right,strand))
        if ExonDict.has_key((chr,left,right,strand)):
            pass
        else:
            ExonDict[(chr,left,right,strand)]=[]
        ExonDict[(chr,left,right,strand)].append((fields,transcriptID,transcriptName,geneID,geneName))

    FinalExonDict={}
    FinalTranscriptDict={}

    print 'finished inputting annotation'

    for (transcriptID,transcriptName,geneID,geneName) in TranscriptDict.keys():
        key = (transcriptID,transcriptName,geneID,geneName)
        TranscriptDict[key].sort()
        if transcriptID == 'ENSMUST00000020408':
             print TranscriptDict[key]
        strand = TranscriptDict[key][0][3]
        TotalLength = 0
        for (chr,left,right,strand) in TranscriptDict[key]:
            TotalLength+=(right-left)
        FinalTranscriptDict[key]=[]
        if TotalLength <= NewLength:
            for (chr,left,right,strand) in TranscriptDict[key]:
                FinalTranscriptDict[key].append((chr,left,right,strand))
                FinalExonDict[(chr,left,right,strand)]=ExonDict[(chr,left,right,strand)]
            TranscriptDict[key].sort()
            continue
        CurrentLength=0
        if strand == '+':
            TranscriptDict[key].reverse()
            for (chr,left,right,strand) in TranscriptDict[key]:
                if CurrentLength > NewLength:
                    TranscriptDict[key].sort()
                    continue
                CurrentLength += (right-left)
                if CurrentLength > NewLength:
                    newleft = left + (CurrentLength - NewLength)
                    FinalTranscriptDict[key].append((chr,newleft,right,strand))
                    FinalExonDict[(chr,newleft,right,strand)]=ExonDict[(chr,left,right,strand)]
                    break
                else:
                    FinalTranscriptDict[key].append((chr,left,right,strand))
                    FinalExonDict[(chr,left,right,strand)]=ExonDict[(chr,left,right,strand)]
        if strand == '-':
            for (chr,left,right,strand) in TranscriptDict[key]:
                if transcriptID == 'ENSMUST00000020408':
                     print transcriptID,chr,left,right,strand
                if CurrentLength > NewLength:
                    TranscriptDict[key].sort()
                    continue
                CurrentLength += (right-left)
                if CurrentLength > NewLength:
                    newright = right - (CurrentLength - NewLength)
                    FinalTranscriptDict[key].append((chr,left,newright,strand))
                    FinalExonDict[(chr,left,newright,strand)]=ExonDict[(chr,left,right,strand)]
                    if transcriptID == 'ENSMUST00000020408':
                        print transcriptID,chr,left,right,strand,chr,left,newright,strand,FinalExonDict[(chr,left,newright,strand)]
                    break
                else:
                    FinalTranscriptDict[key].append((chr,left,right,strand))
                    FinalExonDict[(chr,left,right,strand)]=ExonDict[(chr,left,right,strand)]
                    if transcriptID == 'ENSMUST00000020408':
                        print transcriptID,chr,left,right,strand,chr,left,right,strand,FinalExonDict[(chr,left,right,strand)]
        TranscriptDict[key].sort()
                    
    print 'finished shortening transcripts'

    print len(FinalTranscriptDict.keys()), len(TranscriptDict.keys())

    outfile = open(outfilename, 'w')

    for (transcriptID,transcriptName,geneID,geneName) in FinalTranscriptDict.keys():
        if transcriptID == 'ENSMUST00000020408':
            print FinalTranscriptDict[(transcriptID,transcriptName,geneID,geneName)]
            print FinalExonDict[FinalTranscriptDict[(transcriptID,transcriptName,geneID,geneName)][0]]
        FinalTranscriptDict[(transcriptID,transcriptName,geneID,geneName)].sort()
        for (chr,left,right,strand) in FinalTranscriptDict[(transcriptID,transcriptName,geneID,geneName)]:
            (fields,transcriptID2,transcriptName2,geneID2,geneName2) = FinalExonDict[(chr,left,right,strand)][0]
            outline = chr + '\t' + fields[1] + '\texon\t' + str(left) + '\t' + str(right) + '\t' + fields[5] + '\t' + strand + '\t' + fields[7] + '\t' + 'gene_id "' + geneID + '"; ' + 'transcript_id "' + transcriptID + '"; ' + 'gene_name "' + geneName + '"; ' + 'transcript_name "' + transcriptName + '"; ' 
            outfile.write(outline + '\n')

#    for (transcriptID,transcriptName,geneID,geneName) in FinalTranscriptDict.keys():
#        if FinalTranscriptDict.has_key((transcriptID,transcriptName,geneID,geneName)):
#            pass
#        else:
#            continue
#        for (chr,left,right,strand) in FinalTranscriptDict[(transcriptID,transcriptName,geneID,geneName)]:
#            for (fields,transcriptID2,transcriptName2,geneID2,geneName2) in FinalExonDict[(chr,left,right,strand)]:
#                if FinalTranscriptDict.has_key((transcriptID2,transcriptName2,geneID2,geneName2)):
#                    if transcriptID != transcriptID2 and FinalTranscriptDict[(transcriptID,transcriptName,geneID,geneName)] == FinalTranscriptDict[(transcriptID2,transcriptName2,geneID2,geneName2)]:
#                        print FinalTranscriptDict[(transcriptID,transcriptName,geneID,geneName)]
#                        print FinalTranscriptDict[(transcriptID2,transcriptName2,geneID2,geneName2)]
#                        print '..................'
#                        del FinalTranscriptDict[(transcriptID2,transcriptName2,geneID2,geneName2)]
#
#    print len(FinalTranscriptDict.keys()), 'left after removing duplicates'
#
#    keys = FinalExonDict.keys() 
#    keys.sort()
#
#    for (chr,left,right,strand) in keys:
#        FinalExonDict[(chr,left,right,strand)] = list(Set(FinalExonDict[(chr,left,right,strand)]))
#        for (fields,transcriptID,transcriptName,geneID,geneName) in FinalExonDict[(chr,left,right,strand)]:
#            if transcriptID == 'ENSMUST00000020408':
#                 print TranscriptDict[(transcriptID,transcriptName,geneID,geneName)]
#                 print FinalTranscriptDict[(transcriptID,transcriptName,geneID,geneName)]
#                 print FinalExonDict[(chr,left,right,strand)]
#            if FinalTranscriptDict.has_key((transcriptID,transcriptName,geneID,geneName)):
#                 outline = chr + '\t' + fields[1] + '\texon\t' + str(left) + '\t' + str(right) + '\t' + fields[5] + '\t' + strand + '\t' + fields[7] + '\t' + 'gene_id "' + geneID + '"; ' + 'transcript_id "' + transcriptID + '"; ' + 'gene_name "' + geneName + '"; ' + 'transcript_name "' + transcriptName + '"; ' 
#                 outfile.write(outline + '\n')

    outfile.close()
        
if __name__ == '__main__':
    main(sys.argv)

