##################################
#                                #
# Last modified 10/17/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set

def SubChains(N):
    Subsequences = []
    for i in range(N-1):
        a = [i]
        for j in range(i+1,N):
            b = a + [j]
            Subsequences.append(b)
            a.append(j)

    return Subsequences

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s gtf1 gtf2 outputfilename' % argv[0]
        sys.exit(1)

    GTF1 = argv[1]
    GTF2 = argv[2]
    outfilename = argv[3]

    TranscriptDict1={}
    linelist=open(GTF1)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        chr = fields[0]
        left = int(fields[3])
        right = int(fields[4])
        strand = fields[6]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if TranscriptDict1.has_key(chr):
            pass
        else:
            TranscriptDict1[chr]={}
        if TranscriptDict1[chr].has_key(transcriptID):
            pass
        else:
            TranscriptDict1[chr][transcriptID]=[]
        TranscriptDict1[chr][transcriptID].append((chr,left,right,strand))

    print 'finished inputting', GTF1

    TranscriptDict2={}
    linelist=open(GTF2)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        chr = fields[0]
        left = int(fields[3])
        right = int(fields[4])
        strand = fields[6]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if TranscriptDict2.has_key(chr):
            pass
        else:
            TranscriptDict2[chr]={}
        if TranscriptDict2[chr].has_key(transcriptID):
            pass
        else:
            TranscriptDict2[chr][transcriptID]=[]
        TranscriptDict2[chr][transcriptID].append((chr,left,right,strand))

    print 'finished inputting', GTF2

    ChainDict1={}
    ChainDict2={}

    for chr in TranscriptDict1.keys():
        for transcriptID in TranscriptDict1[chr].keys():
            if len(TranscriptDict1[chr][transcriptID]) == 1:
                continue
            chr = TranscriptDict1[chr][transcriptID][0][0]
            strand = TranscriptDict1[chr][transcriptID][0][3]
            if ChainDict1.has_key((chr,strand)):
                pass
            else:
                ChainDict1[(chr,strand)]={}
            TranscriptDict1[chr][transcriptID].sort()
            i=0
            FullCoordinates=[]
            chain=[]
            for (chr,left,right,strand) in TranscriptDict1[chr][transcriptID]:
                i+=1
                FullCoordinates.append(left)
                FullCoordinates.append(right)
                if i == 1:
                    if len(TranscriptDict1[chr][transcriptID]) == i:
                        chain.append(left)
                    chain.append(right)
                elif (i > 1) and (i == len(TranscriptDict1[chr][transcriptID])):
                    chain.append(left)
                else:
                    chain.append(left)
                    chain.append(right)
            chain=tuple(chain)
            ChainDict1[(chr,strand)][chain]=(transcriptID,FullCoordinates)

    for chr in TranscriptDict2.keys():
        for transcriptID in TranscriptDict2[chr].keys():
            if len(TranscriptDict2[chr][transcriptID]) == 1:
                continue
            chr = TranscriptDict2[chr][transcriptID][0][0]
            strand = TranscriptDict2[chr][transcriptID][0][3]
            if ChainDict2.has_key((chr,strand)):
                pass
            else:
                ChainDict2[(chr,strand)]={}
            TranscriptDict2[chr][transcriptID].sort()
            i=0
            FullCoordinates=[]
            chain=[]
            for (chr,left,right,strand) in TranscriptDict2[chr][transcriptID]:
                i+=1
                FullCoordinates.append(left)
                FullCoordinates.append(right)
                if i == 1:
                    if len(TranscriptDict2[chr][transcriptID]) == i:
                        chain.append(left)
                    chain.append(right)
                elif (i > 1) and (i == len(TranscriptDict2[chr][transcriptID])):
                    chain.append(left)
                else:
                    chain.append(left)
                    chain.append(right)
            chain=tuple(chain)
            ChainDict2[(chr,strand)][chain]=(transcriptID,FullCoordinates)

    keys1 = ChainDict1.keys()
    keys2 = ChainDict2.keys()

    keys= keys1 + keys2
    keys = list(Set(keys))
    keys.sort()

    outfile = open(outfilename, 'w')
    outfile.write("#TranscriptID1\tNumExons\tTranscriptID2\tNumExons\t5'end_distance\t3'end_distance\tPartialOverlapTranscript1\tPartialOverlapExonsTranscript1\tPartialOverlapTranscript2\tPartialOverlapExonsTranscript2" + '\n')

    for (chr,strand) in keys:
        if ChainDict1.has_key((chr,strand)):
            pass
        else:
            for chain in ChainDict2[(chr,strand)].keys():
                outline = '-\t-\t' + ChainDict2[(chr,strand)][chain][0] + '\t' + str(len(ChainDict2[(chr,strand)][chain][1])/2) + '\t-\t-' + '\t-\t-' + '\t-\t-'
                outfile.write(outline + '\n')
            continue
        for chain in ChainDict1[(chr,strand)].keys():
            if ChainDict2.has_key((chr,strand)) and ChainDict2[(chr,strand)].has_key(chain):
                outline = ChainDict1[(chr,strand)][chain][0] + '\t' + str(len(ChainDict1[(chr,strand)][chain][1])/2) + '\t' + ChainDict2[(chr,strand)][chain][0] + '\t' + str(len(ChainDict2[(chr,strand)][chain][1])/2)
                if len(ChainDict1[(chr,strand)][chain][1])/2 != len(ChainDict2[(chr,strand)][chain][1])/2:
                    print ChainDict1[(chr,strand)][chain][0], ChainDict2[(chr,strand)][chain][0]
                    print ChainDict1[(chr,strand)][chain][1], ChainDict2[(chr,strand)][chain][1]
                if strand == '+':
                    dist5 = ChainDict1[(chr,strand)][chain][1][0] - ChainDict2[(chr,strand)][chain][1][0]
                    dist3= ChainDict1[(chr,strand)][chain][1][-1] - ChainDict2[(chr,strand)][chain][1][-1]
                if strand == '-':
                    dist5 = ChainDict2[(chr,strand)][chain][1][-1] - ChainDict1[(chr,strand)][chain][1][-1]
                    dist3 = ChainDict2[(chr,strand)][chain][1][0] - ChainDict1[(chr,strand)][chain][1][0]
                del ChainDict2[(chr,strand)][chain]
                outline = outline + '\t' + str(dist5) + '\t' + str(dist3) + '\t-\t-' + '\t-\t-'
                outfile.write(outline + '\n')
            else:
                HasPartial = False
                partials = []
                SCs = SubChains(len(chain))
                for SC in SCs:
                    subchain = []
                    for i in SC:
                        subchain.append(chain[i])
                    subchain = tuple(subchain)
                    if ChainDict2.has_key((chr,strand)) and ChainDict2[(chr,strand)].has_key(subchain):
                        partials.append(subchain)
                        HasPartial = True
                if HasPartial:
                    for partial in partials:
                        outline = ChainDict1[(chr,strand)][chain][0] + '\t' + str(len(ChainDict1[(chr,strand)][chain][1])/2) + '\t-\t-'  + '\t-\t-' + '\t' + ChainDict2[(chr,strand)][partial][0] + '\t' + str(len(ChainDict2[(chr,strand)][partial][1])/2) + '\t-\t-'
                        outfile.write(outline + '\n')
                else:
                    outline = ChainDict1[(chr,strand)][chain][0] + '\t' + str(len(ChainDict1[(chr,strand)][chain][1])/2) + '\t-\t-\t-\t-\t-\t-\t-\t-'
                    outfile.write(outline + '\n')
        if ChainDict2.has_key((chr,strand)):
            pass
        else:
            continue
        for chain in ChainDict2[(chr,strand)].keys():
            HasPartial = False
            partials = []
            SCs = SubChains(len(chain))
            for SC in SCs:
                subchain = []
                for i in SC:
                    subchain.append(chain[i])
                subchain = tuple(subchain)
                if ChainDict1.has_key((chr,strand)) and ChainDict1[(chr,strand)].has_key(subchain):
                    partials.append(subchain)
                    HasPartial = True
            if HasPartial:
                for partial in partials:
                    outline = '-\t-\t' + ChainDict2[(chr,strand)][chain][0] + '\t' + str(len(ChainDict2[(chr,strand)][chain][1])/2) + '\t-\t-'  + '\t-\t-' + '\t' + ChainDict1[(chr,strand)][partial][0] + '\t' + str(len(ChainDict1[(chr,strand)][partial][1])/2)
                    outfile.write(outline + '\n')
            else:
                outline = '-\t-\t' + ChainDict2[(chr,strand)][chain][0] + '\t' + str(len(ChainDict2[(chr,strand)][chain][1])/2) + '\t-\t-'  + '\t-\t-' + '\t-\t-'
                outfile.write(outline + '\n')
        
    outfile.close()
   
if __name__ == '__main__':
    main(sys.argv)
