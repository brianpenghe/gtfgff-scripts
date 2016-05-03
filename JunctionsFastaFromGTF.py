##################################
#                                #
# Last modified 07/01/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence


def main(argv):

    if len(argv) < 4:
        print 'usage: python %s fasta GTF span outputfilename' % argv[0]
        print '\t this script will take a GTF file and output all circularized junctions, i.e. sequences joining a splice sites to all acceptors upstream of it in the transcript'
        print '\t the span parameter referse to the length of sequence on each side, i.e. if you have 75bp reads, you would want to use a span around 60 for stringency purposes'
        sys.exit(1)

    fasta = argv[1]
    GTF = argv[2]
    span = int(argv[3])
    outfilename = argv[4]

    inputdatafile = open(fasta)
    SequenceDict={}
    sequence = ''
    for line in inputdatafile:
        if line[0]=='>':
            if sequence != '':
                sequence = ''.join(sequence)
                SequenceDict[chr]=sequence
            chr = line.strip().split('>')[1]
            print chr
            sequence=[]
        else:
            sequence.append(line.strip())   
    sequence = ''.join(sequence)
    SequenceDict[chr]=sequence
   
    j=0
    lineslist = open(GTF)
    TranscriptDict={}
    for line in lineslist:
        j+=1
        if j % 100000 == 0:
            print j, 'lines processed'
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        if 'gene_name "' in fields[8]:
            gene=fields[8].split('gene_name "')[1].split('";')[0]
        else:
            gene=fields[8].split('gene_id "')[1].split('";')[0]
        if 'transcript_name "' in fields[8]:
            transcript=fields[8].split('transcript_name "')[1].split('";')[0]
        else:
            transcript=fields[8].split('transcript_id "')[1].split('";')[0]
        ID = (gene,transcript)
        if TranscriptDict.has_key(ID):
            pass
        else:
            TranscriptDict[ID]=[]
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        strand=fields[6]
        TranscriptDict[ID].append((chr,left,right,strand))

    outfile=open(outfilename,'w')

    JunctionSequenceDict={}

    g=0 
    print 'found', len(TranscriptDict.keys()), 'transcripts'
    for (gene,transcript) in TranscriptDict.keys():
        g+=1
        if g % 10000 == 0:
            print g, 'transcripts sequences processed'
        leftEnds=[]
        rightEnds=[] 
        TranscriptDict[(gene,transcript)].sort()
        orientation = TranscriptDict[(gene,transcript)][0][3]
        JunctionPoints=[]
        sequence = ''
        if orientation=='+':
            for (chr,left,right,strand) in TranscriptDict[(gene,transcript)]:
                sequence = sequence + SequenceDict[chr][left:right]
                JunctionPoints.append(len(sequence))
        if orientation=='-':
            TranscriptDict[(gene,transcript)].reverse()
            for (chr,left,right,strand) in TranscriptDict[(gene,transcript)]:
                sequence = sequence + getReverseComplement(SequenceDict[chr][left-1:right])
                JunctionPoints.append(len(sequence))
        NumJunctions = len(JunctionPoints)-1
        if NumJunctions == 0:
            continue
        JunctionPoints.append(0)
        JunctionPoints.sort()
        for j in JunctionPoints[1:-1]:
            JuncSequence = sequence[j-span:j+span]
            outline = '>' + gene + ':' + transcript + ':JunctionPoint:' + str(j)  
            outfile.write(outline + '\n')
            outfile.write(JuncSequence + '\n')

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

