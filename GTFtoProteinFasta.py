##################################
#                                #
# Last modified 10/14/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import math
import random
import string
from cistematic.core.geneinfo import geneinfoDB
from cistematic.genomes import Genome
from sets import Set

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for i in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-i-1]]
    return sequence

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s genome gtf outfilename [-spliced] [-class_code symbol]' % argv[0]
        print '     this script will output the translation of all three possible reading frames; stop codons will be converted to a .'
        sys.exit(1)

    genome = argv[1]
    gtf=argv[2]
    outputfilename = argv[3]

    doSpliced=False
    if '-spliced' in argv:
        doSpliced=True
        print 'will only look at transciprs with more than one exon'

    doClassCode=False
    if '-class_code' in argv:
        doClassCode=True
        class_code=argv[argv.index('-class_code')+1]
        print 'will only look at transciprs if class code', class_code

    CodonDict={'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
               'UUA':'L', 'UUG':'L', 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
               'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
               'AAA':'K', 'AAG':'K',
               'AAU':'N', 'AAC':'N',
               'AUG':'M',
               'GAU':'D', 'GAC':'D',
               'UUU':'F', 'UUC':'F',
               'UGU':'C', 'UGC':'C',
               'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
               'CAA':'Q', 'CAG':'Q',
               'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S', 'AGU':'S', 'AGC':'S',
               'GAA':'E', 'GAG':'E',
               'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
               'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
               'UGG':'W',
               'CAU':'H', 'CAC':'H',
               'UAU':'Y', 'UAC':'Y',
               'AUU':'I', 'AUC':'I', 'AUA':'I',
               'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
               'START':'AUG',
               'UAA':'.',
               'UGA':'.',
               'UAG':'.'}

    outfile = open(outputfilename, 'w')

    hg = Genome(genome)

    j=0
    lineslist = open(gtf)
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
        if doClassCode:
            if 'class_code "' in fields[8]:
                cc = fields[8].split('class_code "')[1].split('";')[0]
                if cc != class_code:
                    continue
            else:
                continue
        if 'transcript_name "' in fields[8]:
            TranscriptID=fields[8].split('transcript_name "')[1].split('";')[0]
        else:
            TranscriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if TranscriptDict.has_key(TranscriptID):
            pass
        else:
            TranscriptDict[TranscriptID]=[]
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        orientation=fields[6]
        TranscriptDict[TranscriptID].append((chr,left,right,orientation))

    g=0 
    print 'Found', len(TranscriptDict.keys()), 'transcripts'
    for transcript in TranscriptDict.keys():
        g+=1
        if g % 1000 == 0:
            print g, 'transcripts sequences processed'
        TranscriptDict[transcript] = list(Set(TranscriptDict[transcript]))
        if doSpliced:
            if len(TranscriptDict[transcript]) == 1:
                del TranscriptDict[transcript]
                continue
        sequence=''
        leftEnds=[]
        rightEnds=[]
        orientation = TranscriptDict[transcript][0][3]
        TranscriptDict[transcript].sort()
        if orientation=='+':
            for (chr,left,right,orientation) in TranscriptDict[transcript]:
                try:
                    sequence=sequence+hg.sequence(chr[3:len(chr)],left,right-left)
                except:
                    print "can't retrieve sequence", chr,left,right,orientation
                    for p in range(left,right-left):
                        try:
                            sequence=sequence+hg.sequence(chr[3:len(chr)],p,1)
                        except:
                            sequence=sequence+'N'
                            missed+=1
            sense='plus_strand'
        if orientation=='-':
            for (chr,left,right,orientation) in reversed(TranscriptDict[transcript]):
                try:
                    exonsequence=hg.sequence(chr[3:len(chr)],left-1,right-left+1)
                    sequence=sequence+getReverseComplement(exonsequence)
                except:
                    print "can not retrieve sequence", chr,left,right,orientation
                    for p in range(left-1,right-left+1):
                        try:
                            sequence=sequence+getReverseComplement(hg.sequence(chr[3:len(chr)],p,1))
                        except:
                            sequence=sequence+'N'
                            missed+=1
            sense='minus_strand'
        if orientation=='.':
            for (chr,left,right,orientation) in TranscriptDict[transcript]:
                try:
                    sequence=sequence+hg.sequence(chr[3:len(chr)],left,right-left)
                except:
                    print "can not retrieve sequence", chr,left,right,orientation
                    for p in range(left,right-left):
                        try:
                            sequence=sequence+hg.sequence(chr[3:len(chr)],p,1)
                        except:
                            sequence=sequence+'N'
                            missed+=1
            sense='unknown_strand'
        LeftEnd=TranscriptDict[transcript][0][1]
        RightEnd=TranscriptDict[transcript][-1][2]
        if orientation == '+' or orientation == '-':
            sequence = sequence.upper().replace('T','U')
            max_protein_length = len(sequence)

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame1'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(0,max_protein_length-3,3):
                if 'N' in sequence[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence[i:i+3]]
            outfile.write(protein+'\n')

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame2'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(1,max_protein_length-4,3):
                if 'N' in sequence[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence[i:i+3]]
            outfile.write(protein+'\n')

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame3'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(2,max_protein_length-5,3):
                if 'N' in sequence[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence[i:i+3]]
            outfile.write(protein+'\n')
        else:
            sequence1 = sequence.upper().replace('T','U')
            sequence2 = getReverseComplement(sequence).upper().replace('T','U')
            max_protein_length = len(sequence1)

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame1'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(0,max_protein_length-3,3):
                if 'N' in sequence1[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence1[i:i+3]]
            outfile.write(protein+'\n')

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame2'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(1,max_protein_length-4,3):
                if 'N' in sequence1[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence1[i:i+3]]
            outfile.write(protein+'\n')

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame3'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(2,max_protein_length-5,3):
                if 'N' in sequence1[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence1[i:i+3]]
            outfile.write(protein+'\n')

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame4'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(0,max_protein_length-3,3):
                if 'N' in sequence2[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence2[i:i+3]]
            outfile.write(protein+'\n')

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame5'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(1,max_protein_length-4,3):
                if 'N' in sequence2[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence2[i:i+3]]
            outfile.write(protein+'\n')

            outline='>'+transcript+':'+chr+':'+str(LeftEnd)+'-'+str(RightEnd)+'-'+sense+'::frame6'
            outfile.write(outline+'\n')
            protein = ''
            for i in range(2,max_protein_length-5,3):
                if 'N' in sequence2[i:i+3]:
                    protein = protein + '.'
                else:
                    protein = protein + CodonDict[sequence2[i:i+3]]
            outfile.write(protein+'\n')

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

