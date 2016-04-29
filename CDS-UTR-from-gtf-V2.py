##################################
#                                #
# Last modified 03/08/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s gtf outfilename' % argv[0]
        print '\tNote: the script relies on the GTF file containing CDS annotations, but explicit UTR annotations are not necessary; also, it will only look at entries in the GTF file which contain the phrase "protein_coding"'
        sys.exit(1)

    inputfilename = argv[1]
    outfilename = argv[2]

    listoflines = open(inputfilename)
    TranscriptDict={}
    i=0
    for line in listoflines:
        if line.startswith('#'):
            continue
        i+=1
        if i % 100000 == 0:
            print i, 'lines processed'
        if 'protein_coding' in line:
            pass
        else:
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'CDS' and fields[2] != 'exon':
            continue
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if 'gene_name' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName = geneID
        if 'trancsript_name' in fields[8]:
            trancsriptName=fields[8].split('trancsript_name "')[1].split('";')[0]
        else:
            transcriptName = transcriptID
        chr=fields[0]
        strand=fields[6]
        left=int(fields[3])
        right=int(fields[4])
        if TranscriptDict.has_key((geneID,geneName,transcriptID,transcriptName)):
            pass
        else:
            TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]={}
            TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['CDS']=[]
            TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['exon']=[]
        if fields[2]=='CDS':
            TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['CDS'].append((chr,left,right,strand))
        if fields[2]=='exon':
            TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['exon'].append((chr,left,right,strand))

    print 'finished unputting annotation'

    outfile = open(outfilename,'w')
    outfile.write('#chr\tleft\tright\tstrand\ttype\tgeneID\tgeneName\ttranscriptID\ttranscriptName\t')

    t=0
    for (geneID,geneName,transcriptID,transcriptName) in TranscriptDict.keys():
        TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['CDS'].sort()
        TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['exon'].sort()
        t+=1
        if t % 1000 == 0:
            print t, 'transcripts processed'
        CDS={}
        exon=[]
        for (chr,left,right,strand) in TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['CDS']:
            for i in range(left,right):
                CDS[i]=0
        for (chr,left,right,strand) in TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['exon']:
            for i in range(left,right):
                exon.append(i)
        UTR5=[]
        UTR3=[]
        CDSList=[]
        exon.sort()
        if strand == '-':
            exon.reverse()
        current = '5UTR'
        for i in exon:
            if current == '5UTR' and CDS.has_key(i):
                current = 'CDS'
            if current == 'CDS':
                if CDS.has_key(i):
                    pass
                else:
                    current = '3UTR'
            if current == '5UTR':
                UTR5.append(i)
            if current == '3UTR':
                UTR3.append(i)
            if current == 'CDS':
                CDSList.append(i)
        CDSList.sort()
        CDS = CDSList
        UTR5.sort()
        UTR3.sort()
        annotation = '\t' + geneID + '\t' + geneName + '\t' + transcriptID + '\t' + transcriptName + '\n'
#        print annotation.strip()
#        print TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['CDS']
#        print TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['exon']
        if len(CDS) == 0:
            for (chr,left,right,strand) in TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]['CDS']:
                outline = chr + '\t' + str(start) + '\t' + str(current + 1) + '\t' + strand + '\tno_CDS_annotated' + annotation
                outfile.write(outline)
            continue
        start = CDS[0]
        current = CDS[0]
        for i in CDS:
            if current != i and current + 1 != i:
                outline = chr + '\t' + str(start) + '\t' + str(current + 1) + '\t' + strand + '\tCDS' + annotation
                outfile.write(outline)
                start = i
                current = i
            else:
                current = i
        outline = chr + '\t' + str(start) + '\t' + str(current + 1) + '\t' + strand + '\tCDS' + annotation
        outfile.write(outline)
        if len(UTR5) == 0:
            pass
        else:
            start = UTR5[0]
            current = UTR5[0]
            for i in UTR5:
                if current != i and current + 1 != i:
                    outline = chr + '\t' + str(start) + '\t' + str(current + 1) + '\t' + strand + '\tUTR5' + annotation
                    outfile.write(outline)
                    start = i
                    current = i
                else:
                    current = i
            outline = chr + '\t' + str(start) + '\t' + str(current + 1) + '\t' + strand + '\tUTR5' + annotation
            outfile.write(outline)
        if len(UTR3) == 0:
            pass
        else:
            start = UTR3[0]
            current = UTR3[0]
            for i in UTR3:
                if current != i and current + 1 != i:
                    outline = chr + '\t' + str(start) + '\t' + str(current + 1) + '\t' + strand + '\tUTR3' + annotation
                    outfile.write(outline)
                    start = i
                    current = i
                else:
                    current = i
            outline = chr + '\t' + str(start) + '\t' + str(current + 1) + '\t' + strand + '\tUTR3' + annotation
            outfile.write(outline)

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

