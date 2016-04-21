##################################
#                                #
# Last modified 01/01/2011       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s GENCODE-gtf outfileprefix' % sys.argv[0]
        print '	Note: the script will not output UTRs not in the last or first exons' % sys.argv[0]
        sys.exit(1)

    inputfilename = sys.argv[1]
    outputfilprefix = sys.argv[2]

    outfileCDS = open(outputfilprefix+'-CDS', 'w')
    outfileUTR5 = open(outputfilprefix+'-5UTR', 'w')
    outfileUTR3 = open(outputfilprefix+'-3UTR', 'w')

    listoflines = open(inputfilename)
    CDS=[]
    UTR5=[]
    UTR3=[]
    GeneDict={}
    i=0
    for line in listoflines:
        if line.startswith('#'):
            continue
        i+=1
        if i % 100000 == 0:
            print i, 'lines processed'
        fields=line.split('\t')
        if fields[2] != 'CDS' and fields[2] != 'UTR':
            continue
        ID=fields[8].split('gene_id "')[1].split('";')[0]
        chr=fields[0]
        strand=fields[6]
        left=fields[3]
        right=fields[4]
        type=fields[2]
        if GeneDict.has_key(ID):
            pass
        else:
            GeneDict[ID]=[]
        GeneDict[ID].append((type,chr,left,right,strand))

    print 'finished unputting annotation'            

    for ID in GeneDict.keys():
        GeneDict[ID].sort()
        i=0
        for (type,chr,left,rigth,strand) in GeneDict[ID]:
            i+=1
            if type=='CDS':
                CDS.append((chr,left,rigth,strand))
            if type=='UTR' and strand == '+':
                if i==1:
                    UTR5.append((chr,left,rigth,strand))
                if i==(len(GeneDict[ID])):
                    UTR3.append((chr,left,rigth,strand))
            if type=='UTR' and strand == '-':
                if i==1:
                    UTR3.append((chr,left,rigth,strand))
                if i==(len(GeneDict[ID])):
                    UTR5.append((chr,left,rigth,strand))

    CDS=list(Set(CDS))
    UTR5=list(Set(UTR5))
    UTR3=list(Set(UTR3))

    CDS.sort()
    for (chr,left,right,strand) in CDS:
        outline=chr+'\t'+left+'\t'+right+'\t'+strand
        outfileCDS.write(outline+'\n')
    outfileCDS.close()

    UTR5.sort()
    for (chr,left,right,strand) in UTR5:
        outline=chr+'\t'+left+'\t'+right+'\t'+strand
        outfileUTR5.write(outline+'\n')
    outfileUTR5.close()

    UTR3.sort()
    for (chr,left,right,strand) in UTR3:
        outline=chr+'\t'+left+'\t'+right+'\t'+strand
        outfileUTR3.write(outline+'\n')
    outfileUTR3.close()

run()

