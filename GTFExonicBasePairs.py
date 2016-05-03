##################################
#                                #
# Last modified 04/26/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
from sets import Set

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s GTF outfilename [-extend bp] [-genetype biotype]' % argv[0]
        sys.exit(1)

    GTF = argv[1]
    outputfilename = argv[2]

    doExtend = False
    if '-extend' in argv:
        doExtend = True
        extension = int(argv[argv.index('-extend')+1])

    doBT = False
    if '-genetype' in argv:
        doBT = True
        BT = argv[argv.index('-genetype')+1]

    ExonDict = {}

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
        if doBT:
            biotype = fields[8].split('gene_type "')[1].split('";')[0]
            if biotype != BT:
                continue
        if fields[2]!='exon':
            continue
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        orientation=fields[6]
        if ExonDict.has_key(chr):
            pass
        else:
            ExonDict[chr] = []
        if doExtend:
            ExonDict[chr].append((chr,left-extension,right+extension,orientation))
        else:
            ExonDict[chr].append((chr,left,right,orientation))

    N = 0

    chromosomes = ExonDict.keys()
    chromosomes.sort()

    for chr in chromosomes:
        print chr
        basePairDict = {}
        for (chr,left,right,orientation) in ExonDict[chr]:
            for i in range(left,right):
                basePairDict[i]=0
        N+=len(basePairDict)
        
    outfile = open(outputfilename, 'w')
    outline='Number exonic basepairs = ' + str(N)
    print outline
    outfile.write(outline + '\n')
            
    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

