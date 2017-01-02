##################################
#                                #
# Last modified 2017/01/02       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s gtf radius outputfilename [-bioType biotype1,biotype2...biotype3]' % sys.argv[0]
        sys.exit(1)
    
    GTF = sys.argv[1]
    radius = int(sys.argv[2])
    outfilename = sys.argv[3]

    doBioType=False
    if '-bioType' in sys.argv:
        doBioType=True
        bioTypeDict={}
        bioTypes = sys.argv[sys.argv.index('-bioType')+1].split(',')
        print 'will only consider', bioTypes
        for b in bioTypes:
            bioTypeDict[b]=''

    TranscriptDict={}
    
    linelist = open(GTF)
    i=0
    for line in linelist:
        i+=1
        if i % 100000 == 0:
            print i, 'lines processed'
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        if doBioType:
            if bioTypeDict.has_key(fields[8].split('transcript_type "')[1].split('";')[0]):
                pass
            else:
                continue
        chr=fields[0]
        start=int(fields[3])
        stop=int(fields[4])
        strand=fields[6]
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        if 'gene_name "' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName=geneID
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        if TranscriptDict.has_key((geneName,geneID,transcriptID)):
            pass
        else:
            TranscriptDict[(geneName,geneID,transcriptID)]=[]
        TranscriptDict[(geneName,geneID,transcriptID)].append((chr,start,stop,strand))

    TSSDict={}

    outfile = open(outfilename,'w')

    outfile.write('#chr\tleft\tright\tstrand\tgeneName(s)\tgeneID(s)\ttranscript(s)\n')

    for (geneName,geneID,transcriptID) in TranscriptDict:
        TranscriptDict[(geneName,geneID,transcriptID)].sort()
        chr=TranscriptDict[(geneName,geneID,transcriptID)][0][0]
        strand=TranscriptDict[(geneName,geneID,transcriptID)][0][3]
        if strand=='+':
            TSS=(chr,TranscriptDict[(geneName,geneID,transcriptID)][0][1],strand)
        if TranscriptDict[(geneName,geneID,transcriptID)][0][3]=='-':
            TSS=(chr,TranscriptDict[(geneName,geneID,transcriptID)][-1][2],strand)
        if TSSDict.has_key(TSS):
            pass
        else:
            TSSDict[TSS]={}
            TSSDict[TSS]['transcripts']=[]
            TSSDict[TSS]['genes']=[]
        TSSDict[TSS]['transcripts'].append(transcriptID)
        TSSDict[TSS]['genes'].append((geneName,geneID))
        
    keys=TSSDict.keys()
    keys.sort()
    for (chr,TSS,strand) in keys:
        TSSDict[(chr,TSS,strand)]['genes'] = list(Set(TSSDict[(chr,TSS,strand)]['genes']))
        TSSDict[(chr,TSS,strand)]['transcripts'] = list(Set(TSSDict[(chr,TSS,strand)]['transcripts']))
        outline=chr+'\t'+str(TSS-radius)+'\t'+str(TSS+radius)+'\t'+strand+'\t'
        for (geneName,geneID) in TSSDict[(chr,TSS,strand)]['genes']:
            outline=outline+geneName+','
        outline=outline[0:-1] + '\t'
        for (geneName,geneID) in TSSDict[(chr,TSS,strand)]['genes']:
            outline=outline+geneID+','
        outline=outline[0:-1] + '\t'
        for transcriptID in TSSDict[(chr,TSS,strand)]['transcripts']:
            outline=outline+transcriptID+','
        outline=outline[0:-1]
        outfile.write(outline+'\n')
   
    outfile.close()
   
run()
