##################################
#                                #
# Last modified 01/29/2015       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
from sets import Set

try:
	import psyco
	psyco.full()
except:
	pass

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s inputfilename outfilename' % sys.argv[0]
        sys.exit(1)

    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]

    outfile = open(outputfilename, 'w')

    TranscriptDict={}
    listoflines = open(inputfilename)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        gID=fields[8].split('gene_id "')[1].split('"')[0]
        if 'gene_name "' in fields[8]:
            gname=fields[8].split('gene_name "')[1].split('"')[0]
        else:
            gname = gID
        tID=fields[8].split('transcript_id "')[1].split('"')[0]
        if 'transcript_name "' in fields[8]:
            tname=fields[8].split('transcript_name "')[1].split('"')[0]
        else:
            tname = tID
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        strand=fields[6]
        if TranscriptDict.has_key((gID,gname,tID,tname)):
            pass
        else:
            TranscriptDict[(gID,gname,tID,tname)]={}
            TranscriptDict[(gID,gname,tID,tname)]['strand']=strand
            TranscriptDict[(gID,gname,tID,tname)]['chr']=chr
            TranscriptDict[(gID,gname,tID,tname)]['coordinates']=[]
        TranscriptDict[(gID,gname,tID,tname)]['coordinates'].append(left)
        TranscriptDict[(gID,gname,tID,tname)]['coordinates'].append(right)

    outline='#chr\tleft\tright\tstrand\tgeneID\tgeneName\ttranscriptID\ttranscriptName'
    outfile.write(outline+'\n')

    keys = TranscriptDict.keys()
    keys.sort()

    for (gID,gname,tID,tname) in keys:
        strand=TranscriptDict[(gID,gname,tID,tname)]['strand']
        left=min(TranscriptDict[(gID,gname,tID,tname)]['coordinates'])
        right=max(TranscriptDict[(gID,gname,tID,tname)]['coordinates'])
        chr=TranscriptDict[(gID,gname,tID,tname)]['chr']
        outline=chr+'\t'+str(left)+'\t'+str(right)+'\t'+strand + '\t' + gID + '\t' + gname + '\t' + tID + '\t' + tname
        outfile.write(outline+'\n')
            
    outfile.close()

run()

