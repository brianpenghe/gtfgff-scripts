##################################
#                                #
# Last modified 07/28/2010       # 
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
        print 'usage: python %s inputfilename outfilename [-only middle | first | last]' % sys.argv[0]
        sys.exit(1)

    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]
    doSkip=False
    doMiddle=False
    doFirst=False
    doLast=False
    if '-only' in sys.argv:
        doSkip=True
        if sys.argv[sys.argv.index('-only') + 1]=='middle':
            doMiddle=True
            print 'Will skip first exons'
        if sys.argv[sys.argv.index('-only') + 1]=='first':
            doFirst=True
            print 'Will only output first exons (skipping exons that are both first and middle)'
        if sys.argv[sys.argv.index('-only') + 1]=='last':
            doLast=True
            print 'Will only output last exons (skipping exons that are both last and middle)'
        ExonDict={}
        TranscriptDict={}
        if not doMiddle and not doFirst and not doLast:
            print 'invalid arguments for exon selection, exiting'
            sys.exit(1)

    outfile = open(outputfilename, 'w')
    outline='#Transcript_ID\tchr\tstart\tstop\tstrand'

    i=0
    listoflines = open(inputfilename)
    for line in listoflines:
        if line.startswith('#'):
            continue
        i+=1
        if i % 1000000 == 0:
            print i, 'lines processed'
        fields=line.strip().split('\t')
        name=fields[8].split('transcript_id "')[1].split('"')[0]
        chr=fields[0]
        start=fields[3]
        end=fields[4]
        strand=fields[6]
        if i==1:
            if 'conf_lo' in line:
                outline=outline+'\tconf_lo_RPKM\tconf_high_RPKM'
            outfile.write(outline+'\n')
        if fields[2]!='exon':
            continue
        if doSkip:
            transcript_name=fields[8].strip().split('transcript_name "')[1].split('";')[0] 
            if TranscriptDict.has_key(transcript_name):
                pass
            else:
                TranscriptDict[transcript_name]=[]
            TranscriptDict[transcript_name].append((chr,start,end,strand))
            ExonDict[(chr,start,end,strand)]=[]
        else:
            outline=name+'\t'+chr+'\t'+start+'\t'+end+'\t'+strand
            if 'conf_lo' in line:
                loRPKM=fields[8].split('conf_lo "')[1].split('"')[0]
                hiRPKM=fields[8].split('conf_hi "')[1].split('"')[0]
                outline=outline+'\t'+loRPKM+'\t'+hiRPKM
            outfile.write(outline+'\n')
        
    if doSkip: 
        for transcript_name in TranscriptDict.keys():
            TranscriptDict[transcript_name].sort()
            if TranscriptDict[transcript_name][0][3]=='+':
                ExonDict[TranscriptDict[transcript_name][0]].append('First')
                ExonDict[TranscriptDict[transcript_name][-1]].append('Last')
            if TranscriptDict[transcript_name][0][3]=='-':
                ExonDict[TranscriptDict[transcript_name][-1]].append('First')
                ExonDict[TranscriptDict[transcript_name][0]].append('Last')
            if len(TranscriptDict[transcript_name])>2:
                for (chr,start,end,strand) in TranscriptDict[transcript_name][1:-1]:
                    ExonDict[(chr,start,end,strand)].append('Middle')
        ExonKeys=ExonDict.keys()
        ExonKeys.sort()
        for (chr,start,end,strand) in ExonKeys:
            ExonDict[(chr,start,end,strand)]=list(Set(ExonDict[(chr,start,end,strand)]))
            if doMiddle and ExonDict[(chr,start,end,strand)]==['Middle']:
                outline=name+'\t'+chr+'\t'+start+'\t'+end+'\t'+strand
                outfile.write(outline+'\n')
            if doFirst and ExonDict[(chr,start,end,strand)]==['First']:
                outline=name+'\t'+chr+'\t'+start+'\t'+end+'\t'+strand
                outfile.write(outline+'\n')
            if doLast and ExonDict[(chr,start,end,strand)]==['Last']:
                outline=name+'\t'+chr+'\t'+start+'\t'+end+'\t'+strand
                outfile.write(outline+'\n')
            
    outfile.close()

run()

