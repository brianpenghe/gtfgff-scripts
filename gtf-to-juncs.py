##################################
#                                #
# Last modified 11/14/2011       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s gtf outputfilename ' % argv[0]
        sys.exit(1)
    
    inputfilename = argv[1]
    outfilename = argv[2]

    TranscriptDict={}
    JunctionsDict={}
    
    linelist = open(inputfilename)
    for line in linelist:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        if fields[2]!='exon':
            continue
        chr=fields[0]
        start=int(fields[3])
        stop=int(fields[4])
        strand=fields[6]
        transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
        geneID=fields[8].split('gene_id "')[1].split('";')[0]
        if 'transcript_name' in fields[8]:
            transcriptName=fields[8].split('transcript_name "')[1].split('";')[0]
        else:
            transcriptName=transcriptID
        if 'gene_name' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName=geneID
        if TranscriptDict.has_key((geneID,geneName,transcriptID,transcriptName)):
            pass
        else:
            TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]=[]
        TranscriptDict[(geneID,geneName,transcriptID,transcriptName)].append((chr,start,stop,strand))
        
    for (geneID,geneName,transcriptID,transcriptName) in TranscriptDict.keys():
        TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]=list(Set(TranscriptDict[(geneID,geneName,transcriptID,transcriptName)]))
        TranscriptDict[(geneID,geneName,transcriptID,transcriptName)].sort()
        chr=TranscriptDict[(geneID,geneName,transcriptID,transcriptName)][0][0]
        strand=TranscriptDict[(geneID,geneName,transcriptID,transcriptName)][0][3]
        for i in range(len(TranscriptDict[(geneID,geneName,transcriptID,transcriptName)])-1):
            junction=(chr,TranscriptDict[(geneID,geneName,transcriptID,transcriptName)][i][2],TranscriptDict[(geneID,geneName,transcriptID,transcriptName)][i+1][1],strand)
            if JunctionsDict.has_key(junction):
                pass
            else:
                JunctionsDict[junction]={}
                JunctionsDict[junction]['geneNames']=[]
                JunctionsDict[junction]['geneIDs']=[]
                JunctionsDict[junction]['transcriptNames']=[]
                JunctionsDict[junction]['transcriptIDs']=[]
            JunctionsDict[junction]['geneNames'].append(geneName)
            JunctionsDict[junction]['geneIDs'].append(geneID)
            JunctionsDict[junction]['transcriptNames'].append(transcriptName)
            JunctionsDict[junction]['transcriptIDs'].append(transcriptID)
            
    print 'junctions found:', len(JunctionsDict.keys())
    
    JunctionsList=JunctionsDict.keys()
    JunctionsList.sort()

    outfile = open(outfilename, 'w')
    outline = '#chr\tleft\tright\tstrand\tGeneID(s)\tGeneName(s)\tTranscriptID(s)\tTranscriptName(s)\t'
    outfile.write(outline+'\n')

    for (chr,left,right,strand) in JunctionsList:
        outline=chr+'\t'+str(left-1)+'\t'+str(right-1)+'\t'+strand + '\t'
        if left > right:
             print 'left position larger than right:', outline
             continue
        JunctionsDict[(chr,left,right,strand)]['geneNames']=list(Set(JunctionsDict[(chr,left,right,strand)]['geneNames']))
        JunctionsDict[(chr,left,right,strand)]['geneIDs']=list(Set(JunctionsDict[(chr,left,right,strand)]['geneIDs']))
        JunctionsDict[(chr,left,right,strand)]['transcriptNames']=list(Set(JunctionsDict[(chr,left,right,strand)]['transcriptNames']))
        JunctionsDict[(chr,left,right,strand)]['transcriptIDs']=list(Set(JunctionsDict[(chr,left,right,strand)]['transcriptIDs']))
        for geneID in JunctionsDict[(chr,left,right,strand)]['geneIDs']:
             outline = outline + geneID + ','
        outline = outline[0:-1] + '\t'
        for geneName in JunctionsDict[(chr,left,right,strand)]['geneNames']:
             outline = outline + geneName + ','
        outline = outline[0:-1] + '\t'
        for transcriptID in JunctionsDict[(chr,left,right,strand)]['transcriptIDs']:
             outline = outline + transcriptID + ','
        outline = outline[0:-1] + '\t'
        for transcriptName in JunctionsDict[(chr,left,right,strand)]['transcriptNames']:
             outline = outline + transcriptName + ','
        outline = outline[0:-1]
        outfile.write(outline+'\n')
   
    outfile.close()
   
if __name__ == '__main__':
    main(sys.argv)
