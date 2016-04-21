##################################
#                                #
# Last modified 11/16/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set
import math

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s gtf outfile_prefix [-field1BT] [-JGI] [-GFF3]' % sys.argv[0]
        print '\tuse the -JGI option for files in the following format:'
        print '\tscaffold_1      JGI     exon    2338    2360    .       -       .       name "estExt_fgenesh1_pg.C_10001"; transcriptId 129999'
        print '\tuse the -GFF3 option for files in the following format:'
        print '\tchr1    .       exon    3704    4702    .       -       .       ID=LmjF.01.0010:exon:1;Parent=LmjF.01.0010:mRNA;constitutive=1;rank=1'
        sys.exit(1)

    GTF = sys.argv[1]
    outprefix = sys.argv[2]

    doField1BT = False
    if '-field1BT' in sys.argv:
        doField1BT = True

    GeneDict={}
    linelist=open(GTF)
    if '-JGI' in sys.argv:
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.strip().split('\t')
            chr=fields[0]
            if fields[2] == 'CDS':
                geneID=fields[8].split('name "')[1].split('";')[0]
                transcriptID=fields[8].split('proteinId ')[1].split(';')[0]
                transcriptName=transcriptID
                geneName=geneID
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT'] = 'protein_coding'
                continue
            if fields[2] != 'exon':
                continue
            left = int(fields[3])
            right = int(fields[4])
            strand = fields[6]
            geneID=fields[8].split('name "')[1].split('";')[0]
            transcriptID=fields[8].split('transcriptId ')[1].split(';')[0]
            transcriptName=transcriptID
            geneName=geneID
            if GeneDict.has_key((chr,geneID,geneName,transcriptID,transcriptName)):
                pass
            else:
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]={}
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'] = []
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT'] = '-'
            GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'].append((left,right,strand))
    elif '-GFF3' in sys.argv:
        GeneTypeDict = {}
        TranscriptParentDict = {}
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.strip().split('\t')
            if fields[2].endswith('_gene') or fields[2] == 'gene' or ((fields[2] == 'pseudogene' or fields[2] == 'pseudogenic_tRNA' or fields[2] == 'RNA') and 'Parent=' not in fields[8]):
                try:
                    gBT = fields[8].split('biotype=')[1].split(';')[0]
                except:
                    gBT = '-'
                geneID = fields[8].split('ID=')[1].split(';')[0]
                if 'external_name=' in fields[8]:
                    geneName = fields[8].split('external_name=')[1].split(';')[0]
                else:
                    geneName = geneID
                GeneTypeDict[geneID] = (geneName,gBT)
                continue
            if fields[2] == 'transcript' or fields[2] == 'miRNA' or fields[2] == 'pseudogene' or fields[2] == 'rRNA' or fields[2] == 'snRNA' or fields[2] == 'snoRNA':
                geneID = fields[8].split('Parent=')[1].split(';')[0]
                transcriptID = fields[8].split('ID=')[1].split(';')[0]
                TranscriptParentDict[transcriptID] = geneID
                continue
        linelist=open(GTF)
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.strip().split('\t')
            if fields[2] != 'exon':
                continue
            chr=fields[0]
            left = int(fields[3])
            right = int(fields[4])
            strand = fields[6]
            transcriptID=fields[8].split('Parent=')[1].split(';')[0]
            transcriptName = transcriptID
            geneID = TranscriptParentDict[transcriptID]
            (geneName,gBT) = GeneTypeDict[geneID]
            if GeneDict.has_key((chr,geneID,geneName,transcriptID,transcriptName)):
                pass
            else:
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]={}
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'] = []
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT'] = gBT
            GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'].append((left,right,strand))
    else:
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.strip().split('\t')
            if fields[2] != 'exon':
                continue
            chr=fields[0]
            left = int(fields[3])
            right = int(fields[4])
            strand = fields[6]
            geneID=fields[8].split('gene_id "')[1].split('";')[0]
            transcriptID=fields[8].split('transcript_id "')[1].split('";')[0]
            if 'transcript_name "' in fields[8]:
                transcriptName=fields[8].split('transcript_name "')[1].split('";')[0]
            else:
                transcriptName=transcriptID
            if 'gene_name "' in fields[8]:
                geneName=fields[8].split('gene_name "')[1].split('";')[0]
            else:
                geneName=geneID
            if GeneDict.has_key((chr,geneID,geneName,transcriptID,transcriptName)):
                pass
            else:
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]={}
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'] = []
            GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'].append((left,right,strand))
            if 'gene_biotype "' in fields[8]:
                gene_biotype = fields[8].split('gene_biotype "')[1].split('";')[0]
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT'] = gene_biotype
            if 'gene_ype "' in fields[8]:
                gene_biotype = fields[8].split('gene_type "')[1].split('";')[0]
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT'] = gene_biotype
            if doField1BT:
                gene_biotype=fields[1]
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT'] = gene_biotype

    outfileITS = open(outprefix + '.individual_transcript_stats', 'w')
    outline = '#GeneID\tGeneName\ttranscriptID\ttranscriptName\tchr\tleft\tright\tstrand\tExons\tlength\tbiotype'
    outfileITS.write(outline + '\n')
    keys = GeneDict.keys()
    keys.sort()

    IntronDict = {}

    for (chr,geneID,geneName,transcriptID,transcriptName) in keys:
        outline = geneID + '\t' + geneName + '\t' + transcriptID + '\t' + transcriptName + '\t' + chr
        GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'] = list(Set(GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons']))
        GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'].sort()
        coordinates = []
        length = 0
        for (left,right,strand) in GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons']:
            coordinates.append(left)
            coordinates.append(right)
            length += math.fabs(right-left)
        outline = outline + '\t' + str(min(coordinates)) + '\t' + str(max(coordinates)) + '\t' + strand + '\t' + str(len(GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'])) + '\t' + str(length)
        if GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)].has_key('gBT'):
            outline = outline + '\t' + GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT']
        else:
            outline = outline + '\t' + '-'
        outfileITS.write(outline + '\n')
        if len(GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons']) > 1:
            for i in range(len(GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'])-1):
                intronLeft = GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'][i][1]
                intronRight = GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'][i+1][0]
                intron = (chr,intronLeft,intronRight,strand)
                if IntronDict.has_key(intron):
                    pass
                else:
                    IntronDict[intron] = {}
                    IntronDict[intron]['geneIDs'] = {}
                    IntronDict[intron]['geneNames'] = {}
                    IntronDict[intron]['transcriptIDs'] = {}
                    IntronDict[intron]['transcriptNames'] = {}
                    IntronDict[intron]['gBTs'] = {}
                IntronDict[intron]['geneIDs'][geneID] = 1
                IntronDict[intron]['geneNames'][geneName] = 1
                IntronDict[intron]['transcriptIDs'][transcriptID] = 1
                IntronDict[intron]['transcriptNames'][transcriptName] = 1
                if GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)].has_key('gBT'):
                    IntronDict[intron]['gBTs'][GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT']] = 1
                else:
                    pass
            
    outfileITS.close()

    outfileIS = open(outprefix + '.intron_stats', 'w')
    outline = '#GeneID(s)\tGeneName(s)\ttranscriptID(s)\ttranscriptName(s)\tchr\tleft\tright\tstrand\tlength\tbiotype(s)'
    outfileIS.write(outline + '\n')

    keys = IntronDict.keys()
    keys.sort()
    for (chr,intronLeft,intronRight,strand) in keys:
        outline = ''
        intron = (chr,intronLeft,intronRight,strand)
        for geneID in IntronDict[intron]['geneIDs'].keys():
            outline = outline + geneID + ','
        outline = outline[0:-1] + '\t'
        for geneName in IntronDict[intron]['geneNames'].keys():
            outline = outline + geneName + ','
        outline = outline[0:-1] + '\t'
        for transcriptID in IntronDict[intron]['transcriptIDs'].keys():
            outline = outline + transcriptID + ','
        outline = outline[0:-1] + '\t'
        for transcriptName in IntronDict[intron]['transcriptNames'].keys():
            outline = outline + transcriptName + ','
        outline = outline[0:-1]
        outline = outline + '\t' + chr + '\t' + str(intronLeft) + '\t' + str(intronRight) + '\t' + strand + '\t' + str(intronRight - intronLeft) + '\t'
        if len(IntronDict[intron]['gBTs'].keys()) > 0:
            for gBT in IntronDict[intron]['gBTs'].keys():
                outline = outline + gBT + ','
            outline = outline[0:-1]
        else:
            outline = outline + '-'
        outfileIS.write(outline + '\n')

    outfileIS.close()
   
run()
