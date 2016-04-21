##################################
#                                #
# Last modified 11/14/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
from sets import Set
import math

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s gff JGI|GFF3|GiardiaDB outfile ' % sys.argv[0]
        print '\tNote: for simplicity, the script will only output exons and CDS/UTR annotations'
        sys.exit(1)

    GFF = sys.argv[1]
    GTF = sys.argv[3]

    GeneDict={}
    outfile = open(GTF, 'w')
    if sys.argv[2] == 'JGI':
        linelist=open(GFF)
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
                GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT'] = 'unknown'
            GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['exons'].append((left,right,strand))
        linelist=open(GFF)
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.strip().split('\t')
            if fields[2] == 'exon' or fields[2] == 'CDS':
                pass
            else:
                continue
            geneID=fields[8].split('name "')[1].split('";')[0]
            if fields[2] == 'exon':
                transcriptID=fields[8].split('transcriptId ')[1].split(';')[0]
            if fields[2] == 'CDS':
                transcriptID=fields[8].split('proteinId ')[1].split(';')[0]
            transcriptName=transcriptID
            geneName=geneID
            chr=fields[0]
            gBT = GeneDict[(chr,geneID,geneName,transcriptID,transcriptName)]['gBT']
            outline = fields[0] + '\t'
            outline = outline + fields[1] + '\t'
            outline = outline + fields[2] + '\t'
            outline = outline + fields[3] + '\t'
            outline = outline + fields[4] + '\t'
            outline = outline + fields[5] + '\t'
            outline = outline + fields[6] + '\t'
            outline = outline + fields[7] + '\t'
            outline = outline + 'gene_id "' + geneID + '"; transcript_id "' + transcriptID + '"; gene_name "' + geneName + '"; transcript_name "' + transcriptName + '"; gene_type "' + gBT + '";'
            outfile.write(outline + '\n')
    if sys.argv[2] == 'GFF3':
        linelist=open(GFF)
        GeneTypeDict = {}
        TranscriptParentDict = {}
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.replace('PARENT=','Parent=').strip().split('\t')
            if fields[2].endswith('_gene') or fields[2] == 'gene' or ((fields[2] == 'pseudogene' or fields[2] == 'pseudogenic_tRNA' or fields[2] == 'RNA') and 'Parent=' not in fields[8]):
                gBT = fields[8].split('biotype=')[1].split(';')[0]
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
        linelist=open(GFF)
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.replace('PARENT=','Parent=').strip().split('\t')
            if fields[2] == 'exon' or fields[2] == 'CDS':
                pass
            else:
                continue
            transcriptID=fields[8].split('Parent=')[1].split(';')[0]
            transcriptName = transcriptID
            geneID = TranscriptParentDict[transcriptID]
            (geneName,gBT) = GeneTypeDict[geneID]
            outline = fields[0] + '\t'
            outline = outline + fields[1] + '\t'
            outline = outline + fields[2] + '\t'
            outline = outline + fields[3] + '\t'
            outline = outline + fields[4] + '\t'
            outline = outline + fields[5] + '\t'
            outline = outline + fields[6] + '\t'
            outline = outline + fields[7] + '\t'
            outline = outline + 'gene_id "' + geneID + '"; transcript_id "' + transcriptID + '"; gene_name "' + geneName + '"; transcript_name "' + transcriptName + '"; gene_type "' + gBT + '";'
            outfile.write(outline + '\n')
    if sys.argv[2] == 'GiardiaDB':
        linelist=open(GFF)
        GeneTypeDict = {}
        TranscriptParentDict = {}
        ExonCDSParentDict = {}
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.strip().split('\t')
            if fields[2] == 'mRNA' or fields[2] == 'rRNA' or fields[2] == 'tRNA' or fields[2] == 'RNase_MRP_RNA' or fields[2] == 'RNase_P_RNA' or fields[2] == 'SRP_RNA_encoding' or fields[2] == 'snRNA':
                transcriptID = fields[8].split('ID=')[1].split(';')[0]
                geneID = fields[8].split('Parent=')[1].split(';')[0]
                gBT = fields[2]
                GeneTypeDict[geneID] = gBT
                TranscriptParentDict[transcriptID] = geneID
                continue
            if fields[2] == 'exon' or fields[2] == 'CDS':
                ID = fields[8].split('ID=')[1].split(';')[0]
                transcriptID = fields[8].split('Parent=')[1].split(';')[0]
                ExonCDSParentDict[ID] = transcriptID
                continue
        linelist=open(GFF)
        for line in linelist:
            if line.startswith('#'):
                continue
            fields=line.strip().split('\t')
            if fields[2] == 'exon' or fields[2] == 'CDS':
                pass
            else:
                continue
            ID=fields[8].split('ID=')[1].split(';')[0]
            transcriptID = ExonCDSParentDict[ID]
            geneID = TranscriptParentDict[transcriptID]
            gBT = GeneTypeDict[geneID]
            geneName = geneID
            transcriptName = transcriptID
            outline = fields[0] + '\t'
            outline = outline + fields[1] + '\t'
            outline = outline + fields[2] + '\t'
            outline = outline + fields[3] + '\t'
            outline = outline + fields[4] + '\t'
            outline = outline + fields[5] + '\t'
            outline = outline + fields[6] + '\t'
            outline = outline + fields[7] + '\t'
            outline = outline + 'gene_id "' + geneID + '"; transcript_id "' + transcriptID + '"; gene_name "' + geneName + '"; transcript_name "' + transcriptName + '"; gene_type "' + gBT + '";'
            outfile.write(outline + '\n')

    outfile.close()
   
run()
