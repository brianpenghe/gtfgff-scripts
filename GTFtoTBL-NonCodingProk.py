##################################
#                                #
# Last modified 01/23/2015       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
from sets import Set
import os
import subprocess

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s GTF locus_tag_prefix start_count_from output [-tRNAScanSE]'
        print '\tonly use this script for prokaryote genomes!!!'
        print '\tthe -tRNAScanSE option tells the script to rename the tRNAs in a way that will make the compliant with tbl2asn'
        sys.exit(1)

FIX ACCORDING TO http://www.ncbi.nlm.nih.gov/genbank/genome_ncrna_information

    input = sys.argv[1]
    LT = sys.argv[2]
    SCF = int(sys.argv[3])
    outfilename = sys.argv[4]

    dotRNAScanSE = False
    if '-tRNAScanSE' in sys.argv:
        dotRNAScanSE = True

    FeatureDict = {}

    FC = 0

    lineslist  = open(input)
    for line in lineslist:
        if line.strip() == '':
            continue
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chr = fields[0].strip()
        type = fields[2].strip()
        start = fields[3].strip()
        end = fields[4].strip()
        strand = fields[6].strip()
        geneID = fields[8].split('gene_id "')[1].split('";')[0]
        if 'gene_name "' in fields[8]:
            geneName = fields[8].split('gene_name "')[1].split('";')[0]
        else:
            geneName = geneID
        if FeatureDict.has_key(chr):
            pass
        else:
            FeatureDict[chr] = {}
        if FeatureDict[chr].has_key(type):
            pass
        else:
            FeatureDict[chr][type] = []
        FeatureDict[chr][type].append((start,end,strand,geneID,geneName))
        FC += 1

    outfile = open(outfilename,'w')

    LTC = SCF

    chromosomes = FeatureDict.keys()
    chromosomes.sort()
    for chr in chromosomes:
        outline = '>Features ' + chr
        outfile.write(outline + '\n')
        for type in FeatureDict[chr].keys():
            for (start,end,strand,geneID,geneName) in FeatureDict[chr][type]:
                if type == 'SRP' or type == 'RNaseP' or type == '6S':
                    if strand == '+':
                        outline = start + '\t' + end + '\t' + 'misc_RNA'
                        outfile.write(outline + '\n')
                        outline = '\t\t\tnote\t' + geneID + ' ' + type
                        outfile.write(outline + '\n')
                        outline = '\t\t\tlocus_tag\t' + LT + '_'
                        for i in range(max(len(str(FC)),5) - len(str(LTC))):
                            outline = outline + '0'
                        outline = outline + str(LTC)
                        outfile.write(outline + '\n')
                        outline = '\t\t\tproduct\t' + geneID
                        outfile.write(outline + '\n')
                    if strand == '-':
                        outline = end + '\t' + start + '\t' + 'misc_RNA'
                        outfile.write(outline + '\n')
                        outline = '\t\t\tnote\t' + geneID + ' ' + type
                        outfile.write(outline + '\n')
                        outline = '\t\t\tlocus_tag\t' + LT + '_'
                        for i in range(max(len(str(FC)),5) - len(str(LTC))):
                            outline = outline + '0'
                        outline = outline + str(LTC)
                        outfile.write(outline + '\n')
                        outline = '\t\t\tproduct\t' + geneID
                        outfile.write(outline + '\n')
                elif type == 'intron' or type == 'riboswitch':
                    if strand == '+':
                        outline = start + '\t' + end + '\t' + 'misc_feature'
                        outfile.write(outline + '\n')
                        outline = '\t\t\tnote\t' + geneID + ' ' + type
                        outfile.write(outline + '\n')
                        outline = '\t\t\tlocus_tag\t' + LT + '_'
                        for i in range(max(len(str(FC)),5) - len(str(LTC))):
                            outline = outline + '0'
                        outline = outline + str(LTC)
                        outfile.write(outline + '\n')
                    if strand == '-':
                        outline = end + '\t' + start + '\t' + 'misc_feature'
                        outfile.write(outline + '\n')
                        outline = '\t\t\tnote\t' + geneID + ' ' + type
                        outfile.write(outline + '\n')
                        outline = '\t\t\tlocus_tag\t' + LT + '_'
                        for i in range(max(len(str(FC)),5) - len(str(LTC))):
                            outline = outline + '0'
                        outline = outline + str(LTC)
                        outfile.write(outline + '\n')
                elif type == 'tRNA':
                    if dotRNAScanSE:
                        newGN = 'tRNA' + '-' + geneName.split('-')[1][0:3]
                    else:
                        newGN = geneName
                    if strand == '+':
                        outline = start + '\t' + end + '\t' + 'gene'
                        outfile.write(outline + '\n')
                        outline = '\t\t\tlocus_tag\t' + LT + '_'
                        for i in range(max(len(str(FC)),5) - len(str(LTC))):
                            outline = outline + '0'
                        outline = outline + str(LTC)
                        outfile.write(outline + '\n')
                        outline = start + '\t' + end + '\t' + type
                        outfile.write(outline + '\n')
                        outline = '\t\t\tproduct\t' + newGN
                        outfile.write(outline + '\n')
                    if strand == '-':
                        outline = end + '\t' + start + '\t' + 'gene'
                        outfile.write(outline + '\n')
                        outline = '\t\t\tlocus_tag\t' + LT + '_'
                        for i in range(max(len(str(FC)),5) - len(str(LTC))):
                            outline = outline + '0'
                        outline = outline + str(LTC)
                        outfile.write(outline + '\n')
                        outline = end + '\t' + start + '\t' + type
                        outfile.write(outline + '\n')
                        outline = '\t\t\tproduct\t' + newGN
                        outfile.write(outline + '\n')
                else:
                    if strand == '+':
                        outline = start + '\t' + end + '\t' + 'gene'
                        outfile.write(outline + '\n')
                        outline = '\t\t\tlocus_tag\t' + LT + '_'
                        for i in range(max(len(str(FC)),5) - len(str(LTC))):
                            outline = outline + '0'
                        outline = outline + str(LTC)
                        outfile.write(outline + '\n')
                        outline = start + '\t' + end + '\t' + type
                        outfile.write(outline + '\n')
                        outline = '\t\t\tproduct\t' + geneName
                        outfile.write(outline + '\n')
                    if strand == '-':
                        outline = end + '\t' + start + '\t' + 'gene'
                        outfile.write(outline + '\n')
                        outline = '\t\t\tlocus_tag\t' + LT + '_'
                        for i in range(max(len(str(FC)),5) - len(str(LTC))):
                            outline = outline + '0'
                        outline = outline + str(LTC)
                        outfile.write(outline + '\n')
                        outline = end + '\t' + start + '\t' + type
                        outfile.write(outline + '\n')
                        outline = '\t\t\tproduct\t' + geneName
                        outfile.write(outline + '\n')
                LTC += 1               

    outfile.close()

run()