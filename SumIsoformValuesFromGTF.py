##################################
#                                #
# Last modified 03/09/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s input GTF outputfilename  ' % argv[0]
        sys.exit(1)
    
    input = argv[1]
    GTF = argv[2]
    outfilename = argv[3]
    outfile = open(outfilename, 'w')

    GeneDict={}
    TranscritpToGeneDict={}

    linelist = open(GTF)
    for line in linelist:
        fields=line.strip().split('\t')
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        geneID=fields[8].split('gene_id "')[1].split('"')[0]
        if 'gene_name "' in fields[8]:
            geneName=fields[8].split('gene_name "')[1].split('"')[0]
        else:
            geneName = geneID
        transcriptID=fields[8].split('transcript_id "')[1].split('"')[0]
        TranscritpToGeneDict[transcriptID]=(geneID,geneName)

    print 'Finished inputting GTF'

    outline = '#geneID\tgeneName\ttot_count\tuniq_count\tpost_counds_mode\teff_counts_mode\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi'
    outfile.write(outline + '\n')

    linelist = open(input)
    for line in linelist:
        fields=line.strip().split('\t')
        if line.startswith('bundle_id\t'):
            TotCountFieldID = fields.index('tot_counts')
            UniqCountID = fields.index('uniq_counts')
            PostCountsModeID = fields.index('est_counts')
            EffCountsModeID = fields.index('eff_counts')
            FPKMID = fields.index('fpkm')
            FPKM_lo_ID = fields.index('fpkm_conf_low')
            FPKM_hi_ID = fields.index('fpkm_conf_high')
            continue
        transcriptID = fields[1]
        TotCount = float(fields[TotCountFieldID])
        UniqCount = float(fields[UniqCountID])
        PostCountsMode = float(fields[PostCountsModeID])
        EffCountsMode = float(fields[EffCountsModeID])
        FPKM = float(fields[FPKMID])
        FPKM_conf_lo = float(fields[FPKM_lo_ID])
        FPKM_conf_hi = float(fields[FPKM_hi_ID])
        try:
            (geneID,geneName) = TranscritpToGeneDict[transcriptID]
        except:
            print transcriptID, fields
        if GeneDict.has_key((geneID,geneName)):
            GeneDict[(geneID,geneName)]['FPKM'] += FPKM
            GeneDict[(geneID,geneName)]['FPKM_conf_lo'] += FPKM_conf_lo
            GeneDict[(geneID,geneName)]['FPKM_conf_hi'] += FPKM_conf_hi
            GeneDict[(geneID,geneName)]['TotCount'] += TotCount
            GeneDict[(geneID,geneName)]['UniqCount'] += UniqCount
            GeneDict[(geneID,geneName)]['PostCountsMode'] += PostCountsMode
            GeneDict[(geneID,geneName)]['EffCountsMode'] += EffCountsMode
        else:
            GeneDict[(geneID,geneName)]={}
            GeneDict[(geneID,geneName)]['FPKM'] = FPKM
            GeneDict[(geneID,geneName)]['FPKM_conf_lo'] = FPKM_conf_lo
            GeneDict[(geneID,geneName)]['FPKM_conf_hi'] = FPKM_conf_hi
            GeneDict[(geneID,geneName)]['TotCount'] = TotCount
            GeneDict[(geneID,geneName)]['UniqCount'] = UniqCount
            GeneDict[(geneID,geneName)]['PostCountsMode'] = PostCountsMode
            GeneDict[(geneID,geneName)]['EffCountsMode'] = EffCountsMode

    genes = GeneDict.keys()
    genes.sort()

    for (geneID,geneName) in genes:
        outline = geneID + '\t' + geneName + '\t' + str(GeneDict[geneName]['TotCount'])  + '\t' + str(GeneDict[geneName]['UniqCount'])  + '\t' + str(GeneDict[geneName]['PostCountsMode'])  + '\t' + str(GeneDict[geneName]['EffCountsMode'])  + '\t' + str(GeneDict[geneName]['FPKM']) + '\t' + str(GeneDict[geneName]['FPKM_conf_lo']) + '\t' + str(GeneDict[geneName]['FPKM_conf_hi'])
        outfile.write(outline + '\n')                        

    outfile.close()
   
if __name__ == '__main__':
    main(sys.argv)







