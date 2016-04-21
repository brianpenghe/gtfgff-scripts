##################################
#                                #
# Last modified 07/13/2013       # 
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

    GeneDict={}
    listoflines = open(inputfilename)
    for line in listoflines:
        if line.startswith('#'):
            continue
        fields=line.strip().split('\t')
        ID=fields[8].split('gene_id "')[1].split('"')[0]
        if 'gene_name "' in fields[8]:
            name=fields[8].split('gene_name "')[1].split('"')[0]
        else:
            name = ID
        chr=fields[0]
        left=int(fields[3])
        right=int(fields[4])
        strand=fields[6]
        if GeneDict.has_key((ID,name)):
            pass
        else:
            GeneDict[(ID,name)]={}
            GeneDict[(ID,name)]['strand']=strand
            GeneDict[(ID,name)]['chr']=chr
            GeneDict[(ID,name)]['coordinates']=[]
        GeneDict[(ID,name)]['coordinates'].append(left)
        GeneDict[(ID,name)]['coordinates'].append(right)

    outline='#chr\tleft\tright\tstrand\tID\tname'
    outfile.write(outline+'\n')

    for (ID,name) in GeneDict.keys():
        strand=GeneDict[(ID,name)]['strand']
        left=min(GeneDict[(ID,name)]['coordinates'])
        right=max(GeneDict[(ID,name)]['coordinates'])
        chr=GeneDict[(ID,name)]['chr']
        outline=chr+'\t'+str(left)+'\t'+str(right)+'\t'+strand+'\t'+ID+'\t'+name
        outfile.write(outline+'\n')
            
    outfile.close()

run()

