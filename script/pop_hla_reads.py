import sys
import os
import re

indir = sys.argv[1]
outdir = sys.argv[2]
extract_sh = '/mnt/d/HLAPro_backup/HLAPro/script/ExtractHLAread.sh'


os.system('find %s/*/*bam>tmp.bam.list'%(indir))
for line in open('tmp.bam.list'):
    bam = line.strip()
    id_re = re.search('(.*?).bam', bam.split('/')[-1])
    id = id_re.group(1)
    order = """
    bash %s -s %s/%s -b %s -r hg38
    """%(extract_sh, outdir, id, bam)
    os.system(order)
    print (id, 'is done.')

