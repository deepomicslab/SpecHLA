import sys
import re

#python3 ../split_hic.py sim.fastq split_test

file=sys.argv[1]
sample=sys.argv[2]
outdir = sys.argv[3]

fwd = open('%s/%s.fwd.fq'%(outdir, sample),'w')
rev = open('%s/%s.rev.fq'%(outdir,sample),'w')

flag = False
i = 0
for line in open(file):
    if line[0] == '@':
        if i % 2 == 0:#line.strip()[-1] == 'R':
            fh = rev
        # elif line.strip()[-1] == 'F':
        else:
            fh = fwd
        i += 1
    print (line, end = '', file = fh)
fwd.close()
rev.close()

