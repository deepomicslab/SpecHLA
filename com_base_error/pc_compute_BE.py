import sys
import os

out = open('select.sample.hla.list', 'w')
sa = open('sample.list', 'w')

#fq_dir = sys.argv[1]
outdir = sys.argv[1]

for i in range(1,51):
#for i in range(id,id+1):
    sample = 'child_%s'%(i)
    print (sample, file = sa)
    true_fa = '/mnt/e/hla_tgs/child_fasta/%s.fasta'%(sample)
    i = 0
    for line in open(true_fa):
        line = line.strip()

        if line[0] == '>':
            allele = line[1:]
            gene = 'HLA_' + allele.split('_')[0]
            if i % 2 == 0:
                print (sample, allele, end = '\t', file = out)
            else:
                print (allele, gene, end = '\n', file = out)
            i += 1

out.close()
sa.close()
#os.system('sh /mnt/d/HLAPro_backup/HLAPro/com_base_error//work.sh %s'%(outdir))
