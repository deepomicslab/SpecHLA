import gzip
import sys
import os
import re

def read_assign(assign_file):
    dict = {}
    for line in open(assign_file, 'r'):
        line=line.strip()
        array = line.split()
        dict[array[0]] = array[1]
    return dict
#print (n)
def filter_fq(file, gene, outdir, index, dict):
    i = 0
    #gene = 'A'
    outfile = outdir + '/%s.R%s.fq'%(gene, index)
    out = open(outfile, 'w')
    flag = False
    for line in gzip.open(file,'rt'):
        line = line.strip()
        if i % 4 == 0:
            if re.search('/1',line) or re.search('/2',line):
                read_name = line.split()[0][1:-2]
            else:
                read_name = line.split()[0][1:]
            if re.search('HSQ1004:134:C0D8DACXX:4:2308:9795:171485', line):
                print (line, read_name)
            if read_name in dict.keys() and dict[read_name] == gene:
                flag = True
                num = 1
                print (line, file = out)
        elif flag:
            print (line, file = out)
            num += 1
            if num == 4:
                flag = False
        i += 1
    out.close()
    os.system('gzip -f %s'%(outfile))

def main(fq1, fq2, outdir, assign_file):
    dict = read_assign(assign_file)
    for gene in ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']:
        filter_fq(fq1, gene, outdir, 1, dict)
        filter_fq(fq2, gene, outdir, 2, dict)


fq1 = sys.argv[1]
fq2 = sys.argv[2]
outdir = sys.argv[3]
assign_file = outdir + '/assign_file.txt'
main(fq1, fq2, outdir, assign_file)



#fq1='/mnt/disk2_workspace/wangmengyao/NeedleHLA/new_WES/fastp/NA18573/NA18573_1.filter.fastq.gz'
#fq2='/mnt/disk2_workspace/wangmengyao/NeedleHLA/new_WES/fastp/NA18573/NA18573_2.filter.fastq.gz'
#fq1='/home/wangmengyao/scripts/NeedleHLA/1000WES/noise_num/NA18573/NA18573.corrected.R1.fq.gz'
#fq2='/home/wangmengyao/scripts/NeedleHLA/1000WES/noise_num/NA18573/NA18573.corrected.R2.fq.gz'
#filter_fq(fq1, 'A.R1.fq')
#filter_fq(fq2, 'A.R2.fq')
