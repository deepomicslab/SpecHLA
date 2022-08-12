"""
In some real samples, the read name is not unique.
This script removes the repeat reads
"""

import sys
import gzip
import os

def read_fq(file, outfile):
    saved_name = {}
    out = gzip.open(outfile, 'wt')
    line_num = 0
    for line in gzip.open(file, 'rt'):
        line = line.strip()
        #if line[0] == '@':
        if line_num % 4 == 0:
            name = line[1:-2]
            if name in saved_name.keys():
                saved_name[name] += 1
                name = name + '_' + str(saved_name[name])
                line = line[0] + name + line[-2:]
            else:
                saved_name[name] = 1
            print (line, file= out)
            #print (line, name)
        #break
        else:
            print (line, file= out)
        line_num += 1
    out.close()
    return line_num

def read_fq_unzip(file, outfile):
    saved_name = {}
    out = open(outfile, 'w')
    line_num = 0
    for line in open(file, 'r'):
        line = line.strip()
        #if line[0] == '@':
        if line_num % 4 == 0:
            name = line[1:-2]
            if name in saved_name.keys():
                saved_name[name] += 1
                name = name + '_' + str(saved_name[name])
                line = line[0] + name + line[-2:]
            else:
                saved_name[name] = 1
            print (line, file= out)
            #print (line, name)
        #break
        else:
            print (line, file= out)
        line_num += 1
    out.close()
    return line_num

if __name__ == "__main__":

    if sys.argv[1][-3:] == ".gz":
        line_num = read_fq(sys.argv[1], sys.argv[2])
    else:
        line_num = read_fq_unzip(sys.argv[1], sys.argv[2][:-3])
        os.system("gzip -f %s"%(sys.argv[2][:-3]))
    reads_num = int(line_num/4)
    if reads_num > 1000000:
        print ("WARNING: %s reads are detected, please check if they are HLA-related reads.\
            If not, please extract HLA-related reads first. Otherwise, the process will be slow."%(reads_num))
