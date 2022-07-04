"""
Assess the read assignment in simulation reads
return the recall and precision

wangshuai July 5, 2022
"""

import gzip
import re
import argparse


def read_truth():
    # file = '/mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20200901/merge.groudtruth.txt'
    file = args["t"]
    f = open(file, 'r')
    print ("# sample    gene    recall  precision")
    for line in f:
        line=line.strip()
        array = line.split()
        if line[0] == 'S':
            name = array
        else:
            print (array[0], end = '\t')
            for i in range(1, len(array), 2):
                gene_name = name[i][:-2]
                gene1 = array[i].replace('*', '_').replace(':', '_')
                gene2 = array[i+1].replace('*', '_').replace(':', '_')
                # true_file = '/mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20200901/simu/fastq/%s/%s.read1.fastq.gz'%(array[0], array[0])
                true_file = args["q"] + '/%s/%s.read1.fastq.gz'%(array[0], array[0])
                true_reads = read_gene_fq(true_file, gene_name)
                # result_file = '/mnt/disk2_workspace/wangshuai/00.strain/08.NeedleHLA/wgs/new_simu/output/%s/%s.R1.fq.gz'%(array[0], gene_name)
                result_file = args["r"] + '/%s/%s.R1.fq.gz'%(array[0], gene_name)
                result_reads = read_fq(result_file)
                recall_rate, precision_rate = eva(true_reads, result_reads)
                print (gene_name, recall_rate, precision_rate, end = '\t' )
            #break
            print ('')
    f.close()

def eva(true_reads, result_reads):
    true_recall_num = 0
    for read in result_reads:
        if read in true_reads:
            true_recall_num += 1
    recall_rate = float(true_recall_num)/len(true_reads)
    precision_rate = float(true_recall_num)/len(result_reads)
    return round(recall_rate, 3), round(precision_rate, 3)

def read_gene_fq(file, gene):
    reads = []
    for line in gzip.open(file, 'rt'):
        line = line.strip()
        if line[0] == '@':
            name = line[1:-2]
            #if name not in reads and re.search(gene, line):
            #    reads.append(name)
            if re.search('@%s_'%(gene), line):
                reads.append(name)
                #print (name, gene)
            #else:
             #   print ('same name.')
    return reads

def read_fq(file):
    reads = []
    for line in gzip.open(file, 'rt'):
        line = line.strip()
        if line[0] == '@':
            name = line[1:-2]
            reads.append(name)
            #if name not in reads:
            #    reads.append(name)
            #else:
            #    print ('same name.')
    return reads
            
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Assess read assignment of SpecHLA", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-t", type=str, help="The file records the allele in each sample. <merge.groudtruth.txt>", metavar="\b")
    required.add_argument("-q", type=str, help="The folder saves fastq file of each allele in simulation. <simu_20200901/simu/fastq>", metavar="\b")
    required.add_argument("-r", type=str, help="The folder save results of SpecHLA. <wgs/new_simu/output/>", metavar="\b")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    read_truth()
