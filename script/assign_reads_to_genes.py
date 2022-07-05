"""
assin reads to original gene loci with the alignment 
"""

import sys
import os
import pysam
import numpy as np
import re
import time
import gzip

B = ['A', 'C', 'G', 'T']
BN = ['A', 'C', 'G', 'T', 'N']

def get_names():
    with open(reads_file, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n

def count_alignment(alignment):
    match_num = 0
    soft_num = 0
    all_num = 0
    for ci in alignment.cigar:
        if ci[0] == 0:
            match_num += ci[1]
        elif ci[0] == 4:
            soft_num += ci[1]
        all_num += ci[1]

    mis_NM = 0
    for ta in alignment.get_tags():
        if ta[0] == 'NM':
            match_num -= ta[1]
            mis_NM += ta[1]
    focus_len = all_num - soft_num

    return mis_NM, soft_num, match_num, focus_len

def check_score(dict, options, name, pair_dict):
    gene_dict = {}
    for align in dict.keys():
        # print (l)
        if pair_dict[align] < 2:
            dict[align] = 0
        gene = align.split('*')[0]
        if gene not in gene_dict.keys():
            gene_dict[gene] = dict[align]
        else:
            if gene_dict[gene] < dict[align]:
                gene_dict[gene] = dict[align]
    new_l = sorted(gene_dict.items(), key=lambda gene_dict:gene_dict[1], reverse = True)
    if float(new_l[0][1]) < 0.1:
        #if new_l[0][0] == 'DRB1':
        #    print (new_l[0][0], name, 'too much mismatch', new_l[0][1])
        return 'REMOVE'
    elif len(new_l) == 1:
        return new_l[0][0]
    elif new_l[0][1] - new_l[1][1] < options.diff_score:
        return 'REMOVE'
    else:
        return new_l[0][0]
    # print (gene_dict)

class Each_read():
    """
    for each pair-end reads
    record all the alignment scores
    assign the pair-end read to gene loci by alignment scores
    """

    def __init__(self):
        self.dict = {}
        self.pair_dict = {}
        self.len_dict = {}    
        self.read_name = ''   

    def add_one_alignment(self, alignment):
        """
        given a alignment record,
        save the mapping situations to the read's dicts.
        """
        flag = True
        if alignment.is_unmapped or alignment.reference_name != alignment.next_reference_name:
            flag = False
        t_name = alignment.reference_name

        mis_NM, soft_num, match_num, focus_len = count_alignment(alignment)
        if soft_num > 0 or mis_NM > options.max_nm:
            flag = False
        if flag:
            self.read_name = alignment.query_name
            if t_name not in self.dict.keys():
                self.dict[t_name] = match_num#round(s,3)
                self.len_dict[t_name] = focus_len
                self.pair_dict[t_name] = 1
            else:
                self.dict[t_name] += match_num#round(s,3)
                self.len_dict[t_name] += focus_len
                self.pair_dict[t_name] += 1
    
    def assign(self):
        """
        determine which locus the read should be assigned
        """
        
        flag = True
        if len(self.dict) == 0:
            flag = False
        for key in self.dict.keys():
            # print (len_dict[key])
            if self.len_dict[key] < 0:  #make sure the reads is paired mapped.
                self.dict[key] = 0
            else:
                self.dict[key] = float(self.dict[key])/self.len_dict[key]
        if flag:
            first_align = check_score(self.dict, options, self.read_name, self.pair_dict)
        else:
            first_align = 'REMOVE'
        return first_align

def assign_fastq(file, gene, index, assign_dict):
    """
    extract gene-specific reads from the raw read file,
    generate gene-specific fastq files
    read name in the original fastq should be
    @<read name> or @<read name>/1
    """
    i = 0
    #gene = 'A'
    outfile = options.outdir + '/%s.R%s.fq'%(gene, index)
    out = open(outfile, 'w')
    flag = False
    if file.split(".")[-1] == "gz":
        f = gzip.open(file,'rt')
    else:
        f = open(file)
    for line in f:
        line = line.strip()
        if i % 4 == 0:
            if re.search('/1',line) or re.search('/2',line):
                read_name = line.split()[0][1:-2]
            else:
                read_name = line.split()[0][1:]
            if read_name in assign_dict.keys() and assign_dict[read_name] == gene:
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

def main():
    print ('start assigning reads...')
    read_dict = {} # record the aligment scores for each read
    assign_dict = {} # record assigned gene for each read
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    t0 = time.time()

    # record alignment scores for all the read
    for alignment in bamfile:
        t_name = alignment.reference_name
        alignment.query_name = alignment.query_name.split("/")[0]
        read_name = alignment.query_name
        # print (read_name)
        if read_name not in read_dict.keys():
            
            new_read = Each_read()
            read_dict[read_name] = new_read
        read_dict[read_name].add_one_alignment(alignment)
    bamfile.close()
    t1 = time.time()
    print ("read bam cost %s"%(t1 - t0))

    # assign genes for each read
    for read_name in read_dict:
        assigned_locus = read_dict[read_name].assign()
        if assigned_locus == 'REMOVE':
            continue
        assign_dict[read_name] = assigned_locus

    # generate gene-specific fastq
    for gene in ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']:
        assign_fastq(options.fq1, gene, 1, assign_dict)
        assign_fastq(options.fq2, gene, 2, assign_dict)
    t2 = time.time()
    print ("read assigment cost %s"%(t2 - t0))

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Read assignment')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-o', '--outdir', help='outdir', required=True)  
    parser.add_argument('-1', '--fq1', help='bin dir', required=True) 
    parser.add_argument('-2', '--fq2', help='bin dir', required=True)  
    parser.add_argument('-n', '--bin_dir', help='bin dir', required=True)  
    parser.add_argument('-nm', '--max_nm', help='MAX NM', required=False, default = 2, type=int)
    parser.add_argument('-d', '--diff_score', help='The score for the best gene must be at least this higher\
         than the second gene', required=False, default = 0.5, type=float)
    options = parser.parse_args()

    main()