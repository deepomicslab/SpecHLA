import sys
import os
import pysam
import numpy as np
import re

B = ['A', 'C', 'G', 'T']
BN = ['A', 'C', 'G', 'T', 'N']


def extract_reads_name(options):
    reads_file = '%s/read_names.txt'%(options.outdir)
    order = 'samtools view %s | cut -f1 | sort | uniq > %s'%(options.bam, reads_file)
    os.system(order)

# def extract_reads_bam(read, bamfile, outdir):
#     order = 'samtools view %s | grep -f %s > %s/read.mapped.info.txt'%(bamfile, read, outdir)

# def iterate_reads(bamfile, outdir):
#     reads_file = '%s/read_names.txt'%(outdir)
#     for line in open(reads_file, 'r'):
#         read = line.strip()

def get_names(names):
    with open(names, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n

def extract_reads(options):
    reads_file = '%s/read_names.txt'%(options.outdir)
    assign_file = '%s/assign_file.txt'%(options.outdir)
    out = open(assign_file, 'w')
    n = get_names(reads_file)
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    error, total, remove = 0, 0, 0
    error_set = []
    for name in n:
        try:
            iterator = name_indexed.find(name)
            dict = {}
            len_dict = {}
            for x in iterator:   
                if x.is_unmapped:
                    continue
                s = 0
                # locus_set = x.get_reference_positions(full_length=True)
                t_name = x.reference_name
                # if name == 'A00132:58:HFL2TDSXX:4:2276:6515:29183':
                #     print ('#', t_name, x.next_reference_name)
                if t_name != x.next_reference_name:
                    continue
                # print (t_name, x.next_reference_name)
                # if not re.search('DRB1', t_name):
                #     continue
                match_num = 0
                soft_num = 0
                all_num = 0
                for ci in x.cigar:
                    if ci[0] == 0:
                        match_num += ci[1]
                    elif ci[0] == 4:
                        soft_num += ci[1]
                    all_num += ci[1]
                if soft_num > 0:
                    continue
      

                mis_NM = 0
                for ta in x.get_tags():
                    if ta[0] == 'NM':
                        match_num -= ta[1]
                        mis_NM += ta[1]
                if mis_NM > options.max_nm:
                    continue
                        # print (ta, match_num)
                # match_rate = match_num/(all_num - soft_num)
                focus_len = all_num - soft_num

                
                # if name == 'A00217:72:HFKHYDSXX:4:2603:32542:29637':
                #     print (x.cigar, match_rate, all_num, soft_num)
                # refSequence = x.get_reference_sequence()
                # pre_ref, pre_read = '', ''
                # j = 0
                # for i in range(len(locus_set)):
                #     read_allele = x.query_sequence[i].upper()
                #     if str(locus_set[i]) != 'None':
                #         ref_allele = refSequence[j].upper()
                #         j += 1
                #     else:
                #         ref_allele = 'NONE'
                #     base_quality = x.query_qualities[i]
                #     p = 10 ** (- base_quality / 10)
                    
                #     alpha = calculate_alpha(read_allele, ref_allele, pre_ref, pre_read, p, options)
                #     beta = calculate_beta(read_allele, options)
                #     # print (p, alpha, beta)
                #     s += alpha
                #     s += beta
                #     pre_ref, pre_read = ref_allele, read_allele
                if t_name not in dict.keys():
                    dict[t_name] = match_num#round(s,3)
                    len_dict[t_name] = focus_len
                else:
                    dict[t_name] += match_num#round(s,3)
                    len_dict[t_name] += focus_len
            #evaluation
            total += 1
            if len(dict) == 0:
                continue

            # first_align = l[0][0].split('*')[0]
            # reads_len = #list(len_dict.values())[0]
            for key in dict.keys():
                # print (len_dict[key])
                if len_dict[key] < 0:  #make sure the reads is paired mapped.
                    dict[key] = 0
                else:
                    dict[key] = float(dict[key])/len_dict[key]
            first_align = check_score(dict, options, name)
            # print (dict, first_align, reads_len)
            # break
            if first_align == 'REMOVE':
                remove += 1
                continue
            print (name, first_align, file = out)
        except KeyError:
            pass
    out.close()

def check_score(dict, options, name):
    gene_dict = {}
    for align in dict.keys():
        # print (l)
        gene = align.split('*')[0]
        if gene not in gene_dict.keys():
            gene_dict[gene] = dict[align]
        else:
            if gene_dict[gene] < dict[align]:
                gene_dict[gene] = dict[align]
    new_l = sorted(gene_dict.items(), key=lambda gene_dict:gene_dict[1], reverse = True)
    # if float(new_l[0][1]) < 0.999:#options.reads_map:
    if float(new_l[0][1]) < 0:
        if new_l[0][0] == 'DRB1':
            print (new_l[0][0], name, 'too much mismatch', new_l[0][1])
        return 'REMOVE'
    elif len(new_l) == 1:
        return new_l[0][0]
    elif new_l[0][1] - new_l[1][1] < 0.5:# or new_l[0][1] < options.theta_pm:
        # print ('remove reads', new_l)
        # print (name.split('_')[0], new_l[0][0], new_l[1][0], new_l)
        if new_l[0][0] == 'DRB1':
            print (new_l[0][0], name, 'too similar with other genes', new_l[1][0], new_l[0][1] - new_l[1][1])
        return 'REMOVE'
    else:
        return new_l[0][0]
    # print (gene_dict)

def extract_reads_cigar(options):
    reads_file = '%s/read_names.txt'%(options.outdir)
    n = get_names(reads_file)
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    header = bamfile.header.copy()
    # out = pysam.Samfile(options.out, 'wb', header=header)
    for name in n:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            n = 0
            for x in iterator:
                n += 1
            multi_align_list = []
            iterator = name_indexed.find(name)
            for x in iterator:   
                #if n == 2 or x.cigar == [(0, 150)]:
                    # print (name, x.reference_name, x, x.cigar,cigar2array(x.cigar))  
                 #   continue
                if x.is_unmapped:
                    continue
                cigararray = cigar2array(x.cigar)
                # print (x.cigar)
                s = 0
                # locus_set = x.get_reference_positions(full_length=True)
                t_name = x.reference_name

                # pre_ref, pre_read = '', ''
                for i in range(len(cigararray)):
                    base_quality = x.query_qualities[i]
                    p = 10 ** (- base_quality / 10)
                    cigar = cigararray[i]
                    print (p, cigar)
                #     alpha = calculate_alpha(read_allele, ref_allele, pre_ref, pre_read, p, options)
                #     beta = calculate_beta(read_allele, options)
                #     s += alpha
                #     s += beta

def cigar2array(cigar):
    cigar_list = []
    for ci in cigar:
        for i in range(ci[1]):
            cigar_list.append(ci[0])
    return cigar_list

def calculate_alpha(read_allele, ref_allele, pre_ref, pre_read, p, options):
    if read_allele in B and ref_allele in B and read_allele != ref_allele:
        #print (np.log(p/3), read_allele, ref_allele)
        return 0#-1 #np.log(p/3)
    elif read_allele == 'NONE' and pre_read != 'NONE':
        return options.alpha_do
    elif read_allele == 'NONE' and pre_read == 'NONE':
        return options.alpha_de
    elif ref_allele == 'NONE' and pre_ref != 'NONE':
        return options.alpha_io
    elif ref_allele == 'NONE' and pre_ref == 'NONE':
        return options.alpha_ie
    elif (read_allele in BN and ref_allele == 'N') or (ref_allele in BN and read_allele == 'N'):
        return options.alpha_N
    else:
        return 0
    
def calculate_beta(read_allele, options):
    if read_allele in BN:
        return options.beta
    else:
        return 0

def read_gene():
    file = '/mnt/disk2_workspace/wangshuai/00.strain/08.NeedleHLA/alphlard/ref/genes.txt'
    gene_set = []
    for line in open(file, 'r'):
        gene = line[1:].strip()
        gene_set.append(gene)
    return gene_set


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Extract reads by read name from bam file')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-o', '--outdir', help='outdir', required=True)   
    #parameters
    parser.add_argument('-do', '--alpha_do', help='alpha_do', required=False, default = 0, type=float)   #-1
    parser.add_argument('-de', '--alpha_de', help='alpha_de', required=False, default = 0, type=float)
    parser.add_argument('-io', '--alpha_io', help='alpha_io', required=False, default = 0, type=float)
    parser.add_argument('-ie', '--alpha_ie', help='alpha_ie', required=False, default = 0, type=float)
    parser.add_argument('-an', '--alpha_N', help='alpha_N', required=False, default = 0, type=float)   #-1
    parser.add_argument('-be', '--beta', help='beta', required=False, default = 1, type=float)
    parser.add_argument('-pd', '--theta_pd', help='theta_pd', required=False, default = 20, type=float)
    parser.add_argument('-pm', '--theta_pm', help='theta_pm', required=False, default = 100, type=float)
    parser.add_argument('-ma', '--reads_map', help='reads map creterion', required=False, default = 0, type=float)
    parser.add_argument('-nm', '--max_nm', help='MAX NM', required=False, default = 2, type=int)
    options = parser.parse_args()
    extract_reads_name(options)
    # gene_set = read_gene()
    # for gene in ['C', 'H']:
    extract_reads(options)
