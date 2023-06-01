"""
Map the HLA sequence to exon database (2,3 for class I and 2 for class II)
Then get the G group resolution HLA type

wangshuai Feb 27, 2023
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
import re
import pandas as pd
import numpy as np
from collections import Counter
import argparse

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

def count_n_ratio(fasta_file):
    total_bases = 0
    n_bases = 0

    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            total_bases += len(record.seq)
            n_bases += record.seq.count('N')

    ratio = n_bases / total_bases
    ratio = round(ratio, 2)
    return ratio

def read_G_annotation():
    g_file = "%s/../../db/HLA/hla_nom_g.txt"%(sys.path[0])
    G_annotation_dict = {}
    i = 0
    for line in open(g_file):
        if line[0] == "#":
            continue
        # line.replace(":","_", 10000)
        line = re.sub(":","_",line)
        array = line.strip().split(";")
        gene = array[0][:-1]
        
        if len(array[-1]) == 0:
            
            g_name = gene + "_" + array[-2]
            # print (g_name)
            G_annotation_dict[g_name] = g_name
        else:
            g_name = gene + "_" + array[-1]
            alleles = array[-2].split("/")
            for each in alleles:
                each = gene + "_" + each
                G_annotation_dict[each] = g_name
        # print (array, g_name)
        # print (G_annotation_dict)
        # print (len(array))
        # if i > 2:
        #     break
        i += 1
    return G_annotation_dict

def convert_G(allele):
    allele = re.sub(":","_",allele)
    allele = re.sub("HLA-","",allele)
    allele = re.sub("\*","_",allele)
    allele = allele.split(";")[0]
    if allele in G_annotation_dict:
        allele = G_annotation_dict[allele]
    elif allele + "_01" in G_annotation_dict:
        allele = G_annotation_dict[allele + "_01"]
    elif allele + "_01_01" in G_annotation_dict:
        allele = G_annotation_dict[allele + "_01_01"]
    G_type = allele
    return G_type

class G_annotation():
    def __init__(self, sample, spechla_dir):
        self.sample = sample
        self.spechla_dir = spechla_dir
    
    def blast(self, infer_hap_file, blast_result_file):
        # command = f"""
        # blastn -query {infer_hap_file} -out {blast_result_file} -subject {exon_database} -outfmt 7 -max_target_seqs 60000
        # """
        command = f"""
        blastn -query {infer_hap_file} -out {blast_result_file} -db {exon_database} -outfmt 7 -max_target_seqs 60000 -num_threads {args["j"]}
        """
        os.system(command)  
           
    def read_blast(self, blast_result_file):
        # blast_result_file = "/mnt/d/HLAPro_backup/haplotype_v2/spechla/HG00096/test.out"
        f = open(blast_result_file, 'r')
        identity_record = {}
        record_hit_exon_times = {}
        for line in f:
            if line[0] == "#":
                continue
            array = line.strip().split()
            exon = array[1]
            identity = float(array[2])
            match_len = int(array[3])
            allele = exon.split("|")[0]
            if not self.check_pop(allele) and args["p"] != "nonuse": # check allele freq in population
                # print ("<<<")
                continue
            # print (allele, identity, match_len)
            if allele not in identity_record:
                identity_record[allele] = [0, 0]
            if exon not in record_hit_exon_times:
                identity_record[allele][0] += identity
                identity_record[allele][1] += match_len
            record_hit_exon_times[exon] = 1
        f.close()

        if len(identity_record) == 0:
            return "no_match"

        # sort the dictionary by its values in ascending order
        # sorted_identity_record = sorted(identity_record.items(), key=lambda x: x[1], reverse = True)
        sorted_identity_record = sorted(identity_record.items(), key=lambda x: (x[1][0], x[1][1]), reverse = True)
        # print (sorted_identity_record[0])
        # print the sorted dictionary
        top_alleles = []
        max_identity = sorted_identity_record[0][1][0]
        max_length = sorted_identity_record[0][1][1]
        for allele, info in sorted_identity_record:
            if info[0] == max_identity and info[1] == max_length:
                top_alleles.append(convert_G(allele))
                # print(allele, info, convert_G(allele))
            else:
                break
        # print (top_alleles)
        most_common_allele = most_common(top_alleles)
        return (most_common_allele)

    def main(self):
        sample_results = {}
        print ("The region with low read depth is masked by N. The cutoff is specified by -k.")
        for gene in gene_list:
            sample_results[gene] = []
            for hap_index in range(1,3):
                infer_hap_file = f"{self.spechla_dir}/hla.allele.{hap_index}.HLA_{gene}.fasta"
                n_ratio = count_n_ratio(infer_hap_file) # cal the ratio of N in the fasta
                print ("The ratio of N (masked) is %s for the allele %s"%(n_ratio, infer_hap_file))
                blast_result_file = infer_hap_file + ".exon.blast"
                self.blast(infer_hap_file, blast_result_file)
                g_group_type = self.read_blast(blast_result_file)
                sample_results[gene].append(g_group_type)
        # print (self.sample, sample_results)
        # return sample_results
        COUT =  open(f"{self.spechla_dir}/hla.result.g.group.txt", "w")
        COUT.write("Sample\tHLA_A_1\tHLA_A_2\tHLA_B_1\tHLA_B_2\tHLA_C_1\tHLA_C_2\tHLA_DPA1_1\tHLA_DPA1_2\tHLA_DPB1_1\tHLA_DPB1_2\tHLA_DQA1_1\tHLA_DQA1_2\tHLA_DQB1_1\tHLA_DQB1_2\tHLA_DRB1_1\tHLA_DRB1_2\n")
        print (self.sample, end = "\t", file = COUT)
        for gene in gene_list:
            print (format_allele(sample_results[gene][0]), format_allele(sample_results[gene][1]), sep = "\t", end = "\t", file = COUT)
        print ('', file = COUT)
        COUT.close()

    def check_pop(self, allele):
        allele = re.sub("HLA-","",allele)
        array = allele.split(":")
        two_field = array[0] + ":" + array[1]
        flag = False
        if two_field in  hashp:
            if hashp[two_field] > 0:
                flag = True
        return flag

def population(pop, wxs):
    hashp = {}
    with open(f"{freq}", "r") as fin:
        next(fin)
        for line in fin:
            gene, c, b, a = line.strip().split()
            if wxs == "exon":
                a = "%.3f" % float(a)
                b = "%.3f" % float(b)
                c = "%.3f" % float(c)
            elif wxs == "whole":
                a = "%.8f" % float(a)
                b = "%.8f" % float(b)
                c = "%.8f" % float(c)
            if pop == "Unknown":
                hashp[gene] = (float(a) + float(b) + float(c)) / 3
            elif pop == "Asian":
                hashp[gene] = float(a)
            elif pop == "Black":
                hashp[gene] = float(b)
            elif pop == "Caucasian":
                hashp[gene] = float(c)
            elif pop == "nonuse":
                hashp[gene] = 1
            else:
                print ("Wrong value for the parameter -p.")
                sys.exit(0)
    return hashp

def most_common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]

def format_allele(allele):
    
    if allele != "no_match":
        array = allele.split("_")
        new_name = array[0] + "*"
        new_name += ":".join(array[1:]) 
        return new_name
    else:
        return allele

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="G group resolution HLAtyping annotation", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    required.add_argument("-s", type=str, help="sample name", metavar="\b")
    required.add_argument("-i", type=str, help="the directory of phased sequence.", metavar="\b", default="./output")
    optional.add_argument("-p", type=str, help="The population of the sample [Asian, Black, Caucasian, Unknown, nonuse] for annotation. Unknown means use mean allele frequency in all populations. nonuse indicates only adopting mapping score and considering zero-frequency alleles.", metavar="\b", default="Unknown")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)


    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    exon_database = "%s/../../db/HLA//hla_exons.fasta"%(sys.path[0])
    freq = "%s/../../db/HLA/HLA_FREQ_HLA_I_II.txt"%(sys.path[0])

    print ("Start G group resolution annotation...")
    G_annotation_dict = read_G_annotation()
    hashp = population(args["p"], "whole")
    g_ann = G_annotation(args["s"], args["i"])
    g_ann.main()