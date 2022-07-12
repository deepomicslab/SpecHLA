"""
To test the full-resolution HLA typing performance of SpecHLA 
after integrating different data types.
We randomly add SNP to the reference allele, and generate different 
datatypes for the sample.

wangshuai July 11, 2022
"""


import os 
import numpy as np
import os
import numpy as np
import random
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re



def add_snp(sequence, rate):
    locus_set = {}
    allele = ['A', 'T', 'C', 'G']
    ref_len = len(sequence)
    num = int(ref_len* rate)
    i = 0
    #'''
    while i < num:
        random_locus = np.random.randint(ref_len - 1 )
        if random_locus not in locus_set.keys():
            ref = sequence[random_locus]
            while True:
                allele_index = np.random.randint(4)
                alt = allele[allele_index]
                if alt != ref:
                    break
            # print (random_locus, ref, alt)
            sequence = sequence[:random_locus] + alt + sequence[random_locus+1:]
            # locus_set.append(random_locus)
            locus_set[random_locus] = 1
            i += 1
    return sequence

def read_fasta(file):
    seq_dict = {}
    fasta_sequences = SeqIO.parse(open(file),'fasta')
    for record in fasta_sequences:
        seq_dict[str(record.id)] = str(str(record.seq))
    #     scaffold_name.append(str(record.id))
    return seq_dict

def generate_fasta():  
    # add snp to the reference allele in each gene
    # to generate the novel allele
    seq_dict = read_fasta(origin)
    for index in range(50):
        sample = 'child_%s'%(index)
        s = open(outdir + 'fasta/%s.fasta'%(sample), 'w')
        for gene in seq_dict.keys():
            gene_seq = seq_dict[gene]
            for h in range(1,3):
                if gene != 'HLA_DRB1':
                    novel_seq = gene_seq[:1050] + add_snp(gene_seq[1050:-1050], snp_rate) + gene_seq[-1050:]
                else:
                    novel_seq = gene_seq[:1050] + add_snp(gene_seq[1050:3800], snp_rate) + add_snp(gene_seq[4500:-1050], snp_rate) + gene_seq[-1050:]
                print('>'+ gene.split('_')[1] + '_' + str(index) + '_' + str(h) , file = s)
                print(novel_seq, file = s)
    s.close()
    os.system('cat %s/fasta/* > %s'%(outdir, all_allele))

def generate_illumina():
    for index in range(sample_num):
        sample = 'child_%s'%(index)
        wgs_sim = """dwgsim -r 0 -1 100 -2 100 -C %s %s/fasta/%s.fasta %s/illumina/%s_illumina
                     gzip -f %s/illumina/%s*illumina.bwa.read*.fastq"""%(depth,outdir,sample,outdir,sample,outdir,sample)
        os.system(wgs_sim)

def generate_trio():
    snp_rate = 0.001
    seq_dict = read_fasta(origin)

    for index in range(50):
        child = 'child_%s'%(index)
        father = 'new_father_%s'%(index)
        mother = 'new_mother_%s'%(index)
        child_fasta = 'fasta/%s.fasta'%(child)
        father_fasta = 'fasta/%s.fasta'%(father)
        mother_fasta = 'fasta/%s.fasta'%(mother)
        os.system('samtools faidx %s'%(child_fasta))

        father_order = 'samtools faidx %s '%(child_fasta)
        for gene in gene_list:
            father_order += ' %s_%s_1 '%(gene, index)
        father_order += '>%s'%(father_fasta)
        os.system(father_order)

        #the other hap
        s = open(father_fasta, 'a')
        for gene in seq_dict.keys():
            gene_seq = seq_dict[gene]
            if gene != 'HLA_DRB1':
                novel_seq = gene_seq[:1050] + add_snp(gene_seq[1050:-1050], snp_rate) + gene_seq[-1050:]
            else:
                novel_seq = gene_seq[:1050] + add_snp(gene_seq[1050:3800], snp_rate) + add_snp(gene_seq[4500:-1050], snp_rate) + gene_seq[-1050:]
            print('>'+ gene.split('_')[1] + '_' + str(index) + '_f' , file = s)
            print(novel_seq, file = s)
        s.close()


        wgs_sim = '/mnt/d/HLAPro_backup/insert/dwgsim -e 0 -E 0 -1 150 -2 150 -C 10 -r 0 -o 1 fasta/%s.fasta ./ngs/%s_cov10\n'%(father, father)
        os.system(wgs_sim)

        mother_order = 'samtools faidx %s '%(child_fasta)
        for gene in gene_list:
            mother_order += ' %s_%s_2 '%(gene, index)
        mother_order += '>%s'%(mother_fasta)
        os.system(mother_order)

        #the other hap
        s = open(mother_fasta, 'a')
        for gene in seq_dict.keys():
            gene_seq = seq_dict[gene]
            if gene != 'HLA_DRB1':
                novel_seq = gene_seq[:1050] + add_snp(gene_seq[1050:-1050], snp_rate) + gene_seq[-1050:]
            else:
                novel_seq = gene_seq[:1050] + add_snp(gene_seq[1050:3800], snp_rate) + add_snp(gene_seq[4500:-1050], snp_rate) + gene_seq[-1050:]
            print('>'+ gene.split('_')[1] + '_' + str(index) + '_m' , file = s)
            print(novel_seq, file = s)
        s.close()


        wgs_sim = '/mnt/d/HLAPro_backup/insert/dwgsim -e 0 -E 0 -1 150 -2 150 -C 10 -r 0 -o 1 fasta/%s.fasta ./ngs/%s_cov10\n'%(mother, mother)
        os.system(wgs_sim)


if __name__ == "__main__":  
    all_allele = '/mnt/d/HLAPro_backup/data_types/simulated.fasta'
    outdir = "/mnt/d/HLAPro_backup/data_types/"

    gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
    origin = '/mnt/d/HLAPro_backup/HLAPro/db/ref/hla.ref.extend.fa'
    snp_rate = 0.005
    depth = 30
    sample_num = 1
    
    generate_fasta()
    generate_illumina()
    # generate_trio()
