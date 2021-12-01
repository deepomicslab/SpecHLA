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

# truth_list = '/mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20201116/merge.groudtruth.txt'
truth_list = './merge.groudtruth.txt'
all_allele = './simulated.fasta'
gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
origin = '/mnt/d/HLAPro_backup/HLAPro/db/ref/hla.ref.extend.fa'

def select_alleles():
    f = open(truth_list, 'w')
    print ('Sample  A_1     A_2     B_1     B_2     C_1     C_2     DPA1_1  DPA1_2  DPB1_1  DPB1_2  DQA1_1  DQA1_2  DQB1_1  DQB1_2  DRB1_1  DRB1_2', file = f)
    alleles_dict = {}
    for line in open(all_allele):
        if line[0] == '>':
            line = line.strip()
            gene = line[1:].split('_')[0]
            if gene not in alleles_dict:
                alleles_dict[gene] = [line[1:]]
            else:
                alleles_dict[gene] += [line[1:]]
    for i in range(11, 51):
        child_allele = []
        father_allele = []
        mother_allele = []
        for gene in gene_list:
            alleles = alleles_dict[gene]
            while True:
                random_index = np.random.randint(len(alleles), size = 4)
                if alleles[random_index[0]].split('_')[1] == alleles[random_index[2]].split('_')[1] and alleles[random_index[0]] != alleles[random_index[2]]:
                    break
            #print (alleles[random_index[0]])
            child_allele += [alleles[random_index[0]], alleles[random_index[2]]]
            father_allele += [alleles[random_index[0]], alleles[random_index[1]]]
            mother_allele += [alleles[random_index[3]], alleles[random_index[2]]]
        # print (child_allele[0], convert_name(child_allele[0]))
        # print (father_allele)
        # print (mother_allele)
        a = [child_allele, father_allele, mother_allele]
        b = ['child', 'father', 'mother']
        for j in range(3):
            cont = b[j] + '_%s\t'%(i)
            for z in range(len(a[j])):
                cont += convert_name(a[j][z]) + '\t'
            print (cont, file = f)
    f.close()

def convert_name(name):
    array = name.split('_')
    new = array[0] + '*'
    new += ':'.join(array[1:])
    return new

def read_samples():
    for line in open(truth_list):
        array = line.strip().split()
        if array[0] == 'Sample':
            continue
        #if array[0] != 'HLA_2_T_50-50':
        #    continue
        sample = array[0]
        #order = 'samtools faidx hla_gen.format.filter.fasta '
        outdir = 'fq/%s'%(sample)
        os.system('mkdir %s'%(outdir))
        order = 'cat '
        #order = 'samtools faidx /home/wangmengyao/scripts/PHLAT/database/ref/hla_gen.format.filter.fasta '
        for arr in array[1:]:
            arr = arr.replace('*','_')
            arr = arr.replace(':','_')
            gene=arr.split('_')[0]
            os.system('sh ./simu.hla.sh %s %s'%(gene, arr))
            order += ' fasta/simu/fasta/%s.fasta '%(arr)
            #order += ' /mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20201116/simu/fasta/%s.fasta '%(arr)
            #order += ' %s '%(arr)
        order += '>%s/%s.fasta\n'%(outdir, sample)
        wgs_sim = '/home/wangmengyao/packages/DWGSIM/dwgsim -e 0 -E 0 -1 150 -2 150 -d 500 -C 10 -r 0 -o 1 %s/%s.fasta %s/%s_cov10\n'%(outdir, sample, outdir, sample)
        # sim = '/home/yuyonghan/PBSIM-PacBio-Simulator/src/pbsim --data-type CLR --seed 88 --accuracy-mean 0.85 --accuracy-min 0.80 --prefix fq/%s/%s --depth 200 --model_qc /home/yuyonghan/PBSIM-PacBio-Simulator/data/model_qc_clr %s/%s.fasta\n'%(sample,sample,outdir,sample)
        #sim = '/home/yuyonghan/PBSIM-PacBio-Simulator/src/pbsim --depth 200 --prefix fq/%s/%s --model_qc /home/yuyonghan/PBSIM-PacBio-Simulator/data/model_qc_ccs %s/%s.fasta\n'%(sample,sample,outdir,sample)
        merge_sim = ''#'cat %s/%s_*.fastq >%s/%s.fastq\ngzip -f %s/%s.fastq'%(outdir, sample, outdir, sample, outdir, sample)
        order = order + wgs_sim  + merge_sim
        #order += '\n sh bbmap.sh %s'%(sample)
        # print (order)
        os.system(order)
        #break

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

    snp_rate = 0.001
    seq_dict = read_fasta(origin)
    for index in range(50):
        sample = 'child_%s'%(index)
        s = open('fasta/%s.fasta'%(sample), 'w')
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
    os.system('cat fasta/* > %s'%(all_allele))

def generate_ngs():
    for index in range(50):
        sample = 'child_%s'%(index)
        wgs_sim = '/mnt/d/HLAPro_backup/insert/dwgsim -e 0 -E 0 -1 150 -2 150 -C 10 -r 0 -o 1 fasta/%s.fasta ./ngs/%s_cov10\n'%(sample, sample)
        os.system(wgs_sim)

def generate_parents():
    for index in range(50):
        child = 'child_%s'%(index)
        father = 'father_%s'%(index)
        mother = 'mother_%s'%(index)
        child_fasta = 'fasta/%s.fasta'%(child)
        father_fasta = 'fasta/%s.fasta'%(father)
        mother_fasta = 'fasta/%s.fasta'%(mother)
        os.system('samtools faidx %s'%(child_fasta))

        father_order = 'samtools faidx %s '%(child_fasta)
        for gene in gene_list:
            father_order += ' %s_%s_1 '%(gene, index)
        father_order += '>%s'%(father_fasta)
        os.system(father_order)
        os.system('cat %s >>%s'%(origin, father_fasta))
        wgs_sim = '/mnt/d/HLAPro_backup/insert/dwgsim -e 0 -E 0 -1 150 -2 150 -C 10 -r 0 -o 1 fasta/%s.fasta ./ngs/%s_cov10\n'%(father, father)
        os.system(wgs_sim)

        mother_order = 'samtools faidx %s '%(child_fasta)
        for gene in gene_list:
            mother_order += ' %s_%s_2 '%(gene, index)
        mother_order += '>%s'%(mother_fasta)
        os.system(mother_order)
        os.system('cat %s >>%s'%(origin, mother_fasta))
        wgs_sim = '/mnt/d/HLAPro_backup/insert/dwgsim -e 0 -E 0 -1 150 -2 150 -C 10 -r 0 -o 1 fasta/%s.fasta ./ngs/%s_cov10\n'%(mother, mother)
        os.system(wgs_sim)

def generate_parents_new():
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



# generate_fasta()
# generate_ngs()
generate_parents_new()
