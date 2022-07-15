"""
simulate samples with allele imbalance
the two haplotypes have different frequency

dependencies: samtools, dwgsim Version: 0.1.13

python sim_imbalance.py <allele database> <outdir> <number of sample>

wangshuai  June 13, 2022
"""

import sys
import os
import numpy as np
import random
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

class Fasta():
    def __init__(self):       
        self.database = database
        self.dir = outdir # the dir to store simulated data
        self.allele_file = "%s/allele_list.txt"%(outdir) # allele name list. in each line, example: >A*01:01
        if not os.path.isfile(self.allele_file):
            os.system(f"cat {self.database} |grep '>' >{self.allele_file}")
        self.allele_dict = {}
        self.selected_alleles = []
        self.replicate_times = sample_num
        self.record_allele_name = {}

    def get_all_allele(self):
        for line in open(self.allele_file):
            allele = line.strip()[1:]
            gene = allele.split("_")[0]
            if gene not in self.allele_dict:
                self.allele_dict[gene] = []
            self.allele_dict[gene].append(allele)

    def select_allele(self, sample):
        self.selected_alleles = []
        for gene in gene_list:
            random.shuffle(self.allele_dict[gene])
            self.selected_alleles += self.allele_dict[gene][:2]
            self.record_allele_name[self.allele_dict[gene][0]] = sample + "." + gene + "_1.fasta"
            self.record_allele_name[self.allele_dict[gene][1]] = sample + "." + gene + "_2.fasta"
        print (self.selected_alleles )
    
    def get_sample_fasta(self, sample):
        self.select_allele(sample)
        sample_fasta = f"{self.dir}/{sample}.fasta"
        command  = f"samtools faidx {self.database} "
        for allele in self.selected_alleles:
            os.system(f"samtools faidx {self.database} {allele} >{self.dir}/truth/{self.record_allele_name[allele]}")
            command += " " + allele
        command += ">" + sample_fasta
        os.system(command)

class Novel():

    def __init__(self):
        self.snp_rate = mutation_rate
        self.origin = origin
        self.seq_dict = self.read_fasta(self.origin)
        self.dir = outdir # the dir to store simulated data

    def add_snp(self, sequence, rate):
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

    def read_fasta(self, file):
        seq_dict = {}
        fasta_sequences = SeqIO.parse(open(file),'fasta')
        for record in fasta_sequences:
            seq_dict[str(record.id)] = str(str(record.seq))
        #     scaffold_name.append(str(record.id))
        return seq_dict

    def get_novel_fasta(self, sample):  
        # add snp to the reference allele in each gene
        # to generate the novel allele
        
        for h in range(1,3):
            s = open(outdir + '/%s.%s.fasta'%(sample, h), 'w')
            for gene in self.seq_dict.keys():
                gene_seq = self.seq_dict[gene]
                record_truth = f"{self.dir}/truth/{sample}.{gene}_{h}.fasta"
                if gene != 'HLA_DRB1':
                    novel_seq = gene_seq[:1050] + self.add_snp(gene_seq[1050:-1050], self.snp_rate) + gene_seq[-1050:]
                else:
                    novel_seq = gene_seq[:1050] + self.add_snp(gene_seq[1050:3800], self.snp_rate)\
                         + self.add_snp(gene_seq[4500:-1050], self.snp_rate) + gene_seq[-1050:]

                print('>novel_'+ gene.split('_')[1] + '_' + str(h) , file = s)
                print(novel_seq, file = s)

                w = open(record_truth, 'w')
                print ('>novel_'+ gene.split('_')[1] + '_' + str(h) , file = w)
                print(novel_seq, file = w)
                w.close()

            s.close()

class Fastq():

    def __init__(self, depth_low, depth_high):
        self.dir = outdir # the dir to store simulated data
        self.depth_low = depth_low
        self.depth_high = depth_high
        self.replicate_times = sample_num
        self.read_len = read_length

    def get_illumina(self, sample):
        command = f"""
        sample={sample}
        mkdir {self.dir}/$sample
        {dwgsim_script} -r 0 -e 0 -E 0 -1 {self.read_len} -2 {self.read_len} -C {self.depth_low}\
             {self.dir}/$sample.1.fasta {self.dir}/$sample/$sample.1.illumina
        {dwgsim_script} -r 0 -e 0 -E 0 -1 {self.read_len} -2 {self.read_len} -C {self.depth_high}\
             {self.dir}/$sample.2.fasta {self.dir}/$sample/$sample.2.illumina
        zcat {self.dir}/$sample/$sample.1.illumina.bwa.read1.fastq.gz {self.dir}/$sample/$sample.2.illumina.bwa.read1.fastq.gz >{self.dir}/$sample/$sample.illumina.bwa.read1.fastq
        zcat {self.dir}/$sample/$sample.1.illumina.bwa.read2.fastq.gz {self.dir}/$sample/$sample.2.illumina.bwa.read2.fastq.gz >{self.dir}/$sample/$sample.illumina.bwa.read2.fastq
        gzip -f {self.dir}/$sample/$sample.illumina.bwa.read*.fastq
        """
        os.system(command)  

def simulate():
    print ("start simulation")
    # fa = Fasta()
    # fa.get_all_allele()
    # fa.get_sample_fasta(sample)
    
    no = Novel()
    
    for depth_set in [[20, 80], [30, 70], [40, 60], [48, 52], [50, 50]]:
        depth_low, depth_high = depth_set[0], depth_set[1]
        fq = Fastq(depth_low, depth_high)
        for i in range(sample_num):
            prefix = "imbalance_" + "%s_%s"%(depth_low, depth_high)
            sample = f"{prefix}_{i}"
            print (sample, depth_set)
            no.get_novel_fasta(sample)
            fq.get_illumina(sample)

if __name__ == "__main__":  
    dwgsim_script = "/mnt/d/HLAPro_backup/insert/dwgsim"
    origin = '/mnt/d/HLAPro_backup/HLAPro/db/ref/hla.ref.extend.fa'

    mutation_rate = 0.002
    read_length = 75

    database = sys.argv[1]
    outdir = sys.argv[2]
    sample_num = int(sys.argv[3])

    
    simulate()

