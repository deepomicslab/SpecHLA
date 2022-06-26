"""
randomly select alleles from the database
two alleles for each gene

the allele is too short for some simulater, add additional seq to the allele

dependencies: 

wangshuai  June 25, 2022
"""

import random
import os

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

def read_fasta(file):
    # return the sequence saved in fasta file
    seq=''
    for line in open(file, 'r'):
        if line[0] == '>':
            continue
        line=line.strip()
        seq+=line
    return seq

class Sim_Tenx():

    def __init__(self):
        self.allele_file = "/mnt/d/HLAPro_backup/pacbio/allele_list.txt" # allele name list. in each line, example: >A*01:01
        self.database = "/mnt/d/HLAPro_backup/pacbio/hla_gen.format.filter.extend.DRB.no26789.v2.fasta" # allele database
        self.addi_seq_dir = "/mnt/d/HLAPro_backup/tenx/prefix_seq/" # folder saves left and right 10kb of each gene in hg19
        self.dir = "/mnt/d/HLAPro_backup/tenx/fasta" # the dir to store simulated data
        self.allele_dict = {}
        self.selected_alleles = []
        self.replicate_times = 1
        self.spechla_comman = "run_spechla.sh"
        self.addition_seq = {}
        self.allele_seq = {}
        self.allele_name = {}

    def read_addition_seq(self):
        for gene in gene_list:
            self.addition_seq[gene] = ['', '']
            for index in [1, 2]:
                seq_file = self.addi_seq_dir + "/HLA_%s_%s.fa"%(gene, index)
                seq = read_fasta(seq_file)
                self.addition_seq[gene][index - 1] = seq

    def get_all_allele(self):
        for line in open(self.allele_file):
            allele = line.strip()[1:]
            gene = allele.split("*")[0]
            if gene not in self.allele_dict:
                self.allele_dict[gene] = []
            self.allele_dict[gene].append(allele)

    def select_allele(self):
        self.selected_alleles = []
        for gene in gene_list:
            random.shuffle(self.allele_dict[gene])
            self.selected_alleles += self.allele_dict[gene][:2]
        print (self.selected_alleles)
    
    def get_sample_fasta(self, sample):
        sample_fasta = f"{self.dir}/{sample}.fasta"
        command  = f"samtools faidx {self.database} "
        for allele in self.selected_alleles:
            command += " " + allele
        command += ">" + sample_fasta
        os.system(command)
        self.get_allele_seq(sample_fasta)
        self.read_addition_seq()
        # print (self.addition_seq)
        self.add_addition(sample_fasta)

    def get_allele_seq(self, sample_fasta):
        for line in open(sample_fasta, 'r'):
            line=line.strip()
            if line[0] == '>':
                allele = line[1:]
                self.allele_seq[allele] = ''
            else:
                self.allele_seq[allele] += line

    def add_addition(self, sample_fasta):
        f = open(sample_fasta, 'w')
        for allele in self.allele_seq.keys():
            gene = allele.split("*")[0]
            self.allele_seq[allele] = self.addition_seq[gene][0] + self.allele_seq[allele] + self.addition_seq[gene][1]
            print (">"+allele, file = f)
            print (self.allele_seq[allele], file = f)
        f.close()

    def get_sample_fastq(self, sample):
        sample_prefix = f"{self.dir}/{sample}"
        command = f"""
        sample={sample}
        mkdir {self.dir}/$sample
        pbsim --data-type CLR --seed 88 --accuracy-mean 0.85 --accuracy-min 0.80 --prefix {self.dir}/$sample/$sample --depth 200 --model_qc model_qc_clr {self.dir}/$sample.fasta
        cat {self.dir}/$sample/{sample}_*fastq>{self.dir}/$sample/$sample.fastq
        rm {self.dir}/$sample/{sample}_*fastq
        gzip {self.dir}/$sample/$sample.fastq
        """
        os.system(command)

    def get_benchmark(self):
        print ("start")
        self.get_all_allele()
        for i in range(self.replicate_times):
            sample = f"sample_{i}"
            self.select_allele()
            self.get_sample_fasta(sample)
            # self.get_sample_fastq(sample)
    


sim = Sim_Tenx()
sim.get_benchmark()

        
