"""
randomly select alleles from the database
two alleles for each gene
simulate pacbio reads with pbsim

dependencies: pbsim, samtools, allele_database

wangshuai  June 21, 2022
"""

import random
import os

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

class Sim_Pac():

    def __init__(self):
        self.allele_file = "/mnt/d/HLAPro_backup/pacbio/allele_list.txt" # allele name list. in each line, example: >A*01:01
        self.database = "/mnt/d/HLAPro_backup/pacbio/hla_gen.format.filter.extend.DRB.no26789.v2.fasta" # allele database
        self.dir = "simulation" # the dir to store simulated data
        self.allele_dict = {}
        self.selected_alleles = []
        self.replicate_times = 10
        self.spechla_comman = "run_spechla.sh"

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
        print (self.selected_alleles )
    
    def get_sample_fasta(self, sample):
        sample_fasta = f"{self.dir}/{sample}.fasta"
        command  = f"samtools faidx {self.database} "
        for allele in self.selected_alleles:
            command += " " + allele
        command += ">" + sample_fasta
        os.system(command)

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
            sample = f"pacbio_{i}"
            self.select_allele()
            self.get_sample_fasta(sample)
            self.get_sample_fastq(sample)
    
    def get_spechla(self):
        # f = open(self.spechla_comman, 'w')
        # for i in range(self.replicate_times):
        for i in range(1, self.replicate_times):
            sample = f"pacbio_{i}"     
            print (f"run {sample}...")   
            command = f"python ../../script/long_read_typing.py -n {sample} -r {self.dir}/{sample}/{sample}.fastq.gz -p unuse -j 10"
            print (command)
            os.system(command)


sim_pacbio = Sim_Pac()
# sim_pacbio.get_benchmark()
sim_pacbio.get_spechla()
        
