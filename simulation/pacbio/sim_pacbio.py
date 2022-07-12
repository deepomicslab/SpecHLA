"""
randomly select alleles from the database
two alleles for each gene
simulate sequencing reads with different platfors, like
Pacbio, Nanopore, Hi-c, Illumina, 10x

dependencies: pbsim, samtools, dwgsim Version: 0.1.13

python sim_pacbio.py /mnt/d/HLAPro_backup/haplotype/hla_gen.format.filter.fasta /mnt/d/HLAPro_backup/haplotype/pac_sim_test 10 2
python sim_pacbio.py <allele database> <outdir> <depth> <number of sample>

wangshuai  June 21, 2022
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
        self.snp_rate = 0.002
        self.origin = '/mnt/d/HLAPro_backup/HLAPro/db/ref/hla.ref.extend.fa'
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
        s = open(outdir + '/%s.fasta'%(sample), 'w')
        for gene in self.seq_dict.keys():
            gene_seq = self.seq_dict[gene]
            for h in range(1,3):
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
        # os.system('cat %s/fasta/* > %s'%(outdir, all_allele))


class Fastq():

    def __init__(self):
        self.dir = outdir # the dir to store simulated data
        self.depth = round(int(depth)/2)
        self.replicate_times = sample_num

    def get_pacbio(self, sample):
        command = f"""
        sample={sample}
        mkdir {self.dir}/$sample
        pbsim --data-type CLR --seed 88 --accuracy-mean 0.85 --accuracy-min 0.80 --prefix {self.dir}/$sample/$sample --depth {self.depth} --model_qc {sys.path[0]}/model_qc_clr {self.dir}/$sample.fasta
        cat {self.dir}/$sample/{sample}_*fastq>{self.dir}/$sample/$sample.fastq
        rm {self.dir}/$sample/{sample}_*fastq
        gzip {self.dir}/$sample/$sample.fastq
        """
        os.system(command)

    def get_nanopore(self, sample):
        command = """
            sample={sample}
            /mnt/e/hla_tgs/nanopore/DeepSimulator/deep_simulator.sh\
                -i {self.dir}/$sample.fasta\
                -l 3000 -B 2 -c 6 -K 10 -P 0\
                -o {self.dir}/nanopore/$sample
        """
        os.system(command)

    def get_illumina(self, sample):
        command = f"""
        sample={sample}
        mkdir {self.dir}/$sample
        /mnt/d/HLAPro_backup/insert/dwgsim -r 0 -e 0 -E 0 -1 100 -2 100 -C {self.depth} {self.dir}/$sample.fasta {self.dir}/$sample/$sample.illumina
        # gzip -f {self.dir}/$sample/$sample.illumina.*fastq
        """
        #-e 0.01 -E 0.01 -d 200 -r 0 -1 70 -2 70
        os.system(command)  
    
    def get_hic(self, sample):
        command = f"""
        sample={sample}
        outdir={self.dir}/hic/
        sim3C --dist uniform -n 10000 -l 150 -e NlaIII -m hic {self.dir}/$sample.fasta $outdir/$sample.fastq --simple-reads --insert-mean 500 --anti-rate 0 --trans-rate 0 --spurious-rate 0 -r 88 
        rm $outdir/profile.tsv
        python3 {sys.path[0]}/../split_hic.py $outdir/$sample.fastq $sample hic
        """
        os.system(command)

def simulate():
    print ("start simulation")
    fa = Fasta()
    fq = Fastq()
    no = Novel()
    fa.get_all_allele()
    for i in range(sample_num):
        sample = f"{prefix}_{i}"
        fa.get_sample_fasta(sample)
        no.get_novel_fasta(sample)
        fq.get_pacbio(sample)  
        fq.get_illumina(sample)

if __name__ == "__main__":  

    database = sys.argv[1]
    outdir = sys.argv[2]
    depth = sys.argv[3]
    sample_num = int(sys.argv[4])

    prefix = "novel"
    simulate()
    # sim_pacbio = Sim_Pac()
    # sim_pacbio.get_benchmark()
    # sim_pacbio.get_illumina()
    # sim_pacbio.get_spechla()
