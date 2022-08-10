"""
randomly select alleles from the database
two alleles for each gene
simulate sequencing reads with different platfors, like
Pacbio, Nanopore, Hi-c, Illumina, 10x

dependencies: pbsim, sim3C, samtools, DeepSimulator, dwgsim Version: 0.1.13, Nanosim

python sim_diff_platforms.py <allele database> <outdir> <depth> <number of sample>

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

def read_G_annotation():
    g_file = "%s/../db/HLA/hla_nom_g.txt"%(sys.path[0])
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
            record_selected_allele = []
            for allele in self.allele_dict[gene]:
                if allele not in G_annotation_dict:
                    continue
                if G_annotation_dict[allele][-1] == "G":
                    record_selected_allele.append(allele)
                if len(record_selected_allele) == 2:
                    break
            self.selected_alleles += record_selected_allele
            self.record_allele_name[record_selected_allele[0]] = sample + "." + gene + "_1.fasta"
            self.record_allele_name[record_selected_allele[1]] = sample + "." + gene + "_2.fasta"
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
    
    def get_parent_fasta(self, parent_name, index):
        father_alleles = []
        for i in range(index, len(self.selected_alleles), 2):
            inherit_allele = self.selected_alleles[i]
            gene = inherit_allele.split("_")[0]
            random.shuffle(self.allele_dict[gene])
            not_inherit_allele = self.allele_dict[gene][0]
            father_alleles += [inherit_allele, not_inherit_allele]
        command  = f"samtools faidx {self.database} "
        for allele in father_alleles:
            # os.system(f"samtools faidx {self.database} {allele} >{self.dir}/truth/{self.record_allele_name[allele]}")
            command += " " + allele
        sample_fasta = f"{self.dir}/{parent_name}.fasta"
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

class Fastq():

    def __init__(self):
        self.dir = outdir # the dir to store simulated data
        self.depth = round(int(depth)/2)
        self.replicate_times = sample_num
        self.read_len = read_length

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

    def get_nanopore_old(self, sample):
        command = f"""
            sample={sample}
            {deep_simulator_script} -i {self.dir}/$sample.fasta -o {self.dir}/nanopore/$sample\
                 -S 1 -l 3000 -B 1 -c 12 -K {self.depth} -P 0 -M 0      
        """
        os.system(command)

    def get_nanopore(self, sample):
        command = f"""
            sample={sample}
            simulator.py genome -rg {self.dir}/$sample.fasta -o {self.dir}/nanopore/$sample \
                -c /mnt/d/HLAPro_backup/pacbio/simulation/human_NA12878_DNA_FAB49712_guppy/training -max 5000 -n 4000  --seed 66
            python {sys.path[0]}/fasta2fastq.py  {self.dir}/nanopore/{sample}_aligned_reads.fasta > {self.dir}/nanopore/{sample}_aligned_reads.fastq   
        """# Nanosim
        print (command)
        os.system(command)

    def get_illumina(self, sample):
        command = f"""
        sample={sample}
        mkdir {self.dir}/$sample
        {dwgsim_script} -r 0 -e {read_error} -E {read_error} -d {frag_size} -1 {self.read_len} -2 {self.read_len} -C {self.depth}\
             {self.dir}/$sample.fasta {self.dir}/$sample/$sample.illumina
        # gzip -f {self.dir}/$sample/$sample.illumina.*fastq
        """
        os.system(command)  
    
    def get_hic(self, sample):
        command = f"""
        sample={sample}
        outdir={self.dir}/hic/
        sim3C --dist uniform -n 10000 -l {self.read_len} -e NlaIII -m hic {self.dir}/$sample.fasta $outdir/$sample.fastq \
            --simple-reads --insert-mean 500 --anti-rate 0 --trans-rate 0 --spurious-rate 0 -r 88 
        rm $outdir/profile.tsv
        python3 {sys.path[0]}/split_hic.py $outdir/$sample.fastq $sample $outdir
        """
        os.system(command)

def simulate():
    print ("start simulation")
    # fa = Fasta()
    # fa.get_all_allele()
    # fa.get_sample_fasta(sample)
    fq = Fastq()
    no = Novel()
    
    for i in range(sample_num):
        sample = f"{prefix}_{i}"
        print (sample)
        # no.get_novel_fasta(sample)
        # fq.get_pacbio(sample)  
        # fq.get_illumina(sample)
        # fq.get_nanopore(sample)
        # fq.get_hic(sample)

def simulate_hybrid():
    print ("start simulation")
    fa = Fasta()
    fa.get_all_allele()
    
    fq = Fastq()
    # no = Novel()
    
    for i in range(sample_num):
        sample = f"{prefix}_{i}"
        print (sample)
        fa.get_sample_fasta(sample)
        # no.get_novel_fasta(sample)
        fq.get_pacbio(sample)  
        fq.get_illumina(sample)
        # fq.get_nanopore(sample)
        # fq.get_hic(sample)

def simulate_trio():
    print ("start simulation")
    fa = Fasta()
    fa.get_all_allele()
    
    fq = Fastq()
    # no = Novel()
    
    for i in range(sample_num):
        sample = f"{prefix}_{i}"
        print (sample)
        fa.get_sample_fasta(sample)
        # fq.get_illumina(sample)
        father_name = "father_%s"%(i)
        fa.get_parent_fasta(father_name, 0)
        mother_name = "mother_%s"%(i)
        fa.get_parent_fasta(mother_name, 1)

        fq.get_illumina(sample)
        fq.get_illumina(father_name)
        fq.get_illumina(mother_name)

if __name__ == "__main__":  

    deep_simulator_script = "/mnt/e/hla_tgs/nanopore/DeepSimulator/deep_simulator.sh"
    dwgsim_script = "/mnt/d/HLAPro_backup/insert/dwgsim"
    origin = '/mnt/d/HLAPro_backup/HLAPro/db/ref/hla.ref.extend.fa'
    G_annotation_dict = read_G_annotation()
    
    mutation_rate = 0.001


    read_length = 75 
    frag_size = 200 
    read_error = 0.01
    prefix = "hybrid"

    # read_length = 75 #150
    # frag_size = 300 #500
    # read_error = 0.01
    # prefix = "child"


    

    database = sys.argv[1]
    outdir = sys.argv[2]
    depth = sys.argv[3]
    sample_num = int(sys.argv[4])

    if not os.path.isdir(outdir + "/truth/"):
        os.system("mkdir %s/truth"%(outdir))
    # read_error = 0
    # prefix = "novel"
    # simulate()

    simulate_hybrid()
    # simulate_trio()


