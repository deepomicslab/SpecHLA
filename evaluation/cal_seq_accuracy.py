"""
calculate base error and gap percentage
step1: map the sequence to the truth with blastn
step2: get mapped intervals
step3: get unique intervals
step4: calculate the mismatch rate in the unique mapped region
step5: calculate gap recall and gap precision

dependency: Blastn 2.6.0+, samtools 1.14

wangshuai July 1, 2022
"""

import os
import sys
import pandas as pd
from Bio import SeqIO

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['DQB1']

class Align(object):

    def __init__(self, mapped_len,infer_hap_len,truth_hap_len,mismatch_num,gap_open_num,true_mapped_len):
        self.mapped_len, self.infer_hap_len, self.truth_hap_len, self.mismatch_num, self.gap_open_num,\
            self.true_mapped_len=mapped_len,infer_hap_len,truth_hap_len,mismatch_num,gap_open_num,true_mapped_len

        self.short_gap_error = self.gap_open_num/self.mapped_len
        self.base_error = self.mismatch_num/self.mapped_len
        self.gap_recall = self.true_mapped_len/self.truth_hap_len
        # self.gap_recall = self.mapped_len/self.truth_hap_len
        self.gap_precision = self.mapped_len/self.infer_hap_len
          

class Seq_error():

    def __init__(self, infer_hap_file, truth_hap_file):
        self.infer_hap_file = infer_hap_file
        self.truth_hap_file = truth_hap_file
        self.blast_file = f"{self.infer_hap_file}.blast"
        self.infer_hap_len = None
        self.truth_hap_len = None
        self.mapped_interval = [] # for inferred hap
        self.true_mapped_interval = [] # for true hap
        self.mapped_len = 0
        self.true_mapped_len = 0
        self.mismatch_num = 0
        self.gap_open_num = 0

    def blast_map(self):
        command = f"""
        blastn -query {self.infer_hap_file} -out {self.blast_file} -subject {self.truth_hap_file} \
            -outfmt 7 
        # cat  {self.blast_file}
        """
        # print (command)
        os.system(command)
    
    def get_fasta_len(self):
        with open(self.infer_hap_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.infer_hap_len = len(record.seq)
        with open(self.truth_hap_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.truth_hap_len = len(record.seq)

    def read_blast(self):
        f = open(self.blast_file, 'r')
        i = 0
        for line in f:
            if line[0] == "#":
                continue
            array = line.strip().split()
            map_s = int(array[6])
            map_e = int(array[7])
            mismatches = int(array[4])
            gap_opens = int(array[5])
            # print (map_s, map_e)
            self.mapped_interval = self.get_uniq_map(map_s, map_e, mismatches, gap_opens, self.mapped_interval)
            true_map_s = int(array[8])
            true_map_e = int(array[9])
            self.true_mapped_interval = self.get_uniq_map(true_map_s, true_map_e, 0, 0, self.true_mapped_interval)
            i += 1
            # if i > 5:
            #     break
        # print (self.mapped_interval)
        # print (self.true_mapped_interval)

    def get_uniq_map(self, map_s, map_e, mismatches, gap_opens, mapped_interval):
        if map_e < map_s:
            a = map_e
            map_e = map_s
            map_s = a
        flag = True
        origin_map_len = map_e - map_s + 1
        for interval in mapped_interval:
            if map_s >= interval[0] and map_s <= interval[1]:
                if map_e > interval[1]:
                    map_s = interval[1] + 1
            if map_e >= interval[0] and map_e <= interval[1]:
                if map_s < interval[0]:
                    map_e = interval[0] - 1
            if map_s >= interval[0] and map_e <= interval[1]:
                flag = False
        if flag and map_e - map_s > 5:
            uniq_map_len = map_e - map_s + 1
            mismatches = mismatches * (uniq_map_len/origin_map_len)
            gap_opens = gap_opens * (uniq_map_len/origin_map_len)
            self.mismatch_num += mismatches
            self.gap_open_num += gap_opens
            mapped_interval.append([map_s, map_e])
        return mapped_interval

    def get_gap_per(self):
        for interval in self.mapped_interval:
            self.mapped_len = self.mapped_len + (interval[1] - interval[0] + 1)
        for interval in self.true_mapped_interval:
            self.true_mapped_len = self.true_mapped_len + (interval[1] - interval[0] + 1)
        # print (self.mapped_len, self.infer_hap_len, self.truth_hap_len)
        # print (self.mismatch_num, self.gap_open_num)
        align = Align(self.mapped_len, self.infer_hap_len, self.truth_hap_len, \
            self.mismatch_num, self.gap_open_num,self.true_mapped_len)
        return align

    def main(self):
        self.get_fasta_len()
        self.blast_map()
        self.read_blast()
        align = self.get_gap_per()
        # print (self.mapped_interval, self.infer_hap_len, self.truth_hap_len, self.mapped_len, align.base_error)
        return align

def eva_HG002_spechla():
    outdir = "/mnt/d/HLAPro_backup/trio/HG002/"
    truth_file1 = "/mnt/d/HLAPro_backup/trio/truth_MHC/H1-asm.fa"
    truth_file2 = "/mnt/d/HLAPro_backup/trio/truth_MHC/H2-asm.fa"
    data = []
    for gene in gene_list:
        infer_file1 = outdir + "hla.allele.1.HLA_%s.fasta"%(gene)
        infer_file2 = outdir + "hla.allele.2.HLA_%s.fasta"%(gene)
        seq = Seq_error(infer_file1, truth_file1)
        align_11 = seq.main()
        seq = Seq_error(infer_file2, truth_file2)
        align_22 = seq.main()
        seq = Seq_error(infer_file1, truth_file2)
        align_12 = seq.main()
        seq = Seq_error(infer_file2, truth_file1)
        align_21 = seq.main()
        if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
            choose_align1 = align_11
            choose_align2 = align_22
        else:
            choose_align1 = align_12
            choose_align2 = align_21
        # print ("#", align_11.base_error, align_22.base_error, align_12.base_error, align_21.base_error)
        base_error = (choose_align1.base_error + choose_align2.base_error)/2
        short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
        gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
        gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
        data.append([base_error, short_gap_error, gap_recall, gap_precision, gene])
        print (gene, base_error, short_gap_error, gap_recall, gap_precision)
        # break
    df = pd.DataFrame(data, columns = ["base_error", "short_gap_error", "gap_recall", "gap_precision", "Gene"])
    df.to_csv('/mnt/d/HLAPro_backup/trio/hg002_haplo_assess.csv', sep=',')

def eva_HG002_hisat():
    outdir = "/mnt/d/HLAPro_backup/trio/Trio/hisat/HG002/"
    hisat_fasta = "/mnt/d/HLAPro_backup/trio/Trio/hisat/HG002/assembly_graph-hla.HG002_GRCh38_2x250_extract_1_fq_gz-hla-extracted-1_fq.fasta"
    contig2gene = "HG002_hisat_contig2gene.txt"
    truth_file1 = "/mnt/d/HLAPro_backup/trio/truth_MHC/H1-asm.fa"
    truth_file2 = "/mnt/d/HLAPro_backup/trio/truth_MHC/H2-asm.fa"

    true_allele = {}
    f = open(contig2gene)
    for line in f:
        array = line.strip().split()
        if array[1] in gene_list:
            if array[1] not in true_allele:
                true_allele[array[1]] = []
            true_allele[array[1]].append(array[0])
    f.close()    

    data = []
    for gene in gene_list:
        infer_file1 = outdir + "hisat.hla.allele.1.HLA_%s.fasta"%(gene)
        infer_file2 = outdir + "hisat.hla.allele.2.HLA_%s.fasta"%(gene)
        alleles = true_allele[gene]

        with open(hisat_fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == alleles[0]:
                    with open(infer_file1, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                    break

        with open(hisat_fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == alleles[1]:
                    with open(infer_file2, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                    break

        seq = Seq_error(infer_file1, truth_file1)
        align_11 = seq.main()
        seq = Seq_error(infer_file2, truth_file2)
        align_22 = seq.main()
        seq = Seq_error(infer_file1, truth_file2)
        align_12 = seq.main()
        seq = Seq_error(infer_file2, truth_file1)
        align_21 = seq.main()
        if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
            choose_align1 = align_11
            choose_align2 = align_22
        else:
            choose_align1 = align_12
            choose_align2 = align_21
        # print ("#", align_11.base_error, align_22.base_error, align_12.base_error, align_21.base_error)
        base_error = (choose_align1.base_error + choose_align2.base_error)/2
        short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
        gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
        gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
        data.append([base_error, short_gap_error, gap_recall, gap_precision, gene])
        print (gene, base_error, short_gap_error, gap_recall, gap_precision)
        # break
    df = pd.DataFrame(data, columns = ["base_error", "short_gap_error", "gap_recall", "gap_precision", "Gene"])
    df.to_csv('/mnt/d/HLAPro_backup/trio/hisat_hg002_haplo_assess.csv', sep=',')

def eva_simu(database, record_true_file, outdir):
    outdir = outdir + "/"
    true_allele = {}
    f = open(record_true_file)
    for line in f:
        array = line.strip().split()
        if array[1] in gene_list:
            true_allele[array[1]] = [array[3], array[4]]
    f.close()
    
    data = []
    for gene in gene_list:
        truth_file1 = outdir + ".true.%s.1.fasta"%(gene)
        truth_file2 = outdir + ".true.%s.2.fasta"%(gene)
        alleles = true_allele[gene]
        command = f"samtools faidx {database} {alleles[0]} >{truth_file1}"
        os.system(command)
        command = f"samtools faidx {database} {alleles[1]} >{truth_file2}"
        os.system(command)

        infer_file1 = outdir + "hla.allele.1.HLA_%s.fasta"%(gene)
        infer_file2 = outdir + "hla.allele.2.HLA_%s.fasta"%(gene)
        seq = Seq_error(infer_file1, truth_file1)
        align_11 = seq.main()
        seq = Seq_error(infer_file2, truth_file2)
        align_22 = seq.main()
        seq = Seq_error(infer_file1, truth_file2)
        align_12 = seq.main()
        seq = Seq_error(infer_file2, truth_file1)
        align_21 = seq.main()
        if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
            choose_align1 = align_11
            choose_align2 = align_22
        else:
            choose_align1 = align_12
            choose_align2 = align_21
        # print (align_11.base_error, align_22.base_error, align_12.base_error, align_21.base_error)
        base_error = (choose_align1.base_error + choose_align2.base_error)/2
        short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
        gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
        gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
        data.append([base_error, short_gap_error, gap_recall, gap_precision, gene])
        print (gene, base_error, short_gap_error, gap_recall, gap_precision)
        # break
    df = pd.DataFrame(data, columns = ["base_error", "short_gap_error", "gap_recall", "gap_precision", "Gene"])
    df.to_csv('%s/haplotype_assessment.csv'%(outdir), sep=',')

if __name__ == "__main__":

    eva_HG002_spechla()
    # eva_HG002_hisat()

    # database = sys.argv[1]
    # record_true_file = sys.argv[2]
    # sample_dir = sys.argv[3]
    # eva_simu(database, record_true_file, sample_dir)
    
