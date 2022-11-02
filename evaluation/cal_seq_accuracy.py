"""
assess reconstructed sequences
step1: map the sequence to the truth with blastn
step2: get mapped intervals
step3: get unique intervals
step4: calculate the mismatch rate in the unique mapped region
step5: calculate sequence recall and sequence precision

python cal_seq_accuracy.py hla_gen.format.filter.fasta simu_D100_L100_E0.truth.txt ./sample10 sample10

dependency: Blastn 2.6.0+, samtools 1.14, MUSCLE v3.8.1551

wangshuai July 1, 2022
"""

import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import re

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['DQA1', 'DQB1', 'DRB1']

def merge_exon(infer_file, gene, index):
    used_exon = ["HLA_A:1504-1773","HLA_A:2015-2290","HLA_B:1486-1755","HLA_B:2001-2276","HLA_C:1699-1968","HLA_C:2215-2490",\
        "HLA_DQA1:5600-5848","HLA_DQB1:3073-3342","HLA_DRB1:6972-7241","HLA_DPA1:5208-5453","HLA_DPB1:6002-6265"]
    merge_exon_file = infer_file + ".merged.exon.fasta"
    exon_seq = ''
    with open(infer_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in used_exon:# or record.id.split(":")[0] in ["HLA_DPA1","HLA_DPB1"]:
                exon_seq += str(record.seq)
                exon_seq += "N"*200
    with open(merge_exon_file, "w") as output_handle:
        print (">HLA_%s_%s"%(gene, index), file = output_handle)  
        print (exon_seq, file = output_handle) 
    return  merge_exon_file

def splitN_for_kourami(infer_file, gene, index):
    splitN_file = infer_file + ".split.exon.fasta"
    with open(infer_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            merged_exon = str(record.seq)
    if len(merged_exon) > 270:
        new_exon = merged_exon[:270] + "N"*200 + merged_exon[270:]
    else:
        new_exon = merged_exon
    # print (new_exon)
    output_handle = open(splitN_file, "w")
    print (">HLA_%s_%s"%(gene, index), file = output_handle)  
    print (new_exon, file = output_handle) 
    output_handle.close()
    return splitN_file

class Align(object):

    def __init__(self, mapped_len,infer_hap_len,truth_hap_len,mismatch_num,gap_open_num,true_mapped_len):
        self.mapped_len, self.infer_hap_len, self.truth_hap_len, self.mismatch_num, self.gap_open_num,\
            self.true_mapped_len=mapped_len,infer_hap_len,truth_hap_len,mismatch_num,gap_open_num,true_mapped_len

        if self.mapped_len > 0:
            self.short_gap_error = self.gap_open_num/self.mapped_len
            self.base_error = self.mismatch_num/self.mapped_len
            self.gap_recall = self.true_mapped_len/self.truth_hap_len
            # self.gap_recall = self.mapped_len/self.truth_hap_len
            if self.mapped_len > self.infer_hap_len:
                self.gap_precision = 1
                print ("larger than 1 caused by short N sequence.", self.mapped_len, self.infer_hap_len)
            else:
                self.gap_precision = self.mapped_len/self.infer_hap_len
        else:
            self.short_gap_error = 1
            self.base_error = 1
            self.gap_recall = 0
            # self.gap_recall = self.mapped_len/self.truth_hap_len
            self.gap_precision = 0
          
class Seq_error():

    def __init__(self, infer_hap_file, truth_hap_file, gene):
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
        self.gene = gene

    def blast_map(self, flag):
        #-penalty -1 -reward 1 -gapopen 4 -gapextend 1 -strand plus
        if flag == "strict":
            command = f"""
            blastn -query {self.infer_hap_file} -out {self.blast_file} -subject {self.truth_hap_file} -outfmt 7
            """
            # if not os.path.isfile(f"{self.truth_hap_file}.nhr"):
            #     os.system(f"makeblastdb -in {self.truth_hap_file} -dbtype nucl -out {self.truth_hap_file}")
            # command = f"""
            # blastn -query {self.infer_hap_file} -out {self.blast_file} -db {self.truth_hap_file} -outfmt 7 -num_threads 10
            # """
        elif flag == "somewhat":
            # Somewhat similar sequences (blastn) in https://blast.ncbi.nlm.nih.gov/Blast.cgi
            command = f"""
            blastn -query {self.infer_hap_file} -out {self.blast_file} -subject {self.truth_hap_file} -outfmt 7 -evalue 0.05 \
                -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 
            """
        # print (flag)
        os.system(command)

    def record_blast(self, index):
        # -penalty -1 -reward 1 -gapopen 4 -gapextend 1 -strand plus
        command = f"""
        blastn -query {self.infer_hap_file} -out {self.blast_file}.{index}.fmt7.record -subject {self.truth_hap_file} -outfmt 7 
        blastn -query {self.infer_hap_file} -out {self.blast_file}.{index}.fmt3.record -subject {self.truth_hap_file} -outfmt 3  
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

    def read_blast(self, blast_file):
        f = open(blast_file, 'r')
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
            if i > 20:
                break

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

    def main(self):
        self.get_fasta_len()
        self.blast_map("strict")
        self.read_blast(self.blast_file)
        self.get_gap_per()
        if self.mapped_len == 0:
            self.blast_map("somewhat")
            self.read_blast(self.blast_file)
            self.get_gap_per()            

        # print (self.mapped_interval, self.infer_hap_len, self.truth_hap_len, self.mapped_len, align.base_error)
        align = Align(self.mapped_len, self.infer_hap_len, self.truth_hap_len, \
            self.mismatch_num, self.gap_open_num,self.true_mapped_len)
        return align

class Seq_error_accelerate():

    def __init__(self, infer_hap_file, truth_hap_file, tag, method):
        self.infer_hap_file = infer_hap_file
        self.truth_hap_file = truth_hap_file
        self.blast_file = f"{self.infer_hap_file}.{tag}.blast"
        self.infer_hap_len_dict = {}
        self.infer_hap_N_ratio_dict = {}
        self.truth_hap_len_dict = {}
        self.mapped_interval_dict = {} # for inferred hap
        self.true_mapped_interval_dict = {} # for true hap
        self.mapped_len_dict = {}
        self.true_mapped_len_dict = {}
        self.mismatch_num_dict = {}
        self.gap_open_num_dict = {}
        self.method = method

    def blast_map(self, flag):
        #-penalty -1 -reward 1 -gapopen 4 -gapextend 1 -strand plus
        if flag == "strict":
            # command = f"""
            # blastn -query {self.infer_hap_file} -out {self.blast_file} -subject {self.truth_hap_file} -outfmt 7
            # """
            if not os.path.isfile(f"{self.truth_hap_file}.nhr"):
                os.system(f"makeblastdb -in {self.truth_hap_file} -dbtype nucl -out {self.truth_hap_file}")
            command = f"""
            blastn -query {self.infer_hap_file} -out {self.blast_file} -db {self.truth_hap_file} -outfmt 7 -num_threads 10
            """
        elif flag == "somewhat":
            # Somewhat similar sequences (blastn) in https://blast.ncbi.nlm.nih.gov/Blast.cgi
            command = f"""
            blastn -query {self.infer_hap_file} -out {self.blast_file} -subject {self.truth_hap_file} -outfmt 7 -evalue 0.05 \
                -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 
            """
        # print (flag)
        os.system(command)

    def get_fasta_len(self):
        with open(self.infer_hap_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                nCount = str(record.seq).lower().count('n')
                total_length = len(record.seq)
                if self.method == "spechla" or self.method == "hisat": 
                    self.infer_hap_len_dict[record.id] = total_length - nCount
                else:
                    self.infer_hap_len_dict[record.id] = total_length
                self.infer_hap_N_ratio_dict[record.id] = round(nCount/total_length, 2)
                # print (record.id, total_length, nCount, self.infer_hap_len_dict[record.id], round(nCount/total_length, 2))
        # with open(self.truth_hap_file) as handle: # too large
        #     for record in SeqIO.parse(handle, "fasta"):
        #         self.truth_hap_len_dict[record.id] = len(record.seq)

    def read_blast(self, blast_file):
        f = open(blast_file, 'r')
        i = 0
        line_count = {}
        for line in f:
            if line[0] == "#":
                continue
            array = line.strip().split()
            seq_name = array[0]
            map_s = int(array[6])
            map_e = int(array[7])
            mismatches = int(array[4])
            gap_opens = int(array[5])
            # print (map_s, map_e)
            if seq_name not in self.mapped_interval_dict:
                self.mapped_interval_dict[seq_name] = []
                self.true_mapped_interval_dict[seq_name] = []
                self.mismatch_num_dict[seq_name] = 0
                self.gap_open_num_dict[seq_name] = 0
                line_count[seq_name] = 0
            line_count[seq_name] += 1
            if line_count[seq_name] > 50:
                continue

            self.mapped_interval_dict[seq_name] = self.get_uniq_map(map_s, map_e, mismatches, gap_opens, self.mapped_interval_dict[seq_name],seq_name)
            true_map_s = int(array[8])
            true_map_e = int(array[9])
            self.true_mapped_interval_dict[seq_name] = self.get_uniq_map(true_map_s, true_map_e, 0, 0, self.true_mapped_interval_dict[seq_name],seq_name)
            i += 1
        # print (self.mapped_interval_dict, self.infer_hap_len_dict)

    def get_uniq_map(self, map_s, map_e, mismatches, gap_opens, mapped_interval,seq_name):
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
            if map_s <= interval[0] and map_e >= interval[1]:
                interval[0] = map_s
                interval[1] = map_e
                flag = False
            if map_s >= interval[0] and map_e <= interval[1]:
                flag = False
        if flag and map_e - map_s > 5:
            uniq_map_len = map_e - map_s + 1
            mismatches = mismatches * (uniq_map_len/origin_map_len)
            gap_opens = gap_opens * (uniq_map_len/origin_map_len)
            self.mismatch_num_dict[seq_name] += mismatches
            self.gap_open_num_dict[seq_name] += gap_opens
            mapped_interval.append([map_s, map_e])
        # print (seq_name, mapped_interval)
        return mapped_interval

    def get_gap_per(self):
        for geno_name in self.mapped_interval_dict.keys():
            self.mapped_len_dict[geno_name] = 0
            self.true_mapped_len_dict[geno_name] = 0
            for interval in self.mapped_interval_dict[geno_name]:
                self.mapped_len_dict[geno_name] += (interval[1] - interval[0] + 1)
            for interval in self.true_mapped_interval_dict[geno_name]:
                self.true_mapped_len_dict[geno_name] += (interval[1] - interval[0] + 1)

    def main(self):
        result_dict = {}
        self.get_fasta_len()
        self.blast_map("strict")
        print ("blast is done")
        self.read_blast(self.blast_file)
        self.get_gap_per()
        # if self.mapped_len == 0:
        #     self.blast_map("somewhat")
        #     self.read_blast(self.blast_file)
        #     self.get_gap_per()            
        # print (self.mapped_interval, self.infer_hap_len, self.truth_hap_len, self.mapped_len, align.base_error)
        for geno_name in self.mapped_interval_dict.keys():
            align = Align(self.mapped_len_dict[geno_name], self.infer_hap_len_dict[geno_name], 1000000, \
                self.mismatch_num_dict[geno_name], self.gap_open_num_dict[geno_name],self.true_mapped_len_dict[geno_name])
            # print (geno_name, align.base_error, align.short_gap_error, self.mapped_len_dict[geno_name])
            # print (geno_name, self.mapped_len_dict[geno_name], self.infer_hap_len_dict[geno_name], self.mapped_interval_dict[geno_name])
            result_dict[geno_name] = [align.base_error, align.short_gap_error, align.gap_precision, self.mapped_len_dict[geno_name], self.infer_hap_N_ratio_dict[geno_name]]
        return result_dict

class Seq_error_accelerate_sim(Seq_error_accelerate):
    def read_blast(self, blast_file):
        focus_alleles = ["A_01_01","B_01_01","C_01_01","DPA1_01_01","DPB1_01_01","DQA1_01_01","DQB1_01_01","DRB1_01_01" ]
        f = open(blast_file, 'r')
        i = 0
        line_count = {}
        for line in f:
            if line[0] == "#":
                continue
            array = line.strip().split()
            seq_name = array[0]
            map_s = int(array[6])
            map_e = int(array[7])
            mismatches = int(array[4])
            gap_opens = int(array[5])
            # print (map_s, map_e)
            if seq_name not in self.mapped_interval_dict:
                mapped_name = array[1]
                self.mapped_interval_dict[seq_name] = []
                self.true_mapped_interval_dict[seq_name] = []
                self.mismatch_num_dict[seq_name] = 0
                self.gap_open_num_dict[seq_name] = 0
                line_count[seq_name] = 0
            else:
                if array[1] != mapped_name:
                    continue
            if array[1] not in focus_alleles:
                continue

            line_count[seq_name] += 1
            if line_count[seq_name] > 50:
                continue

            self.mapped_interval_dict[seq_name] = self.get_uniq_map(map_s, map_e, mismatches, gap_opens, self.mapped_interval_dict[seq_name],seq_name)
            true_map_s = int(array[8])
            true_map_e = int(array[9])
            self.true_mapped_interval_dict[seq_name] = self.get_uniq_map(true_map_s, true_map_e, 0, 0, self.true_mapped_interval_dict[seq_name],seq_name)
            i += 1

class Seq_error_muscle():
    def __init__(self, infer_hap_file, truth_hap_file, gene):
        self.infer_hap_file = infer_hap_file
        self.truth_hap_file = truth_hap_file
        self.muscle_file = f"{self.infer_hap_file}.muscle"
        self.merge_file = f"{self.infer_hap_file}.merge.fasta"

        self.infer_hap_len = None
        self.truth_hap_len = None
        self.mapped_interval = [] # for inferred hap
        self.true_mapped_interval = [] # for true hap
        self.mapped_len = 0
        self.true_mapped_len = 0
        self.mismatch_num = 0
        self.gap_open_num = 0
        self.gene = gene

        self.indel_cutoff = 150 #bp
    
    def run_muscle(self):
        command = f"""
        cat {self.infer_hap_file} {self.truth_hap_file} >{self.merge_file}
        muscle -in {self.merge_file} -out {self.muscle_file} 
        echo {self.muscle_file} 
        """
        os.system(command)
    
    def read_muscle(self):
        record_seq = ['', '']
        with open(self.muscle_file) as handle:
            index = 0
            for line in handle:
                line = line.strip()
                if line[0] == ">":
                    index += 1
                else:
                    record_seq[index-1] += line
        infer_small_variant_len, infer_long_indel_len = self.count_seq(record_seq[0])
        true_small_variant_len, true_long_indel_len = self.count_seq(record_seq[1])
        seq_len = len(record_seq[0])
        base_error = true_small_variant_len/(seq_len-true_long_indel_len)
        gap_recall = (seq_len-true_long_indel_len)/seq_len
        gap_precision = (seq_len-infer_long_indel_len)/seq_len
        print (base_error, gap_recall, gap_precision)
    
    def count_seq(self, seq):
        unmapped_len = 0
        total_unmapped_len = 0
        small_variant_len = 0
        long_indel_len = 0

        seq = str(seq)
        seq_len = len(seq)
        mis_len = 0
        for i in range(seq_len):
            if seq[i] == "-":
                unmapped_len += 1
                mis_len += 1
            else:
                total_unmapped_len += mis_len
                if mis_len < self.indel_cutoff:
                    small_variant_len += mis_len
                else:
                    long_indel_len += mis_len
                mis_len = 0

        total_unmapped_len += mis_len
      
            # print (i, seq[i])
            # break
        
        # small_variant_rate = small_variant_len/(seq_len - long_indel_len)
        # print (unmapped_len, seq_len, total_unmapped_len, small_variant_len, long_indel_len, small_variant_rate)
        return small_variant_len, long_indel_len

    
    def main(self):
        self.run_muscle()
        self.read_muscle()

def eva_HG002_spechla():
    # outdir = "/mnt/d/HLAPro_backup/trio/HG002/"
    # outdir = "/mnt/d/HLAPro_backup/trio/trio_1000/output/HG002_single_pac/"
    # outdir = "/mnt/d/HLAPro_backup/trio/trio_1000/spechla/HG002/"
    outdir = "/mnt/d/HLAPro_backup/trio/HG002_exon/"
    truth_file1 = "/mnt/d/HLAPro_backup/trio/truth_MHC/H1-asm.fa"
    truth_file2 = "/mnt/d/HLAPro_backup/trio/truth_MHC/H2-asm.fa"
    data = []
    for gene in gene_list:
        infer_file1 = outdir + "hla.allele.1.HLA_%s.fasta"%(gene)
        infer_file1 = merge_exon(infer_file1, gene, 0)
        infer_file2 = outdir + "hla.allele.2.HLA_%s.fasta"%(gene)
        infer_file2 = merge_exon(infer_file2, gene, 1)
        seq = Seq_error(infer_file1, truth_file1, gene)
        align_11 = seq.main()
        seq = Seq_error(infer_file2, truth_file2, gene)
        align_22 = seq.main()
        seq = Seq_error(infer_file1, truth_file2, gene)
        align_12 = seq.main()
        seq = Seq_error(infer_file2, truth_file1, gene)
        align_21 = seq.main()
        if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
            seq = Seq_error(infer_file1, truth_file1, gene)
            # align_11 = seq.main()
            seq.record_blast(1)
            seq = Seq_error(infer_file2, truth_file2, gene)
            # align_22 = seq.main()
            seq.record_blast(2)
            choose_align1 = align_11
            choose_align2 = align_22
        else:
            seq = Seq_error(infer_file1, truth_file2, gene)
            # align_12 = seq.main()
            seq.record_blast(1)
            seq = Seq_error(infer_file2, truth_file1, gene)
            # align_21 = seq.main()
            seq.record_blast(2)
            choose_align1 = align_12
            choose_align2 = align_21
        # print ("#", align_11.base_error, align_22.base_error, align_12.base_error, align_21.base_error)
        base_error = (choose_align1.base_error + choose_align2.base_error)/2
        short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
        gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
        gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
        mapped_len = (choose_align1.mapped_len + choose_align2.mapped_len)/2
        data.append([base_error, short_gap_error, gap_recall, gap_precision, gene])
        # print (gene, base_error, short_gap_error, gap_recall, gap_precision)
        print (gene, round(base_error,6), round(short_gap_error,6), round(gap_recall,6), round(gap_precision,6),\
             round(choose_align1.base_error,6), round(choose_align2.base_error,6), mapped_len)
        # break
    df = pd.DataFrame(data, columns = ["base_error", "short_gap_error", "gap_recall", "gap_precision", "Gene"])
    df.to_csv('/mnt/d/HLAPro_backup/trio/hg002_high_dp_haplo_assess.csv', sep=',')

def eva_pedigree_spechla():
    outdir = "/mnt/d/HLAPro_backup/trio/trio_1000/spechla/"

    data = []
    # pedigree_samples_list = [["NA12878", "NA12891", "NA12892"], ["NA19240", "NA19238", "NA19239"],["HG002","HG003","HG004"],["HG005","HG006","HG007"]]
    pedigree_samples_list = [["NA12878", "NA12891", "NA12892"], ["NA19240", "NA19238", "NA19239"]]
    mismatch_list, gap_list = [], []
    # pedigree_samples_list = [["NA12878", "NA12891", "NA12892"]]
    for pedigree_samples in pedigree_samples_list:
        for gene in gene_list:
            for j in range(2):
                infer_file1 = outdir + pedigree_samples[0] + "/hla.allele.1.HLA_%s.fasta"%(gene)
                infer_file2 = outdir + pedigree_samples[0] + "/hla.allele.2.HLA_%s.fasta"%(gene)

                parent_file1 = outdir + pedigree_samples[1+j] + "/hla.allele.1.HLA_%s.fasta"%(gene)
                parent_file2 = outdir + pedigree_samples[1+j] + "/hla.allele.2.HLA_%s.fasta"%(gene)

                record_base_error = []
                seq = Seq_error(infer_file1, parent_file1, gene)
                align = seq.main()
                record_base_error.append(align.base_error)
                choose_align = align
                choose_seq = seq
                

                seq = Seq_error(infer_file2, parent_file2, gene)
                align = seq.main()
                record_base_error.append(align.base_error)
                if align.base_error < choose_align.base_error:
                    choose_align = align
                    choose_seq = seq
                    

                seq = Seq_error(infer_file1, parent_file2, gene)
                align = seq.main()
                record_base_error.append(align.base_error)
                if align.base_error < choose_align.base_error:
                    choose_align = align
                    choose_seq = seq

                seq = Seq_error(infer_file2, parent_file1, gene)
                align = seq.main()
                record_base_error.append(align.base_error)
                if align.base_error < choose_align.base_error:
                    choose_align = align
                    choose_seq = seq
                # print (record_base_error)

                choose_seq.record_blast(j)
                base_error = choose_align.base_error 
                short_gap_error = choose_align.short_gap_error 
                gap_recall = choose_align.gap_recall 
                gap_precision = choose_align.gap_precision
                data.append([pedigree_samples[1+j], base_error, gene, "Mismatch_Rate"])
                data.append([pedigree_samples[1+j], short_gap_error,gene, "Gap_Rate"])
                # print (gene, base_error, short_gap_error, gap_recall, gap_precision)
                print (pedigree_samples[1+j], gene, round(base_error,6), round(short_gap_error,6), round(gap_recall,6), round(gap_precision,6),\
                    choose_align.mismatch_num)
                mismatch_list.append(base_error)
                gap_list.append(short_gap_error)
        # break
    df = pd.DataFrame(data, columns = ["parent","value", "Gene", "group"])
    df.to_csv('/mnt/d/HLAPro_backup/trio/pedigree_haplo_assess.csv', sep=',')
    print ("Mean base error:", np.mean(mismatch_list), "mean gap rate is", np.mean(gap_list))

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

        seq = Seq_error(infer_file1, truth_file1, gene)
        align_11 = seq.main()
        seq = Seq_error(infer_file2, truth_file2, gene)
        align_22 = seq.main()
        seq = Seq_error(infer_file1, truth_file2, gene)
        align_12 = seq.main()
        seq = Seq_error(infer_file2, truth_file1, gene)
        align_21 = seq.main()
        if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
            seq = Seq_error(infer_file1, truth_file1, gene)
            choose_align1 = seq.main()
            seq = Seq_error(infer_file2, truth_file2, gene)
            choose_align2 = seq.main()
            # choose_align1 = align_11
            # choose_align2 = align_22
        else:
            seq = Seq_error(infer_file1, truth_file2, gene)
            choose_align1 = seq.main()
            seq = Seq_error(infer_file2, truth_file1, gene)
            choose_align2 = seq.main()
            # choose_align1 = align_12
            # choose_align2 = align_21
        # print ("#", align_11.base_error, align_22.base_error, align_12.base_error, align_21.base_error)
        base_error = (choose_align1.base_error + choose_align2.base_error)/2
        short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
        gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
        gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
        mapped_len = (choose_align1.mapped_len + choose_align2.mapped_len)/2
        data.append([base_error, short_gap_error, gap_recall, gap_precision, gene])
        print (gene, round(base_error,6), round(short_gap_error,6), round(gap_recall,6), round(gap_precision,6),\
             round(choose_align1.base_error,6), round(choose_align2.base_error,6),round(mapped_len))
        # break
    df = pd.DataFrame(data, columns = ["base_error", "short_gap_error", "gap_recall", "gap_precision", "Gene"])
    df.to_csv('/mnt/d/HLAPro_backup/trio/hisat_hg002_haplo_assess.csv', sep=',')

def eva_HG002_kourami():
    outdir = "/mnt/d/HLAPro_backup/trio/Trio/Kourami/HG002/"
    truth_file1 = "/mnt/d/HLAPro_backup/trio/truth_MHC/H1-asm.fa"
    truth_file2 = "/mnt/d/HLAPro_backup/trio/truth_MHC/H2-asm.fa"
 

    data = []
    for gene in gene_list:
        infer_file1 = outdir + "kourami.hla.allele.1.HLA_%s.fasta"%(gene)
        infer_file2 = outdir + "kourami.hla.allele.2.HLA_%s.fasta"%(gene)

        alleles = ["%s_0:0"%(gene), "%s_1:1"%(gene)]
        if gene == "DQA1":
            alleles = ["%s_0"%(gene), "%s_0"%(gene)]
        fasta_file = outdir + "/HG002.hla_%s.typed.fa"%(gene)
        if not os.path.isfile(fasta_file):
            continue

        with open(fasta_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == alleles[0]:
                    with open(infer_file1, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                    break

        with open(fasta_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == alleles[1]:
                    with open(infer_file2, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                    break
        infer_file1 = splitN_for_kourami(infer_file1)
        infer_file2 = splitN_for_kourami(infer_file2)
        seq = Seq_error(infer_file1, truth_file1, gene)
        align_11 = seq.main()
        seq = Seq_error(infer_file2, truth_file2, gene)
        align_22 = seq.main()
        seq = Seq_error(infer_file1, truth_file2, gene)
        align_12 = seq.main()
        seq = Seq_error(infer_file2, truth_file1, gene)
        align_21 = seq.main()
        if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
            seq = Seq_error(infer_file1, truth_file1, gene)
            choose_align1 = seq.main()
            seq = Seq_error(infer_file2, truth_file2, gene)
            choose_align2 = seq.main()
            # choose_align1 = align_11
            # choose_align2 = align_22
        else:
            seq = Seq_error(infer_file1, truth_file2, gene)
            choose_align1 = seq.main()
            seq = Seq_error(infer_file2, truth_file1, gene)
            choose_align2 = seq.main()
            # choose_align1 = align_12
            # choose_align2 = align_21
        # print ("#", align_11.base_error, align_22.base_error, align_12.base_error, align_21.base_error)
        base_error = (choose_align1.base_error + choose_align2.base_error)/2
        short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
        gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
        gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
        mapped_len = (choose_align1.mapped_len + choose_align2.mapped_len)/2
        data.append([base_error, short_gap_error, gap_recall, gap_precision, gene])
        print (gene, round(base_error,6), round(short_gap_error,6), round(gap_recall,6), round(gap_precision,6),\
             round(choose_align1.base_error,6), round(choose_align2.base_error,6),mapped_len)
        # break
    df = pd.DataFrame(data, columns = ["base_error", "short_gap_error", "gap_recall", "gap_precision", "Gene"])
    df.to_csv('/mnt/d/HLAPro_backup/trio/kourami_hg002_haplo_assess.csv', sep=',')   

def get_sim_true_allele(record_true_file):
    true_allele = {}
    f = open(record_true_file)
    for line in f:
        array = line.strip().split()
        if array[0] != sample_name:
            continue
        if array[1] in gene_list:
            true_allele[array[1]] = [array[3], array[4]]
    f.close()
    return true_allele

def split_hisat_fasta(outdir, sample_name):
    # hisat_fasta = "%s/assembly_graph-hla.%s_read1_fastq_gz-hla-extracted-1_fq.fasta"%(outdir, sample_name)   
    hisat_fasta = "%s/assembly_graph-hla.simu_read1_fastq_gz-hla-extracted-1_fq.fasta"%(outdir)  
    record_dict = {}
    with open(hisat_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            gene_name = record.id.split("_")[0].split(")")[1]
            if gene_name not in record_dict:
                record_dict[gene_name] = 0
            record_dict[gene_name] += 1
            if record_dict[gene_name] > 2:
                continue
            infer_file = outdir + "/hisat.hla.allele.%s.HLA_%s.fasta"%(record_dict[gene_name], gene_name)
            with open(infer_file, "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")

def split_kourami_fasta(outdir, sample_name):
    record_dict = {}
    files = os.listdir(outdir)
    for file in files:
        if file[-3:] != ".fa":
            continue
        with open(outdir + "/" + file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                gene_name = record.id.split("_")[0]
                if gene_name not in record_dict:
                    record_dict[gene_name] = 0
                record_dict[gene_name] += 1
                infer_file = outdir + "/kourami.hla.allele.%s.HLA_%s.fasta"%(record_dict[gene_name], gene_name)
                with open(infer_file, "w") as output_handle:
                    SeqIO.write(record, output_handle, "fasta")
                new_infer_file = splitN_for_kourami(infer_file)
                os.system(f"cp {new_infer_file} {infer_file}")

def eva_simu(database, record_true_file, outdir, sample_name):
    outdir = outdir + "/"
    true_allele = get_sim_true_allele(record_true_file)
  
    data = []
    for gene in gene_list:
        if method == "spechla":
            infer_file1 = outdir + "hla.allele.1.HLA_%s.fasta"%(gene)
            infer_file2 = outdir + "hla.allele.2.HLA_%s.fasta"%(gene)
            if exon_flag == "exon":
                infer_file1 = merge_exon(infer_file1, gene, 0)
                infer_file2 = merge_exon(infer_file2, gene, 1)
        elif method == "hisat":
            infer_file1 = outdir + "hisat.hla.allele.1.HLA_%s.fasta"%(gene)
            infer_file2 = outdir + "hisat.hla.allele.2.HLA_%s.fasta"%(gene)
        elif method == "kourami":
            infer_file1 = outdir + "kourami.hla.allele.1.HLA_%s.fasta"%(gene)
            infer_file2 = outdir + "kourami.hla.allele.2.HLA_%s.fasta"%(gene)
        else:
            print ("pls provide method name, spechla, hisat, or kourami.")
        if not os.path.isfile(infer_file1):
            print ("there is no %s."%(infer_file1))
            continue

        truth_file1 = outdir + ".true.%s.1.fasta"%(gene)
        truth_file2 = outdir + ".true.%s.2.fasta"%(gene)

        alleles = true_allele[gene]
        command = f"samtools faidx {database} {alleles[0]} >{truth_file1}"
        os.system(command)
        command = f"samtools faidx {database} {alleles[1]} >{truth_file2}"
        os.system(command)


        seq = Seq_error(infer_file1, truth_file1, gene)
        align_11 = seq.main()
        seq = Seq_error(infer_file2, truth_file2, gene)
        align_22 = seq.main()
        seq = Seq_error(infer_file1, truth_file2, gene)
        align_12 = seq.main()
        seq = Seq_error(infer_file2, truth_file1, gene)
        align_21 = seq.main()



        if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
            choose_align1 = align_11
            choose_align2 = align_22
            seq = Seq_error(infer_file1, truth_file1, gene)
            seq.record_blast(0)
            seq = Seq_error(infer_file2, truth_file2, gene)
            seq.record_blast(1)

            # seq = Seq_error_muscle(infer_file1, truth_file1, gene)
            # seq.main()
            # seq = Seq_error_muscle(infer_file2, truth_file2, gene)
            # seq.main()

        else:
            choose_align1 = align_12
            choose_align2 = align_21
            seq = Seq_error(infer_file1, truth_file2, gene)
            seq.record_blast(0)
            seq = Seq_error(infer_file2, truth_file1, gene)
            seq.record_blast(1)

            # seq = Seq_error_muscle(infer_file1, truth_file2, gene)
            # seq.main()
            # seq = Seq_error_muscle(infer_file2, truth_file1, gene)
            # seq.main()

        # print (align_11.base_error, align_22.base_error, align_12.base_error, align_21.base_error)
        base_error = (choose_align1.base_error + choose_align2.base_error)/2
        short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
        gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
        gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
        data.append([base_error, short_gap_error, gap_recall, gap_precision, gene])
        print (gene, base_error, short_gap_error, gap_recall, gap_precision, choose_align1.mismatch_num, choose_align2.mismatch_num)
        # break
    df = pd.DataFrame(data, columns = ["base_error", "short_gap_error", "gap_recall", "gap_precision", "Gene"])
    df.to_csv('%s/haplotype_assessment.csv'%(outdir), sep=',')

def eva_hgsvc2_spechla_accelerate():
    # outdir = "/mnt/d/HLAPro_backup/trio/HG002/"
    sample = "HG03009"
    outdir = "/mnt/d/HLAPro_backup/haplotype/hgsvc2/HG03009/"
    truth_file1 = "/mnt/d/HLAPro_backup/haplotype/my_HLA/v12_HG03009_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta"
    truth_file2 = "/mnt/d/HLAPro_backup/haplotype/my_HLA/v12_HG03009_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta"
    infer_file = outdir + "/hla.allele.all.fasta"
    os.system("cat %s/hla.allele.*.HLA_*.fasta>%s"%(outdir, infer_file))
    data = []

    seq = Seq_error_accelerate(infer_file, truth_file1, "hap1")
    result_dict_1 = seq.main()
    seq = Seq_error_accelerate(infer_file, truth_file2, "hap2")
    result_dict_2 = seq.main()
    for geno_name in result_dict_1.keys():
        if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
            result_dict_1[geno_name] = result_dict_2[geno_name]
    gene_dict = {}
    for geno_name in result_dict_1.keys():
        gene = geno_name[:-2]
        if gene not in gene_dict:
            gene_dict[gene] = []
        gene_dict[gene].append(result_dict_1[geno_name])
    for gene in gene_dict:
        base_error = (gene_dict[gene][0][0] + gene_dict[gene][1][0])/2
        short_gap_error = (gene_dict[gene][0][1] + gene_dict[gene][1][1])/2
        gap_precision = (gene_dict[gene][0][2] + gene_dict[gene][1][2])/2
        map_len = (gene_dict[gene][0][3] + gene_dict[gene][1][3])/2
        print (sample, gene, base_error, short_gap_error, map_len, gap_precision)
        data.append([sample, gene, base_error, short_gap_error, map_len, gap_precision])


    df = pd.DataFrame(data, columns = ["sample", "gene", "base_error", "short_gap_error", "map_len", "gap_precision"])
    df.to_csv('/mnt/d/HLAPro_backup/haplotype/hgsvc_haplo_assess.csv', sep=',')

def eva_hgsvc2_spechla():
    # outdir = "/mnt/d/HLAPro_backup/trio/HG002/"
    sample = "HG03009"
    outdir = "/mnt/d/HLAPro_backup/haplotype/hgsvc2/HG03009/"
    truth_file1 = "/mnt/d/HLAPro_backup/haplotype/my_HLA/v12_HG03009_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta"
    truth_file2 = "/mnt/d/HLAPro_backup/haplotype/my_HLA/v12_HG03009_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta"
    data = []
    for gene in gene_list:
        infer_file1 = outdir + "hla.allele.1.HLA_%s.fasta"%(gene)
        infer_file2 = outdir + "hla.allele.2.HLA_%s.fasta"%(gene)
        seq = Seq_error(infer_file1, truth_file1, gene)
        align_11 = seq.main()
        seq = Seq_error(infer_file2, truth_file2, gene)
        align_22 = seq.main()
        seq = Seq_error(infer_file1, truth_file2, gene)
        align_12 = seq.main()
        seq = Seq_error(infer_file2, truth_file1, gene)
        align_21 = seq.main()
        if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
            choose_align1 = align_11
            choose_align2 = align_22
        else:
            choose_align1 = align_12
            choose_align2 = align_21
        base_error = (choose_align1.base_error + choose_align2.base_error)/2
        short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
        gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
        gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
        data.append([sample, gene, base_error, short_gap_error, gap_recall, gap_precision])
        print (gene, round(base_error,6), round(short_gap_error,6), round(gap_recall,6), round(gap_precision,6),\
             round(choose_align1.base_error,6), round(choose_align2.base_error,6))
    df = pd.DataFrame(data, columns = ["sample", "gene", "base_error", "short_gap_error", "gap_recall", "gap_precision"])
    df.to_csv('/mnt/d/HLAPro_backup/haplotype/hgsvc_haplo_assess.csv', sep=',')

def each_simulated_sample(gene, outdir, sample, truth_dir):
    truth_file1 = "%s/%s.HLA_%s_1.fasta"%(truth_dir, sample, gene)
    truth_file2 = "%s/%s.HLA_%s_2.fasta"%(truth_dir, sample, gene)
    if not os.path.isfile(truth_file1):
        truth_file1 = "%s/%s.%s_1.fasta"%(truth_dir, sample, gene)
        truth_file2 = "%s/%s.%s_2.fasta"%(truth_dir, sample, gene)
    infer_file1 = outdir + "hla.allele.1.HLA_%s.fasta"%(gene)
    infer_file2 = outdir + "hla.allele.2.HLA_%s.fasta"%(gene)
    seq = Seq_error(infer_file1, truth_file1, gene)
    align_11 = seq.main()
    seq = Seq_error(infer_file2, truth_file2, gene)
    align_22 = seq.main()
    seq = Seq_error(infer_file1, truth_file2, gene)
    align_12 = seq.main()
    seq = Seq_error(infer_file2, truth_file1, gene)
    align_21 = seq.main()
    if align_11.base_error + align_22.base_error <= align_12.base_error + align_21.base_error:
        choose_align1 = align_11
        choose_align2 = align_22
        seq = Seq_error(infer_file1, truth_file1, gene)
        seq.record_blast(1)
        seq = Seq_error(infer_file2, truth_file2, gene)
        seq.record_blast(2)

    else:
        choose_align1 = align_12
        choose_align2 = align_21
        seq = Seq_error(infer_file1, truth_file2, gene)
        seq.record_blast(1)
        seq = Seq_error(infer_file2, truth_file1, gene)
        seq.record_blast(2)
    # print (align_11.base_error, align_22.base_error,align_12.base_error,align_21.base_error)
    base_error = (choose_align1.base_error + choose_align2.base_error)/2
    short_gap_error = (choose_align1.short_gap_error + choose_align2.short_gap_error)/2
    gap_recall = (choose_align1.gap_recall + choose_align2.gap_recall)/2
    gap_precision = (choose_align1.gap_precision + choose_align2.gap_precision)/2
    print (gene, round(base_error,6), round(short_gap_error,6), round(gap_recall,6), round(gap_precision,6),\
            round(choose_align1.base_error,6), round(choose_align2.base_error,6))
    return base_error, short_gap_error, gap_recall, gap_precision

def eva_data_types_spechla():
    # dict = {"PE":"novel", "+PacBio":"pacbio_illumina","+Hi-C":"hic_illumina","+ONT":"ont_illumina","PacBio":"pac_alone","+10X":"10x_illumina" }
    dict = {"PE":"novel", "+PacBio":"pacbio_illumina","+Hi-C":"hic_illumina","+ONT":"ont_illumina","PacBio":"pac_alone","ONT":"ont_alone","+10X":"10x_illumina"}
    # dict = {"+10X":"10x_illumina" }
    data = []
    for i in range(50):
        sample = "novel_%s"%(i)
        for plat in dict.keys():
            outdir = "/mnt/d/HLAPro_backup/pacbio/" + dict[plat]+ "/" + sample + "/"
            truth_dir = "/mnt/d/HLAPro_backup/pacbio/simulation/truth/"
            
            for gene in gene_list:
                base_error, short_gap_error, gap_recall, gap_precision = each_simulated_sample(gene, outdir, sample, truth_dir)
                data.append([sample, gene, base_error, short_gap_error, plat])
                # print (sample, gene, base_error, short_gap_error, plat)
    df = pd.DataFrame(data, columns = ["sample", "gene", "mismatch_rate", "gap_rate", "Data"])
    df.to_csv('/mnt/d/HLAPro_backup/pacbio/novel/haplo_assess.csv', sep=',')

def eva_allele_imblance():
    # sample = "imbalance_48_52_0"
    # sample = "imbalance_40_60_0"
    # outdir = "/mnt/d/HLAPro_backup/imbalance/output/" + sample + "/"
    truth_dir = "/mnt/d/HLAPro_backup/imbalance/data/truth/"
    data = []
    for i in range(20):
        for group in ["50_50", "48_52", "40_60", "30_70", "20_80"]:
        # for group in ["50_50", "40_60", "30_70", "20_80"]:
            sample = "imbalance_%s_%s"%(group, i)
            outdir = "/mnt/d/HLAPro_backup/imbalance/output/" + sample + "/"
            for weight in [0, 0.25, 0.5, 0.75, 1]:
                # print ("bash /home/wangshuai/softwares/SpecHLA/script/whole/test_SpecHLA.sh -n %s -1 x -2 x -w %s -o /mnt/d/HLAPro_backup/imbalance/output"%(sample,weight))
                os.system("bash /home/wangshuai/softwares/SpecHLA/script/whole/test_SpecHLA.sh -n %s -1 x -2 x -w %s -o /mnt/d/HLAPro_backup/imbalance/output"%(sample,weight))
                weight_base_error = []
                for gene in gene_list:
                    base_error, short_gap_error, gap_recall, gap_precision = each_simulated_sample(gene, outdir, sample, truth_dir)
                    data.append([sample, gene, base_error, short_gap_error, weight, group])
                    weight_base_error.append(base_error)
                print (sample, "#####", weight, np.mean(weight_base_error))

    df = pd.DataFrame(data, columns = ["sample", "gene", "base_error", "short_gap_error", "weight", "group"])
    df.to_csv('/mnt/d/HLAPro_backup/imbalance/haplo_assess.csv', sep=',')

class Assess_hgsvc2():

    def __init__(self):
        self.record_truth_file_dict = {}
        self.work_dir = "/mnt/d/HLAPro_backup/haplotype_v2"
        self.spechla_dir = self.work_dir + "/spechla/"
        self.hisat_dir = self.work_dir + "/shell/Hisat/"
        self.kourami_dir = self.work_dir + "/shell/Kourami/"
        self.get_phased_assemblies()
        self.hisat_gene_count = {}
        self.spechla_gene_count = {}
        self.kourami_gene_count = {}
        self.kourami_exon_accuracy = {}
        self.spechla_exon_accuracy = {}
        self.record_hisat_output = {} # record whether hisat reconstruct the seq
        self.N_ratio_cutoff = Min_N_ratio
       
    def get_phased_assemblies(self):
        record_truth_file_dict = {}
        inpath = "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/"
        for file in os.listdir(inpath):
            if file[-5:] != "fasta":
                continue
            sample = file.split("_")[1]
            if sample not in record_truth_file_dict:
                record_truth_file_dict[sample] = ['', '']
            full_file = inpath + file
            if re.search(".h1-", full_file):
                record_truth_file_dict[sample][0] = full_file
            else:
                record_truth_file_dict[sample][1] = full_file
        # print (record_truth_file_dict)
        self.record_truth_file_dict =  record_truth_file_dict

    def main(self):
        data = []
        record_used_samples = []
        sample_num = 0
        print ("number of samples with phased assembly", len(self.record_truth_file_dict))
        for sample in self.record_truth_file_dict.keys():
            if not os.path.isfile(self.spechla_dir + f"/{sample}/hla.allele.1.HLA_A.fasta"):
                continue
            # if sample != "HG03732":
            #     continue
            record_used_samples.append([sample, self.record_truth_file_dict[sample][0].split("/")[-1], self.record_truth_file_dict[sample][1].split("/")[-1]])
            print (sample)
            sample_num += 1
            self.spechla_dir = self.work_dir + "/spechla/"

            data = self.for_hisat(sample, data)
            data = self.for_spechla(sample, data, "SpecHLA")
            
            # self.spechla_dir = self.work_dir + "/spechla_sv/"
            # data = self.for_spechla(sample, data, "SpecHLA-SV")
        print ("total sample number:", sample_num)
        df = pd.DataFrame(data, columns = ["sample", "gene", "mismatch_rate", "gap_rate", "map_len", "sequence_precision", "Methods"])
        df.to_csv(self.work_dir + '/hgsvc_haplo_assess.csv', sep=',')
        print ("reconstructed gene num of hisat", self.hisat_gene_count)
        print ("gene num of spechla with < %s N"%(self.N_ratio_cutoff), self.spechla_gene_count)
        data = []
        for gene in gene_list:
            gene = "HLA_" + gene
            if gene in self.hisat_gene_count:
                allele_num = self.hisat_gene_count[gene]
            else:
                allele_num = 0
            data.append([allele_num, "HISAT", gene])
            data.append([self.spechla_gene_count[gene], "SpecHLA", gene])
        df = pd.DataFrame(data, columns = ["Allele_num", "Methods", "gene"])
        df.to_csv(self.work_dir + '/hgsvc_haplo_assess_allele_num.csv', sep=',')
        df = pd.DataFrame(record_used_samples, columns = ["Sample", "Haplotype1", "Haplotype2"])
        df.to_csv(self.work_dir + '/hgsvc_used_samples.csv', sep=',')

    def for_spechla(self, sample, data, Method):
        outdir = self.spechla_dir + "/" + sample
        low_depth_dict = self.get_low_depth_bed(sample)
        truth_file1 = self.record_truth_file_dict[sample][0]
        # print (self.record_truth_file_dict[sample], truth_file1)
        truth_file2 = self.record_truth_file_dict[sample][1]
        infer_file = outdir + "/hla.allele.all.fasta"
        os.system("cat %s/hla.allele.*.HLA_*.fasta>%s"%(outdir, infer_file))

        if not os.path.isfile(truth_file1) or not os.path.isfile(truth_file2):
            print ("############",truth_file1)
            return data       

        seq = Seq_error_accelerate(infer_file, truth_file1, "hap1", "spechla")
        result_dict_1 = seq.main()
        seq = Seq_error_accelerate(infer_file, truth_file2, "hap2", "spechla")
        result_dict_2 = seq.main()
        for geno_name in result_dict_1.keys():
            if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
                result_dict_1[geno_name] = result_dict_2[geno_name]
        gene_dict = {}
        for geno_name in result_dict_1.keys():
            gene = geno_name[:-2]
            if gene not in gene_dict:
                gene_dict[gene] = []
            gene_dict[gene].append(result_dict_1[geno_name])
        for gene in gene_dict:
            base_error = (gene_dict[gene][0][0] + gene_dict[gene][1][0])/2
            short_gap_error = (gene_dict[gene][0][1] + gene_dict[gene][1][1])/2
            gap_precision = (gene_dict[gene][0][2] + gene_dict[gene][1][2])/2
            map_len = (gene_dict[gene][0][3] + gene_dict[gene][1][3])/2
            N_ratio = (gene_dict[gene][0][4] + gene_dict[gene][1][4])/2
            # if gene in self.record_hisat_output[sample]:
            # if gene not in low_depth_dict:
            # if True:
            if N_ratio < self.N_ratio_cutoff:
                if gene not in self.spechla_gene_count:
                    self.spechla_gene_count[gene] = 0
                self.spechla_gene_count[gene] += 1
                # print (sample, gene, base_error, short_gap_error, map_len, gap_precision, Method, N_ratio)
                data.append([sample, gene, base_error, short_gap_error, map_len, gap_precision, Method])
            else:
                print (sample, gene, base_error, short_gap_error, map_len, gap_precision, Method, N_ratio)
        return data

    def for_spechla_exon(self, sample, data, Method):
        outdir = self.spechla_dir + "/" + sample + "/"
        truth_file1 = self.record_truth_file_dict[sample][0]
        # print (self.record_truth_file_dict[sample], truth_file1)
        truth_file2 = self.record_truth_file_dict[sample][1]
        infer_file = outdir + "/hla.allele.all.fasta"
        for gene in gene_list:
            infer_file1 = outdir + "hla.allele.1.HLA_%s.fasta"%(gene)
            infer_file1 = merge_exon(infer_file1, gene, 0)
            infer_file2 = outdir + "hla.allele.2.HLA_%s.fasta"%(gene)
            infer_file2 = merge_exon(infer_file2, gene, 1)
        os.system("cat %s/*.merged.exon.fasta>%s"%(outdir, infer_file))

        if not os.path.isfile(truth_file1) or not os.path.isfile(truth_file2):
            print ("############", truth_file1)
            return data       

        seq = Seq_error_accelerate(infer_file, truth_file1, "hap1", "spechla")
        result_dict_1 = seq.main()
        seq = Seq_error_accelerate(infer_file, truth_file2, "hap2", "spechla")
        result_dict_2 = seq.main()
        for geno_name in result_dict_1.keys():
            if geno_name in result_dict_1 and geno_name not in result_dict_2:
                pass
            elif geno_name not in result_dict_1 and geno_name in result_dict_2:
                result_dict_1[geno_name] = result_dict_2[geno_name]
            else:
                if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
                    result_dict_1[geno_name] = result_dict_2[geno_name]
        gene_dict = {}
        for geno_name in result_dict_1.keys():
            gene = geno_name[:-2]
            if gene not in gene_dict:
                gene_dict[gene] = []
            gene_dict[gene].append(result_dict_1[geno_name])
        for gene in gene_dict:
            if gene not in self.spechla_exon_accuracy:
                self.spechla_exon_accuracy[gene] = 0
            for i in range(2):
                if gene_dict[gene][i][0]  + gene_dict[gene][i][1] == 0:
                    self.spechla_exon_accuracy[gene] += 1
            base_error = (gene_dict[gene][0][0] + gene_dict[gene][1][0])/2
            short_gap_error = (gene_dict[gene][0][1] + gene_dict[gene][1][1])/2
            gap_precision = (gene_dict[gene][0][2] + gene_dict[gene][1][2])/2
            map_len = (gene_dict[gene][0][3] + gene_dict[gene][1][3])/2
            print (sample, gene, base_error, short_gap_error, map_len, Method)
            data.append([sample, gene, base_error, short_gap_error, map_len, Method])
        return data

    def for_kourami(self, sample, data, Method):
        truth_file1 = self.record_truth_file_dict[sample][0]
        truth_file2 = self.record_truth_file_dict[sample][1]
        outdir = self.kourami_dir + "/" + sample + "/"
        for file in os.listdir(outdir):    
            mat = re.search("(.*?)_(.*?).typed.fa", file) 
            if not mat:  
                continue
            gene = mat.group(2)
            if gene not in self.kourami_gene_count:
                self.kourami_gene_count[gene] = 0
            self.kourami_gene_count[gene] += 1
            infer_file1 = outdir + "kourami.hla.allele.1.HLA_%s.fasta"%(gene)
            infer_file2 = outdir + "kourami.hla.allele.2.HLA_%s.fasta"%(gene)
            with open(outdir + file) as handle:
                i = 1
                for record in SeqIO.parse(handle, "fasta"):
                    if i == 1:
                        with open(infer_file1, "w") as output_handle:
                            SeqIO.write(record, output_handle, "fasta")
                    else:
                        with open(infer_file2, "w") as output_handle:
                            SeqIO.write(record, output_handle, "fasta")
                    i += 1
            infer_file1 = splitN_for_kourami(infer_file1, gene, 0)
            infer_file2 = splitN_for_kourami(infer_file2, gene, 1)
        infer_file = outdir + "/kourami.hla.allele.all.fasta"
        os.system("cat %s/*.split.exon.fasta>%s"%(outdir, infer_file))


        seq = Seq_error_accelerate(infer_file, truth_file1, "hap1", "kourami")
        result_dict_1 = seq.main()
        seq = Seq_error_accelerate(infer_file, truth_file2, "hap2", "kourami")
        result_dict_2 = seq.main()
        for geno_name in result_dict_1.keys():
            if geno_name in result_dict_1 and geno_name not in result_dict_2:
                pass
            elif geno_name not in result_dict_1 and geno_name in result_dict_2:
                result_dict_1[geno_name] = result_dict_2[geno_name]
            else:
                if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
                    result_dict_1[geno_name] = result_dict_2[geno_name]
        gene_dict = {}
        for geno_name in result_dict_1.keys():
            gene = geno_name[:-2]
            if gene not in gene_dict:
                gene_dict[gene] = []
            gene_dict[gene].append(result_dict_1[geno_name])
        for gene in gene_dict:
            if gene not in self.kourami_exon_accuracy:
                self.kourami_exon_accuracy[gene] = 0
            for i in range(2):
                if gene_dict[gene][i][0]  + gene_dict[gene][i][1] == 0:
                    self.kourami_exon_accuracy[gene] += 1
            base_error = (gene_dict[gene][0][0] + gene_dict[gene][1][0])/2
            short_gap_error = (gene_dict[gene][0][1] + gene_dict[gene][1][1])/2
            gap_precision = (gene_dict[gene][0][2] + gene_dict[gene][1][2])/2
            map_len = (gene_dict[gene][0][3] + gene_dict[gene][1][3])/2
            print (sample, gene, base_error, short_gap_error, map_len,  Method)
            data.append([sample, gene, base_error, short_gap_error, map_len,  Method])
        return data

    def for_hisat(self, sample, data):
        if sample not in self.record_hisat_output:
            self.record_hisat_output[sample] = {}

        outdir = self.hisat_dir + "/" + sample
        report_file = outdir + "/assembly_graph-hla.%s_final_extract_1_fq_gz-hla-extracted-1_fq.report"%(sample)
        hisat_fasta = outdir + "/assembly_graph-hla.%s_final_extract_1_fq_gz-hla-extracted-1_fq.fasta"%(sample)
        if not os.path.isfile(hisat_fasta):
            print ("no", hisat_fasta)
            return data
        contig_gene = self.contig2gene(report_file)
        # print (contig_gene)

        truth_file1 = self.record_truth_file_dict[sample][0]
        truth_file2 = self.record_truth_file_dict[sample][1]
        if not os.path.isfile(truth_file1) or not os.path.isfile(truth_file2):
            print ("############",truth_file1)
            return data

        with open(hisat_fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id[1] == "0":
                    hap_index = int(record.id[3]) + 1
                    gene = contig_gene[record.id]
                    record.id = "HLA_%s_%s"%(gene, hap_index-1)
                    infer_file1 = outdir + "/hisat.hla.allele.%s.HLA_%s.fasta"%(hap_index, gene)
                    with open(infer_file1, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")

        infer_file = outdir + "/hisat.hla.allele.all.fasta"
        os.system("cat %s/hisat.hla.allele.*.HLA_*.fasta>%s"%(outdir, infer_file))

        seq = Seq_error_accelerate(infer_file, truth_file1, "hap1", "hisat")
        result_dict_1 = seq.main()
        seq = Seq_error_accelerate(infer_file, truth_file2, "hap2", "hisat")
        result_dict_2 = seq.main()
        for geno_name in result_dict_1.keys():
            if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
                result_dict_1[geno_name] = result_dict_2[geno_name]
        gene_dict = {}
        for geno_name in result_dict_1.keys():
            gene = geno_name[:-2]
            if gene not in gene_dict:
                gene_dict[gene] = []
            gene_dict[gene].append(result_dict_1[geno_name])

        for gene in gene_dict:
            if gene not in self.record_hisat_output[sample]:
                self.record_hisat_output[sample][gene] = 0
            if gene not in self.hisat_gene_count:
                self.hisat_gene_count[gene] = 0
            self.hisat_gene_count[gene] += 1
            if len(gene_dict[gene]) == 2:
                base_error = (gene_dict[gene][0][0] + gene_dict[gene][1][0])/2
                short_gap_error = (gene_dict[gene][0][1] + gene_dict[gene][1][1])/2
                gap_precision = (gene_dict[gene][0][2] + gene_dict[gene][1][2])/2
                map_len = (gene_dict[gene][0][3] + gene_dict[gene][1][3])/2
            else:
                base_error = gene_dict[gene][0][0] 
                short_gap_error = gene_dict[gene][0][1] 
                gap_precision = gene_dict[gene][0][2] 
                map_len = gene_dict[gene][0][3] 
            print (sample, gene, base_error, short_gap_error, map_len, gap_precision, "HISAT")
            data.append([sample, gene, base_error, short_gap_error, map_len, gap_precision, "HISAT"])
        return data             

    def main_exon(self):
        data = []
        sample_num = 0
        for sample in self.record_truth_file_dict.keys():
            if not os.path.isfile(self.work_dir + f"/spechla//{sample}/hla.allele.1.HLA_A.fasta"):
                continue
            # if sample != "HG02587":
            #     continue
            print (sample)
            sample_num += 1
            self.spechla_dir = self.work_dir + "/spechla_exon/"
            data = self.for_spechla_exon(sample, data, "SpecHLA")
            data = self.for_kourami(sample, data, "Kourami")
            # break
        print ("total sample number:", sample_num)
        df = pd.DataFrame(data, columns = ["sample", "gene", "mismatch_rate", "gap_rate", "map_len", "Methods"])
        df.to_csv(self.work_dir + '/hgsvc_haplo_exon_assess.csv', sep=',')
        print ("reconstructed gene num of kourami", self.kourami_gene_count)
        print ("kourami", self.kourami_exon_accuracy)
        print ("spechla", self.spechla_exon_accuracy)
        for gene in self.spechla_exon_accuracy:
            print ("spechla", gene, self.spechla_exon_accuracy[gene], round(self.spechla_exon_accuracy[gene]/(sample_num*2), 2))
            if gene in self.kourami_exon_accuracy:
                print ("kourami", gene, self.kourami_exon_accuracy[gene], round(self.kourami_exon_accuracy[gene]/(sample_num*2), 2))

    def contig2gene(self, report_file):
        contig_gene = {}
        f = open(report_file, "r")
        for line in f:
            line = line.strip()
            if re.search("Node", line):
                contig_name = line.split()[-1].strip() 
            elif re.search("vs", line):
                array = line.split()
                gene = array[-1].split("*")[0]
                # if contig_name[:5] == "(0-0)" or contig_name[:5] == "(0-1)":
                contig_gene[contig_name] = gene
            else:
                continue
        return contig_gene

    def test(self):
        data = []
        sample_num = 0
        for sample in self.record_truth_file_dict.keys():
            if not os.path.isfile(f"/mnt/d/HLAPro_backup/haplotype/spechla//{sample}/hla.allele.1.HLA_A.fasta"):
                continue
            if sample != "NA19240":
                continue
            print (sample)
            sample_num += 1
            self.spechla_dir = "/mnt/d/HLAPro_backup/haplotype/spechla_pac/"
            data = self.for_spechla(sample, data, "SpecHLA")
            # data = self.for_hisat(sample, data)
            # self.spechla_dir = "/mnt/d/HLAPro_backup/haplotype/spechla_sv/"
            # data = self.for_spechla(sample, data, "SpecHLA-SV")
        print ("total sample number:", sample_num)

    def get_low_depth_bed(self, sample):
        low_depth_dict = {}
        low_bed = self.spechla_dir + sample + "/low_depth.bed" 
        f = open(low_bed, 'r')
        for line in f:
            array = line.strip().split()
            gene = array[0]
            low_depth_dict[gene] = 1
        f.close()
        return low_depth_dict

    def main_HG002(self):
        data = []
        sample="HG002"
        self.record_truth_file_dict[sample] = ["/mnt/d/HLAPro_backup/trio/truth_MHC/H1-asm.fa","/mnt/d/HLAPro_backup/trio/truth_MHC/H2-asm.fa"]
        # self.hisat_dir = "/mnt/d/HLAPro_backup/trio/Trio/hisat/"
        # data = self.for_hisat(sample, data)
        self.spechla_dir = "/mnt/d/HLAPro_backup/trio/"
        data = self.for_spechla(sample, data, "SpecHLA")
        print (data)

class Assess_novel():

    def __init__(self):
        self.record_truth_file_dict = {}
        self.work_dir = "/mnt/e/my_HLA/novel_seq"
        self.spechla_dir = self.work_dir + "/SpecHLA/"
        self.hisat_dir = self.work_dir + "/Hisat/"
        self.get_phased_assemblies()
        self.hisat_gene_count = {}
        self.spechla_gene_count = {}
        self.record_hisat_output = {} # record whether hisat reconstruct the seq
        self.N_ratio_cutoff = Min_N_ratio
       
    def get_phased_assemblies(self):
        record_truth_file_dict = {}
        inpath = "/mnt/e/my_HLA/novel_seq/simu/"
        for file in os.listdir(inpath):
            if file[-5:] != "fasta" or len(file.split(".")) == 8 :
                continue
            sample = file.split(".")[0] + "_T_50-50"
            fasta = inpath + file
            # self.split_hap(fasta)
            record_truth_file_dict[sample] = [fasta + ".0.fasta", fasta + ".1.fasta"]
        # print (record_truth_file_dict)
        self.record_truth_file_dict =  record_truth_file_dict

    def split_hap(self, fasta):
        allele_list = []
        for line in open(fasta):
            line = line.strip()
            if line[0] == ">":
                allele = line[1:]
                allele_list.append(allele)
        for i in range(2):
            file = fasta + "." + str(i) + ".fasta"
            for j in range(i, len(allele_list), 2):
                if j <= 1:
                    os.system(f"samtools faidx {fasta} {allele_list[j]} >{file}")
                else:
                    os.system(f"samtools faidx {fasta} {allele_list[j]} >>{file}")

    def main(self):
        data = []
        sample_num = 0
        print ("number of samples with phased assembly", len(self.record_truth_file_dict))
        for sample in self.record_truth_file_dict.keys():
            print (sample)
            # if sample != "HLA_50_T_50-50":
                # continue
            sample_num += 1
            data = self.for_hisat(sample, data)
            data = self.for_spechla(sample, data, "SpecHLA")
            # break
        print ("total sample number:", sample_num)
        df = pd.DataFrame(data, columns = ["sample", "gene", "mismatch_rate", "gap_rate", "map_len", "sequence_precision", "Methods"])
        df.to_csv('/mnt/d/HLAPro_backup/haplotype_v2/novel_haplo_assess.csv', sep=',')
        print ("reconstructed gene num of hisat", self.hisat_gene_count)
        print ("gene num of spechla with < %s N"%(self.N_ratio_cutoff), self.spechla_gene_count)
        data = []
        for gene in gene_list:
            gene = "HLA_" + gene
            if gene in self.hisat_gene_count:
                allele_num = self.hisat_gene_count[gene]
            else:
                allele_num = 0
            data.append([allele_num, "HISAT", gene])
            data.append([self.spechla_gene_count[gene], "SpecHLA", gene])
        df = pd.DataFrame(data, columns = ["Allele_num", "Methods", "gene"])
        df.to_csv('/mnt/d/HLAPro_backup/haplotype_v2/novel_haplo_assess_allele_num.csv', sep=',')

    def for_spechla(self, sample, data, Method):
        outdir = self.spechla_dir + "/" + sample
        low_depth_dict = self.get_low_depth_bed(sample)
        truth_file1 = self.record_truth_file_dict[sample][0]
        # print (self.record_truth_file_dict[sample], truth_file1)
        truth_file2 = self.record_truth_file_dict[sample][1]
        infer_file = outdir + "/hla.allele.all.fasta"
        os.system("cat %s/hla.allele.*.HLA_*.fasta>%s"%(outdir, infer_file))

        if not os.path.isfile(truth_file1) or not os.path.isfile(truth_file2):
            print ("############", truth_file1)
            return data       

        seq = Seq_error_accelerate_sim(infer_file, truth_file1, "hap1", "spechla")
        result_dict_1 = seq.main()
        seq = Seq_error_accelerate_sim(infer_file, truth_file2, "hap2", "spechla")
        result_dict_2 = seq.main()
        for geno_name in result_dict_1.keys():
            if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
                result_dict_1[geno_name] = result_dict_2[geno_name]
        gene_dict = {}
        for geno_name in result_dict_1.keys():
            gene = geno_name[:-2]
            if gene not in gene_dict:
                gene_dict[gene] = []
            gene_dict[gene].append(result_dict_1[geno_name])
        for gene in gene_dict:
            # print (gene_dict[gene][0][0], gene_dict[gene][1][0], gene_dict[gene][0][3], gene_dict[gene][1][3])
            if gene_dict[gene][0][0] < gene_dict[gene][1][0]:
                base_error = gene_dict[gene][0][0]
                short_gap_error = gene_dict[gene][0][1]
                gap_precision = gene_dict[gene][0][2]
                map_len = gene_dict[gene][0][3]
                N_ratio = gene_dict[gene][0][4]
            else:
                base_error = gene_dict[gene][1][0]
                short_gap_error = gene_dict[gene][1][1]
                gap_precision = gene_dict[gene][1][2]
                map_len = gene_dict[gene][1][3]
                N_ratio = gene_dict[gene][1][4]                

            # if gene in self.record_hisat_output[sample]:
            # if gene not in low_depth_dict:
            # if True:
            if N_ratio < self.N_ratio_cutoff:
                if gene not in self.spechla_gene_count:
                    self.spechla_gene_count[gene] = 0
                self.spechla_gene_count[gene] += 1
                print (sample, gene, base_error, short_gap_error, map_len, gap_precision, Method, N_ratio, gene_dict[gene][0][0], gene_dict[gene][1][0])
                data.append([sample, gene, base_error, short_gap_error, map_len, gap_precision, Method])
            else:
                print (sample, gene, base_error, short_gap_error, map_len, gap_precision, Method, N_ratio)
        return data

    def for_hisat(self, sample, data):
        if sample not in self.record_hisat_output:
            self.record_hisat_output[sample] = {}

        outdir = self.hisat_dir + "/" + sample
        report_file = outdir + "/assembly_graph-hla.%s_read1_fastq_gz-hla-extracted-1_fq.report"%(sample)
        hisat_fasta = outdir + "/assembly_graph-hla.%s_read1_fastq_gz-hla-extracted-1_fq.fasta"%(sample)
        if not os.path.isfile(hisat_fasta):
            print ("no", hisat_fasta)
            return data
        contig_gene = self.contig2gene(report_file)
        # print (contig_gene)

        truth_file1 = self.record_truth_file_dict[sample][0]
        truth_file2 = self.record_truth_file_dict[sample][1]
        if not os.path.isfile(truth_file1) or not os.path.isfile(truth_file2):
            print ("############",truth_file1)
            return data

        with open(hisat_fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id[1] == "0":
                    hap_index = int(record.id[3]) + 1
                    gene = contig_gene[record.id]
                    record.id = "HLA_%s_%s"%(gene, hap_index-1)
                    infer_file1 = outdir + "/hisat.hla.allele.%s.HLA_%s.fasta"%(hap_index, gene)
                    with open(infer_file1, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")

        infer_file = outdir + "/hisat.hla.allele.all.fasta"
        os.system("cat %s/hisat.hla.allele.*.HLA_*.fasta>%s"%(outdir, infer_file))

        seq = Seq_error_accelerate_sim(infer_file, truth_file1, "hap1", "hisat")
        result_dict_1 = seq.main()
        seq = Seq_error_accelerate_sim(infer_file, truth_file2, "hap2", "hisat")
        result_dict_2 = seq.main()
        for geno_name in result_dict_1.keys():
            if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
                result_dict_1[geno_name] = result_dict_2[geno_name]
        gene_dict = {}
        for geno_name in result_dict_1.keys():
            gene = geno_name[:-2]
            if gene not in gene_dict:
                gene_dict[gene] = []
            gene_dict[gene].append(result_dict_1[geno_name])

        for gene in gene_dict:
            if gene not in self.record_hisat_output[sample]:
                self.record_hisat_output[sample][gene] = 0
            if gene not in self.hisat_gene_count:
                self.hisat_gene_count[gene] = 0
            self.hisat_gene_count[gene] += 1
            if len(gene_dict[gene]) == 2:
                if gene_dict[gene][0][0] < gene_dict[gene][1][0]:
                    base_error = gene_dict[gene][0][0]
                    short_gap_error = gene_dict[gene][0][1]
                    gap_precision = gene_dict[gene][0][2]
                    map_len = gene_dict[gene][0][3]
                else:
                    base_error = gene_dict[gene][1][0]
                    short_gap_error = gene_dict[gene][1][1]
                    gap_precision = gene_dict[gene][1][2]
                    map_len = gene_dict[gene][1][3] 
            else:
                base_error = gene_dict[gene][0][0] 
                short_gap_error = gene_dict[gene][0][1] 
                gap_precision = gene_dict[gene][0][2] 
                map_len = gene_dict[gene][0][3] 
            print (sample, gene, base_error, short_gap_error, map_len, gap_precision, "HISAT")
            data.append([sample, gene, base_error, short_gap_error, map_len, gap_precision, "HISAT"])
        return data             

    def contig2gene(self, report_file):
        contig_gene = {}
        f = open(report_file, "r")
        for line in f:
            line = line.strip()
            if re.search("Node", line):
                contig_name = line.split()[-1].strip() 
            elif re.search("vs", line):
                array = line.split()
                gene = array[-1].split("*")[0]
                # if contig_name[:5] == "(0-0)" or contig_name[:5] == "(0-1)":
                contig_gene[contig_name] = gene
            else:
                continue
        return contig_gene

    def test(self):
        data = []
        sample_num = 0
        for sample in self.record_truth_file_dict.keys():
            if not os.path.isfile(f"/mnt/d/HLAPro_backup/haplotype/spechla//{sample}/hla.allele.1.HLA_A.fasta"):
                continue
            if sample != "NA19240":
                continue
            print (sample)
            sample_num += 1
            self.spechla_dir = "/mnt/d/HLAPro_backup/haplotype/spechla_pac/"
            data = self.for_spechla(sample, data, "SpecHLA")
            # data = self.for_hisat(sample, data)
            # self.spechla_dir = "/mnt/d/HLAPro_backup/haplotype/spechla_sv/"
            # data = self.for_spechla(sample, data, "SpecHLA-SV")
        print ("total sample number:", sample_num)

    def get_low_depth_bed(self, sample):
        low_depth_dict = {}
        low_bed = self.spechla_dir + sample + "/low_depth.bed" 
        f = open(low_bed, 'r')
        for line in f:
            array = line.strip().split()
            gene = array[0]
            low_depth_dict[gene] = 1
        f.close()
        return low_depth_dict

def eva_simu_trio():
    truth_dir = "/mnt/d/HLAPro_backup/trio/simu_pedigree/data/truth/"
    data = []
    for i in range(50):
        sample = "child_%s"%(i)
        outdir = "/mnt/d/HLAPro_backup/trio/simu_pedigree/output/" + sample + "/"
        print (sample, "SpecHLA-trio")
        for gene in gene_list:
            base_error, short_gap_error, gap_recall, gap_precision = each_simulated_sample(gene, outdir, sample, truth_dir)
            data.append([sample, gene, base_error, short_gap_error, "SpecHLA-trio"])
        print (sample, "SpecHLA")
        outdir = "/mnt/d/HLAPro_backup/trio/simu_pedigree/output_no_trio/" + sample + "/"
        for gene in gene_list:
            base_error, short_gap_error, gap_recall, gap_precision = each_simulated_sample(gene, outdir, sample, truth_dir)
            data.append([sample, gene, base_error, short_gap_error, "SpecHLA"])
    df = pd.DataFrame(data, columns = ["sample", "gene", "base_error", "short_gap_error", "Methods"])
    df.to_csv('/mnt/d/HLAPro_backup/trio/simu_pedigree/haplo_assess_trio.csv', sep=',')

def eva_real_trio():
    # outdir = "/mnt/d/HLAPro_backup/trio/HG002/"
    dict = {"SpecHLA-trio":"spechla", "SpecHLA":"spechla_no_trio"}
    truth_file1_dict = {"NA12878":"/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta",
                        "NA19240":"/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_NA19240_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta"}
    truth_file2_dict = {"NA12878":"/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta",
                        "NA19240":"/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_NA19240_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta"}
    data = []
    for sample in ["NA12878", "NA19240"]:
        for method in dict.keys():
            # sample = "NA12878"
            # truth_file1 = "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta"
            # truth_file2 = "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta"
            # sample = "NA19240"
            # truth_file1 =  "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_NA19240_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta"
            # truth_file2 =  "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_NA19240_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta"
            # sample = "HG00733"
            # truth_file1 =  "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_HG00733_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta"
            # truth_file2 =  "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_HG00733_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta"
            # sample = "HG00514"
            # truth_file1 =  "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_HG00514_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta"
            # truth_file2 =  "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/v12_HG00514_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta"
            truth_file1, truth_file2 = truth_file1_dict[sample], truth_file2_dict[sample]
            outdir = "/mnt/d/HLAPro_backup/trio/real_trio/%s/"%(dict[method]) + sample + "/"
            infer_file = outdir + "/hla.allele.all.fasta"
            os.system("cat %s/hla.allele.*.HLA_*.fasta>%s"%(outdir, infer_file))
            

            seq = Seq_error_accelerate(infer_file, truth_file1, "hap1", "spechla")
            result_dict_1 = seq.main()
            seq = Seq_error_accelerate(infer_file, truth_file2, "hap2", "spechla")
            result_dict_2 = seq.main()
            for geno_name in result_dict_1.keys():
                if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
                    result_dict_1[geno_name] = result_dict_2[geno_name]
            gene_dict = {}
            for geno_name in result_dict_1.keys():
                gene = geno_name[:-2]
                if gene not in gene_dict:
                    gene_dict[gene] = []
                gene_dict[gene].append(result_dict_1[geno_name])
            for gene in gene_dict:
                base_error = (gene_dict[gene][0][0] + gene_dict[gene][1][0])/2
                short_gap_error = (gene_dict[gene][0][1] + gene_dict[gene][1][1])/2
                gap_precision = (gene_dict[gene][0][2] + gene_dict[gene][1][2])/2
                map_len = (gene_dict[gene][0][3] + gene_dict[gene][1][3])/2
                N_ratio = (gene_dict[gene][0][4] + gene_dict[gene][1][4])/2
                if N_ratio < Min_N_ratio:
                    print (sample, gene, base_error, short_gap_error, map_len, gap_precision, method)
                    data.append([sample, gene, base_error, short_gap_error, map_len, gap_precision, method])


    df = pd.DataFrame(data, columns = ["sample", "gene", "base_error", "short_gap_error", "map_len", "gap_precision", "Methods"])
    df.to_csv('/mnt/d/HLAPro_backup/trio/real_trio//real_trio_haplo_assess.csv', sep=',')

def eva_real_hybrid():
    # outdir = "/mnt/d/HLAPro_backup/trio/HG002/"
    dict = {"SpecHLA-hybrid":"spechla_with_pac", "SpecHLA":"spechla_no_pac"}
    # dict = {"SpecHLA-hybrid":"spechla_with_pac"}
    # dict = {"SpecHLA-hybrid-sv":"spechla"}
    truth_dir = "/mnt/d/HLAPro_backup/haplotype/my_HLA/assembly/"
    truth_file1_dict = {
                        "HG00733":"v12_HG00733_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta",
                        "HG00731":"v12_HG00731_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta",
                        "HG00732":"v12_HG00732_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta",
                        "NA19238":"v12_NA19238_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta",
                        "NA19239":"v12_NA19239_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta",
                        "NA19240":"v12_NA19240_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta",
                        "HG00514":"v12_HG00514_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta"
                        }

    truth_file2_dict = {
                        "HG00733":"v12_HG00733_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta",
                        "HG00731":"v12_HG00731_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta",
                        "HG00732":"v12_HG00732_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta",
                        "NA19238":"v12_NA19238_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta",
                        "NA19239":"v12_NA19239_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta",
                        "NA19240":"v12_NA19240_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta",
                        "HG00514":"v12_HG00514_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta"
                        }
    data = []
    for sample in truth_file1_dict.keys():
    # for sample in ["HG00733"]:
        for method in dict.keys():
            truth_file1, truth_file2 = truth_dir + truth_file1_dict[sample], truth_dir + truth_file2_dict[sample]
            outdir = "/mnt/d/HLAPro_backup/hybrid/real_hybrid/%s/"%(dict[method]) + sample + "/"
            infer_file = outdir + "/hla.allele.all.fasta"
            os.system("cat %s/hla.allele.*.HLA_*.fasta>%s"%(outdir, infer_file))
            

            seq = Seq_error_accelerate(infer_file, truth_file1, "hap1", "spechla")
            result_dict_1 = seq.main()
            seq = Seq_error_accelerate(infer_file, truth_file2, "hap2", "spechla")
            result_dict_2 = seq.main()
            for geno_name in result_dict_1.keys():
                if result_dict_1[geno_name][0] > result_dict_2[geno_name][0]:
                    result_dict_1[geno_name] = result_dict_2[geno_name]
            gene_dict = {}
            for geno_name in result_dict_1.keys():
                gene = geno_name[:-2]
                if gene not in gene_dict:
                    gene_dict[gene] = []
                gene_dict[gene].append(result_dict_1[geno_name])
            for gene in gene_dict:
                base_error = (gene_dict[gene][0][0] + gene_dict[gene][1][0])/2
                short_gap_error = (gene_dict[gene][0][1] + gene_dict[gene][1][1])/2
                gap_precision = (gene_dict[gene][0][2] + gene_dict[gene][1][2])/2
                map_len = (gene_dict[gene][0][3] + gene_dict[gene][1][3])/2
                N_ratio = (gene_dict[gene][0][4] + gene_dict[gene][1][4])/2

                if N_ratio < Min_N_ratio:
                    print (sample, gene, base_error, short_gap_error, map_len, gap_precision, method)
                    data.append([sample, gene, base_error, short_gap_error, map_len, gap_precision, method])


    df = pd.DataFrame(data, columns = ["sample", "gene", "base_error", "short_gap_error", "map_len", "gap_precision", "Methods"])
    df.to_csv('/mnt/d/HLAPro_backup/hybrid/real_hybrid/real_hybrid_haplo_assess.csv', sep=',')

if __name__ == "__main__":
    Min_N_ratio = 0.3

    if len(sys.argv) == 1:
        print ("############")
        # nov = Assess_novel()
        # nov.main()
        # ass = Assess_hgsvc2()
        # ass.main()
        # ass.main_exon()
        # ass.test()
        # eva_simu_trio()
        # eva_real_trio()
        # eva_real_hybrid()
        # eva_data_types_spechla()
        eva_allele_imblance()
        # eva_HG002_kourami()
        # eva_pedigree_spechla()
        # eva_HG002_spechla()
        # eva_hgsvc2_spechla()
        # eva_hgsvc2_spechla_accelerate()
        # eva_HG002_hisat()
        # eva_data_types_spechla()
    else:
        database = sys.argv[1]
        record_true_file = sys.argv[2]
        sample_dir = sys.argv[3]
        sample_name = sys.argv[4]
        method = sys.argv[5] # spechla, hisat, or kourami
        exon_flag = sys.argv[6] # exon or not
                
        if method == "hisat":
            split_hisat_fasta(sample_dir, sample_name)
        elif method == "kourami":
            split_kourami_fasta(sample_dir, sample_name)
        eva_simu(database, record_true_file, sample_dir, sample_name)
        
