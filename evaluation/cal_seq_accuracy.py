"""
calculate base error and gap percentage
step1: map the sequence to the truth with blastn
step2: get mapped intervals
step3: get unique intervals
step4: calculate the mismatch rate in the unique mapped region
step5: calculate gap recall and gap precision

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

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['DQA1', 'DQB1', 'DRB1']

def merge_exon(infer_file):
    used_exon = ["HLA_A:1504-1773","HLA_A:2015-2290","HLA_B:1486-1755","HLA_B:2001-2276","HLA_C:1699-1968","HLA_C:2215-2490",\
        "HLA_DQA1:5600-5848","HLA_DQB1:3073-3342","HLA_DRB1:6972-7241"]
    merge_exon_file = infer_file + ".merged.exon.fasta"
    exon_seq = ''
    with open(infer_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in used_exon or record.id.split(":")[0] in ["HLA_DPA1","HLA_DPB1"]:
                exon_seq += str(record.seq)
                exon_seq += "N"*200
    with open(merge_exon_file, "w") as output_handle:
        print (">merged_exon", file = output_handle)  
        print (exon_seq, file = output_handle) 
    return  merge_exon_file

def splitN_for_kourami(infer_file):
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
    print (">split_exon", file = output_handle)  
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

    def __init__(self, infer_hap_file, truth_hap_file, tag):
        self.infer_hap_file = infer_hap_file
        self.truth_hap_file = truth_hap_file
        self.blast_file = f"{self.infer_hap_file}.{tag}.blast"
        self.infer_hap_len_dict = {}
        self.truth_hap_len_dict = {}
        self.mapped_interval_dict = {} # for inferred hap
        self.true_mapped_interval_dict = {} # for true hap
        self.mapped_len_dict = {}
        self.true_mapped_len_dict = {}
        self.mismatch_num_dict = {}
        self.gap_open_num_dict = {}

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
                self.infer_hap_len_dict[record.id] = len(record.seq)
        with open(self.truth_hap_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.truth_hap_len_dict[record.id] = len(record.seq)

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
            if line_count[seq_name] > 20:
                continue

            self.mapped_interval_dict[seq_name] = self.get_uniq_map(map_s, map_e, mismatches, gap_opens, self.mapped_interval_dict[seq_name],seq_name)
            true_map_s = int(array[8])
            true_map_e = int(array[9])
            self.true_mapped_interval_dict[seq_name] = self.get_uniq_map(true_map_s, true_map_e, 0, 0, self.true_mapped_interval_dict[seq_name],seq_name)
            i += 1

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
            if map_s >= interval[0] and map_e <= interval[1]:
                flag = False
        if flag and map_e - map_s > 5:
            uniq_map_len = map_e - map_s + 1
            mismatches = mismatches * (uniq_map_len/origin_map_len)
            gap_opens = gap_opens * (uniq_map_len/origin_map_len)
            self.mismatch_num_dict[seq_name] += mismatches
            self.gap_open_num_dict[seq_name] += gap_opens
            mapped_interval.append([map_s, map_e])
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
        # self.blast_map("strict")
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
            print (geno_name, align.base_error, align.short_gap_error, self.mapped_len_dict[geno_name])
            result_dict[geno_name] = [align.base_error, align.short_gap_error, align.gap_precision, self.mapped_len_dict[geno_name]]
        return result_dict

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
        infer_file1 = merge_exon(infer_file1)
        infer_file2 = outdir + "hla.allele.2.HLA_%s.fasta"%(gene)
        infer_file2 = merge_exon(infer_file2)
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
                infer_file1 = merge_exon(infer_file1)
                infer_file2 = merge_exon(infer_file2)
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
    # dict = {"PE":"novel", "+PacBio":"pacbio_illumina","+Hi-C":"hic_illumina","+ONT":"ont_illumina","PacBio":"pac_alone","ONT":"ont_alone"}
    dict = {"+10X":"10x_illumina" }
    data = []
    for i in range(3):
        sample = "novel_%s"%(i)
        for plat in dict.keys():
            outdir = "/mnt/d/HLAPro_backup/pacbio/" + dict[plat]+ "/" + sample + "/"
            truth_dir = "/mnt/d/HLAPro_backup/pacbio/simulation/truth/"
            
            for gene in gene_list:
                base_error, short_gap_error, gap_recall, gap_precision = each_simulated_sample(gene, outdir, sample, truth_dir)
                data.append([sample, gene, base_error, short_gap_error, plat])
                # print (sample, gene, base_error, short_gap_error, plat)
    df = pd.DataFrame(data, columns = ["sample", "gene", "base_error", "short_gap_error", "Data"])
    df.to_csv('/mnt/d/HLAPro_backup/pacbio/novel/haplo_assess.csv', sep=',')

def eva_allele_imblance():
    # sample = "imbalance_48_52_0"
    # sample = "imbalance_40_60_0"
    # outdir = "/mnt/d/HLAPro_backup/imbalance/output/" + sample + "/"
    truth_dir = "/mnt/d/HLAPro_backup/imbalance/data/truth/"
    data = []
    for i in range(20):
        for group in ["50_50", "48_52", "40_60", "30_70", "20_80"]:
            sample = "imbalance_%s_%s"%(group, i)
            outdir = "/mnt/d/HLAPro_backup/imbalance/output/" + sample + "/"
            for weight in [0, 0.5, 1]:
                os.system("sh /mnt/d/HLAPro_backup//HLAPro/script/whole/test_SpecHLA.sh -n %s -1 x -2 x -w %s -o /mnt/d/HLAPro_backup/imbalance/output"%(sample,weight))
                weight_base_error = []
                for gene in gene_list:
                    base_error, short_gap_error, gap_recall, gap_precision = each_simulated_sample(gene, outdir, sample, truth_dir)
                    data.append([sample, gene, base_error, short_gap_error, weight, group])
                    weight_base_error.append(base_error)
                print ("################", weight, np.mean(weight_base_error))

    df = pd.DataFrame(data, columns = ["sample", "gene", "base_error", "short_gap_error", "weight", "group"])
    df.to_csv('/mnt/d/HLAPro_backup/imbalance/haplo_assess.csv', sep=',')

if __name__ == "__main__":

    if len(sys.argv) == 1:
        eva_data_types_spechla()
        # eva_allele_imblance()
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
        
