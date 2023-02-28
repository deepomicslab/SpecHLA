"""
Infer G group HLA type from alignment of phased assembly and exon database

wangshuai, Feb 28, 2023
"""

import os, re
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
import re
import pandas as pd
import numpy as np
from collections import Counter
import argparse


from eva_type_accuracy import read_G_annotation
from get_HLA_alleles_from_assembly import get_phased_assemblies


def infer_best_match_exon(exon_sam, hap_index, sample):
    
    contig_allele_perfect_exon_dict = {}
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")
    # Split the SAM record into fields
    record_hit_exon_times = {}
    
    for line in open(exon_sam, "r"):
        # Skip header lines
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        exon_name = fields[0]
        cigar = fields[5]
        contig_name = fields[2]

        nm_tag = [tag for tag in fields[11:] if tag.startswith("NM:i:")]
        if len(nm_tag) == 1:
            num_mismatches = int(nm_tag[0].split(":")[2])
        else:
            num_mismatches = 100000


        allele_name = fields[0].split("|")[0]
        gene = allele_name.split("*")[0].split("-")[1]
        truth_contig = all_sample_true_dict[sample][gene][hap_index]
        if fields[2] != truth_contig:
            continue

        if gene not in contig_allele_perfect_exon_dict:
            contig_allele_perfect_exon_dict[gene] = {}

        if exon_name not in record_hit_exon_times:
            if allele_name not in contig_allele_perfect_exon_dict[gene]:
                contig_allele_perfect_exon_dict[gene][allele_name] = {}
            if contig_name not in contig_allele_perfect_exon_dict[gene][allele_name]:
                contig_allele_perfect_exon_dict[gene][allele_name][contig_name] = num_mismatches
            else:
                contig_allele_perfect_exon_dict[gene][allele_name][contig_name] += num_mismatches
            record_hit_exon_times[exon_name] = 1
    # print (contig_allele_perfect_exon_dict)
    infer_best(contig_allele_perfect_exon_dict)


def most_common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]



def infer_best(contig_allele_perfect_exon_dict):
    allele_perfect_exon_dict = {}  
    for gene in contig_allele_perfect_exon_dict:
        allele_perfect_exon_dict[gene] = {}
        for allele_name in contig_allele_perfect_exon_dict[gene]:
            min_mismatch = 10000000
            for contig_name in contig_allele_perfect_exon_dict[gene][allele_name]:
                if contig_allele_perfect_exon_dict[gene][allele_name][contig_name] < min_mismatch:
                    min_mismatch =  contig_allele_perfect_exon_dict[gene][allele_name][contig_name]
            allele_perfect_exon_dict[gene][allele_name] = min_mismatch
        sorted_identity_record = sorted(allele_perfect_exon_dict[gene].items(), key=lambda x: x[1], reverse = False)

        top_alleles = []
        allele_g_dict = {}
        max_identity = sorted_identity_record[0][1]
        for allele, info in sorted_identity_record:
            if info == max_identity:
                if convert_G(allele) not in allele_g_dict:
                    allele_g_dict[convert_G(allele)] = allele
                top_alleles.append(convert_G(allele))
            else:
                break
        # print (top_alleles)
        most_common_allele = most_common(top_alleles)
        if gene not in record_truth[sample]:
            record_truth[sample][gene] = []
        record_truth[sample][gene] += [allele_g_dict[most_common_allele], most_common_allele]
        # return (most_common_allele)
        # print (most_common_allele)

        # best_match_allele = sorted_identity_record[0][0][4:]
        # best_G_type = convert_G(best_match_allele)
        # # print (best_match_allele, best_G_type, sorted_identity_record[:30])
        # if gene not in record_truth[sample]:
        #     record_truth[sample][gene] = []
        # record_truth[sample][gene] += [best_match_allele, best_G_type]
        # print (record_truth)


def convert_G(allele):
    allele = re.sub(":","_",allele)
    allele = re.sub("HLA-","",allele)
    allele = re.sub("\*","_",allele)
    allele = allele.split(";")[0]
    if allele in G_annotation_dict:
        allele = G_annotation_dict[allele]
    elif allele + "_01" in G_annotation_dict:
        allele = G_annotation_dict[allele + "_01"]
    elif allele + "_01_01" in G_annotation_dict:
        allele = G_annotation_dict[allele + "_01_01"]
    G_type = allele
    return G_type


def get_truth_contig():
    all_sample_true_dict = {}
    
    for line in open(extract_HLA_seq):
        if line[0] == ">":
            array = line[1:].split()
            array[2] = re.sub(":","_",array[2])
            array[2] = re.sub("\*","_",array[2])
            allele = array[2]
            sample = array[0].split(".")[0]
            gene = array[0].split("-")[1]
            hap_index = int(array[0].split(".")[1][-1]) - 1
            contig = array[1].split(":")[0]

            if sample not in all_sample_true_dict:
                all_sample_true_dict[sample] = {}
            if gene not in all_sample_true_dict[sample]:
                all_sample_true_dict[sample][gene] = ['', '']
            all_sample_true_dict[sample][gene][hap_index] = contig
    return all_sample_true_dict


def blast_exon(sample, hap_index):
    # command = f"{minimap_path} {record_truth_file_dict[sample][hap_index]} {HLA_data}  -o {result_path}/{sample}.h{hap_index+1}.paf -t 10"
    command = f"blastn -query {record_truth_file_dict[sample][hap_index]} -db {single_exon_database_fasta} -out {result_path}/{sample}.h{hap_index+1}.exon.blast -outfmt 7 -max_target_seqs 60000 -num_threads 15"
    # print (command)
    os.system(command)

if __name__ == "__main__":

    samples_list = ['HG00096', 'HG00171', 'HG00512', 'HG00513', 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'HG00864', 'HG01114', 'HG01505', 'HG01596', 'HG02011', 'HG02492', 'HG02587', 'HG02818', 'HG03009', 'HG03065', 'HG03125', 'HG03371', 'HG03486', 'HG03683', 'HG03732', 'NA12878', 'NA18534', 'NA18939', 'NA19238', 'NA19239', 'NA19240', 'NA19650', 'NA19983', 'NA20509', 'NA20847', 'NA24385']
    # samples_list = ["HG00096"]
    gene_list = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    extract_HLA_seq = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/extracted_HLA_alleles.fasta"
    single_exon_database_fasta = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/xml/hla_exons.fasta"
    result_path = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/"

    all_sample_true_dict = get_truth_contig()
    all_sample_true_dict["HG01505"]["DQB1"][1] = "cluster9_contig_10"
    all_sample_true_dict["HG01505"]["DRB1"][1] = "cluster9_contig_10"
    all_sample_true_dict["HG03732"]["DRB1"][1] = "cluster21_contig_60"
    all_sample_true_dict["NA19240"]["DRB1"][0] = "cluster16_scaffold_10"
    all_sample_true_dict["NA19240"]["DRB1"][1] = "cluster16_contig_7"
    record_truth_file_dict = get_phased_assemblies()

    G_annotation_dict = read_G_annotation()
    record_truth = {}
    for sample in samples_list:
        print (sample)
        record_truth[sample] = {}
        sample_save_alignments_dict = {}
        for hap_index in range(2):
            # blast_exon(sample, hap_index)

            input_sam_exon = f"/mnt/d/HLAPro_backup/minor_rev/extract_alleles/{sample}.h{hap_index+1}.exon.sam"
            infer_best_match_exon(input_sam_exon, hap_index, sample)
        # break

    file = "/mnt/d/HLAPro_backup/haplotype_v2/hgsvc2_G_group_truth.csv"
    f = open(file, 'w')
    print ("sample,gene,allele1,G-group1,allele2,G-group2", file = f)   

    for sample in samples_list:
        for gene in gene_list:
            print (sample, gene,record_truth[sample][gene])
            print (sample, gene, record_truth[sample][gene][0], record_truth[sample][gene][1], record_truth[sample][gene][2], record_truth[sample][gene][3], sep=",", file = f)
    f.close()