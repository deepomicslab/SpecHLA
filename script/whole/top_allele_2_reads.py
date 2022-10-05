import sys
import os
from Bio import SeqIO

outdir = sys.argv[1]
detail_anno = outdir + "/hla.result.details.txt"
top_allele_fasta = outdir + "/top_allele.fasta"
top_allele_fastq = outdir + "/top_allele.fastq"
samtools = sys.path[0] + "/../../bin/samtools"
db_dir = sys.path[0] + "/../../db/HLA/whole/"
f=open(top_allele_fasta, 'w')
f.close()

def read():
    dict = {"A":'', "B":'', "C":'', "DPA1":'',"DPB1":'',"DQA1":'',"DQB1":'',"DRB1":''}
    top_allele = []
    f = open(detail_anno, 'r')
    for line in f:
        array = line.strip().split()
        if array[2] == "allele":
            continue
        allele_list = array[2].split(";")
        for allele in allele_list:
            if allele == '':
                continue
            allele_array = allele.split("*")
            gene = allele_array[0]
            dict[gene] += " " + allele
            top_allele.append(allele)
    # print (dict)
    for gene in dict.keys():
        fasta = db_dir + "HLA_" + gene + ".fasta"
        command = samtools + " faidx " + fasta + dict[gene] + ">>" + top_allele_fasta
        os.system(command)
        # print (command)

def fasta2fastq():
    f = open(top_allele_fastq, 'w')
    for r in SeqIO.parse(top_allele_fasta, "fasta"):
        r.letter_annotations["solexa_quality"] = [40] * len(r)
        print(r.format("fastq"), end='', file = f)
    f.close()

read()
fasta2fastq()
