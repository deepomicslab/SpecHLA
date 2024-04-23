"""
For each HGT gene sequence and the blast result
assign the suitable types for it

It is perfect choice if the max-match-len allele and the max-identity allele are same.
Otherwise, we need to balance the match length and identity to chose the suitable alleles.

wangshuai, Apr 23, 2024
"""

import os
import re
import sys
import pysam
import argparse


def get_1_element(lst):
    return lst[1]

def get_2_element(lst):
    return lst[2]

def get_3_element(lst):
    return lst[3]

def change_allele_name(raw, new):
    with open(raw, "r") as infile, open(new, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                header = line.strip()[1:]
                contig_name = header.split()[1]
                new_header = f">{contig_name}\n"
                outfile.write(new_header)
            else:
                outfile.write(line)

def resort_list_with_same_alleles(sorted_list, first_index, second_index):
    flag = True
    while flag:
        flag = False
        new_sorted_list = sorted_list.copy()
        for i in range(len(sorted_list) - 1):
            if sorted_list[i][first_index] == sorted_list[i+1][first_index] and sorted_list[i+1][second_index] > sorted_list[i][second_index]:
                new_sorted_list[i] = sorted_list[i+1]
                new_sorted_list[i+1] = sorted_list[i]
                flag = True
        sorted_list = new_sorted_list.copy()
    # print (sorted_list[:5])
    return sorted_list
    
def get_max_alleles(sorted_list, index):
    max_value = sorted_list[0][index]
    max_allele_list = []
    for list in sorted_list:
        if list[index] == max_value:
            # max_allele_list.append(list[0])
            list = [str(x) for x in list]
            max_allele_list.append(">".join(list))
        else:
            break
    return max_allele_list

def extract_four_digits(full_name):
    a = full_name.split("*")[1]
    array = a.split(":")
    return array[0] + ":" + array[1]

def compare_match_len_and_identity(match_sorted_list, identity_sorted_list, truth_alleles):
    max_match_len = match_sorted_list[0][1]
    match_len_with_max_identity = identity_sorted_list[0][1]

    max_identity = identity_sorted_list[0][3]
    identiy_with_max_match_len = match_sorted_list[0][3]

    match_len_diff_ratio = (max_match_len - match_len_with_max_identity) / match_len_with_max_identity
    identity_diff_ratio = (max_identity - identiy_with_max_match_len) / identiy_with_max_match_len

    print ("match_len_diff_ratio", match_len_diff_ratio, "identity_diff_ratio", identity_diff_ratio)
    get_help_from_1000G = False

    if extract_four_digits(match_sorted_list[0][0]) in truth_alleles and extract_four_digits(identity_sorted_list[0][0]) not in truth_alleles:
        select_allele_list = match_sorted_list[0]
        get_help_from_1000G = True
    elif extract_four_digits(match_sorted_list[0][0]) not in truth_alleles and extract_four_digits(identity_sorted_list[0][0]) in truth_alleles:
        select_allele_list = identity_sorted_list[0]
        get_help_from_1000G = True
    elif identiy_with_max_match_len < 0.999:
        select_allele_list = identity_sorted_list[0]
    elif match_len_diff_ratio < identity_diff_ratio:
        select_allele_list = identity_sorted_list[0]
    elif match_len_diff_ratio < 0.3:
        select_allele_list = identity_sorted_list[0]
    # elif identity_diff_ratio < 0.005:
    #     select_allele_list = match_sorted_list[0]
    else:
        print (" no determine")
        
    # if get_help_from_1000G == False:
    print ("check to determine use highest identity or match length in person.")
    for allele_info in match_sorted_list[:5]:
        print(allele_info)
    print ("match bases**************************")

    
    for allele_info in identity_sorted_list[:5]:
        print(allele_info)
    print ("identity **************************")
    for allele_info in identity_sorted_list:
        if allele_info[0] == "DRB1*16:02:01:03":
            print (allele_info)
    
    print ("selected allele is ", select_allele_list[0])
    return select_allele_list
    
def get_IMGT_version():
    g_file = "%s/../db/HLA/hla_nom_g.txt"%(sys.path[0])
    G_annotation_dict = {}
    i = 0
    version_info = "N/A"
    for line in open(g_file):
        if re.search("# version:", line):
            version_info = line.strip()
    return version_info

def select_by_alignment(align_list):
    if len(align_list) == 0:
        return []
    full_result_list = []
    # match_sorted_list = sorted(align_list, key=get_1_element, reverse = True)
    match_sorted_list = sorted(align_list, key=get_2_element, reverse = True)  # match length - mismatch
    match_sorted_list = resort_list_with_same_alleles(match_sorted_list, 1, 3)
    identity_sorted_list = sorted(align_list, key=get_3_element, reverse = True)
    identity_sorted_list = resort_list_with_same_alleles(identity_sorted_list, 3, 1)
    max_match_len_alleles = get_max_alleles(match_sorted_list, 2)
    max_identity_alleles = get_max_alleles(identity_sorted_list, 3)

    # f = open(all_align_result_file, 'w')
    # for mat in match_sorted_list:
    #     print (mat, file = f)
    # f.close()

    # print (identity_sorted_list)
    intersection_alleles = list(set(max_match_len_alleles) & set(max_identity_alleles))   
    # print (">>>>>>>>>", match_sorted_list[:10])

    if len(intersection_alleles) > 0:
        select_allele_list = intersection_alleles[0].split(">")
        select_allele = select_allele_list[0]
        print (">>>>>>>>>>perfect:", select_allele)  
        full_result_list.append(select_allele_list)    
        # return select_allele_list
    else:
        # 
        print (">>>>>>>>>>not perfect, report possible alleles")  
        len_diff_cutoff = 0.01 
        ide_diff_cutoff = 0.001
        max_match_len = match_sorted_list[0][2]
        match_len_with_max_identity = identity_sorted_list[0][1]

        max_identity = identity_sorted_list[0][3]
        identiy_with_max_match_len = match_sorted_list[0][3]

        # for allele_info in match_sorted_list[:10]:
        #     print(allele_info, "length")
        # # print ("match bases**************************")

        
        # for allele_info in identity_sorted_list[:10]:
        #     print(allele_info, "identity")

        # match_len_diff_ratio = (max_match_len - match_len_with_max_identity) / match_len_with_max_identity
        # identity_diff_ratio = (max_identity - identiy_with_max_match_len) / identiy_with_max_match_len

        # if match_len_diff_ratio > 0.3:
        #     full_result_list = match_sorted_list[:5]
        # else:
        #     full_result_list = identity_sorted_list[:5] + match_sorted_list[:5]

        good_length_list = []
        for i in range(len(match_sorted_list)):
            if (max_match_len - match_sorted_list[i][2])/max_match_len <= len_diff_cutoff:
                good_length_list.append(match_sorted_list[i])
        # print (len(good_length_list), "len(good_length_list)")
        identity_sorted_list = sorted(good_length_list, key=get_3_element, reverse = True)
        identity_sorted_list = resort_list_with_same_alleles(identity_sorted_list, 3, 1)
        match_len_with_max_identity = identity_sorted_list[0][3]
        # full_result_list = identity_sorted_list[:20]
        for i in range(len(identity_sorted_list)):
            if (match_len_with_max_identity - identity_sorted_list[i][3])/match_len_with_max_identity <= ide_diff_cutoff:
                full_result_list.append(identity_sorted_list[i])
                # print (match_len_with_max_identity, identity_sorted_list[i][3], (match_len_with_max_identity - identity_sorted_list[i][3])/match_len_with_max_identity)
                if len(full_result_list) >= 30:
                    break
        
        # for allele_info in identity_sorted_list[:20]:
        #     print(allele_info, "length and identity")
            
        # # if get_help_from_1000G == False:
        # print ("check to determine use highest identity or match length in person.")

        # print ("identity **************************")

        # print ("selected allele is ", select_allele_list[0])
        # print (select_allele_list)
        # return select_allele_list
    # print (full_result_list)
    return full_result_list

def output(record_best_match, gene_list, result_file, version_info):
    f = open(result_file, 'w')
    # print ("#", version_info, file = f)
    print ("Locus   Chromosome      Allele", file = f)
    for gene in gene_list:
        for ch in [1, 2]:
            alleles = ''
            for a in record_best_match[gene][ch]:
                # print (a)
                alleles += a[0] + ";"
            out_alleles = alleles[:-1]
            print (gene, ch, out_alleles, sep="\t", file = f)


    f.close()

def get_blast_info(blastfile, tag):
    align_list = []
    for line in open(blastfile):
        field = line.strip().split()
        if field[0] != tag:
            continue
        align_info = [field[1], int(field[3]), int(field[3])-float(field[2]), 1-float(field[2])/int(field[3]), 'x', 0, 0]
        align_list.append(align_info)
    identity_sorted_list = sorted(align_list, key=get_3_element, reverse = True)
    return identity_sorted_list



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="HLA Typing from diploid assemblies.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    # required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    # required.add_argument("-1", type=str, help="Assembly file of the first haplotype in fasta formate", metavar="\b")
    # required.add_argument("-2", type=str, help="Assembly file of the second haplotype in fasta formate", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    HLA_data = "%s/../db/HLA/whole/hla_gen.fasta"%(sys.path[0])
    version_info = get_IMGT_version()

    if not os.path.exists(args["o"]):
        os.system("mkdir %s"%(args["o"]))
    if not os.path.isfile(HLA_data):
        os.system("cat %s/../db/HLA/whole/HLA_*.fasta > %s/../db/HLA/whole/hla_gen.fasta"%(sys.path[0], sys.path[0]))
    
    # inputdir = "/mnt/d/HLAPro_backup/Nanopore_optimize/output/fredhutch-hla-1347-4843/"
    inputdir = args['o']
    result_path = args['o']
    sample = args['n']

    result_file = result_path + "/" + "hlala.like.results.txt"
    all_align_result_file = result_path + "/" + "hla.map.results.txt"
    gene_list = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

    record_best_match = {}
    blastfile = f"{result_path}/hla.blast.summary.txt"

    print (f"Start optimize the typing results by balancing the alignment length and the alignment identity. ")

    for hap_index in range(2):

        for gene in gene_list:
            if gene not in record_best_match:
                record_best_match[gene] = {}
            tag = f"HLA_{gene}_{hap_index+1}"
            align_list = get_blast_info(blastfile, tag)

            full_result_list = select_by_alignment(align_list)
            record_best_match[gene][hap_index+1] = full_result_list
    
    output(record_best_match, gene_list, result_file, version_info)
    print (f"The refined typing results is in {result_file}")



