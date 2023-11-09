"""
Extract HLA allele from phased assemblies
1. map the HLA database to the assembly
2. obtain the matched length and identity of each allele
3. choose the best allele by balancing the matched length and identity
4. extract the assembly sequence that mapped to the best allele

wangshuai, Nov 9, 2023
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

def minimap(sample, hap_index):
    command = f"{minimap_path} {record_truth_file_dict[sample][hap_index]} {HLA_data}  -o {result_path}/{sample}.h{hap_index+1}.sam -a -t {args['j']}"
    print (command)
    os.system(command)

def minimap_exon(sample, hap_index):
    # command = f"{minimap_path} {record_truth_file_dict[sample][hap_index]} {HLA_data}  -o {result_path}/{sample}.h{hap_index+1}.paf -t 10"
    command = f"{minimap_path} {record_truth_file_dict[sample][hap_index]} {single_exon_database_fasta}  -o {result_path}/{sample}.h{hap_index+1}.exon.sam -a -t {args['j']}"
    # print (command)
    os.system(command)

def ana_paf(input_paf, gene, sample):
    # Open the PAF file
    align_list = []
    paf_file =  open(input_paf, "r") 
    # Read all lines into a list
    for line in paf_file:
        if not line.startswith(f"{gene}*"):
            continue

        array = line.split("\t")
        matching_bases = int(array[9])
        Alignment_block_length = int(array[10])
        Target_sequence_name = array[5]
        Target_start_position = array[7]
        Target_end_position = array[8]
        identity = round(float(matching_bases)/Alignment_block_length, 6)
        allele = array[0]
        align_list.append([allele, matching_bases, Alignment_block_length, identity, Target_sequence_name, Target_start_position, Target_end_position])
    paf_file.close()
   
    match_sorted_list = sorted(align_list, key=get_1_element, reverse = True)
    identity_sorted_list = sorted(align_list, key=get_3_element, reverse = True)
    print (sample, gene)
    if match_sorted_list[0][0] == identity_sorted_list[0][0]:
        print ("perfect:", match_sorted_list[0])
        select_allele = match_sorted_list[0]
    
    else:
        print ("check to determine use highest identity or match length in person.")
        for allele_info in match_sorted_list[:5]:
            print(allele_info)
        print ("match bases**************************")

        
        for allele_info in identity_sorted_list[:5]:
            print(allele_info)
        print ("identity **************************")

def ana_sam(input_sam, gene, sample):
    # Open the PAF file
    align_list = []
    # # Open the SAM file
    f = open(input_sam, "r")
    for line in f:
        # Skip header lines
        if line.startswith("@"):
            continue
        if not line.startswith(f"{gene}*"):
            continue
        
        allele_name = line.split("\t")[0]
        exon_allele_name = "HLA-" + allele_name
        
        align_info = read_sam_line(line)
        align_list.append(align_info)
    identity_sorted_list = sorted(align_list, key=get_3_element, reverse = True)
    return identity_sorted_list

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
    # print (sorted_list)
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
    
def read_sam_line(line):
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")
    # Split the SAM record into fields
    fields = line.split("\t")

    # Extract the CIGAR string and sequence from the record using regular expressions
    allele_name = fields[0]
    cigar = fields[5]
    sequence = fields[9]
    match_length = 0
    block_length = 0
    target_start = int(fields[3])
    Target_sequence_name = fields[2]


    for length, op in re.findall(pattern, cigar):
        # print (length, op)
        if op == "M":
            match_length += int(length)
        if op != "S" and op != "H":
            block_length += int(length)

    nm_tag = [tag for tag in fields[11:] if tag.startswith("NM:i:")]
    if len(nm_tag) == 1:
        num_mismatches = int(nm_tag[0].split(":")[2])
    else:
        num_mismatches = 0

    # Calculate the match identity
    match_identity = round(float(match_length-num_mismatches)/block_length, 6)
    target_end = target_start + block_length
    # Print the match length and identity to the console
    # print(cigar, f"{allele_name} Match length: {match_length}, Match identity: {match_identity}", num_mismatches, block_length)
    # break
    return [allele_name, match_length, block_length, match_identity, Target_sequence_name, target_start, target_end]

def extract_seq(select_allele_list, assembly_file, hap_index, sample, gene, out_fasta, in_fasta):

    # define the segment name, start position, and end position
    segment_name = select_allele_list[4]
    start_pos = int(select_allele_list[5]) - 1
    end_pos = int(select_allele_list[6]) - 1

    # extract the sequence for the interval
    sequence = in_fasta.fetch(segment_name, start_pos, end_pos)

    # write the segment name and sequence to the output file
    out_fasta.write(f'>{sample}.h{hap_index+1}.HLA-{gene}\t{segment_name}:{start_pos}-{end_pos}\t{select_allele_list[0]}\t{version_info}\n{sequence}\n')

    # close the input and output files
    
def check_trio_consistency(record_best_match, trio_list):
    for gene in gene_list:
        child_alleles = record_best_match[trio_list[0]][gene]
        parent1_alleles = record_best_match[trio_list[1]][gene]
        parent2_alleles = record_best_match[trio_list[2]][gene]
        if (child_alleles[0] in parent1_alleles and child_alleles[1] in parent2_alleles) or (child_alleles[1] in parent1_alleles and child_alleles[0] in parent2_alleles):
            print (trio_list[0], "consistency", gene)
        else:
            print (trio_list[0], "not consistency", gene, child_alleles, parent1_alleles,  parent2_alleles)

def get_exons_databse(single_exon_database):
    out = open(single_exon_database_fasta, 'w')
    test_file = single_exon_database + "A2.exon.txt"
    for item in os.listdir(single_exon_database):
        if re.search(".exon.txt", item):
            test_file = single_exon_database + "/" + item
            # print (test_file)
            f = open(test_file)
            for line in f:
                line = line.replace('\"', '')

                array = line.split()
                # print (array)
                if len(array) == 1:
                    continue
                allele = array[0] + "|" + array[1]
                seq = array[-1].strip()
                print (f">{allele}\n{seq}", file = out)
            f.close()
    out.close()


class Assign_allele():

    def __init__(self, sample_save_alignments_dict, sample):
        self.sample_save_alignments_dict = sample_save_alignments_dict
        self.sample = sample

    def main(self):
        record_selection = {}
        for gene in gene_list:
            gene_alignments = self.sample_save_alignments_dict[gene]
            truth_alleles = [[]]
            first_hap_selection, second_hap_selection = self.handle_each_gene(gene_alignments, truth_alleles, gene)
            record_selection[gene] = [first_hap_selection, second_hap_selection]
            print (self.sample, gene, "selection", first_hap_selection[0], second_hap_selection[0])
        return record_selection
    
    def handle_each_gene(self, gene_alignments, truth_alleles, gene):
        if len(truth_alleles) > 0 and len(truth_alleles[0]) > 0:
            print (self.sample, gene, "1000G", truth_alleles)
            align_00 = self.filter_by_1000G(truth_alleles[0], gene_alignments[0])
            align_11 = self.filter_by_1000G(truth_alleles[1], gene_alignments[1])
            align_01 = self.filter_by_1000G(truth_alleles[0], gene_alignments[1])
            align_10 = self.filter_by_1000G(truth_alleles[1], gene_alignments[0])
            if len(align_00) == 0 or len(align_11) == 0:
                truth_alleles.reverse()
            elif len(align_01) == 0 or len(align_10) == 0:
                pass
            elif align_11[0][3] + align_00[0][3] < align_01[0][3] + align_10[0][3]:
                truth_alleles.reverse()
            # else:
            # print (align_00, "\n", align_11, "\n",align_01, "\n",align_10)
            my_align_00 = self.filter_by_1000G(truth_alleles[0], gene_alignments[0])
            my_align_11 = self.filter_by_1000G(truth_alleles[1], gene_alignments[1])
            # print (my_align_00, my_align_11)
            # print (truth_alleles[0], my_align_00)
            # print (truth_alleles[1], my_align_11)
            return my_align_00[0], my_align_11[0]
        
        else:
            first_hap_selection = self.select_by_alignment(gene_alignments[0], truth_alleles)
            second_hap_selection = self.select_by_alignment(gene_alignments[1], truth_alleles)
            return first_hap_selection, second_hap_selection

    def select_by_alignment(self, align_list, truth_alleles):

        match_sorted_list = sorted(align_list, key=get_1_element, reverse = True)
        match_sorted_list = resort_list_with_same_alleles(match_sorted_list, 1, 3)
        identity_sorted_list = sorted(align_list, key=get_3_element, reverse = True)
        identity_sorted_list = resort_list_with_same_alleles(identity_sorted_list, 3, 1)
        max_match_len_alleles = get_max_alleles(match_sorted_list, 1)
        max_identity_alleles = get_max_alleles(identity_sorted_list, 3)

        # print (identity_sorted_list)
        intersection_alleles = list(set(max_match_len_alleles) & set(max_identity_alleles))   
        print (">>>>>>>>>", match_sorted_list[:10])

        if len(intersection_alleles) > 0:
            select_allele_list = intersection_alleles[0].split(">")
            select_allele = select_allele_list[0]
            print (">>>>>>>>>>perfect:", select_allele)      
            return select_allele_list

        max_match_len = match_sorted_list[0][1]
        match_len_with_max_identity = identity_sorted_list[0][1]

        max_identity = identity_sorted_list[0][3]
        identiy_with_max_match_len = match_sorted_list[0][3]

        match_len_diff_ratio = (max_match_len - match_len_with_max_identity) / match_len_with_max_identity
        identity_diff_ratio = (max_identity - identiy_with_max_match_len) / identiy_with_max_match_len

        print ("match_len_diff_ratio", match_len_diff_ratio, "identity_diff_ratio", identity_diff_ratio)
        get_help_from_1000G = False
        # select_allele_list = identity_sorted_list[0]
        # if extract_four_digits(match_sorted_list[0][0]) in truth_alleles and extract_four_digits(identity_sorted_list[0][0]) not in truth_alleles:
        #     select_allele_list = match_sorted_list[0]
        #     get_help_from_1000G = True
        # elif extract_four_digits(match_sorted_list[0][0]) not in truth_alleles and extract_four_digits(identity_sorted_list[0][0]) in truth_alleles:
        #     select_allele_list = identity_sorted_list[0]
        #     get_help_from_1000G = True
        if identiy_with_max_match_len < 0.999:
            select_allele_list = identity_sorted_list[0]
        elif match_len_diff_ratio < identity_diff_ratio:
            select_allele_list = identity_sorted_list[0]
        elif identity_diff_ratio < 0.005:
            select_allele_list = match_sorted_list[0]
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

        print ("selected allele is ", select_allele_list[0])
        return select_allele_list
 
    def filter_by_1000G(self, truth, align_list):
        new_align_list = []
        if len(truth) == 5:
            for align in align_list:
                array = align[0].split("*")[1].split(":")
                two_field = array[0] + ":" + array[1]
                if two_field == truth:
                    new_align_list.append(align)
        else:
            truth = truth[:2]
            for align in align_list:
                array = align[0].split("*")[1].split(":")
                one_field = array[0]
                if one_field == truth:
                    new_align_list.append(align)
        return new_align_list

def get_IMGT_version():
    g_file = "%s/../db/HLA/hla_nom_g.txt"%(sys.path[0])
    G_annotation_dict = {}
    i = 0
    version_info = "N/A"
    for line in open(g_file):
        if re.search("# version:", line):
            version_info = line.strip()
    return version_info

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="HLA Typing from diploid assemblies.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    # required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    required.add_argument("-1", type=str, help="Assembly file of the first haplotype in fasta formate", metavar="\b")
    required.add_argument("-2", type=str, help="Assembly file of the second haplotype in fasta formate", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=10)
    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    minimap_path = "minimap2"
    HLA_data = "%s/../db/HLA/whole/hla_gen.fasta"%(sys.path[0])
    version_info = get_IMGT_version()

    # https://github.com/ANHIG/IMGTHLA/blob/Latest/fasta/hla_gen.fasta
    # raw_HLA_data = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/hla_gen.fasta"
    # HLA_data = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/hla_gen.rename.fasta"
    
    # raw_HLA_exon_data = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/hla_nuc.fasta"
    # HLA_exon_data = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/hla_nuc.rename.fasta"
    # truth_1000G_file = "/mnt/d/HLAPro_backup/wgs1000/20181129_HLA_types_full_1000_Genomes_Project_panel.txt"
    # single_exon_database = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/xml/"
    # single_exon_database_fasta = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/xml/hla_exons.fasta"
    # change_allele_name(raw_HLA_data, HLA_data)
    # change_allele_name(raw_HLA_exon_data, HLA_exon_data)
    # get_exons_databse(single_exon_database)

    if not os.path.exists(args["o"]):
        os.system("mkdir %s"%(args["o"]))
    if not os.path.isfile(HLA_data):
        os.system("cat %s/../db/HLA/whole/HLA_*.fasta > %s/../db/HLA/whole/hla_gen.fasta"%(sys.path[0], sys.path[0]))
    

    result_path = args['o']
    sample = args['n']
    samples_list = [sample]
    gene_list = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    record_truth_file_dict = {sample : [args['1'], args['2']]}

    # create an output file for the extracted segment
    out_fasta = open(result_path + f"/{sample}_extracted_HLA_alleles.fasta", 'w')

    record_best_match = {}
    for sample in samples_list:
        print (sample)
        sample_save_alignments_dict = {}
        for hap_index in range(2):
            minimap(sample, hap_index)
            # minimap_exon(sample, hap_index)
        for hap_index in range(2):
            # input_sam_exon = f"/mnt/d/HLAPro_backup/minor_rev/extract_alleles/{sample}.h{hap_index+1}.exon.sam"
            
            # input_paf = f"/mnt/d/HLAPro_backup/minor_rev/extract_alleles/{sample}.h{hap_index+1}.paf"
            input_sam = f"{result_path}/{sample}.h{hap_index+1}.sam"
            assembly_file = record_truth_file_dict[sample][hap_index]
            # open the input FASTA file
            
            for gene in gene_list:
                if gene not in sample_save_alignments_dict:
                    sample_save_alignments_dict[gene] = []
                align_list = ana_sam(input_sam, gene, sample)
                sample_save_alignments_dict[gene].append(align_list)
                # print (assembly_file, input_sam)
        ass = Assign_allele(sample_save_alignments_dict, sample)
        record_selection = ass.main()
        record_best_match[sample] = record_selection
        for hap_index in range(2):
            assembly_file = record_truth_file_dict[sample][hap_index]
            in_fasta = pysam.FastaFile(assembly_file)
            for gene in gene_list: 
                select_allele_list = record_selection[gene][hap_index]
            
                extract_seq(select_allele_list, assembly_file, hap_index, sample, gene, out_fasta, in_fasta)
    #     #         break
            in_fasta.close()
    #     # break
    out_fasta.close()
    # check_trio_consistency(record_best_match, trio_list)
    # """

