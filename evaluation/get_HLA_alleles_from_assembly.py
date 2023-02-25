"""
extract HLA allele from phased assemblies
1. map the HLA database to the assembly
2. obtain the matched length and identity of each allele
3. remove the allele does not fit the 1000G typing results
4. remove the allele that not perfect match in exons (2,3 for class I, and 2 for class II genes)
5. choose the best allele by balancing the matched length and identity
6. extract the assembly sequence that mapped to the best allele

wangshuai, Feb 20, 2023
"""

import os, re
import pysam


def get_1_element(lst):
    return lst[1]

def get_2_element(lst):
    return lst[2]

def get_3_element(lst):
    return lst[3]

def get_phased_assemblies():
    record_truth_file_dict = {}
    inpath = "/mnt/d/my_HLA/assembly/"
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
    return  record_truth_file_dict

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
    # command = f"{minimap_path} {record_truth_file_dict[sample][hap_index]} {HLA_data}  -o {result_path}/{sample}.h{hap_index+1}.paf -t 10"
    command = f"{minimap_path} {record_truth_file_dict[sample][hap_index]} {HLA_data}  -o {result_path}/{sample}.h{hap_index+1}.sam -a -t 15"
    print (command)
    os.system(command)

def minimap_exon(sample, hap_index):
    # command = f"{minimap_path} {record_truth_file_dict[sample][hap_index]} {HLA_data}  -o {result_path}/{sample}.h{hap_index+1}.paf -t 10"
    command = f"{minimap_path} {record_truth_file_dict[sample][hap_index]} {single_exon_database_fasta}  -o {result_path}/{sample}.h{hap_index+1}.exon.sam -a -t 15"
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

def recheck_fit_num(input_sam, gene, allele_perfect_exon_dict):
    have_perfect_exon_allele = False
    f = open(input_sam, "r")
    for line in f:
        # Skip header lines
        if line.startswith("@"):
            continue
        if not line.startswith(f"{gene}*"):
            continue
        allele_name = line.split("\t")[0]
        if "HLA-" + allele_name in allele_perfect_exon_dict[gene]:
            if allele_perfect_exon_dict[gene]["HLA-"+allele_name] == True:
                have_perfect_exon_allele = True
                # print ("perfect exon", allele_name)
        elif "HLA-"+allele_name[:-3] in allele_perfect_exon_dict[gene]:
            if allele_perfect_exon_dict[gene]["HLA-"+allele_name[:-3]] == True:
                have_perfect_exon_allele = True
                # print ("perfect exon", allele_name)
    f.close()
    return have_perfect_exon_allele

def ana_sam(input_sam, gene, sample, allele_perfect_exon_dict, fit_num_each_gene):
    # Open the PAF file
    align_list = []
    have_perfect_exon_allele = recheck_fit_num(input_sam, gene, allele_perfect_exon_dict) 
    # print ("<<<<<<<<<", allele_perfect_exon_dict[gene])
    ### check if the gene has allele with 100% matched exon
    ### if not, skip the criteria
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
        if have_perfect_exon_allele == True:
            flag = False #exon is not 100% matched
            if exon_allele_name in allele_perfect_exon_dict[gene]:
                if allele_perfect_exon_dict[gene][exon_allele_name]:
                    flag = True
            elif exon_allele_name[:-3] in allele_perfect_exon_dict[gene]:
                if allele_perfect_exon_dict[gene][exon_allele_name[:-3]]:
                    flag = True
            if flag == False:
                continue
        
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
    
def read_1000G_truth():
    truth_1000_dict = {}
    i = 0
    for line in open(truth_1000G_file):
        array = line.strip().split("\t")
        if i == 0 :
            header_list = array
        else:
            sample = array[2]
            truth_1000_dict[sample] = {}
            for j in range(3, 13):
                header = header_list[j]
                gene = header.split("_")[1]
                if gene not in truth_1000_dict[sample]:
                    truth_1000_dict[sample][gene] = []
                # print (sample, j, array)
                typed_allele = array[j]
                typed_allele = typed_allele.replace("*", '')

                truth_1000_dict[sample][gene].append(typed_allele)
        i += 1
    
    return truth_1000_dict

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
    out_fasta.write(f'>{sample}.h{hap_index+1}.HLA-{gene}\t{segment_name}:{start_pos}-{end_pos}\t{select_allele_list[0]}\n{sequence}\n')

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

def get_alleles_with_perfect_exon(exon_sam):
    allele_perfect_exon_dict = {}  
    contig_allele_perfect_exon_dict = {}
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")
    # Split the SAM record into fields
    
    for line in open(exon_sam, "r"):
        # Skip header lines
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        allele_name = fields[0]
        cigar = fields[5]
        contig_name = fields[2]

        perfect = True
        for length, op in re.findall(pattern, cigar):
            if op != "M":
                perfect = False
        nm_tag = [tag for tag in fields[11:] if tag.startswith("NM:i:")]
        if len(nm_tag) == 1:
            num_mismatches = int(nm_tag[0].split(":")[2])
        else:
            num_mismatches = 100000
        if num_mismatches != 0:
            perfect = False

        allele_name = allele_name.split("|")[0]
        gene = allele_name.split("*")[0].split("-")[1]
        if gene not in contig_allele_perfect_exon_dict:
            contig_allele_perfect_exon_dict[gene] = {}
        if allele_name not in contig_allele_perfect_exon_dict[gene]:
            contig_allele_perfect_exon_dict[gene][allele_name] = {}
        if contig_name not in contig_allele_perfect_exon_dict[gene][allele_name]:
            contig_allele_perfect_exon_dict[gene][allele_name][contig_name] = perfect
        else:
            contig_allele_perfect_exon_dict[gene][allele_name][contig_name] = contig_allele_perfect_exon_dict[gene][allele_name][contig_name] and perfect
    for gene in contig_allele_perfect_exon_dict:
        allele_perfect_exon_dict[gene] = {}
        for allele_name in contig_allele_perfect_exon_dict[gene]:
            flag = False
            for contig_name in contig_allele_perfect_exon_dict[gene][allele_name]:
                if  contig_allele_perfect_exon_dict[gene][allele_name][contig_name]:
                    flag = True
            allele_perfect_exon_dict[gene][allele_name] = flag

    #     if perfect and re.search("DRB1", allele_name):
    #         print (allele_name)
    # print ("test", allele_perfect_exon_dict["DRB1"]["HLA-DRB1*14:54:01:02"])
    # new_allele_perfect_exon_dict = {}
    # for gene in allele_perfect_exon_dict:  #consider the alleles with same 3-field as equal
    #     new_allele_perfect_exon_dict[gene] = {}
    #     for allele in allele_perfect_exon_dict[gene]:
    #         # print (allele)
    #         new_allele_perfect_exon_dict[gene][allele] = allele_perfect_exon_dict[gene][allele]
    #         if len(allele.split(":")) == 4:
    #             three_field = ":".join(allele.split(":")[:-1])
    #             # print (three_field)
    #             if three_field not in new_allele_perfect_exon_dict[gene]:
    #                 new_allele_perfect_exon_dict[gene][three_field] = allele_perfect_exon_dict[gene][allele]
    #             else:
    #                 new_allele_perfect_exon_dict[gene][three_field] = allele_perfect_exon_dict[gene][allele] or new_allele_perfect_exon_dict[gene][three_field]
    # print (new_allele_perfect_exon_dict["DRB1"])
    return allele_perfect_exon_dict

def count_perferct_exon_num(allele_perfect_exon_dict):
    fit_num_each_gene = {} # the number of allele with 100% matched exon
    for gene in allele_perfect_exon_dict:
        for allele in allele_perfect_exon_dict[gene]:
            if gene not in fit_num_each_gene:
                fit_num_each_gene[gene] = 0
            if allele_perfect_exon_dict[gene][allele] == True:
                # print (allele, allele_perfect_exon_dict[gene][allele])
                fit_num_each_gene[gene] += 1
    return fit_num_each_gene


class Assign_allele():

    def __init__(self, sample_save_alignments_dict, sample):
        self.sample_save_alignments_dict = sample_save_alignments_dict
        self.sample = sample

    def main(self):
        record_selection = {}
        for gene in gene_list:
            gene_alignments = self.sample_save_alignments_dict[gene]
            truth_alleles = self.get_1000G_truth(gene)
            first_hap_selection, second_hap_selection = self.handle_each_gene(gene_alignments, truth_alleles, gene)
            record_selection[gene] = [first_hap_selection, second_hap_selection]
            print (self.sample, gene, "selection", first_hap_selection[0], second_hap_selection[0])
        return record_selection
    
    def handle_each_gene(self, gene_alignments, truth_alleles, gene):
        print (self.sample, gene, "1000G", truth_alleles)
        # print (gene_alignments[0])
        # print (gene_alignments[1])
        if len(truth_alleles) > 0 and len(truth_alleles[0]) > 0:
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
        # print (">>>>>>>>>", match_sorted_list)

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
        select_allele_list = identity_sorted_list[0]
        # if extract_four_digits(match_sorted_list[0][0]) in truth_alleles and extract_four_digits(identity_sorted_list[0][0]) not in truth_alleles:
        #     select_allele_list = match_sorted_list[0]
        #     get_help_from_1000G = True
        # elif extract_four_digits(match_sorted_list[0][0]) not in truth_alleles and extract_four_digits(identity_sorted_list[0][0]) in truth_alleles:
        #     select_allele_list = identity_sorted_list[0]
        #     get_help_from_1000G = True
        # elif identiy_with_max_match_len < 0.999:
        #     select_allele_list = identity_sorted_list[0]
        # elif match_len_diff_ratio < identity_diff_ratio:
        #     select_allele_list = identity_sorted_list[0]
        # elif identity_diff_ratio < 0.005:
        #     select_allele_list = match_sorted_list[0]
        # else:
        #     print (" no determine")
            
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


    def get_1000G_truth(self, gene):
        truth_alleles = []
        if self.sample in truth_1000_dict:
            if gene in truth_1000_dict[self.sample]:
                truth_alleles = truth_1000_dict[self.sample][gene]
        return truth_alleles


if __name__ == "__main__":
    # sample = "HG00096"
    minimap_path = "minimap2"
    # https://github.com/ANHIG/IMGTHLA/blob/Latest/fasta/hla_gen.fasta
    raw_HLA_data = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/hla_gen.fasta"
    HLA_data = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/hla_gen.rename.fasta"
    raw_HLA_exon_data = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/hla_nuc.fasta"
    HLA_exon_data = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/hla_nuc.rename.fasta"
    truth_1000G_file = "/mnt/d/HLAPro_backup/wgs1000/20181129_HLA_types_full_1000_Genomes_Project_panel.txt"
    single_exon_database = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/xml/"
    single_exon_database_fasta = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/xml/hla_exons.fasta"
    # change_allele_name(raw_HLA_data, HLA_data)
    # change_allele_name(raw_HLA_exon_data, HLA_exon_data)
    get_exons_databse(single_exon_database)
    # """
    result_path = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/"
    samples_list = ['HG00096', 'HG00171', 'HG00512', 'HG00513', 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'HG00864', 'HG01114', 'HG01505', 'HG01596', 'HG02011', 'HG02492', 'HG02587', 'HG02818', 'HG03009', 'HG03065', 'HG03125', 'HG03371', 'HG03486', 'HG03683', 'HG03732', 'NA12878', 'NA18534', 'NA18939', 'NA19238', 'NA19239', 'NA19240', 'NA19650', 'NA19983', 'NA20509', 'NA20847', 'NA24385']
    trio_list = ["NA19240", "NA19239", "NA19238"]
    # trio_list = ["HG00733", "HG00731", "HG00732"]
    # trio_list = ["HG00514", "HG00512", "HG00513"]
    gene_list = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    gene_list = ["A"]
    record_truth_file_dict = get_phased_assemblies()
    truth_1000_dict = read_1000G_truth()
    # truth_1000_dict["HG03009"]["C"] = []
    # print (record_truth_file_dict.keys())
    # 


    # create an output file for the extracted segment
    out_fasta = open(result_path + "/extracted_HLA_alleles.fasta", 'w')

    record_best_match = {}
    # for sample in samples_list:
    # for sample in trio_list:
    for sample in ["HG00514"]:
        print (sample)
        sample_save_alignments_dict = {}
        # for hap_index in range(2):
            # minimap(sample, hap_index)
            # minimap_exon(sample, hap_index)
        for hap_index in range(2):
            input_sam_exon = f"/mnt/d/HLAPro_backup/minor_rev/extract_alleles/{sample}.h{hap_index+1}.exon.sam"
            allele_perfect_exon_dict = get_alleles_with_perfect_exon(input_sam_exon)
            fit_num_each_gene = count_perferct_exon_num(allele_perfect_exon_dict)
            # input_paf = f"/mnt/d/HLAPro_backup/minor_rev/extract_alleles/{sample}.h{hap_index+1}.paf"
            input_sam = f"/mnt/d/HLAPro_backup/minor_rev/extract_alleles/{sample}.h{hap_index+1}.sam"
            assembly_file = record_truth_file_dict[sample][hap_index]
            # open the input FASTA file
            
            for gene in gene_list:
                if gene not in sample_save_alignments_dict:
                    sample_save_alignments_dict[gene] = []
                align_list = ana_sam(input_sam, gene, sample, allele_perfect_exon_dict, fit_num_each_gene)
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

