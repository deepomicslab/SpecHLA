"""
Calculate typing accuracy of HLA*LA and SpecHLA at the G group level in HGSVC2 samples.
In the samples with both PE and PacBio data.

wangshuai July 11, 2022
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
import re
import pandas as pd
import numpy as np
from collections import Counter

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['DRB1']

class Eva_typing():

    def __init__(self):
        self.hla_la_result = "/mnt/d/HLAPro_backup/hybrid/HLA-LA.merge.result.pacbio.txt"
        # self.hla_la_result = "/mnt/d/HLAPro_backup/hybrid/HLA-LA.merge.result.ngs.txt"

        # self.spechla_outdir = "/mnt/d/HLAPro_backup/hybrid/test_illumina/"
        # self.spechla_outdir = "/mnt/d/HLAPro_backup/hybrid/illumina/"
        self.spechla_outdir = "/mnt/d/HLAPro_backup/hybrid/pacbio/"
        # self.spechla_outdir = "/mnt/d/HLAPro_backup/hybrid/pacbio_illumina/"
        self.spechla_result = "/mnt/d/HLAPro_backup/hybrid/spechla.merge.result.txt"
        self.true_dir = "/mnt/d/HLAPro_backup/hybrid/data/"
        self.sample_list = []

    def get_spechla_merge_result(self):
        command = f"cat {self.spechla_outdir}/*/hla.result.txt|grep -v Sample >{self.spechla_result}"
        os.system(command)

    def extract_truth(self, simulated_fasta):
        sample_true_dict = {}
        with open(simulated_fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                array = record.id.split("_")
                gene = array[0]
                # if digit == 6 and len(array) > 3:
                #     # print (array, len(array))
                #     type = array[1] + ":" + array[2] + ":" + array[3]
                # else:
                #     type = array[1] + ":" + array[2]
                if gene not in sample_true_dict:
                    sample_true_dict[gene] = []
                type=G_annotation_dict[record.id]
                sample_true_dict[gene].append(type)
        return sample_true_dict
    
    def get_all_truth(self):
        all_sample_true_dict = {}
        for sample in self.sample_list:
            simulated_fasta = self.true_dir + sample + ".fasta"
            sample_true_dict = self.extract_truth(simulated_fasta)
            all_sample_true_dict[sample] = sample_true_dict
        return all_sample_true_dict

    def extract_inferred(self, inferred_result):
        all_sample_infer_dict = {}
        f = open(inferred_result, 'r')
        for line in f:
            array = line.strip().split()
            if array[0] == "Sample":
                continue
            sample = array[0]
            # sample = sample.split("_")[0]
            all_sample_infer_dict[sample] = {}
            for allele in array[1:]:
                gene = allele.split("*")[0]
                # if digit == 4:
                #     type = allele.split("*")[1][:5]
                # elif digit == 6:
                #     type = allele.split("*")[1][:8]
                if gene not in all_sample_infer_dict[sample]:
                    all_sample_infer_dict[sample][gene] = []
                allele = re.sub(":","_",allele)
                allele = re.sub("\*","_",allele)
                allele = allele.split(";")[0]
                if allele in G_annotation_dict:
                    allele = G_annotation_dict[allele]

                elif allele + "_01" in G_annotation_dict:
                    allele = G_annotation_dict[allele + "_01"]
                elif allele + "_01_01" in G_annotation_dict:
                    allele = G_annotation_dict[allele + "_01_01"]
                # else:
                #     print ("no G-group resolution:", allele)
                type = allele
                all_sample_infer_dict[sample][gene].append(type)
        return all_sample_infer_dict

    def main(self): 
        data = []   
        self.sample_list = []
        for i in range(5):
            sample="hybrid_%s"%(i)
            self.sample_list.append(sample)
        all_sample_true_dict = self.get_all_truth()

        # self.hla_la_result = "/mnt/d/HLAPro_backup/hybrid/pacbio.HLA-LA.merge.result.all.txt"
        # hla_la_all_sample_infer_dict = self.extract_inferred(self.hla_la_result)
        # data = self.assess(all_sample_true_dict, hla_la_all_sample_infer_dict, data, "HLA*LA_PB")

        # self.hla_la_result = "/mnt/d/HLAPro_backup/hybrid/ngs.HLA-LA.merge.result.all.txt"
        # hla_la_all_sample_infer_dict = self.extract_inferred(self.hla_la_result)
        # data = self.assess(all_sample_true_dict, hla_la_all_sample_infer_dict, data, "HLA*LA_PE")

        self.spechla_outdir = "/mnt/d/HLAPro_backup/hybrid/illumina/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_PE")

        self.spechla_outdir = "/mnt/d/HLAPro_backup/hybrid/pacbio/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_PB")

        self.spechla_outdir = "/mnt/d/HLAPro_backup/hybrid/pacbio_illumina/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_hybrid")

        # print ("HLA*LA")
        
        df = pd.DataFrame(data, columns = ["Accuracy", "Gene", "Methods"])
        df.to_csv('/mnt/d/HLAPro_backup/hybrid/hybrid_G_assess.csv', sep=',')

    def get_truth_allele(self):
        all_sample_true_dict = {}
        
        for line in open(truth_g_group):
            line = line.strip()
            array = line.split(",")
            if array[0] == "sample":
                continue
            sample = array[0]
            gene = array[1]
            if sample not in all_sample_true_dict:
                all_sample_true_dict[sample] = {}
            all_sample_true_dict[sample][gene] = [array[3], array[5]]

        return all_sample_true_dict

    def check_spechla_accuracy(self):
        # check results in all 32 samples
        spechla_dir = "/mnt/d/HLAPro_backup/haplotype_v2/spechla/"
        spechla_result = "/mnt/d/HLAPro_backup/haplotype_v2/spechla_merge_results.txt"
        command = f"cat {spechla_dir}/*/hla.result.g.group.txt|grep -v Sample >{spechla_result}"
        os.system(command)
        all_sample_true_dict = self.get_truth_allele()
        data  = []
        spechla_32_sample_infer_dict = self.extract_inferred(spechla_result)
        self.sample_list = []
        for sample in all_sample_true_dict.keys():
            if sample in spechla_32_sample_infer_dict:
                self.sample_list.append(sample)
        print (len(spechla_32_sample_infer_dict), len(all_sample_true_dict))
        self.assess(all_sample_true_dict, spechla_32_sample_infer_dict, data, "spechla in 32 samples")

    def check_new_annotation_spechla_accuracy(self):
        # check results in all 32 samples
        spechla_dir = "/mnt/d/HLAPro_backup/haplotype_v2/spechla/"

        all_sample_true_dict = self.get_truth_allele()
        self.sample_list = []
        for sample in all_sample_true_dict.keys():
            if os.path.isdir(f"{spechla_dir}/{sample}"):
                self.sample_list.append(sample)
        spechla_32_sample_infer_dict = {}
        for sample in self.sample_list:
            # print (sample)
            g_ann = G_annotation(sample, spechla_dir)
            sample_result = g_ann.main()
            spechla_32_sample_infer_dict[sample] = sample_result
        # print (spechla_32_sample_infer_dict)
        print (all_sample_true_dict)
        self.assess(all_sample_true_dict, spechla_32_sample_infer_dict, [], "spechla with new G annotation")

    def main_real(self): 
        data = []   
        type_data = []
        self.sample_list = ["NA12878", 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240', 'NA12878']
        # self.sample_list = ["NA12878", 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240', 'NA12878_nanopore', 'NA12878_pacbio']
        all_sample_true_dict = self.get_truth_allele()
        new_all_sample_true_dict = {}
        for key in all_sample_true_dict:
            if key in self.sample_list:
                new_all_sample_true_dict[key] = all_sample_true_dict[key] 
        all_sample_true_dict = new_all_sample_true_dict

        all_sample_true_dict["NA12878"] = {'A': ['A_01_01_01G', 'A_11_01_01G'], 'B': ['B_56_01_01G', 'B_08_01_01G'], 'C': ['C_01_02_01G', 'C_07_01_01G'], \
        'DPA1': ['DPA1_01_03_01G', 'DPA1_02_01_01G'], 'DPB1': ['DPB1_04_01_01G', 'DPB1_14_01_01G'], 'DQA1': ['DQA1_01_01_01G', 'DQA1_05_01_01G'], \
        'DQB1': ['DQB1_02_01_01G', 'DQB1_05_01_01G'], 'DRB1': ['DRB1_01_01_01G', 'DRB1_03_01_01G']}
        # print (all_sample_true_dict["HG00514"])

        self.hla_la_result = "/mnt/d/HLAPro_backup/compare_hlala/hifi_hlala/HLA-LA.merge.result.txt"
        # self.hla_la_result = "/mnt/d/HLAPro_backup/compare_hlala/hifi_hlala/original_bam/HLA-LA.merge.result.txt"
        hla_la_all_sample_infer_dict = self.extract_inferred(self.hla_la_result)
        # print (hla_la_all_sample_infer_dict)
        data, type_data = self.assess(all_sample_true_dict, hla_la_all_sample_infer_dict, data, "HLA*LA_PacBio", type_data)
        # # print (hla_la_all_sample_infer_dict)

        self.hla_la_result = "/mnt/d/HLAPro_backup/compare_hlala/ngs_hlala/HLA-LA.merge.result.txt"
        hla_la_all_sample_infer_dict = self.extract_inferred(self.hla_la_result)
        data, type_data = self.assess(all_sample_true_dict, hla_la_all_sample_infer_dict, data, "HLA*LA_PE", type_data)

        self.spechla_outdir = "/mnt/d/HLAPro_backup/compare_hlala/pacbio/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        data, type_data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_PacBio", type_data)
        
        # # print ("<<<", spechla_all_sample_infer_dict.keys())
        self.spechla_outdir = "/mnt/d/HLAPro_backup/compare_hlala/spechla_no_pac/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        data, type_data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_PE", type_data)

        self.spechla_outdir = "/mnt/d/HLAPro_backup/compare_hlala/spechla_with_pac/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        # print(spechla_all_sample_infer_dict)
        data, type_data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_hybrid", type_data)

        # print ("HLA*LA")
        
        df = pd.DataFrame(data, columns = ["Accuracy", "Gene", "Methods"])
        df.to_csv('/mnt/d/HLAPro_backup/hybrid/hybrid_G_assess_real.csv', sep=',')
        df = pd.DataFrame(type_data, columns = ["Methods", "Sample", "Gene", "Truth_allele_1", "Truth_allele_2", "Typed_allele_1", "Typed_allele_2"])
        df.to_csv('/mnt/d/HLAPro_backup/hybrid/hybrid_type_results.csv', sep=',')
    
    
    # def save_type_result(self):

    def assess(self, all_sample_true_dict, infer_all_sample_infer_dict, data, method, type_data):
        # print (infer_all_sample_infer_dict)
        gene_count = {}
        for gene in gene_list:
            gene_count[gene] = {"right":0, "all":0}
        for sample in infer_all_sample_infer_dict.keys():
            pure_sample = sample.split("_")[0]
            if pure_sample not in all_sample_true_dict:
                continue
            for gene in all_sample_true_dict[pure_sample]:
                true_alleles = all_sample_true_dict[pure_sample][gene]
                
                if gene in infer_all_sample_infer_dict[sample]:
                    infer_alleles = infer_all_sample_infer_dict[sample][gene]
                else:
                    infer_alleles = ["no result", "no result"]
                # print (sample, gene, "true:", true_alleles, infer_alleles)
                type_data.append([method, sample, gene] + true_alleles + infer_alleles)
                right_num = self.compare_allele(true_alleles, infer_alleles)
                gene_count[gene]["right"] += right_num
                gene_count[gene]["all"] += 2
                if right_num != 2:
                    print (sample, gene, true_alleles, infer_alleles)
        accuracy_list = []
        total_num = 0
        right_num = 0
        for gene in gene_list:
            print (method, gene, gene_count[gene]["right"], gene_count[gene]["all"], gene_count[gene]["right"]/gene_count[gene]["all"])
            accuracy = gene_count[gene]["right"]/gene_count[gene]["all"]
            accuracy_list.append(accuracy)
            data.append([accuracy, gene, method])
            total_num += gene_count[gene]["all"]
            right_num += gene_count[gene]["right"]
        print (right_num/total_num, total_num, right_num)
        print ("\n")
        return data, type_data
    
    def compare_allele(self, my_true_alleles, my_infer_alleles):
        true_alleles = my_true_alleles.copy()
        
        if len(my_infer_alleles) == 1:
            my_infer_alleles = my_infer_alleles + my_infer_alleles
        infer_alleles = my_infer_alleles.copy()
        # print (my_infer_alleles)
        for i in range(2):
            true_alleles[i] = re.sub("G","",my_true_alleles[i])
            infer_alleles[i] = re.sub("G","",my_infer_alleles[i])
        right_num, test_1, test_2 = 0, 0, 0
        for i in range(2):
            if true_alleles[i] == infer_alleles[i]:# or re.search(infer_alleles[i], true_alleles[i]):
                test_1 += 1
        true_alleles = true_alleles[::-1]
        for i in range(2):
            if true_alleles[i] == infer_alleles[i]:# or re.search(infer_alleles[i], true_alleles[i]):
                test_2 += 1
        right_num = max(test_1, test_2)
        return right_num

    def single(self, gene):
        error_num = 0
        allele_num =0
        for line in open('compare4.log','r'):
            line = line.strip()
            array = line.split()
            if array[1] != gene:
                continue
            truth = array[4:]
            allele_num += 2
            origin_num = error_num
            if array[2] in truth:
                truth.remove(array[2])
            else:
                error_num += 1
            if array[3] not in truth:
                error_num += 1
            if error_num > origin_num:
                print (line, '\t',  error_num)
        error_rate = error_num/allele_num
        print (gene, allele_num, 1 - error_rate)

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

class G_annotation():
    def __init__(self, sample, spechla_dir):
        self.sample = sample
        self.spechla_dir = spechla_dir
    
    def blast(self, infer_hap_file, blast_result_file):
        command = f"""
        blastn -query {infer_hap_file} -out {blast_result_file} -subject {exon_database} -outfmt 7 -max_target_seqs 60000
        """
        # os.system(command)  
           
    def read_blast(self, blast_result_file):
        # blast_result_file = "/mnt/d/HLAPro_backup/haplotype_v2/spechla/HG00096/test.out"
        f = open(blast_result_file, 'r')
        identity_record = {}
        record_hit_exon_times = {}
        for line in f:
            if line[0] == "#":
                continue
            array = line.strip().split()
            exon = array[1]
            identity = float(array[2])
            match_len = int(array[3])
            allele = exon.split("|")[0]
            if not self.check_pop(allele):
                continue

            if allele not in identity_record:
                identity_record[allele] = [0, 0]
            if exon not in record_hit_exon_times:
                identity_record[allele][0] += identity
                identity_record[allele][1] += match_len
            record_hit_exon_times[exon] = 1
        f.close()

        if len(identity_record) == 0:
            return "no_match"

        # sort the dictionary by its values in ascending order
        # sorted_identity_record = sorted(identity_record.items(), key=lambda x: x[1], reverse = True)
        sorted_identity_record = sorted(identity_record.items(), key=lambda x: (x[1][0], x[1][1]), reverse = True)
        # print (sorted_identity_record[0])
        # print the sorted dictionary
        top_alleles = []
        max_identity = sorted_identity_record[0][1][0]
        max_length = sorted_identity_record[0][1][1]
        for allele, info in sorted_identity_record:
            if info[0] == max_identity and info[1] == max_length:
                top_alleles.append(convert_G(allele))
                # print(allele, info, convert_G(allele))
            else:
                break
        most_common_allele = most_common(top_alleles)
        return (most_common_allele)

    def main(self):
        sample_results = {}
        for gene in gene_list:
            sample_results[gene] = []
            for hap_index in range(1,3):
                infer_hap_file = f"{self.spechla_dir}/{self.sample}/hla.allele.{hap_index}.HLA_{gene}.fasta"
                blast_result_file = infer_hap_file + ".exon.blast"
                self.blast(infer_hap_file, blast_result_file)
                g_group_type = self.read_blast(blast_result_file)
                sample_results[gene].append(g_group_type)
        print (self.sample, sample_results)
        return sample_results

    def check_pop(self, allele):
        allele = re.sub("HLA-","",allele)
        array = allele.split(":")
        two_field = array[0] + ":" + array[1]
        flag = False
        if two_field in  hashp:
            if hashp[two_field] > 0:
                flag = True
        return flag

def population(pop, wxs):
    hashp = {}
    freq = "/home/wangshuai/softwares/SpecHLA/db/HLA/HLA_FREQ_HLA_I_II.txt"
    with open(f"{freq}", "r") as fin:
        next(fin)
        for line in fin:
            gene, c, b, a = line.strip().split()
            if wxs == "exon":
                a = "%.3f" % float(a)
                b = "%.3f" % float(b)
                c = "%.3f" % float(c)
            elif wxs == "whole":
                a = "%.8f" % float(a)
                b = "%.8f" % float(b)
                c = "%.8f" % float(c)
            if pop == "Unknown":
                hashp[gene] = (float(a) + float(b) + float(c)) / 3
            elif pop == "Asian":
                hashp[gene] = float(a)
            elif pop == "Black":
                hashp[gene] = float(b)
            elif pop == "Caucasian":
                hashp[gene] = float(c)
            elif pop == "nonuse":
                hashp[gene] = 0
    return hashp

def most_common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]

class All_HGSCV2(Eva_typing):

    def get_spechla_merge_result(self):
        command = f"cat {self.spechla_outdir}/*/hla.result.g.group.txt|grep -v Sample >{self.spechla_result}"
        os.system(command)

    def main_real(self): 
        data = []   
        type_data = []
        self.sample_list = ["NA12878", 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240', 'NA12878']
        # self.sample_list = ["NA12878", 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240', 'NA12878_nanopore', 'NA12878_pacbio']
        all_sample_true_dict = self.get_truth_allele()
        # new_all_sample_true_dict = {}
        # for key in all_sample_true_dict:
        #     if key in self.sample_list:
        #         new_all_sample_true_dict[key] = all_sample_true_dict[key] 
        # all_sample_true_dict = new_all_sample_true_dict

        all_sample_true_dict["NA12878"] = {'A': ['A_01_01_01G', 'A_11_01_01G'], 'B': ['B_56_01_01G', 'B_08_01_01G'], 'C': ['C_01_02_01G', 'C_07_01_01G'], \
        'DPA1': ['DPA1_01_03_01G', 'DPA1_02_01_01G'], 'DPB1': ['DPB1_04_01_01G', 'DPB1_14_01_01G'], 'DQA1': ['DQA1_01_01_01G', 'DQA1_05_01_01G'], \
        'DQB1': ['DQB1_02_01_01G', 'DQB1_05_01_01G'], 'DRB1': ['DRB1_01_01_01G', 'DRB1_03_01_01G']}
        # print (all_sample_true_dict["HG00514"])


        
        # # print ("<<<", spechla_all_sample_infer_dict.keys())
        self.spechla_outdir = "/mnt/d/HLAPro_backup/haplotype_v2/spechla/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        data, type_data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA", type_data)

        
        # df = pd.DataFrame(data, columns = ["Accuracy", "Gene", "Methods"])
        # df.to_csv('/mnt/d/HLAPro_backup/hybrid/hybrid_G_assess_real.csv', sep=',')
        # df = pd.DataFrame(type_data, columns = ["Methods", "Sample", "Gene", "Truth_allele_1", "Truth_allele_2", "Typed_allele_1", "Typed_allele_2"])
        # df.to_csv('/mnt/d/HLAPro_backup/hybrid/hybrid_type_results.csv', sep=',')





if __name__ == "__main__":
    # for gene in ['A', 'B', 'C', 'DQB1','DRB1']:
    #     single(gene)
    # matched_allele_file = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/extracted_HLA_alleles.fasta"
    truth_g_group = "/mnt/d/HLAPro_backup/haplotype_v2/hgsvc2_G_group_truth.csv"
    exon_database = "/mnt/d/HLAPro_backup/minor_rev/extract_alleles/xml/hla_exons.fasta"
    G_annotation_dict = read_G_annotation()
    hashp = population("Unknown", "whole")
    digit = 6
    # typ = Eva_typing()
    # typ.main_real()

    typ = All_HGSCV2()
    typ.main_real()


    # typ.print_truth() # save g-group truth in a table
    # typ.check_spechla_accuracy()
    # typ.check_new_annotation_spechla_accuracy()
    # g_ann = G_annotation("NA19239", "/mnt/d/HLAPro_backup/haplotype_v2/spechla/")
    # g_ann.main()
    
