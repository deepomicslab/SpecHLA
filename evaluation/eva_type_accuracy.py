"""
Calculate typing accuracy in simulated data
at the four digits

wangshuai July 11, 2022
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
import re
import pandas as pd

gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

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
            sample = sample.split("_")[0]
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

    def main_real(self): 
        data = []   
        self.sample_list = []
        all_sample_true_dict = {}
        for line in open("/mnt/d/HLAPro_backup/compare_hlala/merge.anno.out", "r"):
            array = line.strip().split()
            sample = array[0][:-3]
            if sample not in self.sample_list:
                self.sample_list.append(sample)
            gene = array[1]
            if sample not in all_sample_true_dict:
                all_sample_true_dict[sample] = {}
            array[2] = re.sub(":","_",array[2])
            array[2] = re.sub("\*","_",array[2])
            # print (array)
            if gene not in all_sample_true_dict[sample]:
                all_sample_true_dict[sample][gene] = []
            all_sample_true_dict[sample][gene] += [array[2]]


        self.hla_la_result = "/mnt/d/HLAPro_backup/compare_hlala/hifi_hlala/HLA-LA.merge.result.txt"
        hla_la_all_sample_infer_dict = self.extract_inferred(self.hla_la_result)
        data = self.assess(all_sample_true_dict, hla_la_all_sample_infer_dict, data, "HLA*LA_PacBio")
        # print (hla_la_all_sample_infer_dict)

        self.hla_la_result = "/mnt/d/HLAPro_backup/compare_hlala/ngs_hlala/HLA-LA.merge.result.txt1"
        hla_la_all_sample_infer_dict = self.extract_inferred(self.hla_la_result)
        data = self.assess(all_sample_true_dict, hla_la_all_sample_infer_dict, data, "HLA*LA_PE")

        self.spechla_outdir = "/mnt/d/HLAPro_backup/compare_hlala/pacbio/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_PacBio")

        self.spechla_outdir = "/mnt/d/HLAPro_backup/compare_hlala/spechla_no_pac/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_PE")

        self.spechla_outdir = "/mnt/d/HLAPro_backup/compare_hlala/spechla_with_pac/"
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        # print(spechla_all_sample_infer_dict)
        data = self.assess(all_sample_true_dict, spechla_all_sample_infer_dict, data, "SpecHLA_hybrid")

        # print ("HLA*LA")
        
        df = pd.DataFrame(data, columns = ["Accuracy", "Gene", "Methods"])
        df.to_csv('/mnt/d/HLAPro_backup/hybrid/hybrid_G_assess_real.csv', sep=',')

    def assess(self, all_sample_true_dict, infer_all_sample_infer_dict, data, method):
        # print (infer_all_sample_infer_dict)
        gene_count = {}
        for gene in gene_list:
            gene_count[gene] = {"right":0, "all":0}
        for sample in self.sample_list:
            for gene in all_sample_true_dict[sample]:
                true_alleles = all_sample_true_dict[sample][gene]
                infer_alleles = infer_all_sample_infer_dict[sample][gene]
                right_num = self.compare_allele(true_alleles, infer_alleles)
                gene_count[gene]["right"] += right_num
                gene_count[gene]["all"] += 2
                if right_num != 2:
                    print (sample, gene, true_alleles, infer_alleles)
        for gene in gene_list:
            print (gene, gene_count[gene]["right"], gene_count[gene]["all"], gene_count[gene]["right"]/gene_count[gene]["all"])
            accuracy = gene_count[gene]["right"]/gene_count[gene]["all"]
            data.append([accuracy, gene, method])
        return data
    
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



if __name__ == "__main__":
    # for gene in ['A', 'B', 'C', 'DQB1','DRB1']:
    #     single(gene)
    G_annotation_dict = read_G_annotation()
    digit = 6
    typ = Eva_typing()
    typ.main_real()
    
