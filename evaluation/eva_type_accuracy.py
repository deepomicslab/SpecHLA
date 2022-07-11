"""
Calculate typing accuracy in simulated data
at the four digits

wangshuai July 11, 2022
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os


gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

class Eva_typing():

    def __init__(self):
        self.hla_la_result = "/mnt/d/HLAPro_backup/pacbio/HLA-LA.merge.result.txt"
        # self.spechla_outdir = "/mnt/d/HLAPro_backup/pacbio/output/"
        self.spechla_outdir = "/mnt/d/HLAPro_backup/pacbio/hybrid/"
        self.spechla_result = "/mnt/d/HLAPro_backup/pacbio/spechla.merge.result.txt"
        self.true_dir = "/mnt/d/HLAPro_backup/pacbio/simulation/"
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
                type = array[1] + ":" + array[2]
                if gene not in sample_true_dict:
                    sample_true_dict[gene] = []
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
            all_sample_infer_dict[sample] = {}
            for allele in array[1:]:
                gene = allele.split("*")[0]
                type = allele.split("*")[1][:5]
                if gene not in all_sample_infer_dict[sample]:
                    all_sample_infer_dict[sample][gene] = []
                all_sample_infer_dict[sample][gene].append(type)
        return all_sample_infer_dict

    def main(self):    
        self.get_spechla_merge_result()
        spechla_all_sample_infer_dict = self.extract_inferred(self.spechla_result)
        hla_la_all_sample_infer_dict = self.extract_inferred(self.hla_la_result)
        self.sample_list = list(hla_la_all_sample_infer_dict.keys())
        all_sample_true_dict = self.get_all_truth()
        self.assess(all_sample_true_dict, spechla_all_sample_infer_dict)
        self.assess(all_sample_true_dict, hla_la_all_sample_infer_dict)
    
    def assess(self, all_sample_true_dict, infer_all_sample_infer_dict):
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
                    print (sample,gene, true_alleles,infer_alleles )
        for gene in gene_list:
            print (gene, gene_count[gene]["right"], gene_count[gene]["all"], gene_count[gene]["right"]/gene_count[gene]["all"])
    
    def compare_allele(self, true_alleles, infer_alleles):
        right_num, test_1, test_2 = 0, 0, 0
        for i in range(2):
            if true_alleles[i] == infer_alleles[i]:
                test_1 += 1
        true_alleles = true_alleles[::-1]
        for i in range(2):
            if true_alleles[i] == infer_alleles[i]:
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

if __name__ == "__main__":
    # for gene in ['A', 'B', 'C', 'DQB1','DRB1']:
    #     single(gene)
    typ = Eva_typing()
    typ.main()
