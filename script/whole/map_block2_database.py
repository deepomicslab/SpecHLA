"""
to obtain the linkage between blocks,
get block seq,
map seq to the database
get highest score of two blocks mapped to a same allele

wangshuai July 8, 2022
"""
import sys
import os




def blast_map(fragment1):
    command = f"""
    bin={sys.path[0]}/../../bin/
    $bin/samtools faidx {hla_ref} {fragment1}  | $bin/bcftools consensus -H 1 {vcf} >{outdir}/{fragment1}_hap1.fasta
    $bin/samtools faidx {hla_ref} {fragment1}  | $bin/bcftools consensus -H 2 {vcf} >{outdir}/{fragment1}_hap2.fasta

    $bin/blastn -query {outdir}/{fragment1}_hap1.fasta -out {outdir}/{fragment1}_hap1.fasta.out -subject {gene_db} -outfmt 6 -max_target_seqs 10000 -strand plus
    $bin/blastn -query {outdir}/{fragment1}_hap2.fasta -out {outdir}/{fragment1}_hap2.fasta.out -subject {gene_db} -outfmt 6 -max_target_seqs 10000 -strand plus

    """
    os.system(command)

class Construct_Graph():
    # get the linkage of blocks from database

    def __init__(self):
        self.fragments = []
        self.noise = 150
        self.break_point_list = [1000]
        self.link_type = [[[1,1],[2,2]],  [[2,1],[1,2]] ]  
        self.dup_start = 3950 #3898
        self.dup_end = 4300 #4400

    def split_fragments(self):
        f = open(break_point_file)
        record_breakpoint_num = 0
        for line in f:
            if line[0] == "#":
                continue
            record_breakpoint_num += 1
            array = line.strip().split()
            # break_point = int(array[1])
            break_point = round((int(array[1]) + int(array[2]))/2)
            if break_point > self.dup_start - self.noise and break_point < self.dup_end + self.noise:
                continue
            self.break_point_list.append(break_point)
        self.break_point_list.append(gene_length)
        if gene == "HLA_DRB1":
            self.break_point_list += [self.dup_start, self.dup_end]
        self.break_point_list = sorted(self.break_point_list)
        self.save_fragments()

        if record_breakpoint_num == 0: # no breakpoint, thus no fragment splitted
            self.fragments = []
    
    def save_fragments(self):
        for i in range(len(self.break_point_list) - 1):
            start = self.break_point_list[i]+1
            end = self.break_point_list[i+1]
            fragment = "%s:%s-%s"%(gene, str(start), str(end))
            if fragment == "HLA_DRB1:%s-%s"%(self.dup_start+1, self.dup_end):
                continue
            self.fragments.append(fragment)
        # print (self.fragments)
    
    def get_edge(self):
        score_out = open(score_file, "w")
        print ("#frag1 frag2 00_edge_score;01_edge_score  00_allele;frag1_map_score;frag1_map_len;frag2_map_score;frag2_map_len \
        01_allele;frag1_map_score;frag1_map_len;frag2_map_score;frag2_map_len", file = score_out)
        if len(self.fragments) <= 1:
            print ("No need to phase block.")
            return 0
        for i in range(len(self.fragments)):
            blast_map(self.fragments[i])

        for i in range(len(self.fragments)):
            for j in range(i+1, len(self.fragments)):
                fragment1 = self.fragments[i]
                fragment2 = self.fragments[j]
                # for two fragments, get the weight of their linkage from blast file
                # the weight is the edge weight in the graph
                egde_info = []

                # choose a higher score from 00 and 11, or 01 and 10 as the edge weight
                # for x in range(2):
                #     edge_score = 0
                #     edge_link = None
                #     for y in range(2):
                #         # print (self.link_type[x][y])
                #         blast_file_1 = f"{outdir}/{fragment1}_hap{self.link_type[x][y][0]}.fasta.out"
                #         blast_file_2 = f"{outdir}/{fragment2}_hap{self.link_type[x][y][1]}.fasta.out"
                #         analyze = Analyze_map()
                #         link = analyze.main(blast_file_1, blast_file_2)
                #         if link.high_score >= edge_score:
                #             edge_link = link
                #             edge_score = link.high_score
                #     egde_info += [edge_score, edge_link.support_allele]

                # sum the score of 00 and 11, or 01 and 10 as the edge weight
                egde_info = [0, [], 0, []]
                for x in range(2):
                    for y in range(2):
                        # print (self.link_type[x][y])
                        blast_file_1 = f"{outdir}/{fragment1}_hap{self.link_type[x][y][0]}.fasta.out"
                        blast_file_2 = f"{outdir}/{fragment2}_hap{self.link_type[x][y][1]}.fasta.out"
                        analyze = Analyze_map()
                        link = analyze.main(blast_file_1, blast_file_2)
                        if x == y:
                            egde_info[0] += link.high_score
                            egde_info[1] += link.support_allele
                        else:
                            egde_info[2] += link.high_score
                            egde_info[3] += link.support_allele        
                
                support_1 = ""
                for allele in egde_info[1]:  
                    for ele in allele:
                        support_1 = support_1 + str(ele) + ";" 
                support_2 = ""
                for allele in egde_info[3]:  
                    for ele in allele:
                        support_2 = support_2 + str(ele) + ";" 
                print (fragment1, fragment2, "%s;%s"%(egde_info[0], egde_info[2]), support_1, support_2, file = score_out)
        score_out.close()
   
class Analyze_map():
    # def __init__(self, gene):
    #     self.gene = gene
    def read_blast(self, blast_file): # in outfmt 6
        score_dict = {}
        f = open(blast_file)
        for line in f:
            array = line.strip().split()
            allele = array[1]
            score = round(float(array[2]),2)
            map_len = int(array[3])
            score_dict[allele] = [score, map_len]
        return score_dict
    
    def merge_score(self, score_dict_1, score_dict_2):
        merged_score_dict = {}
        for allele in score_dict_1:
            if allele in score_dict_2:
                merged_score_dict[allele] = score_dict_1[allele] + score_dict_2[allele]
        return merged_score_dict
    
    def main(self, blast_file_1, blast_file_2):
        # blast_file_1 = "/mnt/d/HLAPro_backup/haplotype/sample20/frag1_1.out"
        # blast_file_2 = "/mnt/d/HLAPro_backup/haplotype/sample20/frag1_1.out"
        score_dict_1 = self.read_blast(blast_file_1)
        score_dict_2 = self.read_blast(blast_file_2)
        merged_score_dict = self.merge_score(score_dict_1, score_dict_2)
        link = Linkage(merged_score_dict)
        # print (link.high_score, link.support_allele)
        return link

class Linkage():
    # get the edge score, and the support allele
    def __init__(self, merged_score_dict):
        high_score = 0
        support_allele = []
        if len(merged_score_dict) > 0:
            i = 0
            for allele in merged_score_dict:
                score_list = merged_score_dict[allele]
                weight_score = score_list[0] * score_list[1] + score_list[2] * score_list[3]
                # score is frag_1_mapped_length*frag_1_identity + frag_2_mapped_length*frag_2_identity
                if i == 0:
                    high_score = weight_score
                    support_allele.append([allele] + score_list )
                else:
                    if weight_score == high_score:
                        support_allele.append([allele] + score_list )
                # if i < 5:
                #     print (allele, merged_score_dict[allele], weight_score)

                i += 1
        self.high_score = high_score
        self.support_allele = support_allele
        # print (support_allele, high_score)


# analyze = Analyze_map()
# analyze.main()
# map = Map_database()
# map.blast_map("HLA_DRB1:1001-3950", "HLA_DRB1:4300-6378")
if __name__ == "__main__":

    # gene = "HLA_DRB1"
    # outdir = "/mnt/d/HLAPro_backup/haplotype/sample20/"
    gene = sys.argv[1]
    outdir = sys.argv[2] + "/"

    hla_ref = '%s/../../db/ref/hla.ref.extend.fa'%(sys.path[0])
    gene_db = '%s/../../db/HLA/whole/%s.fasta'%(sys.path[0], gene)
    vcf = "%s/%s.spechap.vcf.gz"%(outdir, gene)

    break_point_file = outdir + "/%s_break_points_spechap.txt"%(gene)
    score_file = outdir + "/%s_break_points_score.txt"%(gene)

    gene_length_dict = {'HLA_A':[1000,4503],'HLA_B':[1000,5081],'HLA_C':[1000,5304],'HLA_DPA1':[1000,10775],\
        'HLA_DPB1':[1000,12468],'HLA_DQA1':[1000,7492],'HLA_DQB1':[1000,8480],'HLA_DRB1':[1000,12229]}

    gene_length = gene_length_dict[gene][1]


    cons = Construct_Graph()
    cons.split_fragments()
    cons.get_edge()
