"""
Choose best-mapped allele as reference

cmd: python select_reference.py -n fredhutch-hla-FH1 -o /mnt/d/HLAPro_backup/Nanopore_optimize/output

wangshuai, wshuai294@gmail.com
"""


import os
import numpy as np
import pickle
import sys
import argparse
from collections import defaultdict


gene_list = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# gene_list = ['A']

def run_depth(args):
    cmd = f"""
    echo search  {args["o"]}/{args["n"]}/{args["n"]}.db.sam ... 
    samtools view -bS -F 0x800  {args["o"]}/{args["n"]}/{args["n"]}.db.sam | samtools sort - >{args["o"]}/{args["n"]}/{args["n"]}.db.bam
    samtools depth -aa {args["o"]}/{args["n"]}/{args["n"]}.db.bam>{args["o"]}/{args["n"]}/{args["n"]}.db.depth
    """
    os.system(cmd)
    return """%s/%s/%s.db.depth"""%(args["o"], args["n"], args["n"])

def mapping_p(MAPQ):
    ## MAPQ: MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality is not available.
    return 1 - 10**(-(0.1+MAPQ)/10)

class Get_depth():

    def __init__(self, depth_file):
        self.depth_file = depth_file
        self.depth_dict = {}

    def record_depth(self):
        f = open(self.depth_file)
        for line in f:
            array = line.strip().split()

            allele = array[0]

            gene = allele.split("*")[0]

            depth = int(array[2])
            if gene not in self.depth_dict:
                self.depth_dict[gene] = {}
            if allele not in self.depth_dict[gene]:
                self.depth_dict[gene][allele] = []

            self.depth_dict[gene][allele].append(depth)
    
    def select(self, output):
        f = open(output, 'w')
        print ("Gene\tAllele\tDepth\tAllele_length", file = f)

        record_candidate_alleles = defaultdict(set)
        record_allele_length = {}
        for gene in self.depth_dict:
            if gene not in gene_list:
                continue
            record_allele_depth = {}
            
            record_allele_info = {}
            for allele in self.depth_dict[gene]:
                mean_depth = np.mean(self.depth_dict[gene][allele])
                median_depth = np.median(self.depth_dict[gene][allele])

                over_0 = 0
                for e in self.depth_dict[gene][allele]:
                    if e > 0:
                       over_0 += 1
                coverage =  float(over_0)/len(self.depth_dict[gene][allele])
                record_allele_length[allele] = len(self.depth_dict[gene][allele])

                record_allele_depth[allele] = mean_depth
                record_allele_info[allele] = [mean_depth, median_depth, coverage]

                # print (allele, mean_depth, median_depth, coverage)
            sorted_dict = sorted(record_allele_depth.items(), key=lambda x: x[1], reverse=True)
            for i in range(len(sorted_dict)):
                print (gene, sorted_dict[i][0], sorted_dict[i][1], record_allele_length[sorted_dict[i][0]], file = f)
        #         print (sorted_dict[i], record_allele_length[sorted_dict[i][0]], record_allele_length[sorted_dict[i][0]]*sorted_dict[i][1] )
        #         record_candidate_alleles[gene].add(sorted_dict[i][0])
        # return record_candidate_alleles
        f.close()

def map2db(args, gene):

    minimap_para = ''
    if args["y"] == "pacbio":
        minimap_para = " -x map-pb "
    elif args["y"] == "nanopore":
        minimap_para = " -x map-ont "

    outdir = args["o"] + "/" + args["n"]
    sam = outdir + "/" + args["n"] + "." + gene + ".db.sam"
    bam = outdir + "/" + args["n"] + "." + gene + ".db.bam"
    depth_file = outdir + "/" + args["n"] + "." + gene + ".db.depth"

    # map raw reads to database
    alignDB_order = f"""
    fq={args["r"] }
    ref={args["f"] }
    outdir={args["o"]}/{args["n"] }
    sample={args["n"] }

    fq=/mnt/d/HLAPro_backup/Nanopore_optimize/output/$sample/{gene}.long_read.fq.gz
    ref=/mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/db/HLA/whole/HLA_{gene}.fasta

    minimap2 -t {args["j"] } {minimap_para} -p 0.1 -N 100000 -a $ref $fq > {sam}
    samtools view -bS -F 0x800  {sam} | samtools sort - >{bam}
    samtools depth -aa {bam}>{depth_file}
    echo alignment done.
    """
    os.system(alignDB_order)
    return sam, depth_file


if __name__ == "__main__":
    # depth_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/output/fredhutch-hla-1408-1012/fredhutch-hla-1408-1012.db.depth"
    # get_depth = Get_depth(depth_file)
    # get_depth.record_depth()
    # record_candidate_alleles = get_depth.select()

    # model()
    # construct_matrix()
    # sam = "/mnt/d/HLAPro_backup/Nanopore_optimize/output/fredhutch-hla-1408-1012/fredhutch-hla-1408-1012.db.bam"
    # main(sam)

    parser = argparse.ArgumentParser(description="Sort allele by depth.", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    # required.add_argument("-f", type=str, help="IMGT reference.", metavar="\b")
    # required.add_argument("-r", type=str, help="Long-read fastq file. PacBio or Nanopore.", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="The output folder to store the typing results.", metavar="\b", default="./output")

    # optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    # optional.add_argument("-m", type=int, help="Maintain this number of alleles for ILP step.", metavar="\b", default=10)
    # optional.add_argument("-b", type=float, help="The match length increase ratio higher than this value is homo [0-1].", metavar="\b", default=0.9)

    # optional.add_argument("-g", type=int, help="Whether use G group resolution annotation [0|1].", metavar="\b", default=0)
    # optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    # optional.add_argument("-y", type=str, help="Read type, [nanopore|pacbio].", metavar="\b", default="pacbio")
    # optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    output =  """%s/%s/%s.map.txt"""%(args["o"], args["n"], args["n"])   

    depth_file = run_depth(args)
    get_depth = Get_depth(depth_file)
    get_depth.record_depth()
    record_candidate_alleles = get_depth.select(output)

    print ("result is in", output)




    
    # ref = "/mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"

    # main(args)