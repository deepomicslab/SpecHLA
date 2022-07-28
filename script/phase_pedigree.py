#!/usr/bin/env python3
from argparse import ArgumentParser
from argparse import ArgumentTypeError
import sys
import time
import re
from itertools import combinations, permutations
import os
from pysam import VariantFile
import numpy as np
import pysam

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Please give right flag (True or False).')

Usage = \
"""
python3 phase_pedigree.py [options] 

Example: python3 phase_pedigree.py --trio NA12878:NA12891:NA12892 --workdir output/

Help information can be found by python3 phase_pedigree.py -h/--help, additional information can be found in \
README.MD or https://github.com/deepomicslab/SpecHLA.
"""
scripts_dir=sys.path[0]+'/'
parser = ArgumentParser(description="SpecHLA_pedigree.",prog='python3 phase_pedigree.py',usage=Usage)
optional=parser._action_groups.pop()
required=parser.add_argument_group('required arguments')
flag_parser = parser.add_mutually_exclusive_group(required=False)
flag_data = parser.add_mutually_exclusive_group(required=False)
#necessary parameter
required.add_argument("--trio",help="The trio infromation; give sample names in the order of child:mother:father.\
 Example: NA12878:NA12891:NA12892. The order of mother and father can be ambiguous.",dest='trio',metavar='', type=str)
required.add_argument("-o", "--workdir",help="The directory consists of all samples' results.\
 The workdir running previous scipts.",dest='workdir',metavar='')

parser._action_groups.append(optional)
args = parser.parse_args()

def generate_ped_file():
    sample_list = args.trio.split(':')
    ped = '%s/%s/trio.ped'%(workdir, sample_list[0])
    f = open(ped, 'w')
    if len(sample_list) == 3:
        content = """0 %s %s %s 0 1\n0 %s 0 0 2 1\n0 %s 0 0 1 1"""%(sample_list[0], sample_list[2], sample_list[1], sample_list[1], sample_list[2])
    elif len(sample_list) == 2:
        content = """0 %s %s 0 0 1\n0 %s 0 0 0 1"""%(sample_list[0], sample_list[1], sample_list[1])
    else:
        print ('The trio info may be incorrect.')
    print (content, end='',file = f)
    f.close()
    return sample_list

def remove_cofict_allele(vcffile, newvcffile):
    #remove confict, for bcftools consensus
    in_vcf = VariantFile(vcffile)
    out = VariantFile(newvcffile,'w',header=in_vcf.header)
    sample = list(in_vcf.header.samples)[0]
    snp_list, beta_set, allele_set = [], [], []
    for record in in_vcf.fetch():
        # record.samples[sample]['AD'] = tuple(list(depth)[:2])
        # record.alts = tuple([record.alts[0]])
        # print (record.ref, record.alts, record.samples[sample]['AD'])

        
        geno = sorted(list(record.samples[sample]['GT']))
        
        single = True
        if geno == [0, 1] or geno == [1, 1]:
            index = 0
        elif geno == [0, 2] or geno == [2, 2] :
            index = 1
        elif geno == [0, 3] or geno == [3, 3]:
            index = 2
        else:
            single = False
        # print (record.pos, record.ref, record.alts, geno, index, single)
        if single:
            same = 0
            for i in range(len(record.ref)):
                if record.ref[-1-i] !=  record.alts[index][-1-i]:
                    break
                same += 1
            if same >= len(record.ref):
                same = len(record.ref) - 1
            if same > 0:
                # print (same, record.alts[0])
                record.ref = record.ref[:-same]
                alts = list(record.alts)
                alts[index] = alts[index][:-same]  
                record.alts = tuple(alts)
        elif geno == [1, 2]:
            same = 0
            for i in range(len(record.ref)):
                if record.ref[-1-i] !=  record.alts[0][-1-i] or record.ref[-1-i] !=  record.alts[1][-1-i]:
                    break
                same += 1
            if same >= len(record.ref):
                same = len(record.ref) - 1
            if same > 0:
                # print (same, record.alts[0])
                record.ref = record.ref[:-same]
                alts = list(record.alts)
                alts[0] = alts[0][:-same]  
                alts[1] = alts[1][:-same] 
                record.alts = tuple(alts)
        #     for i in range(len(record.alts[0])):
        #         if i >= len(record.alts[1]):
        #             break
        #         if record.alts[0][-1-i] !=  record.alts[1][-1-i]:
        #             break
        #         same += 1
        #     if same >= len(record.alts[0]):
        #         same = len(record.alts[0]) - 1
        #     if same >= len(record.alts[1]):
        #         same = len(record.alts[1]) - 1
        #     if same > 0:
        #         # print (same, record.alts[0])
        #         alts = list(record.alts)
        #         alts[0] = alts[0][:-same]  
        #         alts[1] = alts[1][:-same]
        #         record.alts = tuple(alts)
        # else:
        # print (vcffile, record.pos, record.ref, record.alts, geno, index, single)


        out.write(record)

    out.close()
    in_vcf.close()

if __name__ == "__main__":  
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})

    workdir = args.workdir
    sample_list = generate_ped_file()
    gene_list = ['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']
    # gene_list = ['HLA_DQB1']
    for sample in sample_list:
        outdir =  workdir + '/' + sample
        os.system('mkdir %s/trio/'%(outdir))
        for gene in gene_list:
            zipvcf = """
            tabix -f %s/%s.rephase.vcf.gz
            """%(outdir, gene)
            os.system(zipvcf)


        
    #--threshold1 0.6 --threshold2 0 
    for gene in gene_list:
        # for sample in sample_list:
        if len(sample_list) == 3:
            order = """
            bin=%s/../bin/
            workdir=%s
            gene=%s
            child=%s
            mother=%s
            father=%s
            $bin/bcftools merge $workdir/$child/$gene.rephase.vcf.gz $workdir/$mother/$gene.rephase.vcf.gz $workdir/$father/$gene.rephase.vcf.gz -o $workdir/$child/$gene.trio.merge.vcf.gz -Oz -0
            tabix -f $workdir/$child/$gene.trio.merge.vcf.gz
            python3 %s/pedhap/main.py --threshold1 0.6 --threshold2 0 -v $workdir/$child/$gene.trio.merge.vcf.gz -p $workdir/$child/trio.ped -o $workdir/$child/$gene.trio.rephase.vcf.gz
            file=$workdir/$child/$gene.trio.rephase.vcf.gz
            tabix -f $file
            for sample in `$bin/bcftools query -l $file`; do
                $bin/bcftools view -c1 -Oz -s $sample -o $workdir/$sample/trio/$sample.$gene.trio.vcf.gz $file
                tabix -f $workdir/$sample/trio/$sample.$gene.trio.vcf.gz
            done        
            """%(sys.path[0],workdir,gene, sample_list[0], sample_list[1], sample_list[2], sys.path[0])
        # elif len(sample_list) == 2:
        #     order = """
        #     bin=%s/../bin/
        #     workdir=%s
        #     gene=%s
        #     child=%s
        #     mother=%s
        #     $bin/bcftools merge $workdir/$child/$gene.specHap.phased.vcf.gz $workdir/$mother/$gene.specHap.phased.vcf.gz -o $workdir/$child/$gene.trio.merge.vcf.gz -Oz -0
        #     $bin/tabix -f $workdir/$child/$gene.trio.merge.vcf.gz
        #     $bin/python3 %s/pedhap/main.py --threshold1 0.6 --threshold2 0 -v $workdir/$child/$gene.trio.merge.vcf.gz -p $workdir/$child/trio.ped -o $workdir/$child/$gene.trio.rephase.vcf.gz
        #     file=$workdir/$child/$gene.trio.rephase.vcf.gz
        #     $bin/tabix -f $file
        #     for sample in `$bin/bcftools query -l $file`; do
        #         $bin/bcftools view -c1 -Oz -s $sample -o $workdir/$sample/trio/$sample.$gene.trio.vcf.gz $file
        #     done        
        #     """%(sys.path[0],workdir,gene, sample_list[0], sample_list[1], sys.path[0])  
        os.system(order)     
        # print (order)    




