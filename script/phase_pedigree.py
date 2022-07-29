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

def add_PS_tag(raw, new):
    # raw = "/mnt/d/HLAPro_backup/trio/trio_1000/spechla/HG005/test.vcf"
    # new = "/mnt/d/HLAPro_backup/trio/trio_1000/spechla/HG005/test.new.vcf"
    unzip_tmp = raw + ".1"
    unzip_tmp_ps = raw + ".2"
    os.system("zcat %s >%s"%(raw, unzip_tmp))
    out_f = open(unzip_tmp_ps, 'w')
    for line in open(unzip_tmp):
        line = line.strip()
        if line[0] == "#":
            if re.search("INFO=<ID=DP", line):
                line = line.replace("Number=1", "Number=A")
            if re.search("FORMAT=<ID=AD", line):
                line = line.replace("Number=R", "Number=.")
            if re.search("FORMAT=<ID=AO", line):
                line = line.replace("Number=A", "Number=.")
            if re.search("FORMAT=<ID=QA", line):
                line = line.replace("Number=A", "Number=.")
            if re.search("INFO=<ID=CIGAR", line):
                line = line.replace("Number=A", "Number=.")
            # line= re.sub("Number=(.),", "", line)
            
            # line = line.replace("Type=Integer", "Type=String")
            # line = line.replace("Type=Float", "Type=String")
            # print (line)
            print (line, file = out_f)
            
        if re.search("##FORMAT=<ID=GT", line):
            print ("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing Block No.\">", file = out_f)
        if line[0] != "#":
            array = line.split()
            array[-2] += ":PS"
            array[-1] += ":1"
            new_line = "\t".join(array)
            print (new_line, file = out_f)
    out_f.close()
    os.system("%s -kg %s >%s"%(vcf_split_tool, unzip_tmp_ps, new))
    in_vcf = VariantFile(new)
    out = VariantFile(new+".gz",'w',header=in_vcf.header)
    sample = list(in_vcf.header.samples)[0]
    for record in in_vcf.fetch():
        # record.samples[sample]['PS'] = 1
        # record.info[] [('PS', 1)]
        # record.samples[sample].info['PS'] = 1
        record.info['CIGAR'] = "2X"
        # print (record.info['CIGAR'])
        out.write(record)
    out.close()

    os.system("tabix -f %s.gz"%(new))

if __name__ == "__main__":  
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})

    workdir = args.workdir
    sample_list = generate_ped_file()
    gene_list = ['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']
    # gene_list = ['HLA_DRB1']

    vcf_split_tool = "%s/../bin/vcfallelicprimitives"%(sys.path[0])
    
    for sample in sample_list:
        outdir =  workdir + '/' + sample
        os.system('mkdir %s/trio/'%(outdir))
        
        for gene in gene_list:
            mnp_vcf = "%s/%s.rephase.vcf.gz"%(outdir, gene)
            
            snp_vcf = "%s/%s.rephase.ps.vcf"%(outdir, gene)
            add_PS_tag(mnp_vcf, snp_vcf)

    #--threshold1 0.6 --threshold2 0 
    for gene in gene_list:
        print (gene)
        if len(sample_list) == 3:
            order = """
            bin=%s/../bin/
            workdir=%s
            gene=%s
            child=%s
            mother=%s
            father=%s
            $bin/bcftools merge $workdir/$child/$gene.rephase.ps.vcf.gz $workdir/$mother/$gene.rephase.ps.vcf.gz $workdir/$father/$gene.rephase.ps.vcf.gz -o $workdir/$child/$gene.trio.merge.vcf.gz -Oz -0
            tabix -f $workdir/$child/$gene.trio.merge.vcf.gz
            python3 %s/pedhap/main.py --threshold1 0.6 --threshold2 0 -v $workdir/$child/$gene.trio.merge.vcf.gz -p $workdir/$child/trio.ped -o $workdir/$child/$gene.trio.rephase.vcf.gz
            file=$workdir/$child/$gene.trio.rephase.vcf.gz
            tabix -f $file
            for sample in `$bin/bcftools query -l $file`; do
                $bin/bcftools view -c1 -Oz -s $sample -o $workdir/$sample/trio/$sample.$gene.trio.vcf.gz $file
                tabix -f $workdir/$sample/trio/$sample.$gene.trio.vcf.gz
                $bin/samtools faidx $bin/../db/ref/hla.ref.extend.fa $gene |$bin/bcftools consensus -H 1 $workdir/$sample/trio/$sample.$gene.trio.vcf.gz>$workdir/$sample/trio/hla.allele.1.$gene.fasta
                $bin/samtools faidx $bin/../db/ref/hla.ref.extend.fa $gene |$bin/bcftools consensus -H 2 $workdir/$sample/trio/$sample.$gene.trio.vcf.gz>$workdir/$sample/trio/hla.allele.2.$gene.fasta
            done        
            """%(sys.path[0],workdir,gene, sample_list[0], sample_list[1], sample_list[2], sys.path[0])

        os.system(order)     
        # print (order)    

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


