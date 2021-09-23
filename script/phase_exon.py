#!/usr/bin/env python3

from argparse import ArgumentParser
from argparse import ArgumentTypeError
from my_imports import *
import time
import re
from itertools import combinations, permutations
import realign_and_sv_break
#from algorithm_test import Workflow

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Please give right flag (True or False).')

def population_bool(x):
    if x.lower() == 'asian': 
        return 'Asian'
    elif x.lower() == 'black':
        return 'Black'
    elif x.lower() == 'caucasian':
        return 'Caucasian'
    else:
        raise ArgumentTypeError('Please give right population type, you should choose from Asian, Black or Caucasian.')

Usage = \
"""
python3 phase_exon.py [options] 

Help information can be found by phase_exon.py -h/--help, additional information can be found in \
README.MD or https://github.com/deepomicslab/HLAPro.
"""
scripts_dir=sys.path[0]+'/'
parser = ArgumentParser(description="SpecHLA",prog='python3 phase_exon.py',usage=Usage)
optional=parser._action_groups.pop()
required=parser.add_argument_group('required arguments')
flag_parser = parser.add_mutually_exclusive_group(required=False)
flag_data = parser.add_mutually_exclusive_group(required=False)
#necessary parameter
required.add_argument("--ref",help="The hla reference file used in alignment",dest='ref',metavar='', type=str)
required.add_argument("-b", "--bam",help="The bam file of the input samples.",dest='bamfile',metavar='')
required.add_argument("-v", "--vcf",help="The vcf file of the input samples.",dest='vcf',metavar='')
required.add_argument("-a", "--phase",help="Choose phasing method, SpecHap if True, otherwise use PStrain.\
     Default is SpecHap.",dest='phase_flag',metavar='',default=True,type=str2bool)
# required.add_argument("-e", "--anno",help="Choose annotation, allele with population freq>0 if True, otherwise use all.\
#      Default is True.",dest='anno_flag',metavar='',default=True,type=str2bool)
required.add_argument("-o", "--outdir",help="The output directory.",dest='outdir',metavar='')
#alternative parameter
optional.add_argument("--freq_bias",help="freq_bias (default is 0.05)",dest='freq_bias',\
    metavar='',default=0.05, type=float)
flag_parser.add_argument("-g", "--pstrain_bk",help="Use Pstrain bk if True, else use spechap."
     ,dest='pstrain_bk',metavar='', default=True, type=str2bool)  
optional.add_argument("--snp_dp",help="The minimum depth of SNPs to be considered in HLAtyping\
     step (default is 5).",dest='snp_dp',metavar='',default=5, type=int)
optional.add_argument("--allele_support",help="The minimum depth of allele to be considered as hete in HLAtyping\
     step (default is 2).",dest='allele_support',metavar='',default=2, type=int)     
optional.add_argument("--indel_len",help="The maximum length for indel to be considered in HLAtyping\
     step (default is 150).",dest='indel_len',metavar='',default=150, type=int)
optional.add_argument("--block_len",help="The minimum length for block to be considered in final\
     result (default is 300).",dest='block_len',metavar='',default=50, type=int)
optional.add_argument("--points_num",help="The minimum hete loci number for block to be considered\
     in final result (default is 2).",dest='points_num',metavar='',default=3, type=int)
optional.add_argument("--prefix",help="sample name",dest='prefix',metavar='',default='sample', type=str)
#break points
optional.add_argument("--reads_num",help="The number of supporting reads between two adjcent loci\
     lower than this value will be regard as break points.(default is 10)",dest='reads_num',\
     metavar='',default=10, type=int)
optional.add_argument("--noise_num",help="If the haplotype number is 2, there will be at most two \
    types of linked reads. If the third type of reads number is over this value, then these two \
    loci will be regarded as break points.(default is 5)",dest='noise_num',metavar='',default=5, \
    type=int)
optional.add_argument( "--lambda1",help="The weight of prior knowledge while rectifying genotype\
 frequencies. The value is between 0~1. (default is 0.0)",dest='lambda1',metavar='',default=0,\
  type=float)
optional.add_argument( "--lambda2",help="The weight of prior estimation while rectifying second\
 order genotype frequencies. The value is between 0~1. (default is 0.0)",dest='lambda2',\
 metavar='',default=0, type=float)
optional.add_argument("--elbow",help="The cutoff of elbow method while identifying HLAs number. \
If the loss reduction ratio is less than the cutoff, then the HLAs number is determined.",\
    dest='elbow',metavar='',default=0.24, type=float)
optional.add_argument("-w", "--weight",help="The weight of genotype frequencies while computing\
     loss, then the weight of linked read type frequencies is 1-w. The value is between 0~1.\
      (default is 1)",dest='weight',metavar='',default=1, type=float)
parser._action_groups.append(optional)
args = parser.parse_args()

def exon_region():
    file = '%s/exon/exon_extent.bed'%(sys.path[0])
    exon_region_list = []
    for line in open(file, 'r'):
        line = line.strip()
        array = line.split()
        exon_region_list.append(array)
    return exon_region_list

def if_in_exon(chrom, locus, exon_region_list):
    in_flag = False
    for exon in exon_region_list:
        if chrom == exon[0] and locus >= float(exon[1]) and locus <= float(exon[2]):
            in_flag = True
    return in_flag

def update_locus_for_exon(chrom, locus, exon_region_accumulate_list, ref='', alt_tuple=()):
    # print (chrom, locus, exon_region_accumulate_list)
    new_alt_tuple = alt_tuple
    for exon in exon_region_accumulate_list:
        if chrom == exon[0] and float(locus) >= float(exon[1]) and float(locus) <= float(exon[2]):
            # locus = locus - float(exon[1]) + 1 + float(exon[3])
            if locus + len(ref) - 1 > float(exon[2]):
                ref = ref[:int(float(exon[2])-locus) + 1]
            new_alt_tuple = []
            for alt in alt_tuple:
                if locus + len(alt) - 1 > float(exon[2]):
                    alt = alt[:int(float(exon[2])-locus) + 1]
                    new_alt_tuple.append(alt)
                else:
                    new_alt_tuple.append(alt)
            new_alt_tuple = tuple(new_alt_tuple)
            locus = locus - float(exon[1]) + 1 + float(exon[3])
        elif chrom == exon[0] and float(locus) + len(ref) >= float(exon[1]) and float(locus) + len(ref) <= float(exon[2]):
            # print (chrom, locus, exon_region_accumulate_list)
            over_len = int(float(locus) + len(ref) - float(exon[1]) - 1)
            ref = ref[over_len:]
            new_alt_tuple = []            
            for alt in alt_tuple:
                alt = alt[over_len:]
                new_alt_tuple.append(alt)
            new_alt_tuple = tuple(new_alt_tuple)
            locus = 1 + float(exon[3]) #float(exon[1])
            # print ('special', chrom, str(int(locus)), ref, new_alt_tuple)
        #else:
        #    print ('The locus is not in the exon region.', chrom, locus)
    if ref == '' and alt_tuple == ():
        return chrom, str(int(locus))
    else:
        # print ('locus', ref, alt_tuple)
        # print (chrom, str(int(locus)))
        return chrom, str(int(locus)), ref, new_alt_tuple

def exon_accumulate_locus(exon_region_list):   
    gene = ''
    accumulate_locus = 0
    for exon in exon_region_list:
        if exon[0] != gene:
            gene = exon[0]
            accumulate_locus = 0
            exon.append(accumulate_locus)
            accumulate_locus = accumulate_locus + float(exon[2]) - float(exon[1]) + 1
        else:
            exon.append(accumulate_locus)
            accumulate_locus = accumulate_locus + float(exon[2]) - float(exon[1]) + 1
    return exon_region_list
    # print (exon_region_list)

def read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,chrom_name,freq_bias,strainsNum):
    exon_region_list = exon_region()
    exon_region_accumulate_list = exon_accumulate_locus(exon_region_list)
    pysam.index(bamfile)
    samfile = pysam.AlignmentFile(bamfile, "rb")
    if not os.path.exists(outdir):
        os.system('mkdir '+outdir)
    in_vcf = VariantFile(vcffile)
    md_vcf = VariantFile('%s/middle.vcf.gz'%(outdir),'w',header=in_vcf.header)
    sample = list(in_vcf.header.samples)[0]
    snp_list, beta_set, allele_set = [], [], []
    for record in in_vcf.fetch():
        if record.chrom != chrom_name:
            continue
        if 'DP' not in record.info.keys() or record.info['DP'] < snp_dp:
            continue
        geno = record.samples[sample]['GT']    
        depth = record.samples[sample]['AD']
        # if record.chrom == 'HLA_DRB1' and record.pos > 9000:
        #     continue

        if geno == (1,2,3) or geno == (0, 1, 2):
            norm_depth = np.array(depth)/sum(depth)
            dp_sort=np.argsort(norm_depth)
            max_two_freq = norm_depth[dp_sort[-1]] + norm_depth[dp_sort[-2]]
            alt_genos = []
            fresh_depth = []
            sum_freq = 0
            ref_flag = False
            for i in range(strainsNum):
                sum_freq += norm_depth[dp_sort[-(1+i)]]
                
                if dp_sort[-(1+i)] == 0:
                    ref_flag = True
                if dp_sort[-(1+i)] != 0:
                    alt_genos.append(record.alts[dp_sort[-(1+i)] - 1])
                    fresh_depth.append(depth[dp_sort[-(1+i)]])
            
            if ref_flag == False :
                star_geno = 1
                fresh_depth = [depth[0]] + fresh_depth
            else:
                star_geno = 0
                fresh_depth = [depth[0]] + fresh_depth
            new_geno = []
            for i in range(strainsNum):
                new_geno.append(i + star_geno)
            #there is a bug for three haps
            new_geno.append(new_geno[-1])

            if sum_freq > 0.7:
                # print ('before', record)
                record.alts = alt_genos
                record.samples[sample]['GT'] = new_geno
                record.samples[sample]['AD'] = fresh_depth
                geno = tuple(new_geno)
                
            else:
                continue


        if record.qual == None:
            print ('WARNING: no vcf quality value.')

        if len(record.ref) > indel_len:
            continue
        if len(record.alts) > 2:
            continue
        alt_too_long = False
        for alt_allele in record.alts:
            if len(alt_allele) > indel_len:
                alt_too_long = True
        if alt_too_long:
            continue

        if geno != (1,1,1):
            if geno == (1,1,2) or geno == (1,2,2):                
                snp = [record.chrom,record.pos,record.alts[0],record.alts[1],record.ref]
                # beta=depth[2]/(depth[1] + depth[2])
            else:
                snp=[record.chrom,record.pos,record.ref,record.alts[0],record.ref]
                # beta=depth[1]/(depth[0] + depth[1])

            reads_list = reads_support(samfile, snp)
            allele_dp = [len(reads_list[0]), len(reads_list[1])]
            new_dp=sum(allele_dp)
            if new_dp == 0:
                print ('%s is removed due to the low depth %sx'%(snp,new_dp))
                continue
            beta=float(allele_dp[1])/new_dp
            if beta <= freq_bias or len(reads_list[1]) < args.allele_support:
                if geno == (1,1,2) or geno == (1,2,2):
                    record.samples[sample]['GT']= generate_geno(strainsNum, 1)
                else:
                    record.samples[sample]['GT']= generate_geno(strainsNum, 0)
                record.samples[sample].phased=True
            elif beta >= 1 - freq_bias or len(reads_list[0]) < args.allele_support:
                if geno == (1,1,2) or geno == (1,2,2):
                    record.samples[sample]['GT']= generate_geno(strainsNum, 2)
                else:
                    record.samples[sample]['GT']= generate_geno(strainsNum, 1)
                record.samples[sample].phased=True
            else:
                snp_list.append(snp)
                beta_set.append(allele_dp)
                allele_set.append(2)
        else:
            record.samples[sample]['GT']= generate_geno(strainsNum, 1)
            record.samples[sample].phased=True
        md_vcf.write(record)
    in_vcf.close()
    md_vcf.close()
    # os.system('%s/../bin/tabix -f %s/middle.vcf.gz'%(sys.path[0],outdir))
    print ("The number of hete loci is %s."%(len(snp_list)))
    return snp_list, beta_set, allele_set

def generate_geno(strainsNum, geno):
    genotype = []
    for i in range(strainsNum):
        genotype.append(geno)
    return tuple(genotype)

def delta(outdir,extractHAIRS,bamfile,beta_set,reffile,vcffile):
    hapcut_order='%s --bam %s --VCF %s --indels 1 --ref %s --out %s/vcf.conn'%(extractHAIRS,bamfile,vcffile,reffile,outdir)
    os.system(hapcut_order)
    snp_num=len(beta_set)
    snp_num =500
    delta_set=[]
    for i in range(snp_num-1):
        # delta_set.append([0]*(len(beta_set[i])*len(beta_set[i+1])))
        delta_set.append([0]*4)
    for line in open('%s/vcf.conn'%(outdir)):
        line=line.strip()
        array=line.split()
        if array[0] == '1': 
            delta_index=int(array[2])-1
            geno_type=array[3]
            for i in range(len(geno_type)-1):
                # fol_allele=len(beta_set[delta_index+i+1])
                fol_allele = 2
                array_index=int(geno_type[i])*fol_allele+int(geno_type[i+1])
                # print (array_index, len(delta_set), delta_index + 1)
                delta_set[delta_index+i][array_index]+=1
                # delta_set[delta_index+i][int(geno_type[i])][int(geno_type[i+1])]+=1
    for i in range(len(delta_set)):
        delta=delta_set[i]
    # for delta in delta_set:
        delta=np.array(delta)
        sum_dp=sum(delta)
        if sum_dp>4:
            delta=delta/sum_dp
            delta=np.round(delta,6).tolist()
            # delta_set[i]=delta
            # print (delta)
        else:
            delta=[0]*len(delta)
        delta_set[i]=delta
    return delta_set 

def output(outdir,final_alpha,seq_list,snp_list,gene):
    exon_region_list = exon_region()
    exon_region_accumulate_list = exon_accumulate_locus(exon_region_list)
    ra_file=open(outdir+'/%s_freq.txt'%(gene),'w')    
    print ('# HLA\tFrequency',file=ra_file)
    seq=np.array(seq_list)
    seq=np.transpose(seq)
    # print (len(snp_list),len(seq),'points num')
    snp=snp_list  
    alpha=final_alpha   
    k=len(alpha) 
    for j in range(k):
        print ('str-'+str(j+1),alpha[j],file=ra_file)
    ra_file.close()
    ##############################
    he=0
    m = VariantFile('%s/middle.vcf.gz'%(outdir))
    out = VariantFile('%s/%s.vcf.gz'%(outdir,gene),'w',header=m.header)
    sample = list(m.header.samples)[0]
    for record in m.fetch():
        geno = record.samples[sample]['GT']    
        depth = record.samples[sample]['AD']

        if record.samples[sample].phased != True:
            if geno == (1,1,2) or geno == (1,2,2):
                phased_locus=seq[he]
                for i in range(len(phased_locus)):
                    phased_locus[i]+=1
            else:
                phased_locus=seq[he]

            record.samples[sample]['GT']= tuple(phased_locus)
            record.samples[sample].phased=True
            he+=1
        # chrom, locus, record.ref, record.alts = update_locus_for_exon(record.chrom, float(record.pos), exon_region_accumulate_list, record.ref, record.alts)
        # record.chrom, record.pos = chrom, int(locus)
        out.write(record)
    m.close()
    out.close()

def rectify(snp_list,beta_set,delta_set,lambda1,lambda2,germline_flag):
    nucleotide={'A':0,'T':1,'C':2,'G':3}
    first_beta=[]
    for i in range(len(beta_set)):
        #prior freq
        snp=snp_list[i]
        pos=int(snp[1])
        ref=snp[2]
        alt=snp[3][0]
        # print (snp_list[i])
        prior_f=np.array([1.0/len(beta_set[i])]*len(beta_set[i]))
        dp=sum(beta_set[i])
        hat_beta=np.array(beta_set[i])/dp
        beta=hat_beta*(1-lambda1/(1+dp))+prior_f*(lambda1/(1+dp))
        beta=beta/sum(beta)
        beta=np.round(beta,6)
        first_beta.append(beta.tolist())

    # rectify beta^2
    sec_beta=[]
    for i in range(len(delta_set)):
        c=sum(delta_set[i])
        if c==0:
            hat_delta=np.array(delta_set[i])
        else:
            hat_delta=np.array(delta_set[i])/sum(delta_set[i])
        inde_delta=[]
        for m in range(len(first_beta[i])):
            for n in range(len(first_beta[i+1])):
                inde_delta.append(first_beta[i][m]*first_beta[i+1][n])
        inde_delta=np.array(inde_delta)
        delta=inde_delta*(lambda2/(1+c)) + hat_delta*(1-lambda2/(1+c))  #retify

        delta=np.round(delta,6)
        
        #turn the two biggest values to 0.5/0.5, used to handle normal samples
        index_sort=np.argsort(delta)
        if germline_flag and c >= 3 and index_sort[0] + index_sort[1] == 3:
            index_sort=np.argsort(delta)
            delta[index_sort[0]],delta[index_sort[1]]=0,0
            delta[index_sort[2]],delta[index_sort[3]]=0.5,0.5
            # delta = [0, 0, 0, 0]
            # delta[index_sort[3]], delta[3 - index_sort[3]] = 0.5, 0.5
        elif germline_flag:
            if index_sort[0] + index_sort[2] == 3:
                delta[index_sort[0]],delta[index_sort[2]]=0,0
                delta[index_sort[1]],delta[index_sort[3]]=0.5,0.5
            elif index_sort[1] + index_sort[2] == 3:
                delta[index_sort[1]],delta[index_sort[2]]=0,0
                delta[index_sort[0]],delta[index_sort[3]]=0.5,0.5
            else:
                delta = np.array([0, 0, 0, 0])

        sec_beta.append(delta.tolist())
    return first_beta,sec_beta

def extract(first,second,file): #the index of first and second should be 0-index
    allele_index={'A':0,'T':1,'C':2,'G':3}
    dict={}
    i=0
    hla_name=[]
    for line in open(file,'r'):
        i+=1
        line=line.strip()
        array=line.split()
        if i==1 or array[1] in hla_name:
            continue
        hla_name.append(array[1])        
        freq=0
        for j in range(2,len(array)):
            if array[j] != '-' and isfloat(array[j]) :
                freq+=float(array[j])
        if i == 2:
            ref_seq=array[22]
            ref_index=[]
            for z in range(len(ref_seq)):
                if ref_seq[z] != '.' and ref_seq[z] != '|':
                    ref_index.append(z)
        new_first=ref_index[first]
        new_second=ref_index[second]
        seq=array[22]
        first_allele=seq[new_first]
        second_allele=seq[new_second]
        if first_allele == '-':
            first_allele=ref_seq[new_first]
        if second_allele == '-':
            second_allele=ref_seq[new_second]
        key=first_allele+second_allele
        if key in dict.keys():
            dict[key]+=freq
        else:
            dict[key]=freq
    return dict
        
def isfloat(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def second_prior(snp1,snp2):
    allele_index={'A':0,'T':1,'C':2,'G':3}
    if str(snp1[0]) !=  str(snp2[0]) or str(snp1[0]) not in ['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1']:
        return np.array([0])
    first=int(snp1[1])-1-1000
    second=int(snp2[1])-1-1000
    snp1_allele=[snp1[2][0]]+[snp1[3][0]]
    snp2_allele=[snp2[2][0]]+[snp2[3][0]]

    file=sys.path[0]+'/'+'database/%s.fre.alignments.mat.txt'%(str(snp1[0]))
    dict=extract(first,second,file)
    prior_sec=[]
    for i in range(len(snp1_allele)):
        for j in range(len(snp2_allele)):           
            link_type=snp1_allele[i]+snp2_allele[j]
            #print (link_type,dict)
            if link_type in dict.keys():
                prior_sec.append(float(dict[link_type]))
            else:
                prior_sec.append(0)
    #print (snp1,snp2,snp1_allele,snp2_allele,prior_sec)
    prior_sec=np.array(prior_sec)
    if sum(prior_sec) != 0:
        prior_sec=prior_sec/sum(prior_sec)
    else:
        prior_sec=np.array([0])
    return prior_sec

def reads_support(samfile,first):   
    reads_list=[]
    allele_num=len(first[3])+1
    for i in range(allele_num):
        reads_list.append([])
    num=0
    for read in samfile.fetch(str(first[0]),int(first[1])-1,int(first[1])):
        
        if int(first[1])-1 in read.get_reference_positions(full_length=True) and read.mapping_quality >30:   
            
            reads_index=read.get_reference_positions(full_length=True).index(int(first[1])-1)
            if first[2][0] != first[3][0] and 1 > 2:  #if there are more than two allele, there will be noise, if we only consider the first .
                #if the first allele is not same for indel alleles, we can just focus on the first locus 
                if read.query_sequence[reads_index] == first[2][0]:
                    reads_list[0].append(read.query_name)
                elif read.query_sequence[reads_index] == first[3][0]:
                    reads_list[1].append(read.query_name)
            else:  
                index_list=[]
                true_ref=first[4]
                for i in range(len(true_ref)):
                    position=first[1]+i
                    point_flag=isin(position-1,read.get_reference_positions(full_length=True))
                    if point_flag:
                        position_index=read.get_reference_positions(full_length=True).index(position-1)
                        index_list.append(position_index)
                allele_list=read.query_sequence[index_list[0]:index_list[-1]+1].upper()
                if allele_list == first[2]:
                    reads_list[0].append(read.query_name)
                elif allele_list == first[3]:
                    reads_list[1].append(read.query_name)

            # if first[0] == 'HLA_A' and first[1] == 2886:
            #     print (read.mapq)
                    # print (first,allele_list,index_list)

    #if int(first[1]) == 2065:
    #    print (len(reads_list[0]),len(reads_list[1]),first)
    return reads_list

def share_reads(samfile,left,right,new_left):
    left_num=2#len(left[3])+1
    right_num=2#len(right[3])+1
    # left_reads=reads_support(samfile,left)
    left_reads=new_left
    right_reads=reads_support(samfile,right)
    delta_count=[]
    for i in range(left_num):
        for j in range(right_num):
            left_set=left_reads[i]
            right_set=right_reads[j]
            same_num=len(set(left_set).intersection(set(right_set)))
            delta_count.append(same_num)
    return delta_count,right_reads

def second_beta(bamfile,snp_list):
    delta_set=[]
    samfile = pysam.AlignmentFile(bamfile, "rb")
    new_left=''
    for i in range(len(snp_list)-1):  
        left=snp_list[i]
        right=snp_list[i+1]  
        if new_left=='':   
            new_left=reads_support(samfile,left)
        delta_count,right_reads=share_reads(samfile,left,right,new_left)
        delta_set.append(delta_count)
        new_left=right_reads
    return delta_set

def isin(x,seq):
    try:
        seq.index(x)
        return True
    except :
        return False

def relate_order(file):
    block_dict={}
    dict={}
    previous_gene=''
    for line in open(file,'r'):
        if line[0] == '#':
            continue
        line=line.strip()
        array=line.split()
        gene_name=array[0]
        locus=array[1]
        array=array[6:]
        hla_num=len(array)
        if gene_name != previous_gene:
            block_dict[previous_gene] = dict
            dict={}
            previous_gene=gene_name            
            ref_order=[]
            for j in range(hla_num):
                ref_order.append(j)
        if len(array)<len(ref_order):
            continue
        new_order=[]
        # print (hla_num,ref_order,array)
        for e in range(hla_num):
            new_order.append(array[int(ref_order[e])])
        ref_order=new_order[:]
        dict[locus]=ref_order
        # print (gene_name,ref_order)
    block_dict[previous_gene] = dict
    # print (block_dict)
    return block_dict

def newphase(outdir,final_alpha,seq_list,snp_list,vcffile,gene):
    exon_region_list = exon_region()
    exon_region_accumulate_list = exon_accumulate_locus(exon_region_list)
    file=outdir+'/%s_break_points_phased.txt'%(gene)
    block_dict=relate_order(file)
    seq=np.array(seq_list)
    seq=np.transpose(seq)
    snp=snp_list  
    update_seqlist=[]
    alpha=final_alpha   
    k=len(alpha) 
    ##############################

    if gene in block_dict.keys():
        gene_dict=block_dict[gene]
    else:
        gene_dict={}
    ref_order=[]
    for orde in range(k):
        ref_order.append(orde)
    he=0
    m = VariantFile('%s/middle.vcf.gz'%(outdir))
    rephase_file = '%s/%s.rephase.vcf.gz'%(outdir,gene)
    if os.path.isfile(rephase_file):
        os.system('rm %s'%(rephase_file))
    out = VariantFile(rephase_file,'w',header=m.header)
    sample = list(m.header.samples)[0]
    for record in m.fetch():
        geno = record.samples[sample]['GT']    
        depth = record.samples[sample]['AD']
        if record.samples[sample].phased != True:
            if geno == (1,1,2) or geno == (1,2,2):
                phased_locus = seq[he]
                for i in range(len(phased_locus)):
                    phased_locus[i] += 1
            else:
                phased_locus=seq[he]

            update_phased_locus=[]
            for pp in range(k):
                update_phased_locus.append(phased_locus[int(ref_order[pp])])
            phased_locus=update_phased_locus[:]
            record.samples[sample]['GT']= tuple(phased_locus)
            record.samples[sample].phased=True

            if phased_locus[0] > 1 or phased_locus[1]>1:
                phased_locus[0]-=1
                phased_locus[1]-=1
            update_seqlist.append(phased_locus)
            he+=1
        chrom, locus, record.ref, record.alts = update_locus_for_exon(record.chrom, float(record.pos), exon_region_accumulate_list, record.ref, record.alts)
        record.chrom, record.pos = chrom, int(locus)
        out.write(record)
        # print ('new order',gene_dict)
        if str(locus) in gene_dict.keys():  #rephase new order
            ref_order=gene_dict[str(locus)]
            # print ('yes', ref_order)
            # print ('new order',record.pos, ref_order)

    m.close()
    out.close()
    os.system('%s/../bin/tabix %s/%s.rephase.vcf.gz'%(sys.path[0],outdir,gene))
    return update_seqlist

def gene_phased(update_seqlist,snp_list):
    gene_profile={}
    gene_name=['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']
    for gene in gene_name:
        gene_snp=[]
        gene_seq=[]
        for i in range(len(snp_list)):
            if snp_list[i][0] == gene:
                gene_snp.append(snp_list[i])
                gene_seq.append(update_seqlist[i])
        gene_seq = np.array(gene_seq) 
        gene_seq = np.transpose(gene_seq)
        gene_profile[gene] = [gene_snp,gene_seq]
    return gene_profile
    # print (gene_profile)

def no_snv_gene_phased(vcffile, outdir, gene, strainsNum):
    exon_region_list = exon_region()
    exon_region_accumulate_list = exon_accumulate_locus(exon_region_list)
    # in_vcf = VariantFile(vcffile)
    in_vcf = VariantFile('%s/middle.vcf.gz'%(outdir))
    out_vcf = VariantFile('%s/%s.rephase.vcf.gz'%(outdir, gene),'w',header=in_vcf.header)
    for_spec_vcf = VariantFile('%s/%s.vcf.gz'%(outdir, gene),'w',header=in_vcf.header)
    sample = list(in_vcf.header.samples)[0]
    for record in in_vcf.fetch():
        if record.chrom == gene:
            chrom, pos = update_locus_for_exon(record.chrom, float(record.pos), exon_region_accumulate_list)
            record.chrom, record.pos = chrom, int(pos)
            out_vcf.write(record)
            for_spec_vcf.write(record)
    in_vcf.close()
    out_vcf.close()
    for_spec_vcf.close()
    os.system('%s/../bin/tabix %s/%s.rephase.vcf.gz'%(sys.path[0],outdir,gene))
    ra_file=open(outdir+'/%s_freq.txt'%(gene),'w')    
    print ('# HLA\tFrequency',file=ra_file)
    print ('str-'+str(1), 1, file=ra_file)
    for j in range(1, strainsNum):
        print ('str-'+str(j+1), 0, file=ra_file)
    ra_file.close()

def read_dup():
    dup_dict={}
    for line in open(sys.path[0]+'/complex_region_wes.txt','r'):
        line=line.strip()
        array=line.split()
        dup_dict[array[0]]=array[1:]
    return dup_dict
########align the segs from SV haps result

def generate_break_points(outdir,snp_list,delta_set,strainsNum):   
    exon_region_list = exon_region()
    exon_region_accumulate_list = exon_accumulate_locus(exon_region_list) 
    #break points
    bp=open(outdir+'/%s_break_points.txt'%(gene),'w')
    print ('#gene\tlocus\t00\t01\t10\t11\tpoints_num\tnext_locus',file=bp) 
    points_number=0
    rephase_break_point_list = []
    break_points_index = []
    b_hard = False
    for i in range(len(snp_list)-1):

        points_number += 1
        sort_index=np.argsort(np.array(delta_set[i]))
        max_index=sort_index[-1]
        sec_index=sort_index[-2]
        
        if points_number < args.points_num:
            continue

        # if delta_set[i][sort_index[-3]] > 5 or abs(max_index - sec_index) == 2:
        #     print (snp_list[i],delta_set[i])
        chrom, locus = snp_list[i][0], snp_list[i][1]
        chrom, locus = update_locus_for_exon(chrom, locus, exon_region_accumulate_list)
        next_chrom, next_locus = snp_list[i+1][0], snp_list[i+1][1]
        next_chrom, next_locus = update_locus_for_exon(next_chrom, next_locus, exon_region_accumulate_list)
        if chrom == 'HLA_B':
            if int(float(locus)) in [250, 290, 306]:  
                b_hard = True
        # b_hard = False
        if sum(delta_set[i]) <args.reads_num or (strainsNum == 2 and max_index + sec_index != 3) or (strainsNum == 2 and delta_set[i][sort_index[-3]] > args.noise_num) or b_hard or (strainsNum == 2 and delta_set[i][sort_index[-2]] < 2):
        #if sum(delta_set[i]) <args.reads_num or (strainsNum == 2 and max_index + sec_index != 3 and sort_index[-3] + sec_index != 3 and sort_index[-3] + max_index != 3) or b_hard: #or (strainsNum == 2 and delta_set[i][sort_index[-3]] > args.noise_num) 
            print ('breakpoints while phasing', chrom, int(float(locus)), delta_set[i][0], delta_set[i][1],delta_set[i][2],delta_set[i][3], b_hard)
            print (chrom, int(float(locus)), delta_set[i][0], delta_set[i][1],delta_set[i][2],delta_set[i][3],points_number,int(float(next_locus)),file=bp)
            b_hard = False
            rephase_break_point_list.append([chrom, int(float(locus))])
            # break_points_index.append(int(float(locus)))
            break_points_index.append(i)
            points_number=0
    bp.close()
    return break_points_index

def all_possibilities(break_points_index, strainsNum, snp_list, seq_list, final_alpha, gene, use_pstrain_break_point_flag): #need to update alpha finally
    my_table = all_table( len(break_points_index) + 1, 2)
    locus_seq = '' 
    if not use_pstrain_break_point_flag: #use spechap breakpoints
        for i in range(len(snp_list)):
            if int(snp_list[i][1]) in break_points_index:
                locus_seq += '|' + str(i)
                locus_seq += '_'
            else:
                locus_seq += '|' + str(i)
    elif use_pstrain_break_point_flag: #use pstrain breakpoints
        for i in range(len(snp_list)):
            if i in break_points_index:
                locus_seq += '|' + str(i)
                locus_seq += '_'
            else:
                locus_seq += '|' + str(i)
    else:
        print ('Wrong bk flag, P or S.')
    locus_set = []
    for aa in locus_seq.split('_'):
        new_aa = []
        for a in aa.split('|'):
            if a != '':
                new_aa.append(int(float(a)))
        locus_set.append(new_aa)
    index_set = []
    for i in range(len(locus_set)):
        for j in range(len(locus_set[i])):
            index_set.append(i)
    new_seq_list_set = []
    for i in range(len(my_table)):
        combination = my_table[i]
        new_seq_list = []
        for j in range(strainsNum):
            reverse_order = combination
            seq = seq_list[j]
            new_seq = []
            for u in range(len(seq)):

                order = reverse_order[index_set[u]]
                if order == 0:
                    new_seq.append(seq[u])
                else:
                    new_seq.append(1 - seq[u])
            new_seq_list.append(new_seq)
        rephase_output(outdir,final_alpha,new_seq_list,snp_list,vcffile,gene, i)
        new_seq_list_set.append(new_seq_list)
               
def rephase_output(outdir,final_alpha,seq_list,snp_list,vcffile,gene, index):
    exon_region_list = exon_region()
    exon_region_accumulate_list = exon_accumulate_locus(exon_region_list)
    seq=np.array(seq_list)
    seq=np.transpose(seq)
    snp=snp_list  
    update_seqlist=[]
    alpha=final_alpha   
    k=len(alpha) 
    ##############################
    he=0
    m = VariantFile('%s/middle.vcf.gz'%(outdir))
    rephase_file = '%s/%s.%s.rephase.vcf.gz'%(outdir,gene,index)
    if os.path.isfile(rephase_file):
        os.system('rm %s'%(rephase_file))
    out = VariantFile(rephase_file,'w',header=m.header)
    sample = list(m.header.samples)[0]
    for record in m.fetch():
        geno = record.samples[sample]['GT']    
        depth = record.samples[sample]['AD']
        if record.samples[sample].phased != True:
            if geno == (1,1,2) or geno == (1,2,2):
                phased_locus = seq[he]
                for i in range(len(phased_locus)):
                    phased_locus[i] += 1
            else:
                phased_locus=seq[he]

            record.samples[sample]['GT']= tuple(phased_locus)
            record.samples[sample].phased=True
            
            if phased_locus[0] > 1 or phased_locus[1]>1:
                phased_locus[0]-=1
                phased_locus[1]-=1
            update_seqlist.append(phased_locus)
            he+=1
        chrom, locus, record.ref, record.alts = update_locus_for_exon(record.chrom, float(record.pos), exon_region_accumulate_list, record.ref, record.alts)
        record.chrom, record.pos = chrom, int(locus)
        out.write(record)

    m.close()
    out.close()
    os.system('%s/../bin/tabix %s/%s.%s.rephase.vcf.gz'%(sys.path[0], outdir, gene, index))
    exon_ref = '%s/../db/ref/hla_exon_ref.fasta'%(sys.path[0])
    for i in range(1, k+1):
        fastq2 = '%s/../bin/samtools faidx %s %s | \
        %s/../bin/bcftools consensus -H %s %s/%s.%s.rephase.vcf.gz\
         >%s/%s.%s.%s.fasta'%(sys.path[0], exon_ref, gene, sys.path[0], i, outdir, gene, index, outdir, gene, index,i)
        os.system(fastq2)

def all_table(k,allele_num):
    mytable=[]
    for j in range(allele_num):
        mytable.append([j])
    if k>1:
        for i in range(k-1):
            double_table=[]
            for j in range(allele_num):
                add_allele=[]
                for list in mytable:
                    newlist=list[:]
                    newlist.append(j)
                    # newlist = np.array(newlist)
                    # newlist = [newlist, 1- newlist]
                    add_allele.append(newlist)
                double_table+=add_allele
            mytable=double_table
    num = len(mytable)
    return mytable[:int(num/2)] 

def convert(outdir, gene, invcf):
    exon_region_list = exon_region()
    exon_region_accumulate_list = exon_accumulate_locus(exon_region_list)
    break_points_index = []
    formate_vcf = outdir + '/%s.vcf.gz'%(gene)
    bp = open(outdir + '/%s_break_points.txt'%(gene), 'w')
    print ('#gene   locus   00      01      10      11      points_num      next_locus', file = bp)
    m = VariantFile(invcf)
    out = VariantFile(formate_vcf,'w',header=m.header)
    sample = list(m.header.samples)[0]
    add_block = 1
    for record in m.fetch():
        # if record.qual < 1
        if record.chrom != gene:
            continue
        geno = record.samples[sample]['GT']    
        depth = record.samples[sample]['AD']
        if record.samples[sample].phased != True:
            #record.samples[sample]['GT']= (1,1)
            record.samples[sample].phased = True
            record.samples[sample]['PS'] = add_block
        if record.samples[sample]['PS'] != add_block:
            print (record.chrom, record.pos, '- - - -', 20, record.pos + 100, file = bp)
            add_block = record.samples[sample]['PS']
        # if record.samples[sample].phased != True:
        #     record.samples[sample].phased = True
        #     if record.samples[sample]['GT'] != (1,1) and record.samples[sample]['GT'] != (0,0):
        #         print (record.chrom, record.pos, '- - - -', 20, record.pos + 100, file = bp)
            #record.samples[sample]['GT']= (1,1)
        #     record.samples[sample].phased = True
        #     record.samples[sample]['PS'] = add_block
        # if record.samples[sample]['PS'] != add_block:
        #     break_points_index.append(int(record.pos))
        #     print (record.chrom, record.pos, '- - - -', 20, record.pos + 100, file = bp)
        #     add_block = record.samples[sample]['PS']
        # print (record.alts)
        #convert pos
        
        chrom, locus, record.ref, record.alts = update_locus_for_exon(record.chrom, float(record.pos), exon_region_accumulate_list, record.ref, record.alts)
        record.chrom, record.pos = chrom, int(locus)
        

        out.write(record)
    m.close()
    out.close()
    bp.close()
    return  break_points_index

def read_spechap_seq(vcf, snp_list):
    exon_region_list = exon_region()
    exon_region_accumulate_list = exon_accumulate_locus(exon_region_list)
    locus_list = []
    for snp in snp_list:
        chrom, locus, ref, alts = update_locus_for_exon(snp[0], float(snp[1]), exon_region_accumulate_list, 'T', ('A'))
        locus_list.append(int(locus))

    seq_list = [[],[]]
    in_vcf = VariantFile(vcf)
    sample = list(in_vcf.header.samples)[0]
    for record in in_vcf.fetch():
        geno = record.samples[sample]['GT']
        #print (record)
        # if geno == (1,1):
        #     continue
        # print (locus_list, record.pos)
        if record.pos not in locus_list:
            continue
        if sum(geno) == 1:
            for i in range(2):
                seq_list[i].append(geno[i])
        else: 
            for i in range(2):
                seq_list[i].append(geno[i] - 1)

    return seq_list 

if __name__ == "__main__":
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})
    else: 
        bamfile,outdir,snp_dp,use_pstrain_break_point_flag,indel_len,freq_bias,vcffile=\
            args.bamfile,args.outdir,args.snp_dp,args.pstrain_bk,args.indel_len,args.freq_bias,args.vcf           
        if not os.path.exists(outdir):
            os.system('mkdir '+ outdir)        
        os.system('rm -f %s/*.rephase.vcf.gz*'%(outdir))
        os.system('rm -f %s/*.*.fasta'%(outdir))
        strainsNum = 2


        for gene in ['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']:
        # for gene in ['HLA_B']:

            ######PStrain-filter-vcf
            snp_list,beta_set,allele_set = read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,gene,\
                freq_bias,strainsNum)     
            if len(snp_list)==0:
                #same seq for two haps
                print ('No hete locus')
                no_snv_gene_phased(vcffile, outdir, gene, strainsNum)
                delta_set, seq_list, final_alpha = [], [[],[]], [0.5, 0.5] 
            else:  
                delta_set=second_beta(bamfile,snp_list)   

                # for i in range(len(snp_list)-1):
                #     print (snp_list[i], beta_set[i], delta_set[i])
                fir_beta,sec_beta=rectify(snp_list,beta_set,delta_set,args.lambda1,args.lambda2,\
                    False)
                if strainsNum==0:
                    wo=Workflow(fir_beta,sec_beta,delta_set,args.weight,args.elbow,allele_set)
                    final_alpha,seq_list,loss=wo.choose_k()
                    strainsNum = len(final_alpha)
                else:
                    wo=Workflow(fir_beta,sec_beta,delta_set,args.weight,args.elbow,allele_set)
                    final_alpha,seq_list,loss = wo.given_k(strainsNum)  
                output(outdir,final_alpha,seq_list,snp_list,gene)

            #####Phase-with-SpecHap
            if args.phase_flag == True:
                my_new_vcf = '%s/%s.vcf.gz'%(outdir, gene)
                os.system('gzip -f -d %s'%(my_new_vcf))
                my_new_vcf = '%s/%s.vcf'%(outdir, gene)
                hla_ref = '%s/../db/ref/hla.ref.extend.fa'%(sys.path[0])
                order = '%s/../bin/ExtractHAIRs --triallelic 1 --indels 1 --ref %s --bam %s --VCF %s --out %s/frament.file'%(sys.path[0], hla_ref, bamfile, my_new_vcf, outdir)
                os.system(order)
                os.system('sort -k3 %s/frament.file >%s/frament.sorted.file'%(outdir, outdir))
                os.system('bgzip %s'%(my_new_vcf))
                my_new_vcf = '%s/%s.vcf.gz'%(outdir, gene)
                os.system('%s/../bin/tabix -f %s'%(sys.path[0],my_new_vcf))
                os.system('cp %s/%s.vcf.gz %s/%s.vcf.pstain.gz'%(outdir, gene, outdir, gene))
                order = '%s/../bin/SpecHap --window_size 15000 --vcf %s --frag %s/frament.sorted.file --out %s/%s.specHap.phased.vcf'%(sys.path[0], my_new_vcf, outdir, outdir,gene)
                os.system(order)
                break_points_index = convert(outdir, gene, '%s/%s.specHap.phased.vcf'%(outdir,gene))
                os.system('%s/../bin/tabix -f %s'%(sys.path[0],my_new_vcf))
                seq_list = read_spechap_seq('%s/%s.vcf.gz'%(outdir, gene), snp_list)
            else:
                use_pstrain_break_point_flag = True

            # use_pstrain_break_point_flag = True
            if use_pstrain_break_point_flag:
                break_points_index = generate_break_points(outdir,snp_list,delta_set,strainsNum)
            print ("The Number of break points is %s."%(len(break_points_index)))
            if len(break_points_index) > 10:  #max break points number is 10.
                break_points_index = break_points_index[:10]

            all_possibilities(break_points_index, strainsNum, snp_list, seq_list, final_alpha, gene, use_pstrain_break_point_flag)   
            print ('Phasing for %s is Done.'%(gene))

        # if args.anno_flag:
        #     if args.phase_flag == True:
        #         anno_parameter = 'spechap'
        #     else:
        #         anno_parameter = 'pstrain'
        # else:
        #     anno_parameter = 'all'
        # annotation = 'perl %s/exon/anno_HLA_pop.pl %s %s %s %s '%(sys.path[0], args.prefix, outdir, args.popu, anno_parameter)
        # os.system(annotation)

        # os.system('cat %s/hla.result.txt'%(outdir))
        # os.system('rm %s/*.rephase.vcf.gz*'%(outdir))
        # os.system('rm %s/HLA_*.*.fasta'%(outdir))




