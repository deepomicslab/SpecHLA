#!/usr/bin/env python3

from argparse import ArgumentParser
from argparse import ArgumentTypeError
from my_imports import *
import time
import re
from itertools import combinations, permutations
import realign_and_sv_break

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Please give right flag (True or False).')

Usage = \
"""
python3 PHLA_multi_haps.py [options] -k/--strainsNum <strains number> -b/--bam <bam file> --ref <reference file>\
 -o/--outdir <output directory>

Before using this script, you should map the reads to your reference, and call SNPs. \
Then you can provide the reference and bam file. Once you have the prior knowledge \
of the strains number, you can provie it to the script by -k/--strainsNum parameter. Otherwise,\
it will choose the strains number itself. 

Help information can be found by PHLA_multi_haps.py -h/--help, additional information can be found in \
README.MD or www.xxx.com.
"""
scripts_dir=sys.path[0]+'/'
parser = ArgumentParser(description="PHLA: type HLA.",prog='python3 PHLA_multi_haps.py',usage=Usage)
optional=parser._action_groups.pop()
required=parser.add_argument_group('required arguments')
flag_parser = parser.add_mutually_exclusive_group(required=False)
flag_data = parser.add_mutually_exclusive_group(required=False)
#necessary parameter
required.add_argument("--ref",help="The hla reference file used in alignment",dest='ref',metavar='', type=str)
required.add_argument("-b", "--bam",help="The bam file of the input samples.",dest='bamfile',metavar='')
required.add_argument("-v", "--vcf",help="The vcf file of the input samples.",dest='vcf',metavar='')
required.add_argument("-k", "--strainsNum",help="The number of haplotypes in the sample (default is chosen by elbow method)."\
    ,dest='k',metavar='', type=int)
required.add_argument("-s", "--sv",help="sv after scanindel",dest='sv',metavar='')
required.add_argument("--gene",help="gene",dest='gene',metavar='', type=str)
required.add_argument("--fq1",help="fq1",dest='fq1',metavar='', type=str)
required.add_argument("--fq2",help="fq2",dest='fq2',metavar='', type=str)
required.add_argument("-o", "--outdir",help="The output directory.",dest='outdir',metavar='')
#alternative parameter
flag_parser.add_argument("-g", "--balance",help="The switch to determine if the sample is balance.\
     True for balance samples and False for imbalance samples(default is False)."\
     ,dest='germline',metavar='', default=False,type=str2bool)
flag_data.add_argument("-d", "--database",help="the switch to determine whether use the information\
     of prior database. true for using (default is false)."\
    ,dest='database_flag',metavar='', default=False,type=str2bool)    
optional.add_argument("-w", "--weight",help="The weight of genotype frequencies while computing\
     loss, then the weight of linked read type frequencies is 1-w. The value is between 0~1.\
      (default is 0.15)",dest='weight',metavar='',default=0, type=float)
optional.add_argument("--freq_bias",help="freq_bias (default is 0.1)",dest='freq_bias',\
    metavar='',default=0.1, type=float)
optional.add_argument( "--lambda1",help="The weight of prior knowledge while rectifying genotype\
 frequencies. The value is between 0~1. (default is 0.0)",dest='lambda1',metavar='',default=0,\
  type=float)
optional.add_argument( "--lambda2",help="The weight of prior estimation while rectifying second\
 order genotype frequencies. The value is between 0~1. (default is 0.0)",dest='lambda2',\
 metavar='',default=0, type=float)
optional.add_argument("--elbow",help="The cutoff of elbow method while identifying HLAs number. \
If the loss reduction ratio is less than the cutoff, then the HLAs number is determined.",\
    dest='elbow',metavar='',default=0.24, type=float)
optional.add_argument("--snp_dp",help="The minimum depth of SNPs to be considered in HLAtyping\
     step (default is 5).",dest='snp_dp',metavar='',default=5, type=int)
optional.add_argument("--snp_qual",help="The minimum quality of SNPs to be considered in HLAtyping\
     step (default is 0.01).",dest='snp_qual',metavar='',default=0.01, type=float)
optional.add_argument("--indel_len",help="The maximum length for indel to be considered in HLAtyping\
     step (default is 150).",dest='indel_len',metavar='',default=150, type=int)
optional.add_argument("--block_len",help="The minimum length for block to be considered in final\
     result (default is 300).",dest='block_len',metavar='',default=300, type=int)
optional.add_argument("--points_num",help="The minimum hete loci number for block to be considered\
     in final result (default is 2).",dest='points_num',metavar='',default=2, type=int)
optional.add_argument("--rlen",help="rlen",dest='rlen',metavar='',default=150, type=int)
#break points
optional.add_argument("--reads_num",help="The number of supporting reads between two adjcent loci\
     lower than this value will be regard as break points.(default is 10)",dest='reads_num',\
     metavar='',default=10, type=int)
optional.add_argument("--noise_num",help="If the haplotype number is 2, there will be at most two \
    types of linked reads. If the third type of reads number is over this value, then these two \
    loci will be regarded as break points.(default is 5)",dest='noise_num',metavar='',default=5, \
    type=int)
parser._action_groups.append(optional)
args = parser.parse_args()

def if_in_deletion(locus, deletion_region):
    dele_flag = False
    for deletion in deletion_region:
        if locus >= deletion[0] and locus < deletion[1]:
            dele_flag = True
    return dele_flag

def read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,chrom_name,freq_bias,strainsNum,deletion_region,snp_qual):
    pysam.index(bamfile)
    samfile = pysam.AlignmentFile(bamfile, "rb")
    if not os.path.exists(outdir):
        os.system('mkdir '+outdir)
    in_vcf = VariantFile(vcffile)
    md_vcf = VariantFile('%s/middle.vcf.gz'%(outdir),'w',header=in_vcf.header)
    sample = list(in_vcf.header.samples)[0]
    snp_list, beta_set, allele_set = [], [], []
    for record in in_vcf.fetch():
        if 'DP' not in record.info.keys() or record.info['DP'] <1:
            continue
        geno = record.samples[sample]['GT']    
        depth = record.samples[sample]['AD']
        dp=sum(depth)  
        if record.qual == None:
            print ('WARNING: no vcf quality value.')
            continue
        if record.qual < snp_qual:
            continue
        if record.chrom != chrom_name:
            continue
        if record.chrom == 'HLA_DRB1' and record.pos >= 3898 and record.pos <= 4400:
            continue
        if dp < snp_dp:
            continue
        if geno == (0,1,2):
            continue
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
        if if_in_deletion(record.pos, deletion_region) and geno != (1,1,1):
            print ('SNV in deletion region, need edit', record)
            if geno == (1,1,2) or geno == (1,2,2):                
                snp = [record.chrom,record.pos,record.alts[0],record.alts[1],record.ref]
            else:
                snp=[record.chrom,record.pos,record.ref,record.alts[0],record.ref]

            reads_list = reads_support(samfile, snp)
            allele_dp = [len(reads_list[0]), len(reads_list[1])]
            new_dp=sum(allele_dp)
            beta=float(allele_dp[1])/new_dp

            if beta <= 1-beta:
                if geno == (1,1,2) or geno == (1,2,2):
                    record.samples[sample]['GT']= generate_geno(strainsNum, 1)
                else:
                    record.samples[sample]['GT']= generate_geno(strainsNum, 0)
                record.samples[sample].phased=True
            elif beta >= 1 - beta:
                if geno == (1,1,2) or geno == (1,2,2):
                    record.samples[sample]['GT']= generate_geno(strainsNum, 2)
                else:
                    record.samples[sample]['GT']= generate_geno(strainsNum, 1)
                record.samples[sample].phased=True
            md_vcf.write(record)
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
            # print (record, allele_dp)
            beta=float(allele_dp[1])/new_dp

            if beta <= freq_bias:
                if geno == (1,1,2) or geno == (1,2,2):
                    record.samples[sample]['GT']= generate_geno(strainsNum, 1)
                else:
                    record.samples[sample]['GT']= generate_geno(strainsNum, 0)
                # record.samples[sample]['GT']= generate_geno(strainsNum, 0)
                #(0, 0, 0) #tuple(phased_locus)
                record.samples[sample].phased=True
            elif beta >= 1 - freq_bias:
                if geno == (1,1,2) or geno == (1,2,2):
                    record.samples[sample]['GT']= generate_geno(strainsNum, 2)
                else:
                    record.samples[sample]['GT']= generate_geno(strainsNum, 1)
                # record.samples[sample]['GT']= generate_geno(strainsNum, 1)
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
    print ("hete loci number",len(snp_list))
    return snp_list, beta_set, allele_set

def generate_geno(strainsNum, geno):
    genotype = []
    for i in range(strainsNum):
        genotype.append(geno)
    return tuple(genotype)

def delta(outdir,extractHAIRS,bamfile,beta_set,reffile):
    hapcut_order='%s --bam %s --VCF %s/filter.vcf --indels 1 --ref %s --out %s/vcf.conn'%(extractHAIRS,bamfile,outdir,reffile,outdir)
    os.system(hapcut_order)
    snp_num=len(beta_set)
    delta_set=[]
    for i in range(snp_num-1):
        delta_set.append([0]*(len(beta_set[i])*len(beta_set[i+1])))
    for line in open('%s/vcf.conn'%(outdir)):
        line=line.strip()
        array=line.split()
        if array[0] == '1': 
            delta_index=int(array[2])-1
            geno_type=array[3]
            for i in range(len(geno_type)-1):
                fol_allele=len(beta_set[delta_index+i+1])
                array_index=int(geno_type[i])*fol_allele+int(geno_type[i+1])
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

def freq_output(outdir, gene, final_alpha, germline_flag):
    if germline_flag:
        final_alpha = [0.5, 0.5]
    ra_file=open(outdir+'/%s_freq.txt'%(gene),'w')    
    print ('# HLA\tFrequency',file=ra_file)
    for j in range(len(final_alpha)):
        print ('str-'+str(j+1),final_alpha[j],file=ra_file)
    ra_file.close()

def output(outdir,final_alpha,seq_list,snp_list,gene,freq_bias):
    # ra_file=open(outdir+'/%s_freq.txt'%(gene),'w')    
    # print ('# HLA\tFrequency',file=ra_file)
    seq=np.array(seq_list)
    seq=np.transpose(seq)
    print (len(snp_list),len(seq),'points num')
    snp=snp_list  
    alpha=final_alpha   
    k=len(alpha) 
    # for j in range(k):
    #     print ('str-'+str(j+1),alpha[j],file=ra_file)
    # ra_file.close()
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
        out.write(record)
    m.close()
    out.close()

def rectify(snp_list,beta_set,delta_set,lambda1,lambda2,germline_flag,database_flag):
    nucleotide={'A':0,'T':1,'C':2,'G':3}
    if database_flag:
        prior_first=first_database()
    has_prior=['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']
    # rectify beta
    first_beta=[]
    for i in range(len(beta_set)):
        #prior freq
        snp=snp_list[i]
        pos=int(snp[1])
        ref=snp[2]
        alt=snp[3][0]
        # print (snp_list[i])
        if len(ref) !=1 or len(alt) != 1 or len(snp[3]) >1 or str(snp[0]) not in has_prior or database_flag==False:
            prior_f=np.array([1.0/len(beta_set[i])]*len(beta_set[i]))
        else:
            tag=str(snp[0])+'_'+str(snp[1]-1000)
            prior=prior_first[tag]
            ref_freq=prior[nucleotide[ref]]
            alt_freq=prior[nucleotide[alt]]
            prior_f=np.array([float(ref_freq),float(alt_freq)])
            if sum(prior_f) > 0:
                prior_f=prior_f/sum(prior_f)
            else:
                prior_f=np.array([1.0/len(beta_set[i])]*len(beta_set[i]))
            # print (snp,prior,ref_freq,alt_freq,prior_f)
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
        if  database_flag==False or str(snp_list[i][0]) not in has_prior or c>0:
            # print (len(hat_delta),len(inde_delta),snp_list[i+1])
            delta=inde_delta*(lambda2/(1+c)) + hat_delta*(1-lambda2/(1+c))  #retify
        else:  #the distance of snp is too far thus we use the link info from database            
            delta=second_prior(snp_list[i],snp_list[i+1])
            if sum(delta) == 0: #no database or loci in different hla gene
                delta=inde_delta            
        # delta=hat_delta  #retify gate
        delta=np.round(delta,6)
        
        #turn the two biggest values to 0.5/0.5, used to handle normal samples
        index_sort=np.argsort(delta)
        if germline_flag :
            index_sort=np.argsort(delta)
            delta[index_sort[0]],delta[index_sort[1]]=0,0
            delta[index_sort[2]],delta[index_sort[3]]=0.5,0.5


        sec_beta.append(delta.tolist())
    return first_beta,sec_beta

def first_database():
    prior_first={}
    for line in open(sys.path[0]+'/'+'database/first_order.database'):
        line=line.strip()
        array=line.split()
        tag=str(array[0])+'_'+str(array[1])
        prior_first[tag]=array[2:]
    return prior_first

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
        
        if int(first[1])-1 in read.get_reference_positions(full_length=True) and read.mapping_quality >10:   
            
            reads_index=read.get_reference_positions(full_length=True).index(int(first[1])-1)
            if first[2][0] != first[3][0]:
                #if the first allele is not same for indel alleles, we can just focus on the first locus 
                if read.query_sequence[reads_index] == first[2][0]:
                    reads_list[0].append(read.query_name)
                elif read.query_sequence[reads_index] == first[3][0]:
                    reads_list[1].append(read.query_name)
            else:  
                index_list=[]
                true_ref=first[4]
                for i in range(len(true_ref)):
                #for i in range(len(first[3])):
                    position=first[1]+i
                    point_flag=isin(position-1,read.get_reference_positions(full_length=True))
                    if point_flag:
                        position_index=read.get_reference_positions(full_length=True).index(position-1)
                        index_list.append(position_index)
                allele_list=read.query_sequence[index_list[0]:index_list[-1]+1].upper()
                ##########for the case that ref is short than alt.#########
                if index_list[-1]+1 < len(read.get_reference_positions(full_length=True)) and len(true_ref) < len(first[2])\
                     and len(true_ref) < len(first[3]) and read.get_reference_positions(full_length=True)[index_list[-1]+1] == None:
                    j = 0
                    # print ('#####', first, allele_list)
                    while index_list[-1]+1+j <  len(read.get_reference_positions(full_length=True)):
                        tag = read.get_reference_positions(full_length=True)[index_list[-1]+1+j]
                        # print (first, tag, type(tag))
                        if read.get_reference_positions(full_length=True)[index_list[-1]+1+j] == None:
                            allele_list+=read.query_sequence[index_list[-1]+1+j]
                        else:
                            break
                        j += 1   
                    # print (first, allele_list)
                ##########for the case that ref is short than alt.#########             
                if allele_list == first[2]:
                    reads_list[0].append(read.query_name)
                elif allele_list == first[3]:
                    reads_list[1].append(read.query_name)
                # if first[1] == 5415:
                #     if 5471-1 in read.get_reference_positions(full_length=True):
                #         new_reads_index=read.get_reference_positions(full_length=True).index(5471-1)
                #         # print (read.query_sequence[new_reads_index])
                #         # print (read.query_name, read.get_reference_positions(full_length=True))
                #         # print (first,allele_list,index_list)
                #     else:
                #         print (first,allele_list,index_list)
                #         print (read.get_reference_positions(full_length=True))
                #         print (read.get_reference_positions(full_length=False))

    # if int(first[1]) == 5415 or int(first[1]) == 5471:
    #    print (len(reads_list[0]),len(reads_list[1]),first)
    #    print (reads_list[0],reads_list[1],first)
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
            # if str(left[1]) == '6983' and i == 0 and j == 1:
            #     print (left, len(set(left_set).intersection(set(right_set))), set(left_set).intersection(set(right_set)))
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

def relate_order(file, snp_list):
    sv_locus = []
    locus_list = []
    n = 0
    block_dict={}
    dict={}
    previous_gene=''
    
    for line in open(file,'r'):
        line=line.strip()
        array=line.split()
        if array[0] == 'gene':
            past_order = []
            for j in range(len(array[6:])):
                past_order.append(j)
            continue
        old_array = array[:]
        gene_name=array[0]
        locus=array[1]
        array=array[6:]

        hla_num=len(array)
        previous_gene=gene_name            
        ref_order=[]
        for j in range(hla_num):
            ref_order.append(j)
        if len(array)<len(ref_order):
            continue

        if old_array[2] == '+':
            sv_locus.append(n)
        # print (locus)
        locus_list.append(int(locus))
        new_order=[]
        # print (hla_num,ref_order,array)
        ##################
        # for e in range(hla_num):
        #     new_order.append(array[int(past_order[e])])        
        # ref_order=new_order[:]
        ref_order = array
        ##################
        past_order = ref_order
        dict[locus]=ref_order
        n += 1
    #update the locus
    replace_dict = {}
    for sv in sv_locus:
        if sv + 1 < len(locus_list):
            gap = [locus_list[sv], locus_list[sv+1]]
        else:
            gap = [locus_list[sv], 'end']
        replace_dict[str(locus_list[sv])] = str(first_snp_in_the_region(gap, snp_list))
    # print (dict)
    # print ('replace_dict', replace_dict)
    for key in replace_dict.keys():
        if key in dict.keys():
            dict[replace_dict[key]] = dict[key]
            del dict[key]
            # print ()

    block_dict[previous_gene] = dict
    # print ('final dict', block_dict)
    return block_dict

def first_snp_in_the_region(gap, snp_list):
    for snp in snp_list:
        snp[1] = int(snp[1])
        if gap[1] == 'end':
            if snp[1] >= gap[0]:
                return snp[1]
        else:
            if snp[1] >= gap[0] and snp[1] < gap[1]:
                return snp[1]
    return gap[0]

def relate_order_other(file, snp_list):
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
    file=outdir+'/%s_break_points_phased.txt'%(gene)
    if os.path.isfile(file):
        if gene == 'HLA_DRB1':
            block_dict=relate_order(file, snp_list)
        else:
            block_dict=relate_order_other(file, snp_list)
    else:
        block_dict={gene:{}}
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
        out.write(record)
        if str(record.pos) in gene_dict.keys():  #rephase new order
            ref_order=gene_dict[str(record.pos)]
            print ('new order',record.pos, ref_order)

    m.close()
    out.close()
    os.system('tabix %s/%s.rephase.vcf.gz'%(outdir,gene))
    return update_seqlist

def newalpha(update_seqlist, sec_beta, strainsNum, allele_set, weight, fir_beta):
    return alpha_step(sec_beta, update_seqlist, strainsNum, allele_set, weight, fir_beta)

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
    in_vcf = VariantFile(vcffile)
    out_vcf = VariantFile('%s/%s.rephase.vcf.gz'%(outdir, gene),'w',header=in_vcf.header)
    sample = list(in_vcf.header.samples)[0]
    for record in in_vcf.fetch():
        if record.chrom == gene and record.samples[sample]['GT'] == (1, 1, 1):
            # print (gene, record)
            phased_locus=[1] * strainsNum
            record.samples[sample]['GT']= tuple(phased_locus)
            record.samples[sample].phased=True
            out_vcf.write(record)
    in_vcf.close()
    out_vcf.close()
    os.system('tabix %s/%s.rephase.vcf.gz'%(outdir,gene))
    ra_file=open(outdir+'/%s_freq.txt'%(gene),'w')    
    print ('# HLA\tFrequency',file=ra_file)
    print ('str-'+str(1), 1, file=ra_file)
    for j in range(1, strainsNum):
        print ('str-'+str(j+1), 0, file=ra_file)
    ra_file.close()

    ####
    gene_profile={}
    gene_name=['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']
    for gene in gene_name:
        gene_snp=[]
        gene_seq=[]
        gene_profile[gene] = [gene_snp,gene_seq]
    return gene_profile

def read_dup():
    dup_dict={}
    for line in open(sys.path[0]+'/complex_region.txt','r'):
        line=line.strip()
        array=line.split()
        dup_dict[array[0]]=array[1:]
    return dup_dict

########align the segs from SV haps result
def isOut(index, myset):
    if index < len(myset):
        return myset[index]
    else:
        return 'NA'

def exists_index(seg_set):
    exists_dict = {}
    for i in range(len(seg_set)):
        for j in range(1000):
            if isOut(j, seg_set[i]) != 'NA':
                if seg_set[i][j] not in exists_dict.keys():
                    exists_dict[seg_set[i][j]] = [j]
                else:
                    exists_dict[seg_set[i][j]].append(j)
    for key in exists_dict.keys():
        exists_dict[key] = sorted(exists_dict[key])
    #to ensure the following segs index is larger than previous segs
    newflag = True
    while newflag:
        newflag = False
        for i in range(len(seg_set)):
            for j in range(1000):
                if isOut(j, seg_set[i]) != 'NA' and isOut(j+1, seg_set[i]) != 'NA':
                    if max(exists_dict[seg_set[i][j+1]]) <= max(exists_dict[seg_set[i][j]]) or\
                        min(exists_dict[seg_set[i][j+1]]) <= min(exists_dict[seg_set[i][j]]) :
                            newflag = True
                            for w in range(len(exists_dict[seg_set[i][j+1]])):
                                exists_dict[seg_set[i][j+1]][w] = exists_dict[seg_set[i][j+1]][w] + 1
    return exists_dict

def back_up_sort_segs(seg_set):
    exists_dict = exists_index(seg_set)
    segs_list = list(exists_dict.keys())
    flag = True
    while flag:
        flag = False
        for i in range(len(segs_list)-1):
            front_seg = segs_list[i]
            after_seg = segs_list[i + 1]
            if min(exists_dict[after_seg]) <= min(exists_dict[front_seg]) and \
                max(exists_dict[after_seg]) <= max(exists_dict[front_seg]):
                segs_list[i] = after_seg
                segs_list[i + 1] = front_seg
                flag = True
    return segs_list

def sort_segs(seg_set):
    exists_dict = exists_index(seg_set)
    segs_list = list(exists_dict.keys())
    new_segs_list = []
    for seg in segs_list:
        new_segs_list.append(int(float(seg[:-1])))
    segs_list = []
    for seg in sorted(new_segs_list):
        segs_list.append(str(seg) + '+')

    return segs_list

def alignments_seg(seg_set):
    strainNum = len(seg_set)
    segs_list = sort_segs(seg_set)
    new_seg_set = []
    for i in range(strainNum):
        new_seg_set.append([])
    for j in range(len(segs_list)):
        for i in range(strainNum):
            if segs_list[j] in seg_set[i]:
                new_seg_set[i].append(segs_list[j])
            else:
                new_seg_set[i].append('*')

    # print (segs_list, new_seg_set)
    return new_seg_set
########align the segs from SV haps result

class AddSV():
    def __init__(self,balance_lh,hap,vcffile,strainsNum,ins):
        self.balance_lh,self.hap,self.vcffile = balance_lh,hap,vcffile 
        self.ins = ins 
        self.strainsNum=strainsNum
        self.copy_num,self.index_locus,self.source,self.sink=self.read_balanced()
        self.segs=list(self.copy_num.keys())  
        self.ins_seq = self.read_ins()

    def read_balanced(self):
        #read the copy number of each seg, and record multicopy segs and corresponding copy number.
        copy_num={}
        index_locus={}
        for line in open(self.balance_lh,'r'):
            line=line.strip()
            array=line.split()
            if str(array[0]) == 'SOURCE':
                locus=array[1].split(':')
                source = locus[1]
            elif str(array[0]) == 'SINK':
                locus=array[1].split(':')
                sink = locus[1]
            elif str(array[0]) == 'SEG':
                locus=array[1].split(':')
                index_locus[locus[1]] = [int(locus[3]),int(locus[4])]
                copy_num[locus[1]] = float(array[3])
        return copy_num,index_locus,source,sink

    def if_normal(self):  #regard the seg with one copy in each hap as normal segs, abnormal otherwise.
        line=open(self.hap,'r').readline().strip()
        seg_order=[]
        array = line.split()
        hap_num = 0
        for i in range(len(array)):
            if array[i] == self.source + '+':
                new_hap = array[i] + ' '
                if self.source == self.sink:
                    seg_order.append(new_hap)
                    hap_num += 1
            elif array[i] == self.sink + '+':
                new_hap = new_hap + array[i] + ' '
                seg_order.append(new_hap)
                hap_num += 1
            else:
                new_hap = new_hap + array[i] + ' '
            if hap_num == self.strainsNum:
                break
        normal_seg=[]
        abnormal_seg=[]
        for seg in self.segs:
            num=0
            fresh_copy=0  #count the copy number once again.
            for arr in seg_order:
                if str(seg)+'+' in arr.split() or str(seg)+'-' in arr.split():
                    num+=1
                for mysegs in arr.split():
                    if str(seg)+'+' == mysegs or str(seg)+'-' == mysegs:
                        fresh_copy += 1
            # if num == self.strainsNum and self.copy_num[seg] == self.strainsNum:
            #     normal_seg.append(seg)
            if num == self.strainsNum and fresh_copy == self.strainsNum:
                normal_seg.append(seg)
            else:
                abnormal_seg.append(seg)
        # print ('if normal', normal_seg, abnormal_seg)

        #insertion locus
        # ins_locus = {}
        # for i in range(len(array)):
        #     if array[i][:-1] in self.ins_seq.keys() and array[i][:-1] not in ins_locus.keys():
        #         ins_locus[array[i][:-1]] = self.index_locus[array[i-1][:-1]][1]
        # print ('insert',ins_locus)

        return seg_order,normal_seg,abnormal_seg,self.copy_num,self.index_locus,self.ins_seq

    def read_seg(self):
        index_locus={}       
        i=0
        for line in open(self.seg_file,'r'):
            if line[0] == 'I':
                continue
            line=line.strip()
            array=line.split()
            if i == 0:
                fir_extend=array[0]
            sec_extend=array[0]
            index_locus[array[0]] = [array[2],array[3]]
            i+=1
        # del index_locus[fir_extend]
        # del index_locus[sec_extend]
        return index_locus

    def read_ins(self):
        ins_seq={}
        if os.path.isfile(self.ins):
            for line in open(self.ins,'r'):
                if line[0] == 'I':
                    continue
                array=line.strip().split()
                seg_ID = array[1].split('_')
                ins_seq[seg_ID[1]] = array[-1]
        else:
            print ('no insertion for this gene')
        return ins_seq

    def deletion_locus(self):
        deletion_region = []
        for seg_ID in self.index_locus.keys():
            if float(self.copy_num[seg_ID]) < self.strainsNum and seg_ID not in self.ins_seq.keys():
                deletion_region.append(self.index_locus[seg_ID])
                # print ('deletion', seg_ID, self.index_locus[seg_ID])    
        return deletion_region

class Seg_phase():
    def __init__(self,strainsNum,vcffile,outdir,snp_dp,bamfile,indel_len,gene_profile,gene,dup_file,freq_bias,ins_seq):
        self.gene, self.freq_bias = gene, freq_bias
        self.outdir,self.snp_dp=outdir,snp_dp
        self.bamfile,self.indel_len=bamfile,indel_len 
        self.vcffile,self.strainsNum=vcffile,strainsNum
        self.normal_sequence=gene_profile[self.gene]
        self.reads_support=self.normal_reads() 
        if len(ins_seq) == 0 :
            self.new_bam = bamfile
            self.new_vcf = vcffile
        else:
            self.new_bam = self.outdir + '/newref_insertion.bam' #ins_bam
            self.new_vcf = self.outdir + '/newref_insertion.freebayes.vcf' #ins_vcf
        self.new_reads_support = self.new_normal_reads()
        self.dup_file = dup_file

    def normal_reads(self):  
        # for the reads support normal segs, divide into two categaries, the reads support the first 
        #hap, and the reads support the second hap  
        reads_support = []
        for i in range(self.strainsNum):
            reads_support.append([])
        samfile = pysam.AlignmentFile(self.bamfile, "rb")
        for seg in self.normal_seg:
            for read in samfile.fetch(self.gene,int(self.index_locus[seg][0])-1,\
                int(self.index_locus[seg][1])):
                rivet_points=False
                support_alleles=[]
                support_loci=[]
                for i in range(len(self.normal_sequence[0])):
                    snv=self.normal_sequence[0][i]
                    if len(snv[2]) != 1 or len(snv[3]) != 1:
                        continue
                    if int(snv[1])-1 in read.get_reference_positions(full_length=True):
                        reads_index=read.get_reference_positions(full_length=True).index(int(snv[1])-1)
                        support_alleles.append(read.query_sequence[reads_index])
                        # print (read.query_name,snv,support_alleles)
                        rivet_points=True
                        support_loci.append(i)
                if rivet_points==True:
                    #if the reads has information, check which hap it belongs to.
                    hap_belong=self.check_hap(support_alleles,support_loci)[0]
                    if hap_belong != 'NA':
                        reads_support[hap_belong].append(read.query_name)
        return reads_support

    def new_normal_reads(self):  
        # for the reads support normal segs, divide into two categaries, the reads support the first 
        #hap, and the reads support the second hap  
        reads_support = []
        for i in range(self.strainsNum):
            reads_support.append([])
        samfile = pysam.AlignmentFile(self.new_bam, "rb")
        for seg in self.normal_seg:
            for read in samfile.fetch(self.gene,int(self.index_locus[seg][0])-1,\
                int(self.index_locus[seg][1])):
                rivet_points=False
                support_alleles=[]
                support_loci=[]
                for i in range(len(self.normal_sequence[0])):
                    snv=self.normal_sequence[0][i]
                    if len(snv[2]) != 1 or len(snv[3]) != 1:
                        continue
                    if int(snv[1])-1 in read.get_reference_positions(full_length=True):
                        reads_index=read.get_reference_positions(full_length=True).index(int(snv[1])-1)
                        support_alleles.append(read.query_sequence[reads_index])
                        # print (read.query_name,snv,support_alleles)
                        rivet_points=True
                        support_loci.append(i)
                if rivet_points==True:
                    # if len(support_loci) < 3:
                    #     continue
                    #if the reads has information, check which hap it belongs to.
                    hap_belong=self.check_hap(support_alleles,support_loci)[0]
                    if hap_belong != 'NA':
                        reads_support[hap_belong].append(read.query_name)
        return reads_support

    def most_support(self,support_num):
        hap_order = np.argsort(np.array(support_num))[::-1]
        return hap_order

    def sort_support(self,support_num):       
        return sorted(support_num,reverse = True)

    def check_hap(self,support_alleles,support_loci):
        support_num=[]  #check the allele num that support the each hap respectively.
        for i in range(self.strainsNum):
            support_num.append(0)
        for i in range(len(support_loci)):
            locus_index=support_loci[i]
            allele=support_alleles[i]
            snv=self.normal_sequence[0][locus_index]
            for j in range(self.strainsNum):
                hap_allele=snv[self.normal_sequence[1][j][locus_index]+2]
                if hap_allele == allele:
                    support_num[j] += 1
        #check which hap has most same alleles support_num.index(max(support_num))
        return self.most_support(support_num)

    def ins_check_hap(self,support_alleles,support_loci,main_snp,main_strain):
        support_num=[]  #check the allele num that support the each hap respectively.
        for i in range(self.strainsNum):
            support_num.append(0)
        for i in range(len(support_loci)):
            locus_index=support_loci[i]
            allele=support_alleles[i]
            snv=main_snp[locus_index]
            for j in range(self.strainsNum):
                hap_allele=snv[main_strain[j][locus_index]+2]
                if hap_allele == allele:
                    support_num[j] += 1
        return self.most_support(support_num)

    def abnormal_reads(self,seg_ID):
        link_reads=[]  #the reads number that shared with different haps, the copy number may be 0
        for i in range(self.strainsNum):
            link_reads.append(0)
        samfile = pysam.AlignmentFile(self.bamfile, "rb")
        for read in samfile.fetch(self.gene,int(self.index_locus[seg_ID][0])+1,\
            int(self.index_locus[seg_ID][1])-1):
            for i in range(self.strainsNum):
                if read.query_name in self.reads_support[i]:
                    link_reads[i] += 1
        return self.most_support(link_reads)

    def insertion_single_reads(self, seg_ID):
        link_reads=[]
        for i in range(self.strainsNum):
            link_reads.append(0)
        samfile = pysam.AlignmentFile(self.new_bam, "rb")
        seg_name = '%s_%s'%(self.gene, seg_ID)
        # seg_name = '%s_%s_%s'%(self.gene, self.index_locus[seg_ID][0], self.index_locus[seg_ID][1])
        for read in samfile.fetch(seg_name, 0, 10000):
            for i in range(self.strainsNum):
                # if read.query_name in self.reads_support[i]:
                if read.query_name in self.new_reads_support[i]:
                    link_reads[i] += 1
        return self.most_support(link_reads)

    def consensus(self):
        consensus_order = [[0, 0], [0, 0]]   #snv <==> sv
        consensus_order = []
        for i in range(self.strainsNum):
            consensus_order.append([0] * self.strainsNum)
        print ('consensus_order',consensus_order)
        for seg_ID in self.abnormal_seg:
            if float(self.copy_num[seg_ID]) < self.strainsNum and seg_ID not in self.ins_seq.keys():
                #undelete2snv indicate which snv hap to link for undeleted seg 
                #the undeleted seg link to the x-th hap        
                undelete2snv_set=self.abnormal_reads(seg_ID) 
                # if undelete2snv != 'NA':
                copy_number=int(self.copy_num[seg_ID])
                for j in range(copy_number):
                    undelete2snv = undelete2snv_set[j]
                    for i in range(self.strainsNum):
                        if seg_ID+'+' in self.seg_order[i]:
                            consensus_order[undelete2snv][i] += 1
                    print ('undeletion', seg_ID, undelete2snv)
        for seg_ID in self.abnormal_seg:
            copy_number=int(self.copy_num[seg_ID])
            if float(copy_number) < self.strainsNum and seg_ID in self.ins_seq.keys():
                #ins2snv indicate which snv hap to link for insertion seg  
                ins2snv_set = self.insertion_single_reads(seg_ID)               
                for j in range(int(self.copy_num[seg_ID])):
                    ins2snv = ins2snv_set[j]
                    for i in range(self.strainsNum):
                        if seg_ID+'+' in self.seg_order[i]:
                            consensus_order[ins2snv][i] += 1
        print ('consensus_order',consensus_order,correlation(consensus_order))
        return correlation(consensus_order)

    def insertion_phase(self, seg_ID):
        # for a insertion existes in both haps, phase it
        samfile = pysam.AlignmentFile(self.new_bam, "rb")
        seg_name = '%s_%s'%(self.gene, seg_ID)
        # seg_name = '%s_%s_%s'%(self.gene, self.index_locus[seg_ID][0], self.index_locus[seg_ID][0])
        print ('insertion segname',seg_name)
        seg_snp_list, beta_set, allele_set = read_vcf(self.new_vcf,self.outdir,self.snp_dp,\
            self.new_bam,self.indel_len,seg_name,self.freq_bias,self.strainsNum,[])        
        if len(seg_snp_list) == 0:
            print ('no snv for this insertion')
            insertion_seq_double = []
            for i in range(self.strainsNum):
                insertion_seq_double.append(self.ins_seq[seg_ID])                
            return insertion_seq_double

        delta_set = second_beta(self.new_bam, seg_snp_list)
        # for i in range(len(delta_set)):
        #     print ('insertion phase', seg_snp_list[i], delta_set[i])
        fir_beta, sec_beta = rectify(seg_snp_list, beta_set, delta_set, args.lambda1,\
            args.lambda2, germline_flag, database_flag)
        wo = Workflow(fir_beta, sec_beta, delta_set, args.weight, args.elbow, allele_set)
        main_alpha, main_strain, main_loss = wo.given_k(int(self.strainsNum))  
        main_snp = seg_snp_list
        print ('insertion phase',len(main_snp),main_alpha)
        #transform genotypes to sequence
        phased_seqs = self.transform(main_snp, main_strain, seg_ID, main_alpha)
        #for the phased haps, split the reads that support each hap.  
        ins_reads_support=[]
        for i in range(self.strainsNum):
            ins_reads_support.append([])
        for read in samfile.fetch(seg_name, 0, 10000):
            rivet_points=False
            support_alleles=[]
            support_loci=[]
            for i in range(len(main_strain[0])):
                snv=main_snp[i]
                if len(snv[2]) != 1 or len(snv[3]) != 1:
                    continue
                if int(snv[1])-1 in read.get_reference_positions(full_length=True):
                    reads_index=read.get_reference_positions(full_length=True).index(int(snv[1])-1)
                    support_alleles.append(read.query_sequence[reads_index])
                    rivet_points=True
                    support_loci.append(i)
            if rivet_points==True:
                #if the reads has information, check which hap it belongs to.
                if len(support_loci) < 4:   #why ??? the reason may be the phase result of the insertion
                    continue                #is too worse
                hap_belong=self.ins_check_hap(support_alleles,support_loci,main_snp,main_strain)[0]
                if hap_belong != 'NA':
                    ins_reads_support[hap_belong].append(read.query_name)
        share_num = []
        for i in range(self.strainsNum):
            num_set = []
            for j in range(self.strainsNum):
                num = 0
                for re in ins_reads_support[i]:
                    # if re in self.reads_support[j]:
                    if re in self.new_reads_support[j]:
                        num+=1
                num_set.append(num)
            share_num.append(num_set)
        snv_ins_order = correlation(share_num)
        # snv_ins_order = (0, 1)
        print ('double insertion!', share_num, snv_ins_order)
        new_phased_seqs = []
        for i in range(self.strainsNum):
            new_phased_seqs.append(phased_seqs[snv_ins_order[i]])
        return new_phased_seqs

    def insertion_double(self):
        dou_add_dict = {}
        for seg_ID in self.ins_seq.keys():
            if seg_ID not in self.abnormal_seg:
                dou_add_dict[seg_ID] = self.insertion_phase(seg_ID)
        return dou_add_dict

    def dup_assign(self):       
        #uniq reads for each dup type
        drb1_complex_seq, uniq_drb1_complex_reads = DRB1_complex_region(self.dup_file)
        #the relation between dup type and snv-haplotype
        share_num = []
        new_drb1_complex_seq = []
        for i in range(self.strainsNum):
            max_num = 0
            max_seq = ''
            for j in range(len(drb1_complex_seq)):
                num = 0
                for re in uniq_drb1_complex_reads[j]:
                    if re in self.reads_support[i]:
                        num+=1
                if num >= max_num:
                    max_num = num
                    max_seq = drb1_complex_seq[j]
            new_drb1_complex_seq.append(max_seq)
        return new_drb1_complex_seq
   
    def transform(self, main_snp, main_strain, seg_ID, main_alpha):
        #transform genotype to sequence 
        seg_name = '%s_%s'%(self.gene, seg_ID)
        output(self.outdir,main_alpha,main_strain,main_snp,seg_name,self.freq_bias)
        out_vcf = '%s/%s.vcf'%(outdir,seg_name)
        os.system('bgzip -f %s\ntabix %s.gz'%(out_vcf,out_vcf))
        insertion_seqs = []
        for i in range(self.strainsNum):
            order = 'samtools faidx %s/newref_insertion.fa\
            %s_%s |/home/wangmengyao/miniconda2/bin/bcftools\
            consensus -H %s %s.gz  >%s/insert_%s.fa'%(self.outdir,\
            self.gene, seg_ID, i+1, out_vcf, self.outdir, i)
            os.system(order)
            seq = read_fasta('%s/insert_%s.fa'%(self.outdir, i))
            insertion_seqs.append(seq)
        return insertion_seqs

    def split_seg(self):
        print ('start split seg.')
        cons = self.consensus()
        print ('consensus is over.', cons)
        dou_add_dict = self.insertion_double()
        if self.gene == 'HLA_DRB1':
            my_drb1_complex_seq = self.dup_assign()
        gap = ''
        id_name = {}
        for seg in self.index_locus.keys():
            print ('seg', seg)
            if seg not in self.ins_seq.keys():
                #to split the dup region and replace the dup to phased one.
                if self.gene == 'HLA_DRB1' and self.index_locus[seg][0] < 3898 and \
                self.index_locus[seg][1] > 4400:
                    contain_dup_seg = seg
                    seg_region = ' %s:%s-3898 %s:3898-4400 %s:4400-%s '%(self.gene,\
                    self.index_locus[seg][0],self.gene,self.gene,self.index_locus[seg][1])
                    seg_region_front = ' %s:%s-3898 '%(self.gene, self.index_locus[seg][0])
                    seg_region_dup = ' %s:3898-4400 '%(self.gene)
                    seg_region_behind = ' %s:4400-%s '%(self.gene,self.index_locus[seg][1])
                    id_name[seg_region_front.strip()] = str(seg)+'>'
                    id_name[seg_region_dup.strip()] = str(seg)+'='
                    id_name[seg_region_behind.strip()] = str(seg)+'<'
                else:
                    seg_region = ' %s:%s-%s '%(self.gene,self.index_locus[seg][0],\
                    self.index_locus[seg][1])
                    id_name[seg_region.strip()] = seg
                gap += seg_region
        for i in range(self.strainsNum):
            # order='samtools faidx /home/wangmengyao/scripts/NeedleHLA/ref/hla.ref.extend.fa\
            order='samtools faidx /home/wangmengyao/scripts/PHLAT/database/ref/extend/hla.ref.extend.fa\
                %s |/home/wangmengyao/miniconda2/bin/bcftools\
                consensus -H %s %s/%s.rephase.vcf.gz  >%s/%s_%s_seg.fa'%(gap,i+1,self.outdir,\
                self.gene,self.outdir,self.gene,i)
            print (order)
            os.system(order)
            fa_file = '%s/%s_%s_seg.fa'%(self.outdir,self.gene,i)
            seg_sequence = chrom_seq(fa_file)
            #covert seg name
            new_seg_sequence = {}
            for segseq in seg_sequence.keys():
                new_seg_sequence[id_name[segseq]] = seg_sequence[segseq]
            # print (new_seg_sequence)
            sv_result = self.align_seq[cons[i]]
            hap_seq = ''
            #hap_seq = '>%s_%s\n'%(self.gene, i)
            # print (new_seg_sequence)
            for seg_ID in sv_result:
                print ('The seg %s in the %s strain! Merging!'%(seg_ID, i))
                if seg_ID == '*':
                    continue
                seg_ID = seg_ID[:-1]
                if seg_ID in new_seg_sequence.keys():
                    hap_seq += new_seg_sequence[seg_ID]
                elif seg_ID in dou_add_dict.keys():
                    hap_seq += dou_add_dict[seg_ID][i]
                    print ('??', i, len(dou_add_dict[seg_ID][i]))
                    # hap_seq += dou_add_dict[seg_ID][1]
                    # continue
                # elif float(self.copy_num[seg_ID]) < self.strainsNum and seg_ID in self.ins_seq.keys():
                elif seg_ID in self.ins_seq.keys():
                    hap_seq += self.ins_seq[seg_ID]
                else:
                    print ('the seg contain dup region')
                    # the seg that contain the dup region
                    hap_seq += new_seg_sequence[str(seg_ID)+'>']
                    #add the seq of chosen dup type
                    hap_seq += my_drb1_complex_seq[i]
                    # the seg that contain the dup region
                    hap_seq += new_seg_sequence[str(seg_ID)+'<']
            #delete front and behind 1000bp
            hap_seq =  '>%s_%s\n'%(self.gene, i) + hap_seq[1001:-1000]
            out = open('%s/hla.allele.%s.%s.fasta'%(self.outdir,i+1,self.gene), 'w')
            print (hap_seq, file = out)
            out.close()

def focus_region():
    return {'HLA_A':[1000,4503],'HLA_B':[1000,5081],'HLA_C':[1000,5304],'HLA_DPA1':[1000,10775],\
        'HLA_DPB1':[1000,12468],'HLA_DQA1':[1000,7492],'HLA_DQB1':[1000,8480],'HLA_DRB1':[1000,12229]}

class Share_reads():

    def __init__(self, deletion_region, outdir, strainsNum, gene, gene_profile, ins_seq, dup_file):
        self.deletion_region = deletion_region
        self.bamfile = outdir + '/newref_insertion.bam'
        self.vcf = outdir + '/newref_insertion.freebayes.vcf'
        self.strainsNum = strainsNum
        self.gene = gene
        self.normal_sequence=gene_profile[self.gene]
        self.reads_support = self.normal_reads()
        self.outdir = outdir
        self.ins_seq = ins_seq
        self.dup_file = dup_file
        
    def generate_normal_region(self):
        gene_area = focus_region()[self.gene]
        normal_region = []
        segs = []
        gene_area[0] = int(float(gene_area[0])) + 1
        gene_area[1] = int(float(gene_area[1]))
        print ('gene_area', gene_area)
        start = gene_area[0]
        print (self.deletion_region)
        for i in range(len(self.deletion_region)):
            if self.deletion_region[i][1] > gene_area[1]:
                self.deletion_region[i][1] = gene_area[1]
            if start < self.deletion_region[i][0]:
                segs.append([start, self.deletion_region[i][0] - 1, 'normal', '.'])
            if self.deletion_region[i][0] == self.deletion_region[i][1]:
                segs.append([self.deletion_region[i][0], self.deletion_region[i][1], 'insertion', i])
                # continue
            else:
                segs.append([self.deletion_region[i][0], self.deletion_region[i][1], 'deletion', i])
            normal_region.append([start, self.deletion_region[i][0]])
            start = self.deletion_region[i][1] + 1
        if int(start) < gene_area[1]:
            normal_region.append([start, gene_area[1]])
            segs.append([start, gene_area[1], 'normal', '.'])
        #print (self.deletion_region)
        #print (segs)
        print ('############normalregion', normal_region, segs)
        return normal_region, segs

    def normal_reads(self):
        normal_region, segs = self.generate_normal_region()
        print (normal_region)
        reads_support = []
        for i in range(self.strainsNum):
            reads_support.append([])
        samfile = pysam.AlignmentFile(self.bamfile, "rb")
        for region in normal_region:
            if abs(region[0] - region[1]) < 10:
                continue
            for read in samfile.fetch(self.gene, region[0] - 1, region[1]):
                rivet_points=False
                support_alleles=[]
                support_loci=[]
                for i in range(len(self.normal_sequence[0])):
                    snv=self.normal_sequence[0][i]
                    if len(snv[2]) != 1 or len(snv[3]) != 1:
                        continue
                    if int(snv[1])-1 in read.get_reference_positions(full_length=True):
                        reads_index=read.get_reference_positions(full_length=True).index(int(snv[1])-1)
                        support_alleles.append(read.query_sequence[reads_index])
                        # print (read.query_name,snv,support_alleles)
                        rivet_points=True
                        support_loci.append(i)
                if rivet_points==True:
                    #if the reads has information, check which hap it belongs to.
                    hap_belong=self.check_hap(support_alleles,support_loci)[0]
                    if hap_belong != 'NA':
                        reads_support[hap_belong].append(read.query_name)
        print (len(reads_support[0]), len(reads_support[1]))
        return reads_support

    def check_hap(self,support_alleles,support_loci):
        support_num=[]  #check the allele num that support the each hap respectively.
        for i in range(self.strainsNum):
            support_num.append(0)
        for i in range(len(support_loci)):
            locus_index=support_loci[i]
            allele=support_alleles[i]
            snv=self.normal_sequence[0][locus_index]
            for j in range(self.strainsNum):
                hap_allele=snv[self.normal_sequence[1][j][locus_index]+2]
                if hap_allele == allele:
                    support_num[j] += 1
        #check which hap has most same alleles support_num.index(max(support_num))
        return self.most_support(support_num)

    def most_support(self,support_num):
        hap_order = np.argsort(np.array(support_num))[::-1]
        return hap_order

    def deletion_reads(self,deletion_index):
        link_reads=[]  #the reads number that shared with different haps, the copy number may be 0
        for i in range(self.strainsNum):
            link_reads.append(0)
        samfile = pysam.AlignmentFile(self.bamfile, "rb")
        print ('###', link_reads)
        for read in samfile.fetch(self.gene,float(self.deletion_region[deletion_index][0])-1,float(self.deletion_region[deletion_index][1])):
            #print (deletion_index, read.query_name, float(self.deletion_region[deletion_index][0])-1, float(self.deletion_region[deletion_index][1]))
            for i in range(self.strainsNum):
                if read.query_name in self.reads_support[i]:
                    #print('####', deletion_index, read.query_name, i, 'link_reads')
                    link_reads[i] += 1
        print ('###', link_reads)
        print ('hihihi', deletion_index, link_reads)
        return self.most_support(link_reads)

    def insertion_reads(self,deletion_index):
        link_reads=[]  #the reads number that shared with different haps, the copy number may be 0
        for i in range(self.strainsNum):
            link_reads.append(0)
        samfile = pysam.AlignmentFile(self.bamfile, "rb")
        for read in samfile.fetch('%s_%s'%(self.gene,self.deletion_region[deletion_index][0])):
            for i in range(self.strainsNum):
                if read.query_name in self.reads_support[i]:
                    link_reads[i] += 1
        print ('insertion', self.most_support(link_reads))
        return self.most_support(link_reads)

    def deletion_phase(self):
        for deletion_index in range(len(self.deletion_region)):
            print (self.deletion_region[deletion_index])
            self.deletion_reads(deletion_index)

    def dup_assign(self):       
        #uniq reads for each dup type
        drb1_complex_seq, uniq_drb1_complex_reads = DRB1_complex_region(self.dup_file)
        #the relation between dup type and snv-haplotype
        share_num = []
        new_drb1_complex_seq = []
        for i in range(self.strainsNum):
            max_num = 0
            max_seq = ''
            for j in range(len(drb1_complex_seq)):
                num = 0
                for re in uniq_drb1_complex_reads[j]:
                    if re in self.reads_support[i]:
                        num+=1
                if num >= max_num:
                    max_num = num
                    max_seq = drb1_complex_seq[j]
            new_drb1_complex_seq.append(max_seq)
        return new_drb1_complex_seq

    def split_seg(self):
        print ('start split seg.')
        if self.gene == 'HLA_DRB1':
            my_drb1_complex_seq = self.dup_assign()
        normal_region, segs = self.generate_normal_region()
        id_name = {}
        gap = ''
        contain_dup_seg = ''
        for seg in segs:
            print ('seg', seg)
            if seg[2] == 'insertion':
                continue
            if self.gene == 'HLA_DRB1' and seg[0] < 3898 and seg[1] > 4400:
                contain_dup_seg = seg
                seg[2] = 'dup'
                seg_region = ' %s:%s-3898 %s:3898-4400 %s:4400-%s '%(self.gene,\
                seg[0],self.gene,self.gene,seg[1])
                seg_region_front = ' %s:%s-3898 '%(self.gene, seg[0])
                seg_region_dup = ' %s:3898-4400 '%(self.gene)
                seg_region_behind = ' %s:4400-%s '%(self.gene,seg[1])
                id_name[seg_region_front.strip()] = str(seg[0])+'>'
                id_name[seg_region_dup.strip()] = str(seg[0])+'='
                id_name[seg_region_behind.strip()] = str(seg[1])+'<'
            else:
                seg_region = ' %s:%s-%s '%(self.gene,seg[0],seg[1])
                id_name[seg_region.strip()] = str(seg[0]) + '_' + str(seg[1]) 
            gap += seg_region
        for i in range(self.strainsNum):
            # order='samtools faidx /home/wangmengyao/scripts/NeedleHLA/ref/hla.ref.extend.fa\
            order='samtools faidx /home/wangmengyao/scripts/PHLAT/database/ref/extend/hla.ref.extend.fa\
                %s |/home/wangmengyao/miniconda2/bin/bcftools\
                consensus -H %s %s/%s.rephase.vcf.gz  >%s/%s_%s_seg.fa'%(gap,i+1,self.outdir,\
                self.gene,self.outdir,self.gene,i)
            print (order)
            os.system(order)
            fa_file = '%s/%s_%s_seg.fa'%(self.outdir,self.gene,i)
            seg_sequence = chrom_seq(fa_file)
            new_seg_sequence = {}
            for segseq in seg_sequence.keys():
                new_seg_sequence[id_name[segseq]] = seg_sequence[segseq]
            print (new_seg_sequence.keys())
            hap_seq = ''
            for seg in segs:
                print (seg, str(seg[0]) + '_' + str(seg[1]))
               # print (str(seg[0]) + '_' + str(seg[1]) )
                if seg[2] == 'normal':
                    hap_seq += new_seg_sequence[str(seg[0]) + '_' + str(seg[1]) ]
                elif seg[2] == 'deletion':
                    deletion_index = seg[3]
                    assign_index = self.deletion_reads(deletion_index)[0]
                    print ('deletion assign_index', assign_index)
                    if  i == assign_index and self.deletion_region[deletion_index][2] == 1:
                        hap_seq += new_seg_sequence[str(seg[0]) + '_' + str(seg[1]) ]
                elif seg[2] == 'insertion':
                    # continue
                    deletion_index = seg[3]
                    assign_index = self.insertion_reads(deletion_index)[0]
                    print ('###########insertion', seg, assign_index)
                    if self.deletion_region[deletion_index][2] == 2:
                        print ('###########insertion two haplotype', seg)
                        hap_seq += self.ins_seq[seg[0]]
                    elif  i == assign_index and self.deletion_region[deletion_index][2] == 1:
                        hap_seq += self.ins_seq[seg[0]]
                elif seg[2] == 'dup':
                    print ('the seg contain dup region')
                    # the seg that contain the dup region
                    hap_seq += new_seg_sequence[str(seg[0])+'>']
                    #add the seq of chosen dup type
                    hap_seq += my_drb1_complex_seq[i]
                    # the seg that contain the dup region
                    hap_seq += new_seg_sequence[str(seg[1])+'<']                    
            hap_seq =  '>%s_%s\n'%(self.gene, i) + hap_seq[:]
            out = open('%s/hla.allele.%s.%s.fasta'%(self.outdir,i+1,self.gene), 'w')
            print (hap_seq, file = out)
            out.close()
            #hap_seq = '>%s_%s\n'%(self.gene, i)
            # print (new_seg_sequence)

def correlation(consensus_order):
    strainNum = len(consensus_order)
    raw_order = []
    for i in range(strainNum):
        raw_order.append(i)
    possible_cases = list(permutations(raw_order, strainNum))
    max_share = 0
    max_case = ''
    for case in possible_cases:
        total_share = 0
        for i in range(strainNum):
            total_share += consensus_order[raw_order[i]][case[i]]
        if total_share >= max_share:
            max_case=case
            max_share = total_share
    return max_case

def chrom_seq(file):
    f = open(file, 'r')
    seg_sequence = {}
    name = 'start'
    seq = ''
    for line in f:
        line = line.strip()
        if line[0] == '>':
            seg_sequence[name] = seq
            name = line[1:]
            seq = ''
        else:
            seq += line
    seg_sequence[name] = seq
    del seg_sequence['start']
    return seg_sequence

def read_fasta(file):
    seq=''
    for line in open(file, 'r'):
        if line[0] == '>':
            continue
        line=line.strip()
        seq+=line
    return seq
        
def segment_mapping(fq1, fq2, ins_seq, outdir, gene, gene_ref):
    newref=outdir+'/newref_insertion.fa'
    os.system('cp %s %s'%(gene_ref, newref))
    for seg in ins_seq.keys():
        #build ref with each insertion seg
        # ref = outdir+'/ref_%s_%s.fa'%(gene, seg)
        
        # print ('cp %s %s'%(gene_ref, ref))
        f = open(newref, 'a')
        print ('>%s_%s\n%s'%(gene, int(seg), ins_seq[seg]), file = f)
        f.close()
        # index the ref
    # os.system('samtools faidx %s'%(newref))
    print ('new mapping start')
    map_call = """\
        outdir=%s/ 
        sample='newref_insertion' 
        samtools faidx %s 
        bwa index %s
        group='@RG\\tID:sample\\tSM:sample' 
        bwa mem -R $group -Y %s %s %s | samtools view -q 1 -F 4 -Sb | samtools sort > $outdir/$sample.sort.bam
        java -jar  /home/BIOINFO_TOOLS/alignment_tools/Picard/picard-tools-2.1.0/picard.jar MarkDuplicates INPUT=$outdir/$sample.sort.bam OUTPUT=$outdir/$sample.bam METRICS_FILE=$outdir/metrics.txt
        rm -rf $outdir/$sample.sort.bam 
        samtools index $outdir/$sample.bam 
        /home/wangmengyao/packages/freebayes/bin/freebayes -f %s -p 3 $outdir/$sample.bam > $outdir/$sample.freebayes.vcf 
        """%(outdir, newref, newref, newref, fq1, fq2, newref)
    # print (map_call)
    os.system(map_call)

def sv_copy_number(outdir, deletion_region, gene, ins_seq):
    os.system('samtools depth -a %s/newref_insertion.bam >%s/newref_insertion.depth'%(outdir, outdir))
    normal_depth = []
    deletions_depth_list = []
    for i in range(len(deletion_region)):
        deletions_depth_list.append([])

    for line in open('%s/newref_insertion.depth'%(outdir)):
        array = line.strip().split()
        if array[0] != gene:
            continue
        deletion_flag = False
        for i in range(len(deletion_region)):
            if float(array[1]) >= float(deletion_region[i][0]) and float(array[1]) < float(deletion_region[i][1]):
                deletions_depth_list[i].append(float(array[2]))
                deletion_flag = True
        if deletion_flag == False:
            normal_depth.append(float(array[2]))

    insertions_depth_list = {}
    for seg in ins_seq.keys():
        chrom_name = '%s_%s'%(gene, int(seg))
        insertions_depth_list[chrom_name] = []
    # print (insertions_depth_list.keys())
    
    for line in open('%s/newref_insertion.depth'%(outdir)):
        array = line.strip().split()
        if array[0] in insertions_depth_list.keys():
            insertions_depth_list[array[0]].append(float(array[2]))


    for i in range(len(deletion_region)):
        print (deletion_region[i])
        # deletion_region[i] += [np.mean(deletions_depth_list[i])/np.mean(normal_depth), zero_per(deletions_depth_list[i])]
        if float(deletion_region[i][0]) == float(deletion_region[i][1]):
            chrom_name = '%s_%s'%(gene, int(deletion_region[i][0]))
            print ('compute assign index', np.mean(insertions_depth_list[chrom_name]),np.mean(normal_depth))
            if np.mean(insertions_depth_list[chrom_name])/np.mean(normal_depth) > 0.5:
                deletion_region[i] += [2]
            else:
                deletion_region[i] += [1]
        else:
            if zero_per(deletions_depth_list[i]) > 0.5:
                deletion_region[i] += [0]
            else:
                deletion_region[i] += [1]
    
    #     print (deletion_region[i], np.mean(deletions_depth_list[i]), zero_per(deletions_depth_list[i]))
    # print ('normal', np.mean(normal_depth), zero_per(normal_depth), np.std(normal_depth))
    print ('####################################copy number', deletion_region)
    return deletion_region

def zero_per(list):
    zero_num = 0
    for li in list:
        if li < 3:
            zero_num  += 1
    return zero_num/len(list)

def uniq_reads(raw_reads_set):
    dup_type_Num = len(raw_reads_set)
    uniq_reads_set = []
    for i in range(dup_type_Num):
        uniq_set = []
        for ele in raw_reads_set[i]:
            uniq_flag = True
            for j in range(dup_type_Num):
                if i == j:
                    continue
                if ele in raw_reads_set[j]:
                    uniq_flag = False
            if uniq_flag == True:
                uniq_set.append(ele)
        uniq_reads_set.append(uniq_set)
    # print ('uniq reads', len(uniq_reads_set[0]), len(uniq_reads_set[1]))
    return uniq_reads_set

def DRB1_complex_region(drb1_complex_file):
    drb1_complex_seq = []
    drb1_complex_reads = []
    for line in open(drb1_complex_file, 'r'):
        line = line.strip()
        array = line.split()
        drb1_complex_seq.append(array[5])
        reads_list = []
        raw_reads_list = array[6].split(';')
        for reads in raw_reads_list[:-1]:
            reads_list.append(reads[:-3])
        drb1_complex_reads.append(reads_list)
    print ('There are %s dup types.'%(len(drb1_complex_seq)))
    uniq_drb1_complex_reads = uniq_reads(drb1_complex_reads)
    return drb1_complex_seq, uniq_drb1_complex_reads

def sv2fasta(ref, seg_order, index_locus, ins_seq, outdir):
    ref_seq = read_fasta(ref)
    normal_dict = {}
    for seg in index_locus.keys():
        if seg in ins_seq.keys():
            continue
        normal_dict[seg] = ref_seq[index_locus[seg][0]-1 : index_locus[seg][1]-1]
    for i in range(len(seg_order)):
        seg_array = seg_order[i].strip().split()
        hap_seq = '>sv_hap_%s\n'%(i)
        out = open(outdir + '/sv_hap_%s.fa'%(i), 'w')
        for arr in seg_array:
            seg_name = arr[:-1]
            if seg_name in ins_seq.keys():
                my_seq = ins_seq[seg_name]
            else:
                my_seq = normal_dict[seg_name]
            hap_seq += my_seq
        print (hap_seq, file = out)
        out.close()

def generate_break_points(outdir,snp_list,delta_set,strainsNum):    
    dup_dict=read_dup() #for dup region 
    #break points
    bp=open(outdir+'/%s_break_points.txt'%(gene),'w')
    print ('#gene\tlocus\t00\t01\t10\t11\tpoints_num\tnext_locus',file=bp) 
    points_number=0
    for i in range(len(snp_list)-1):
        dqa1_break=False

        points_number += 1
        sort_index=np.argsort(np.array(delta_set[i]))
        max_index=sort_index[-1]
        sec_index=sort_index[-2]
        # if delta_set[i][sort_index[-3]] > 5 or abs(max_index - sec_index) == 2:
        #     print (snp_list[i],delta_set[i])
        if points_number < args.points_num:
            continue
        if snp_list[i][0] in dup_dict.keys():
            for du in range(len(dup_dict[snp_list[i][0]])):
                dp_locus=int(dup_dict[snp_list[i][0]][du])
                if dp_locus != 0 and int(snp_list[i][1]) > dp_locus - 30 and int(snp_list[i][1]) < dp_locus + 30:
                    dqa1_break=True
                    dup_dict[snp_list[i][0]][du]=0
        # if sum(delta_set[i]) <args.reads_num or (strainsNum == 2 and abs(max_index - sec_index) == 2) or (strainsNum == 2 and delta_set[i][sort_index[-3]] > args.noise_num) or dqa1_break == True:
        if sum(delta_set[i]) <args.reads_num or (strainsNum == 2 and max_index + sec_index != 3) or (strainsNum == 2 and delta_set[i][sort_index[-3]] > args.noise_num) or dqa1_break == True:
            print (snp_list[i][0],snp_list[i][1], delta_set[i][0], delta_set[i][1],delta_set[i][2],delta_set[i][3],points_number,snp_list[i+1][1],file=bp)
            print ('############break points:', snp_list[i][0],snp_list[i][1], delta_set[i][0], delta_set[i][1],delta_set[i][2],delta_set[i][3],points_number,snp_list[i+1][1], dqa1_break, max_index + sec_index, delta_set[i][sort_index[-3]])
            points_number=0
            # print (snp_list[i],delta_set[i], dqa1_break)
    bp.close()

def dup_region_type(outdir, strainsNum, bamfile):
    order = r"""
        bam=%s
        outdir=%s
        k=%s
        pos=HLA_DRB1:3898-4400        
        ref=/mnt/disk2_workspace/wangmengyao/NeedleHLA/GA_rich/DRB1/bwa/DRB1_dup_extract_ref.fasta
        samtools view -f 64 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##1 "$2,$3}' > $outdir/extract.fa
        samtools view -f 128 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##2 "$2,$3}' >> $outdir/extract.fa
        #/home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/makeblastdb -in DRB1_dup_extract_ref.fasta -dbtype nucl -parse_seqids -out DRB1_dup_extract_ref.fasta
        /home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/blastn -query $outdir/extract.fa -out $outdir/extract.read.blast -db $ref -outfmt 6 -strand plus  -penalty -1 -reward 1 -gapopen 4 -gapextend 1
        perl /mnt/disk2_workspace/wangmengyao/NeedleHLA/GA_rich/DRB1/bwa/count.read.pl $outdir
        less $outdir/DRB1.hla.count| sort -k3,3nr -k4,4nr | head -n $k |awk '$3>0.8'|awk '$4>5' >$outdir/select.DRB1.seq.txt
        """%(bamfile, outdir, strainsNum)
    os.system(order)

def coverage_filter(bamfile, outdir):
    f = open(outdir + '/unreliable.genes.txt', 'w')
    d = open(outdir + '/mean.depth.genes.txt', 'w')
    print ('#gene mean median std', file = d)
    depth_file = '%s/bam.depth'%(outdir)
    depth_order = 'samtools depth %s > %s'%(bamfile, depth_file)
    os.system(depth_order)
    for gene in ['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']:
        locus_list = []
        depth_list = []
        low_depth_window = []
        for line in open(depth_file, 'r'):
            line = line.strip()
            array = line.split()
            array[1] = float(array[1])
            array[2] = float(array[2])
            if array[0] != gene:
                continue
            depth_list.append(array[2])
            if float(array[2]) < 5:
                if len(low_depth_window) == 0:
                    low_depth_window.append(array[1])
                elif float(array[1]) - low_depth_window[-1] == 1:
                    low_depth_window.append(array[1])
                else:
                    low_depth_window = [array[1]]
            if len(low_depth_window) > 10:
                print (gene, file = f)
                # print ('%s coverage dirty!'%(gene), low_depth_window)
                low_depth_window = []
                # break
        print (gene, round(np.mean(depth_list),1), np.median(depth_list), np.std(depth_list), file = d)
    f.close()
    d.close()
            
def edge_score(delta_set):
    my_delta = delta_set[:]
    score = 0  #only for two haplotypes
    edge_num = len(my_delta)
    for i in range(edge_num):
        delta = my_delta[i]
        if sum(delta) == 0:
            score += 1 
        else:
            delta = np.array(delta)
            delta = delta/sum(delta)
            index_sort=np.argsort(delta)
            noise = abs(delta[index_sort[3]] - 0.5) + abs(delta[3 - index_sort[3]] - 0.5)
            # print (delta, index_sort, noise)
            score += noise
    # print (score)
    if edge_num == 0:
        r_score = 0
    else:
        r_score = float(score)/edge_num
    return r_score

def breakpoints(bfile):
    sv_dict = {}
    print (bfile)
    f = open(bfile, 'r')
    for line in f:
        line = line.strip()
        array = line.split()
        if array[0] != array[3]:
            continue
        sv = [array[1], array[4], array[6]]
        if array[0] not in sv_dict:
            sv_dict[array[0]] = [sv]
        else:
            sv_dict[array[0]].append(sv)
    return sv_dict

def backup_find_deletion_region(sv_list):
    deletion_region = []
    ins_seq = {}
    for sv in sv_list:
        if sv[0] != sv[1]:
            deletion_region.append([int(float(sv[0])), int(float(sv[1]))])
        else:
            deletion_region.append([int(float(sv[0])), int(float(sv[1]))])
            seg = float(sv[0])
            ins_seq[seg] = sv[2]
    while True:
        flag = True
        for i in range(len(deletion_region) - 1):
            if deletion_region[i+1][0] < deletion_region[i][0]:
                flag = False
                a = deletion_region[i]
                deletion_region[i] = deletion_region[i+1]
                deletion_region[i+1] = a
        if flag:
            break
    # deletion_file = open('%s/deletion_region.txt'%(outdir), 'w')
    # for region in deletion_region:
    #     print (gene, region[0], region[1], 'deletion', file = deletion_file)
    # # print ('deletion_region', deletion_region)
    # deletion_file.close()
    return deletion_region, ins_seq

def find_deletion_region(sv_list):
    print (sv_list)
    deletion_region = []
    ins_seq = {}
    insertion = []
    points = []
    new_deletion_region = []
    for sv in sv_list:
        if sv[0] != sv[1]:
            points.append(int(float(sv[0])))
            points.append(int(float(sv[1])))
            deletion_region.append([int(float(sv[0])), int(float(sv[1]))])
        else:
            # points.append(sv[0])
            new_deletion_region.append([int(float(sv[0])), int(float(sv[1]))])
            seg = float(sv[0])
            ins_seq[seg] = sv[2]
    start = 1
    split_segs = []
    for p in sorted(points):
        if start == p:
            continue
        split_segs.append([start, p])
        start = p
    print (split_segs)
    
    for segs in split_segs:
        delete_flag = False
        for re in deletion_region:
            if segs[0] >= re[0] and segs[1] <= re[1]:
                delete_flag = True
        if delete_flag and segs[1] - segs[0] > 4:
            new_deletion_region.append(segs)
    deletion_region = new_deletion_region
    while True:
        flag = True
        new_deletion_region = deletion_region[:]
        for i in range(len(deletion_region) - 1):
            if deletion_region[i+1][0] < deletion_region[i][0]:
                flag = False
                new_deletion_region[i] = deletion_region[i+1]
                new_deletion_region[i+1] = deletion_region[i]
                # a = deletion_region[i]
                # deletion_region[i] = deletion_region[i+1]
                # deletion_region[i+1] = a
            elif deletion_region[i][1] > deletion_region[i+1][1]:
                new_deletion_region[i] = [deletion_region[i][0], deletion_region[i+1][0]]
                new_deletion_region[i+1] = [deletion_region[i+1][0], deletion_region[i+1][1]]
                new_deletion_region.append([deletion_region[i+1][1], deletion_region[i][1]])
        deletion_region = new_deletion_region[:]
        print (deletion_region)
        if flag:
            break
    print ('ordered deletion region:', deletion_region)
    return deletion_region, ins_seq

def split_vcf(gene, outdir, deletion_region):
    vcf = '%s/%s.vcf.gz'%(outdir,gene)
    os.system('tabix %s'%(vcf))
    # if len(deletion_region) == 0:
    #     print ('no sv!')
    #     return 0
    vcf_gap = []
    start = 1001
    break_points_list = [3950]
    for dele in deletion_region:
        if abs(start - dele[0]) < 1:
            continue
        break_points_list.append(dele[0])
    break_points_list = sorted(break_points_list)
    start = 1001
    for b_point in break_points_list: 
        if b_point - start < 500:
            continue
        if b_point >6000 and b_point < 7000:
            continue
        vcf_gap.append([start, b_point])
        start = b_point
    if focus_region()[gene][1] - start >= 500:
        vcf_gap.append([start, focus_region()[gene][1]])
    else:
        vcf_gap[-1][1] = focus_region()[gene][1] 
    os.system('rm %s/%s_part_*_*_*.vcf'%(outdir, gene))
    i = 0
    for gap in vcf_gap:
        order = "/home/wangmengyao/miniconda2/bin/bcftools filter -t %s:%s-%s %s -o %s/%s_part_%s_%s_%s.vcf"%(gene, gap[0], gap[1], vcf, outdir, gene, i, gap[0], gap[1])
        os.system(order)
        i+=1
    print (vcf_gap)
    return break_points_list


if __name__ == "__main__":   
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})
    else:       
        raw_bamfile,outdir,strainsNum,snp_dp,germline_flag,database_flag,indel_len,freq_bias=\
            args.bamfile,args.outdir,args.k,args.snp_dp,args.germline,\
            args.database_flag,args.indel_len,args.freq_bias           
        snp_qual = args.snp_qual
        print (outdir)
        if not os.path.exists(outdir):
            os.system('mkdir '+ outdir) 
        #######realign and find break points for SV
        bamfile = raw_bamfile
        vcffile = args.vcf
        # realign_and_sv_break.main_no_realign(raw_bamfile, outdir, args.rlen, 'WGS')
        dup_region_type(outdir, strainsNum, bamfile)
        # sv_break_points = outdir +'/sample.breakpoint.txt'
        fq1 = args.fq1
        fq2 = args.fq2
        dup_file = outdir +'/select.DRB1.seq.txt'
        # to calculate score that determine whether the result is reliable
        sco = open(outdir + '/edge.score.txt', 'w')
        coverage_filter(bamfile, outdir)
        gene = args.gene

        # gene = 'HLA_A'
        # bfile = '/mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20200901/HLA_11_T_50-50/HLA_11_T_50-50.breakpoint.txt'
        sv_dict = breakpoints(args.sv)
        if gene in sv_dict.keys():
            sv_result = sv_dict[gene]
        else:
            sv_result = []
        deletion_region, ins_seq = find_deletion_region(sv_result)
        # print (deletion_region, ins_seq)

        ######for SNV
        snp_list,beta_set,allele_set = read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,gene,\
            freq_bias,strainsNum,deletion_region, snp_qual)      
        if germline_flag == True: 
            #For germline, the first-order genotype frequency is 0.5 for all locus, 
            #thus will not provide any info.
            args.weight = 0.0
        if len(snp_list)==0:
            print ('No hete locus')
            gene_profile = no_snv_gene_phased(vcffile, outdir, gene, strainsNum)
            print (gene, 0, 0, 0, 0, file = sco)
        else:  
            delta_set=second_beta(bamfile,snp_list)   
            r_score = edge_score(delta_set)  #just for double haps
            print (gene, r_score, len(delta_set), end= '\t', file = sco)

            for i in range(len(delta_set)):
                print (snp_list[i], delta_set[i])

            fir_beta,sec_beta=rectify(snp_list,beta_set,delta_set,args.lambda1,args.lambda2,\
                germline_flag,database_flag)

            print ("Rectifying is over! Next start phasing the main sequence.")

            if strainsNum==0:
                wo=Workflow(fir_beta,sec_beta,delta_set,args.weight,args.elbow,allele_set)
                final_alpha,seq_list,loss=wo.choose_k()
                strainsNum = len(final_alpha)
            else:
                # print (fir_beta, sec_beta)
                wo=Workflow(fir_beta,sec_beta,delta_set,args.weight,args.elbow,allele_set)
                final_alpha,seq_list,loss=wo.given_k(strainsNum)  
            print ('The main phasing is over!')

            print (loss, end= '\t', file = sco)
            if len(snp_list)==0:
                print (0, file = sco)
            else:
                print (float(loss)/len(snp_list), file = sco)

            output(outdir,final_alpha,seq_list,snp_list,gene,freq_bias)
            generate_break_points(outdir,snp_list,delta_set,strainsNum)
            if gene == 'HLA_DRB1':
                split_vcf(gene, outdir, deletion_region)
# '''
            print ('Start rephase!')
            if gene == 'HLA_DRB1':
                reph='perl /home/wangmengyao/scripts/PHLAT/script/rephase.DRB1.pl %s/%s_break_points.txt\
                    %s %s %s/%s_break_points_phased.txt %s %s'%(outdir,gene,outdir,strainsNum,outdir,\
                    gene,args.block_len,args.points_num)
            else:
                reph='perl /home/wangmengyao/scripts/NeedleHLA/script/rephaseV1.pl %s/%s_break_points.txt\
                    %s %s %s/%s_break_points_phased.txt %s %s'%(outdir,gene,outdir,strainsNum,outdir,\
                    gene,args.block_len,args.points_num)
            print ('before rephase order: %s'%(reph))
            os.system(str(reph))
            print ('after rephase order: %s'%(reph))
            print ('Rephase is over!')
            update_seqlist=newphase(outdir,final_alpha,seq_list,snp_list,vcffile,gene)   #need to refresh alpha with new result.
            fresh_alpha = newalpha(update_seqlist, sec_beta, strainsNum, allele_set, args.weight, fir_beta)
            print ('fresh_alpha', fresh_alpha)
            freq_output(outdir, gene, final_alpha, germline_flag)
            gene_profile=gene_phased(update_seqlist,snp_list)
            print ('snv of %s is done!'%(gene))
        print ('segment mapping', ins_seq)

        if len(ins_seq) > 0:
            segment_mapping(fq1, fq2, ins_seq, outdir, gene, args.ref)
        else:
            os.system('cp %s/%s.bam %s/newref_insertion.bam'%(outdir, gene.split('_')[-1], outdir))
            os.system('samtools index %s/newref_insertion.bam'%(outdir))
            os.system('cp %s %s/newref_insertion.freebayes.vcf'%(vcffile, outdir))
        deletion_region = sv_copy_number(outdir, deletion_region, gene, ins_seq)
        sh = Share_reads(deletion_region, outdir, strainsNum, gene, gene_profile, ins_seq, dup_file)
        sh.split_seg()
        if gene == 'HLA_DRB1':
            for t in range(2):
                drb1_exon = 'samtools faidx /home/wangmengyao/scripts/PHLAT/database/ref/extend/hla.ref.extend.fa\
                    HLA_DRB1:6972-7241 |/home/wangmengyao/miniconda2/bin/bcftools\
                    consensus -H %s %s/%s.rephase.vcf.gz  >%s/%s_%s_exon.fa'%(t+1, outdir, gene, outdir, gene, t+1)
                print (drb1_exon)
                os.system(drb1_exon)
# '''

