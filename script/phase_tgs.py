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
python3 phase_tgs.py [options] 

Help information can be found by python3 phase_tgs.py -h/--help, additional information can be found in \
README.MD or https://github.com/deepomicslab/HLAPro.
"""
scripts_dir=sys.path[0]+'/'
parser = ArgumentParser(description="SpecHLA.",prog='python3 phase_tgs.py',usage=Usage)
optional=parser._action_groups.pop()
required=parser.add_argument_group('required arguments')
flag_parser = parser.add_mutually_exclusive_group(required=False)
flag_data = parser.add_mutually_exclusive_group(required=False)
#necessary parameter
required.add_argument("--ref",help="The hla reference file used in alignment",dest='ref',metavar='', type=str)
required.add_argument("-b", "--bam",help="The bam file of the input samples.",dest='bamfile',metavar='')
required.add_argument("-v", "--vcf",help="The vcf file of the input samples.",dest='vcf',metavar='')
required.add_argument("--sa",help="Sample ID",dest='sample_id',metavar='', type=str)
required.add_argument("-s", "--sv",help="Long Indel file after scanindel, we will not consider long InDel \
    if not afforded.",dest='sv',metavar='')
required.add_argument("-a", "--phase",help="Choose phasing method, SpecHap if True, otherwise use PStrain.\
     Default is SpecHap.",dest='phase_flag',metavar='',default=True,type=str2bool)
required.add_argument("--gene",help="gene",dest='gene',metavar='', type=str)
required.add_argument("--fq1",help="fq1",dest='fq1',metavar='', type=str)
required.add_argument("--fq2",help="fq2",dest='fq2',metavar='', type=str)
required.add_argument("--tgs",help="PACBIO TGS fastq",dest='tgs',metavar='', type=str)
required.add_argument("--nanopore",help="NANOPORE TGS fastq",dest='nanopore',metavar='', type=str)
required.add_argument("--hic_fwd",help="fwd_hic.fastq",dest='hic_fwd',metavar='', type=str)
required.add_argument("--hic_rev",help="rev_hic.fastq",dest='hic_rev',metavar='', type=str)
required.add_argument("--tenx",help="10X data",dest='tenx',metavar='', type=str)
required.add_argument("-o", "--outdir",help="The output directory.",dest='outdir',metavar='')
#alternative parameter
optional.add_argument("--freq_bias",help="freq_bias (default is 0.05)",dest='freq_bias',\
    metavar='',default=0.05, type=float)
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

def if_in_deletion(locus, deletion_region):
    dele_flag = False
    for deletion in deletion_region:
        if locus >= deletion[0] and locus < deletion[1]:
            dele_flag = True
    return dele_flag

def read_spechap_seq(vcf, snp_list):
    snp_locus_list = []
    for snp in snp_list:
        snp_locus_list.append(int(snp[1]))
    
    seq_list = [[],[]]
    in_vcf = VariantFile(vcf)
    sample = list(in_vcf.header.samples)[0]
    for record in in_vcf.fetch():
        geno = record.samples[sample]['GT']
        #print (record)
        # if geno == (1,1):
        #     continue
        if record.pos not in snp_locus_list:
            continue
        if sum(geno) == 1 or sum(geno) == 0:
            for i in range(2):
                seq_list[i].append(geno[i])
        else: 
            # print (record, geno)
            for i in range(2):
                seq_list[i].append(geno[i] - 1)

    return seq_list 

def read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,chrom_name,freq_bias,strainsNum,deletion_region,snp_qual):
    snp_index = 1
    snp_index_dict = {}
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
        ##########after filter snp################
        snp_index_dict[record.pos] = snp_index
        snp_index += 1
        if if_in_deletion(record.pos, deletion_region) and geno != (1,1,1):
            #print ('SNV in deletion region, need edit', record)
            if geno == (1,1,2) or geno == (1,2,2):                
                snp = [record.chrom,record.pos,record.alts[0],record.alts[1],record.ref]
            else:
                snp=[record.chrom,record.pos,record.ref,record.alts[0],record.ref]

            reads_list = reads_support(samfile, snp)
            allele_dp = [len(reads_list[0]), len(reads_list[1])]
            new_dp=sum(allele_dp)
            if new_dp == 0:
                print ('WARNING: the depth of the locus obtained by pysam is zero!',snp)
                continue
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
            if new_dp == 0:
                print ('WARNING: the depth of the locus obtained by pysam is zero!',snp)
                continue
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
    print ("The number of short hete loci is %s."%(len(snp_list)))
    return snp_list, beta_set, allele_set, snp_index_dict

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
    # if germline_flag:
    #     final_alpha = [0.5, 0.5]
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

def rectify(snp_list,beta_set,delta_set,lambda1,lambda2,germline_flag):
    nucleotide={'A':0,'T':1,'C':2,'G':3}
    # rectify beta
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
        
        if int(first[1])-1 in read.get_reference_positions(full_length=True) and read.mapping_quality >1:   
            
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
    return reads_list

def link_reads(samfile,left,right,new_left,snp_index_dict,f):
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
            reads_name = set(left_set).intersection(set(right_set))
            for name in reads_name:
                if len(left[2]) > 1 or len(left[3]) > 1 or len(right[2]) > 1 or len(right[3]) > 1:
                    left_index = snp_index_dict[int(left[1])]
                    right_index = snp_index_dict[int(right[1])]
                    left_geno = i
                    right_geno = j
                    # if left[4] != left[2]:
                    #     left_geno += 1
                    # if right[4] != right[2]:
                    #     right_geno += 1
                    # print (left,right,'2 %s %s %s %s %s II 60'%(name, left_index, left_geno, right_index, right_geno))
                    if new_formate:
                        print('2 %s 1 -1 -1 %s %s %s %s II 60'%(name, left_index, left_geno, right_index, right_geno), file=f)
                    else:
                        print('2 %s %s %s %s %s II 60'%(name, left_index, left_geno, right_index, right_geno), file=f)
            same_num=len(reads_name)
            delta_count.append(same_num)
    return delta_count,right_reads

def second_beta(bamfile,snp_list,snp_index_dict,outdir):
    f = open(outdir + '/fragment.add.file', 'w')
    delta_set=[]
    samfile = pysam.AlignmentFile(bamfile, "rb")
    new_left=''
    for i in range(len(snp_list)-1):  
        left=snp_list[i]
        right=snp_list[i+1]  
        if new_left=='':   
            new_left=reads_support(samfile,left)
        delta_count,right_reads=link_reads(samfile,left,right,new_left,snp_index_dict, f)
        delta_set.append(delta_count)
        new_left=right_reads
    f.close()
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
            #change name when the break point locus is not same as snp locus.
            if replace_dict[key] != key: 
                dict[replace_dict[key]] = dict[key]
                del dict[key]

    block_dict[previous_gene] = dict
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
        # if gene == 'HLA_DRB1':
        #     block_dict=relate_order(file, snp_list)
        # else:
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
    #m = VariantFile('%s/%s.vcf.gz'%(outdir, gene))
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
            # print (phased_locus)
            update_phased_locus=[]
            for pp in range(k):
                update_phased_locus.append(phased_locus[int(ref_order[pp])])
            phased_locus=update_phased_locus[:]
            # print (record, phased_locus)
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

    m.close()
    out.close()
    os.system('%s/../bin/tabix -f %s/%s.rephase.vcf.gz'%(sys.path[0],outdir,gene))
    return update_seqlist

def newalpha(update_seqlist, sec_beta, strainsNum, allele_set, weight, fir_beta):
    return alpha_step(sec_beta, update_seqlist, strainsNum, allele_set, weight, fir_beta)

def gene_phased(update_seqlist,snp_list, gene):
    gene_profile={}
    # gene_name=['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']
    # for gene in gene_name:
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
    os.system('%s/../bin/tabix -f %s/%s.rephase.vcf.gz'%(sys.path[0],outdir,gene))
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

def focus_region():
    return {'HLA_A':[1000,4503],'HLA_B':[1000,5081],'HLA_C':[1000,5304],'HLA_DPA1':[1000,10775],\
        'HLA_DPB1':[1000,12468],'HLA_DQA1':[1000,7492],'HLA_DQB1':[1000,8480],'HLA_DRB1':[1000,12229]}

class Share_reads():

    def __init__(self, deletion_region, outdir, strainsNum, gene, gene_profile, ins_seq):
        self.deletion_region = deletion_region
        print ('initial', self.deletion_region)
        self.bamfile = outdir + '/newref_insertion.bam'
        self.vcf = outdir + '/%s.insertion.phased.vcf.gz'%(gene)
        self.strainsNum = strainsNum
        self.gene = gene
        self.normal_sequence=gene_profile[self.gene]
        self.reads_support = self.normal_reads()
        self.outdir = outdir
        self.ins_seq = ins_seq
        self.dup_file = outdir +'/select.DRB1.seq.txt'
        
    def generate_normal_region(self):
        gene_area = focus_region()[self.gene]
        normal_region = []
        segs = []
        gene_area[0] = int(float(gene_area[0])) + 1
        gene_area[1] = int(float(gene_area[1]))
        start = gene_area[0]
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
        return normal_region, segs

    def normal_reads(self):
        normal_region, segs = self.generate_normal_region()
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
        for read in samfile.fetch(self.gene,float(self.deletion_region[deletion_index][0])-1,float(self.deletion_region[deletion_index][1])):
            #print (deletion_index, read.query_name, float(self.deletion_region[deletion_index][0])-1, float(self.deletion_region[deletion_index][1]))
            for i in range(self.strainsNum):
                if read.query_name in self.reads_support[i]:
                    link_reads[i] += 1
        print (deletion_index, link_reads, self.most_support(link_reads))
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
        return self.most_support(link_reads)

    def deletion_phase(self):
        for deletion_index in range(len(self.deletion_region)):
            # print (self.deletion_region[deletion_index])
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
                    print ('DRB1 assignment', i, j, max_num, len(drb1_complex_seq))
            new_drb1_complex_seq.append(max_seq)
        return new_drb1_complex_seq

    def split_seg(self):
        # print ('start split seg.')
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
            order='%s/../bin/samtools faidx %s/../db/ref/hla.ref.extend.fa\
                %s |%s/../bin/bcftools\
                consensus -H %s %s/%s.rephase.vcf.gz  >%s/%s_%s_seg.fa'%(sys.path[0],sys.path[0],gap,sys.path[0],i+1,self.outdir,\
                self.gene,self.outdir,self.gene,i)
            os.system(order)
            fa_file = '%s/%s_%s_seg.fa'%(self.outdir,self.gene,i)
            seg_sequence = chrom_seq(fa_file)
            new_seg_sequence = {}
            for segseq in seg_sequence.keys():
                new_seg_sequence[id_name[segseq]] = seg_sequence[segseq]
            hap_seq = ''
            for seg in segs:
                # print (seg, str(seg[0]) + '_' + str(seg[1]))
               # print (str(seg[0]) + '_' + str(seg[1]) )
                if seg[2] == 'normal':
                    hap_seq += new_seg_sequence[str(seg[0]) + '_' + str(seg[1]) ]
                elif seg[2] == 'deletion':
                    deletion_index = seg[3]
                    assign_index = self.deletion_reads(deletion_index)[0]
                    print ('deletion assign_index', assign_index, self.deletion_region[deletion_index])
                    if  i == assign_index and self.deletion_region[deletion_index][2] == 1:
                        hap_seq += new_seg_sequence[str(seg[0]) + '_' + str(seg[1]) ]
                elif seg[2] == 'insertion':
                    # continue
                    deletion_index = seg[3]
                    assign_index = self.insertion_reads(deletion_index)[0]
                    print ('###########insertion', seg, assign_index, self.deletion_region[deletion_index][2])
                    insertion_seg = '%s_%s'%(self.gene, seg[0])
                    if self.deletion_region[deletion_index][2] == 2:                        
                        insert_seq = self.link_diploid_insertion(insertion_seg)
                        hap_seq += insert_seq[i]
                        # print ('#2 copy insertion', seg, assign_index,self.deletion_region[deletion_index][2])
                        # hap_seq += self.ins_seq[seg[0]]
                        
                    elif  i == assign_index and self.deletion_region[deletion_index][2] == 1:
                        # hap_seq += self.ins_seq[seg[0]]
                        # print ('#1 copy insertion', seg, assign_index,self.deletion_region[deletion_index][2])
                        hap_seq += self.consensus_insertion(insertion_seg)

                elif seg[2] == 'dup':
                    # print ('the seg contain dup region')
                    # the seg that contain the dup region
                    hap_seq += new_seg_sequence[str(seg[0])+'>']
                    #add the seq of chosen dup type
                    hap_seq += my_drb1_complex_seq[i]
                    # the seg that contain the dup region
                    hap_seq += new_seg_sequence[str(seg[1])+'<']                    
            hap_seq =  '>%s_%s\n'%(self.gene, i) + hap_seq[:]
            print ('________________', len(hap_seq))
            out = open('%s/hla.allele.%s.%s.fasta'%(self.outdir,i+1,self.gene), 'w')
            print (hap_seq, file = out)
            out.close()

    def link_diploid_insertion(self, insertion_seg):
        # insertion_seg = 'HLA_DRB1_6355'
        insert_reads_support = []
        for i in range(self.strainsNum):
            insert_reads_support.append([])

        in_vcf = VariantFile(self.vcf)
        sample = list(in_vcf.header.samples)[0]
        insert_phase_result = [[],[]]
        for record in in_vcf.fetch():
            geno = record.samples[sample]['GT']    
            depth = record.samples[sample]['AD']
            if record.chrom !=  insertion_seg:
                continue
            if geno == (1, 1):
                continue
            if geno == (0, 1):
                insert_phase_result[0].append([record.pos, record.ref])
                insert_phase_result[1].append([record.pos, record.alts[0]])
            elif geno == (1, 0):
                insert_phase_result[1].append([record.pos, record.ref])
                insert_phase_result[0].append([record.pos, record.alts[0]])  
            if len(insert_phase_result[0]) > 0:
                if len(insert_phase_result[1][-1][1]) != 1 or  len(insert_phase_result[0][-1][1]) != 1:
                    continue          

        samfile = pysam.AlignmentFile(self.bamfile, "rb")
        for read in samfile.fetch(insertion_seg):
            rivet_points=False
            support_alleles=[]
            support_loci=[]
            for i in range(len(insert_phase_result[0])):
                snv1 =insert_phase_result[0][i]
                snv2 =insert_phase_result[1][i]
                if int(snv1[0])-1 in read.get_reference_positions(full_length=True):
                    reads_index=read.get_reference_positions(full_length=True).index(int(snv1[0])-1)
                    if read.query_sequence[reads_index] == snv1[1]:
                        insert_reads_support[0].append(read.query_name)
                    elif read.query_sequence[reads_index] == snv2[1]:
                        insert_reads_support[1].append(read.query_name)
        r00 = 0
        for i in range(2):
            for j in range(len(insert_reads_support[i])):
                if insert_reads_support[i][j] in self.reads_support[i]:
                    r00 += 1
        r01 = 0
        for i in range(2):
            for j in range(len(insert_reads_support[i])):
                if insert_reads_support[i][j] in self.reads_support[1-i]:
                    r01 += 1
        fastq_seq = []
        for i in range(2):
            order = """
            %s/../bin/samtools faidx %s/newref_insertion.fa %s|%s/../bin/bcftools consensus -H %s %s  >%s/seq_%s_%s.fa
            """%(sys.path[0], self.outdir, insertion_seg, sys.path[0], i+1, self.vcf, self.outdir, i, insertion_seg)
            os.system(order)
            fastq_seq.append(read_fasta('%s/seq_%s_%s.fa'%(self.outdir, i, insertion_seg)))
        print ('link long indel supporting reads', r00, r01)
        if  r01 > r00:
            return [fastq_seq[1], fastq_seq[0]]
        else:
            return fastq_seq  

    def consensus_insertion(self, insertion_seg):
        order = """
        %s/../bin/samtools faidx %s/newref_insertion.fa %s|%s/../bin/bcftools consensus -H %s %s  >%s/seq
        """%(sys.path[0], self.outdir, insertion_seg, sys.path[0], 1, self.vcf, self.outdir)
        os.system(order)
        cons_seq = read_fasta('%s/seq'%(self.outdir))
        return cons_seq

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
        
def segment_mapping_pre(fq1, fq2, ins_seq, outdir, gene, gene_ref):
    newref=outdir+'/newref_insertion.fa'
    os.system('cp %s %s'%(gene_ref, newref))
    for seg in ins_seq.keys():

        f = open(newref, 'a')
        print ('>%s_%s\n%s'%(gene, int(seg), ins_seq[seg]), file = f)
        f.close()
        # index the ref
    print ('New mapping starts to link long InDels.')
    map_call = """\
        bindir=%s/../bin/
        outdir=%s/ 
        sample='newref_insertion' 
        $bindir/samtools faidx %s 
        $bindir/bwa index %s
        group='@RG\\tID:sample\\tSM:sample'  #only -B 1
        $bindir/bwa mem -B 1 -O 1,1 -L 1,1 -U 1 -R $group -Y %s %s %s | $bindir/samtools view -q 1 -F 4 -Sb | $bindir/samtools sort > $outdir/$sample.sort.bam
        java -jar  $bindir/picard.jar MarkDuplicates INPUT=$outdir/$sample.sort.bam OUTPUT=$outdir/$sample.bam METRICS_FILE=$outdir/metrics.txt
        rm -rf $outdir/$sample.sort.bam 
        $bindir/samtools index $outdir/$sample.bam 
        $bindir/freebayes -f %s -p 2 $outdir/$sample.bam > $outdir/$sample.freebayes.1.vcf 
        cat $outdir/$sample.freebayes.1.vcf| sed -e 's/\//\|/g'>$outdir/$sample.freebayes.vcf 
        bgzip -f $outdir/$sample.freebayes.vcf 
        tabix -f $outdir/$sample.freebayes.vcf.gz
        """%(sys.path[0], outdir, newref, newref, newref, fq1, fq2, newref)
    os.system(map_call)
    # print (ins_seq)
    for ins in ins_seq.keys():
        ins_call = """%s/../bin/samtools faidx %s/newref_insertion.fa %s_%s |%s/../bin/bcftools consensus -H 1 %s/newref_insertion.freebayes.vcf.gz  >%s/fresh_ins.fa
        """%(sys.path[0],outdir,gene,int(ins),sys.path[0],outdir,outdir)
        # print ('#####################', ins, ins_call)
        os.system(ins_call)
        ins_seq[ins] = read_fasta('%s/fresh_ins.fa'%(outdir))
    # print (ins_seq)
    return ins_seq

    # map_call = """\
    #     bindir=%s/../bin/
    #     outdir=%s/ 
    #     sample='newref_insertion' 
    #     $bindir/samtools faidx %s 
    #     $bindir/bwa index %s
    #     group='@RG\\tID:sample\\tSM:sample'  #only -B 1
    #     /home/wangmengyao/packages/Novoalign/novocraft/novoindex %s.ndx %s
    #     /home/wangmengyao/packages/Novoalign/novocraft/novoalign -g 10 -x 1 -F STDFQ -o SAM -o FullNW  -d %s.ndx -f %s %s| $bindir/samtools view -q 1 -F 4 -Sb | $bindir/samtools sort > $outdir/$sample.sort.bam
    #     java -jar  $bindir/picard.jar MarkDuplicates INPUT=$outdir/$sample.sort.bam OUTPUT=$outdir/$sample.bam METRICS_FILE=$outdir/metrics.txt
    #     rm -rf $outdir/$sample.sort.bam 
    #     $bindir/samtools index $outdir/$sample.bam 
    #     $bindir/freebayes -f %s -p 2 $outdir/$sample.bam > $outdir/$sample.freebayes.vcf 
    #     """%(sys.path[0], outdir, newref, newref, newref, newref, newref, fq1, fq2, newref)
        # -B 1 -O 1,1 -L 1,1 -U 1 
    # print (map_call)
    
def segment_mapping(fq1, fq2, ins_seq, outdir, gene, gene_ref):
    newref=outdir+'/newref_insertion.fa'
    os.system('cp %s %s'%(gene_ref, newref))
    for seg in ins_seq.keys():

        f = open(newref, 'a')
        print ('>%s_%s\n%s'%(gene, int(seg), ins_seq[seg]), file = f)
        f.close()
        # index the ref
    print ('New mapping starts to link long InDels.')
    map_call = """\
        bindir=%s/../bin/
        outdir=%s/ 
        sample='newref_insertion' 
        $bindir/samtools faidx %s 
        $bindir/bwa index %s
        group='@RG\\tID:sample\\tSM:sample'  #only -B 1
        $bindir/bwa mem -B 1 -O 1,1 -L 1,1 -U 1 -R $group -Y %s %s %s | $bindir/samtools view -q 1 -F 4 -Sb | $bindir/samtools sort > $outdir/$sample.sort.bam
        java -jar  $bindir/picard.jar MarkDuplicates INPUT=$outdir/$sample.sort.bam OUTPUT=$outdir/$sample.bam METRICS_FILE=$outdir/metrics.txt
        rm -rf $outdir/$sample.sort.bam 
        $bindir/samtools index $outdir/$sample.bam 
        $bindir/freebayes -f %s -p 2 $outdir/$sample.bam > $outdir/$sample.freebayes.vcf 
        """%(sys.path[0], outdir, newref, newref, newref, fq1, fq2, newref)
    os.system(map_call)

def sv_copy_number_old(outdir, deletion_region, gene, ins_seq):
    os.system('%s/../bin/samtools depth -a %s/newref_insertion.bam >%s/newref_insertion.depth'%(sys.path[0],outdir, outdir))
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
        # deletion_region[i] += [np.mean(deletions_depth_list[i])/np.mean(normal_depth), zero_per(deletions_depth_list[i])]
        if float(deletion_region[i][0]) == float(deletion_region[i][1]):
            chrom_name = '%s_%s'%(gene, int(deletion_region[i][0]))
            # print ('compute assign index', np.mean(insertions_depth_list[chrom_name]),np.mean(normal_depth))
            if np.mean(insertions_depth_list[chrom_name])/np.mean(normal_depth) > 0.5:
                deletion_region[i] += [2]
            else:
                deletion_region[i] += [1]
        else:
            print ('compute copy number for deletion', deletion_region[i], zero_per(deletions_depth_list[i]))
            if zero_per(deletions_depth_list[i]) > 0.2:
                deletion_region[i] += [0]
            else:
                deletion_region[i] += [1]
        print ('### copy number', deletion_region[i])
    
    return deletion_region

def sv_copy_number(deletion_region, sv_list):
    for i in range(len(deletion_region)):
        deletion_region[i] += [0]
    for i in range(len(deletion_region)):
        for sv in sv_list:
            if deletion_region[i][0] >= int(sv[0]) and deletion_region[i][1] <= int(sv[1]):
                deletion_region[i][2] += sv[3]
                # break
        if deletion_region[i][2] > 2:
            deletion_region[i][2] = 2
        print ('### copy number', deletion_region[i])
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
    # print ('There are %s dup types.'%(len(drb1_complex_seq)))
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
    bp=open(outdir+'/%s_break_points_pstrain.txt'%(gene),'w')
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
        # if points_number < args.points_num:
        #     continue
        if snp_list[i][0] in dup_dict.keys():
            for du in range(len(dup_dict[snp_list[i][0]])):
                dp_locus=int(dup_dict[snp_list[i][0]][du])
                if dp_locus != 0 and int(snp_list[i][1]) > dp_locus - 30 and int(snp_list[i][1]) < dp_locus + 30:
                    dqa1_break=True
                    dup_dict[snp_list[i][0]][du]=0
        # if sum(delta_set[i]) <args.reads_num or (strainsNum == 2 and abs(max_index - sec_index) == 2) or (strainsNum == 2 and delta_set[i][sort_index[-3]] > args.noise_num) or dqa1_break == True:
        if sum(delta_set[i]) <args.reads_num or (strainsNum == 2 and max_index + sec_index != 3) or (strainsNum == 2 and delta_set[i][sort_index[-3]] > args.noise_num) or dqa1_break == True:
            print (snp_list[i][0],snp_list[i][1], delta_set[i][0], delta_set[i][1],delta_set[i][2],delta_set[i][3],points_number,snp_list[i+1][1],file=bp)
            # print ('############break points:', snp_list[i][0],snp_list[i][1], delta_set[i][0], delta_set[i][1],delta_set[i][2],delta_set[i][3],points_number,snp_list[i+1][1], dqa1_break, max_index + sec_index, delta_set[i][sort_index[-3]])
            points_number=0
            # print (snp_list[i],delta_set[i], dqa1_break)
    bp.close()

def dup_region_type(outdir, strainsNum, bamfile):
    order = r"""
        bam=%s
        outdir=%s
        k=%s
        pos=HLA_DRB1:3898-4400        
        ref=%s/../db/ref/DRB1_dup_extract_ref.fasta
        %s/../bin/samtools view -f 64 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##1 "$2,$3}' > $outdir/extract.fa
        %s/../bin/samtools view -f 128 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##2 "$2,$3}' >> $outdir/extract.fa
        %s/../bin/blastn -query $outdir/extract.fa -out $outdir/extract.read.blast -db $ref -outfmt 6 -strand plus  -penalty -1 -reward 1 -gapopen 4 -gapextend 1
        perl %s/count.read.pl $outdir
        less $outdir/DRB1.hla.count| sort -k3,3nr -k4,4nr | head -n $k |awk '$3>0.7'|awk '$4>5' >$outdir/select.DRB1.seq.txt
        """%(bamfile, outdir, strainsNum, sys.path[0], sys.path[0], sys.path[0],sys.path[0], sys.path[0])
    os.system(order)

def long_InDel_breakpoints(bfile):
    sv_dict = {}
    if not os.path.isfile(bfile):
        return sv_dict
    f = open(bfile, 'r')
    for line in f:
        line = line.strip()
        array = line.split()
        if array[0] != array[3]:
            continue
        if array[0] == 'HLA_DRB1':
            if int(array[1]) > 3800 and int(array[1]) < 4500:
                continue
            if int(array[4]) > 3800 and int(array[4]) < 4500:
                continue
        sv = [array[1], array[4], array[6]]
        #sv = [array[1], array[4], array[6], int(array[7])]
        if array[0] not in sv_dict:
            sv_dict[array[0]] = [sv]
        else:
            sv_dict[array[0]].append(sv)
    return sv_dict

def find_deletion_region_BK(sv_list):
    # print (sv_list)
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
            # deletion_region.append([int(float(sv[0])), int(float(sv[1])), int(sv[3])])
        else:

            # new_deletion_region.append([int(float(sv[0])), int(float(sv[1])), int(sv[3])])
            # seg = float(sv[0])
            # ins_seq[seg] = sv[2]

            #remove redundant insertions
            uniq_ins = True
            for region in new_deletion_region:
                if region[0] != region[1]:
                    continue
                if abs(int(sv[0]) - region[0]) < 50:
                    uniq_ins = False
            if uniq_ins:
                # new_deletion_region.append([int(float(sv[0])), int(float(sv[1])), int(sv[3])])
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
    # print (split_segs)
    # print (new_deletion_region)
    
    
    for segs in split_segs:
        delete_flag = False
        for re in deletion_region:
            if segs[0] >= re[0] and segs[1] <= re[1]:
                delete_flag = True
        if delete_flag and segs[1] - segs[0] > 4:
            new_deletion_region.append(segs)
    deletion_region = new_deletion_region
    # print ('new_deletion_region',deletion_region)
    while True:
        flag = True
        new_deletion_region = deletion_region[:]
        for i in range(len(deletion_region) - 1):
            if deletion_region[i+1][0] < deletion_region[i][0]:
                flag = False
                new_deletion_region[i] = deletion_region[i+1]
                new_deletion_region[i+1] = deletion_region[i]
                # print ('ite', i, new_deletion_region)
                break
            elif deletion_region[i][1] > deletion_region[i+1][1]:
                new_deletion_region[i] = [deletion_region[i][0], deletion_region[i+1][0]]
                new_deletion_region[i+1] = [deletion_region[i+1][0], deletion_region[i+1][1]]
                new_deletion_region.append([deletion_region[i+1][1], deletion_region[i][1]])
                break
        deletion_region = new_deletion_region[:]
        # print (deletion_region)
        if flag:
            break
    print ('#ordered deletion region:', deletion_region)
    return deletion_region, ins_seq

def find_deletion_region(sv_list):
    # print (sv_list)
    deletion_region = []
    ins_seq = {}
    insertion = []
    points = []

    for sv in sv_list:
        if sv[0] != sv[1]:
            points.append(int(float(sv[0])))
            points.append(int(float(sv[1])))
            deletion_region.append([int(float(sv[0])), int(float(sv[1])), int(sv[3])])
        else:
            # points.append(sv[0])
            deletion_region.append([int(float(sv[0])), int(float(sv[1])), int(sv[3])])
            seg = float(sv[0])
            ins_seq[seg] = sv[2]


    print ('ordered deletion region:', deletion_region)
    return deletion_region, ins_seq

def split_vcf(gene, outdir, deletion_region):
    vcf = '%s/%s.vcf.gz'%(outdir,gene)
    os.system('%s/../bin/tabix -f %s'%(sys.path[0], vcf))
    # if len(deletion_region) == 0:
    #     print ('no sv!')
    #     return 0
    vcf_gap = []
    start = 1001
    break_points_list = [3950]
    # for dele in deletion_region:
    #     if abs(start - dele[0]) < 1:
    #         continue
    #     break_points_list.append(dele[0])
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
        order = "%s/../bin/bcftools filter -t %s:%s-%s %s -o %s/%s_part_%s_%s_%s.vcf"%(sys.path[0],gene, gap[0], gap[1], vcf, outdir, gene, i, gap[0], gap[1])
        os.system(order)
        i+=1
    print (vcf_gap)
    return break_points_list

def convert(outdir, gene, invcf, seq_list, snp_list):
    dict = {}
    for i in range(1, len(snp_list)):
        if seq_list[0][i-1] == seq_list[0][i]:
            dict[int(snp_list[i][1])] = True
        else:
            dict[int(snp_list[i][1])] = False
    # print (seq_list)
    # print (dict,snp_list)

    formate_vcf = outdir + '/%s.vcf.gz'%(gene)
    bp = open(outdir + '/%s_break_points_spechap.txt'%(gene), 'w')
    print ('#gene   locus   00      01      10      11      points_num      next_locus', file = bp)
    m = VariantFile(invcf)
    out = VariantFile(formate_vcf,'w',header=m.header)
    sample = list(m.header.samples)[0]
    add_block = 1
    i = 0
    used_locus = []
    for record in m.fetch():
        # if record.qual < 1
        if record.chrom != gene:
            continue
        if record.samples[sample].phased != True:
            #record.samples[sample]['GT']= (1,1)
            record.samples[sample]['PS'] = add_block
            if record.samples[sample]['GT'] != (1,1) and record.samples[sample]['GT'] != (0,0) and record.samples[sample]['GT'] != (2,2):
                if i > 0 and past_record.pos not in used_locus:
                    used_locus.append(past_record.pos)
                    print (past_record.chrom, past_record.pos, '- - - -', 20, past_record.pos + 100, file = bp)
            record.samples[sample].phased = True
        if record.samples[sample]['PS'] != add_block:
            if i > 0 and past_record.pos not in used_locus:
                used_locus.append(past_record.pos)
                print (past_record.chrom, past_record.pos, '- - - -', 20, past_record.pos + 100, file = bp)
            add_block = record.samples[sample]['PS']
        if record.samples[sample]['GT'] != (1,1) and record.samples[sample]['GT'] != (0,0)  and record.samples[sample]['GT'] != (2,2):
            i += 1
            past_record = record
        # print (record.pos, break_points)
        out.write(record)
    m.close()
    out.close()
    bp.close()

def merge_break_points(outdir, gene, snp_list,use_pstrain_bp):
    bk_list = []
    if not use_pstrain_bp:
        ps = open(outdir + '/%s_break_points_spechap.txt'%(gene), 'r')
        for line in ps:
            if line[0] == '#':
                continue
            array = line.strip().split()
            bk_list.append(int(array[1]))
        ps.close()
    if use_pstrain_bp:
        ps = open(outdir + '/%s_break_points_pstrain.txt'%(gene), 'r')
        for line in ps:
            if line[0] == '#':
                continue
            array = line.strip().split()
            bk_list.append(int(array[1]))
        ps.close()
    bk_list = sorted(bk_list)
    bp = open(outdir + '/%s_break_points.txt'%(gene), 'w')
    print ('#gene   locus   00      01      10      11      points_num      next_locus', file = bp)
    points_num = 0
    for i in range(len(snp_list)):
    # for i in range(len(bk_list)):
        points_num += 1
        if int(snp_list[i][1]) in bk_list:
            # print (snp_list[i])
            if i < len(snp_list) -1:
                print (gene, int(snp_list[i][1]), '- - - -', points_num, snp_list[i+1][1], file = bp)
            else:
                print (gene, int(snp_list[i][1]), '- - - -', points_num, snp_list[i][1], file = bp)
            points_num = 0          
    bp.close()

def phase_insertion(gene, outdir, hla_ref, shdir):
    order = """
    sample=%s
    outdir=%s
    ref=%s
    cat $outdir/newref_insertion.freebayes.vcf|grep '#'>$outdir/filter_newref_insertion.freebayes.vcf
    awk -F'\t' '{if($6>5) print $0}' $outdir/newref_insertion.freebayes.vcf|grep -v '#' >>$outdir/filter_newref_insertion.freebayes.vcf
    %s/../bin/ExtractHAIRs --triallelic 1 --mbq 4 --mmq 0 --indels 1 \
    --ref $ref --bam $outdir/newref_insertion.bam --VCF $outdir/filter_newref_insertion.freebayes.vcf --out $outdir/$sample.fragment.file > spec.log 2>&1
    sort -n -k3 $outdir/$sample.fragment.file >$outdir/$sample.fragment.sorted.file
    bgzip -f $outdir/filter_newref_insertion.freebayes.vcf
    tabix -f $outdir/filter_newref_insertion.freebayes.vcf.gz
    %s/../bin/SpecHap --window_size 15000 -N --vcf $outdir/filter_newref_insertion.freebayes.vcf.gz --frag $outdir/$sample.fragment.sorted.file --out $outdir/$sample.insertion.phased.raw.vcf
    cat $outdir/$sample.insertion.phased.raw.vcf| sed -e 's/1\/1/1\|1/g'>$outdir/$sample.insertion.phased.vcf
    bgzip -f $outdir/$sample.insertion.phased.vcf
    tabix -f $outdir/$sample.insertion.phased.vcf.gz
    """%(gene, outdir, hla_ref, sys.path[0], shdir)
    os.system(order)
    print ('insertion phasing done.')

def clean(outdir, gene):


    formate_vcf = outdir + '/%s.vcf.gz'%(gene)
    new_vcf = outdir + '/%s.new.vcf.gz'%(gene)
    m = VariantFile(formate_vcf)
    out = VariantFile(new_vcf,'w',header=m.header)
    sample = list(m.header.samples)[0]

    for record in m.fetch():
        # if record.qual < 1
        if record.chrom != gene:
            continue
        record.samples[sample].phased = False
        out.write(record)
    m.close()
    out.close()
    return new_vcf


if __name__ == "__main__":   
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})
    else:       
        bamfile,outdir,snp_dp,indel_len,freq_bias=\
            args.bamfile,args.outdir,args.snp_dp,args.indel_len,args.freq_bias           
        snp_qual,gene,fq1,fq2,vcffile = args.snp_qual,args.gene,args.fq1,args.fq2,args.vcf
        strainsNum = 2
        germline_flag = False
        if not os.path.exists(outdir):
            os.system('mkdir '+ outdir) 
        new_formate = False
        if args.hic_fwd != 'NA' or args.tenx != 'NA':
            new_formate = True


        sv_dict = long_InDel_breakpoints(args.sv)
        if gene in sv_dict.keys():
            sv_result = sv_dict[gene]
        else:
            sv_result = []
        # deletion_region, ins_seq = find_deletion_region(sv_result)
        deletion_region, ins_seq = find_deletion_region_BK(sv_result)

        
        ######PStrain-filter-SNV
        snp_list,beta_set,allele_set,snp_index_dict = read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,gene,\
            freq_bias,strainsNum,deletion_region, snp_qual)   
        os.system('%s/../bin/tabix -f %s/middle.vcf.gz'%(sys.path[0], outdir))   
        if germline_flag == True: 
            args.weight = 0.0
        if len(snp_list)==0:
            print ('No heterozygous locus, no need to phase.')
            gene_profile = no_snv_gene_phased(vcffile, outdir, gene, strainsNum)
        else:  
            ######PStrain-phasing
            delta_set=second_beta(bamfile,snp_list,snp_index_dict,outdir)   
            # for i in range(len(delta_set)):
            #     print (snp_list[i], delta_set[i])
            fir_beta,sec_beta=rectify(snp_list,beta_set,delta_set,args.lambda1,args.lambda2,\
                germline_flag)
            if strainsNum==0:
                wo=Workflow(fir_beta,sec_beta,delta_set,args.weight,args.elbow,allele_set)
                final_alpha,seq_list,loss=wo.choose_k()
                strainsNum = len(final_alpha)
            else:
                wo=Workflow(fir_beta,sec_beta,delta_set,args.weight,args.elbow,allele_set)
                final_alpha,seq_list,loss=wo.given_k(strainsNum)  
            generate_break_points(outdir,snp_list,delta_set,strainsNum)
            output(outdir,final_alpha,seq_list,snp_list,gene,freq_bias)
            use_pstrain_bp = True


            ######SpecHap-phasing
            hla_ref = '%s/../db/ref/hla.ref.extend.fa'%(sys.path[0])

            if args.phase_flag == True:
                my_new_vcf = '%s/%s.vcf.gz'%(outdir, gene)
                os.system('%s/../bin/tabix -f %s'%(sys.path[0], my_new_vcf))
                # os.system('gzip -f -d %s'%(my_new_vcf))
                # my_new_vcf = '%s/%s.vcf'%(outdir, gene)
                
                ###############NGS data
                if new_formate:
                    order = '%s/../bin/ExtractHAIRs --new_format 1 --triallelic 1 --indels 1 --ref %s --bam %s --VCF %s --out %s/fragment.file'%(sys.path[0], hla_ref, bamfile, my_new_vcf, outdir)
                else:
                    order = '%s/../bin/ExtractHAIRs --triallelic 1 --indels 1 --ref %s --bam %s --VCF %s --out %s/fragment.file'%(sys.path[0], hla_ref, bamfile, my_new_vcf, outdir)

                # print (order)
                os.system(order)

                
                os.system('cat %s/fragment.file %s/fragment.add.file>%s/fragment.all.file'%(outdir, outdir, outdir))
                # os.system('cat %s/fragment.file >%s/fragment.all.file'%(outdir, outdir))
                
                order='%s/../bin/SpecHap --window_size 15000 --vcf %s --frag %s/fragment.sorted.file --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene)
                ###############PACBIO data
                if args.tgs != 'NA':
                    tgs = """
                    fq=%s
                    ref=%s
                    outdir=%s
                    bin=%s/../bin
                    sample=my
                    $bin/minimap2 -a $ref $fq > $outdir/$sample.tgs.sam
                    $bin/samtools view -F 2308 -b -T $ref $outdir/$sample.tgs.sam > $outdir/$sample.tgs.bam
                    $bin/samtools sort $outdir/$sample.tgs.bam -o $outdir/$sample.tgs.sort.bam
                    $bin/ExtractHAIRs --triallelic 1 --pacbio 1 --indels 1 --ref $ref --bam $outdir/$sample.tgs.sort.bam --VCF %s --out $outdir/fragment.tgs.file
                    """%(args.tgs, hla_ref, outdir, sys.path[0], my_new_vcf)
                    print ('extract linkage info from pacbio TGS data.')
                    os.system(tgs)
                    os.system('cat %s/fragment.tgs.file >> %s/fragment.all.file'%(outdir, outdir))
                    order = '%s/../bin/SpecHap -P --window_size 15000 --vcf %s --frag %s/fragment.sorted.file --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene)
                    # print (order)

                if args.nanopore != 'NA':
                    tgs = """
                    fq=%s
                    ref=%s
                    outdir=%s
                    bin=%s/../bin
                    sample=my
                    $bin/minimap2 -a $ref $fq > $outdir/$sample.tgs.sam
                    $bin/samtools view -F 2308 -b -T $ref $outdir/$sample.tgs.sam > $outdir/$sample.tgs.bam
                    $bin/samtools sort $outdir/$sample.tgs.bam -o $outdir/$sample.tgs.sort.bam
                    $bin/ExtractHAIRs --triallelic 1 --ONT 1 --indels 1 --ref $ref --bam $outdir/$sample.tgs.sort.bam --VCF %s --out $outdir/fragment.nanopore.file
                    # python %s/whole/edit_linkage_value.py $outdir/fragment.raw.nanopore.file 0 $outdir/fragment.nanopore.file
                    # rm $outdir/fragment.nanopore.file
                    # touch $outdir/fragment.nanopore.file

                    """%(args.nanopore, hla_ref, outdir, sys.path[0], my_new_vcf, sys.path[0])
                    print ('extract linkage info from nanopore TGS data.')
                    os.system(tgs)
                    os.system('cat %s/fragment.nanopore.file >> %s/fragment.all.file'%(outdir, outdir))
                    order = '%s/../bin/SpecHap -N --window_size 15000 --vcf %s --frag %s/fragment.sorted.file --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene)
                    

                if args.hic_fwd != 'NA' and args.hic_rev != 'NA':
                    tgs = """
                    fwd_hic=%s
                    rev_hic=%s
                    ref=%s
                    outdir=%s
                    bin=%s/../bin
                    sample=my
                    group='@RG\tID:'$sample'\tSM:'$sample
                    $bin/bwa mem -5SP -Y -U 10000 -L 10000,10000 -O 7,7 -E 2,2 $ref $fwd_hic $rev_hic >$outdir/$sample.tgs.raw.sam
                    cat $outdir/$sample.tgs.raw.sam|grep -v 'XA:'|grep -v 'SA:'>$outdir/$sample.tgs.sam
                    $bin/samtools view -F 2308 -b -T $ref $outdir/$sample.tgs.sam > $outdir/$sample.tgs.bam
                    $bin/samtools sort $outdir/$sample.tgs.bam -o $outdir/$sample.tgs.sort.bam
                    $bin/ExtractHAIRs --new_format 1 --triallelic 1 --hic 1 --indels 1 --ref $ref --bam $outdir/$sample.tgs.sort.bam --VCF %s --out $outdir/fragment.hic.file
                    # python %s/whole/edit_linkage_value.py $outdir/fragment.raw.hic.file 10 $outdir/fragment.hic.file
                    # rm $outdir/fragment.hic.file
                    # touch $outdir/fragment.hic.file
                    """%(args.hic_fwd, args.hic_rev, hla_ref, outdir, sys.path[0], my_new_vcf, sys.path[0])
                    print ('extract linkage info from HiC data.')
                    os.system(tgs)
                    os.system('cat %s/fragment.hic.file >> %s/fragment.all.file'%(outdir, outdir))
                    # my_new_vcf = '%s/%s.new.vcf.gz'%(outdir, gene)
                    order = '%s/../bin/SpecHap -H --new_format --window_size 15000 --vcf %s --frag %s/fragment.sorted.file --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene)
                    print (order)

                if args.tenx != 'NA':
                    # tgs = """
                    # fq=%s
                    # ref=%s
                    # outdir=%s
                    # bin=%s/../bin
                    # sample=my
                    # gene=%s
                    # echo $gene
                    # if [ $gene == "HLA_A" ];
                    # then
                    #     $bin/longranger wgs --id=1 --fastqs=$fq --reference=%s/../db/ref/refdata-hla.ref.extend --sample=%s --sex m --localcores=8 --localmem=32 --jobmode=local --vconly
                    # fi
                    # $bin/ExtractHAIRs --new_format 1 --triallelic 1 --10X 1 --indels 1 --ref $ref --bam ./1/outs/phased_possorted_bam.bam --VCF %s --out $outdir/fragment.tenx.file
                    # if [ $gene == "HLA_DRB1" ];
                    # then
                    # rm -r ./1
                    # fi
                    # """%(args.tenx, hla_ref, outdir, sys.path[0], gene, sys.path[0], args.sample_id, my_new_vcf)
                    # bam=./bam/%s/outs/phased_possorted_bam.bam

                    tgs = """
                    fq=%s
                    ref=%s
                    outdir=%s
                    bin=%s/../bin
                    sample=%s
                    gene=%s
                    echo $gene
                    
                    if [ $gene == "HLA_A" ];
                        then
                            $bin/longranger wgs --id=1 --fastqs=$fq --reference=%s/../db/ref/refdata-hla.ref.extend --sample=$sample --sex m --localcores=8 --localmem=32 --jobmode=local --vconly
                    fi
                    if [ $gene == "HLA_DRB1" ];
                    then
                    rm -r ./1
                    fi
                    bam=./1/outs/phased_possorted_bam.bam

                    $bin/extractHAIRS  --new_format 1 --triallelic 1 --10X 1 --indels 1 --ref $ref --bam $bam --VCF %s --out $outdir/fragment.tenx.file
                    $bin/BarcodeExtract $bam $outdir/barcode_spanning.bed
                    bgzip -f -c $outdir/barcode_spanning.bed > $outdir/barcode_spanning.bed.gz
                    $bin/tabix -f -p bed $outdir/barcode_spanning.bed.gz
                    
                    """%(args.tenx, hla_ref, outdir, sys.path[0], args.sample_id, gene, sys.path[0], my_new_vcf)
                    print ('extract linkage info from 10 X data.')
                    # print (tgs)
                    os.system(tgs)

                    os.system('cat %s/fragment.tenx.file > %s/fragment.all.file'%(outdir, outdir))
                    order = '%s/../bin/SpecHap -T --frag_stat %s/barcode_spanning.bed.gz --new_format --window_size 15000 --vcf %s --frag %s/fragment.sorted.file\
                     --out %s/%s.specHap.phased.vcf'%(sys.path[0],outdir,my_new_vcf, outdir, outdir,gene)

                    # os.system('cat %s/fragment.tenx.file.hapcut.3 >> %s/fragment.all.file'%(outdir, outdir))
                    # order = '%s/../bin/SpecHap -H --window_size 15000 --vcf %s --frag %s/fragment.sorted.file --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene)

                    print (order)
                if new_formate:
                    os.system('sort -n -k6 %s/fragment.all.file >%s/fragment.sorted.file'%(outdir, outdir))
                else:
                    os.system('sort -n -k3 %s/fragment.all.file >%s/fragment.sorted.file'%(outdir, outdir))

            

                os.system('%s/../bin/tabix -f %s'%(sys.path[0], my_new_vcf))
                os.system(order)

                convert(outdir, gene, '%s/%s.specHap.phased.vcf'%(outdir,gene), seq_list, snp_list)
                os.system('%s/../bin/tabix -f %s'%(sys.path[0], my_new_vcf))

                os.system('cat %s/%s_break_points_spechap.txt'%(outdir, gene))

                

                use_pstrain_bp = False



            
            ######Split blocks
            merge_break_points(outdir, gene, snp_list, use_pstrain_bp)
            if gene == 'HLA_DRB1':
                split_vcf(gene, outdir, deletion_region)

            ######Link blocks    
            print ('Start link blocks with database...')
            if gene == 'HLA_DRB1':
                reph='perl %s/whole/rephase.DRB1.pl %s/%s_break_points.txt\
                    %s %s %s/%s_break_points_phased.txt %s %s'%(sys.path[0], outdir,gene,outdir,strainsNum,outdir,\
                    gene,args.block_len,args.points_num)
            else:
                reph='perl %s/whole/rephaseV1.pl %s/%s_break_points.txt\
                    %s %s %s/%s_break_points_phased.txt %s %s'%(sys.path[0],outdir,gene,outdir,strainsNum,outdir,\
                    gene,args.block_len,args.points_num)
            os.system(str(reph))


            if args.phase_flag == True:
                seq_list = read_spechap_seq('%s/%s.vcf.gz'%(outdir, gene), snp_list)

            update_seqlist=newphase(outdir,final_alpha,seq_list,snp_list,vcffile,gene)   #need to refresh alpha with new result.
            fresh_alpha = newalpha(update_seqlist, sec_beta, strainsNum, allele_set, args.weight, fir_beta)
            freq_output(outdir, gene, fresh_alpha, germline_flag)
            gene_profile=gene_phased(update_seqlist,snp_list,gene)
            print ('Phasing of %s is done! Haplotype ratio is %s:%s'%(gene, fresh_alpha[0], fresh_alpha[1]))

        ######link long indels
        if len(ins_seq) > 0:
            ins_seq = segment_mapping_pre(fq1, fq2, ins_seq, outdir, gene, hla_ref)
            # print (ins_seq)
            segment_mapping(fq1, fq2, ins_seq, outdir, gene, hla_ref)
        else:
            os.system('cp %s/%s.bam %s/newref_insertion.bam'%(outdir, gene.split('_')[-1], outdir))
            os.system('%s/../bin/samtools index %s/newref_insertion.bam'%(sys.path[0], outdir))
            os.system('cp %s %s/newref_insertion.freebayes.vcf'%(vcffile, outdir))
        deletion_region = sv_copy_number_old(outdir, deletion_region, gene, ins_seq)
        # deletion_region = sv_copy_number(deletion_region, sv_result)

        if gene == 'HLA_DRB1':
            dup_region_type(outdir, strainsNum, bamfile)
            dup_file = outdir +'/select.DRB1.seq.txt'
        if len(ins_seq) > 0:
            phase_insertion(gene, outdir, args.ref, sys.path[0])

        sh = Share_reads(deletion_region, outdir, strainsNum, gene, gene_profile, ins_seq)
        # print (deletion_region)
        sh.split_seg()


