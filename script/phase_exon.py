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
    ,dest='k',metavar='', default=0, type=int)
required.add_argument("-o", "--outdir",help="The output directory.",dest='outdir',metavar='')
#alternative parameter
flag_parser.add_argument("-g", "--balance",help="The switch to determine if the sample is balance.\
     True for balance samples and False for imbalance samples(default is False)."\
     ,dest='germline',metavar='', default=False, type=str2bool)
flag_data.add_argument("-d", "--database",help="the switch to determine whether use the information\
     of prior database. true for using (default is false)."\
    ,dest='database_flag',metavar='', default=False,type=str2bool)    
optional.add_argument("-p", "--population",help="The population type of the given sample.\
     The type should be one of Asian, Black or Caucasian. Default is Asian."\
    ,dest='popu',metavar='',default='Asian', type=population_bool)
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
optional.add_argument("--allele_support",help="The minimum depth of allele to be considered as hete in HLAtyping\
     step (default is 2).",dest='allele_support',metavar='',default=2, type=int)     
optional.add_argument("--indel_len",help="The maximum length for indel to be considered in HLAtyping\
     step (default is 150).",dest='indel_len',metavar='',default=150, type=int)
optional.add_argument("--block_len",help="The minimum length for block to be considered in final\
     result (default is 300).",dest='block_len',metavar='',default=50, type=int)
optional.add_argument("--points_num",help="The minimum hete loci number for block to be considered\
     in final result (default is 2).",dest='points_num',metavar='',default=3, type=int)
optional.add_argument("--rlen",help="rlen",dest='rlen',metavar='',default=150, type=int)
optional.add_argument("--prefix",help="sample name",dest='prefix',metavar='',default='sample', type=str)
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

def exon_region():
    file = '/mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20200318/exon/exon_extent.bed'
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

def read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,chrom_name,freq_bias,strainsNum,mean_depth, std_depth):
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
    print ("hete loci number",len(snp_list))
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
    print ('##############################')
    ra_file=open(outdir+'/%s_freq.txt'%(gene),'w')    
    print ('# HLA\tFrequency',file=ra_file)
    seq=np.array(seq_list)
    seq=np.transpose(seq)
    print (len(snp_list),len(seq),'points num')
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
        chrom, locus, record.ref, record.alts = update_locus_for_exon(record.chrom, float(record.pos), exon_region_accumulate_list, record.ref, record.alts)
        record.chrom, record.pos = chrom, int(locus)
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
    os.system('tabix %s/%s.rephase.vcf.gz'%(outdir,gene))
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
    sample = list(in_vcf.header.samples)[0]
    for record in in_vcf.fetch():
        if record.chrom == gene:
            chrom, pos = update_locus_for_exon(record.chrom, float(record.pos), exon_region_accumulate_list)
            record.chrom, record.pos = chrom, int(pos)
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
            break_points_index.append(i)
            points_number=0
    bp.close()
    return break_points_index

def all_possibilities(break_points_index, strainsNum, snp_list, seq_list, final_alpha, gene): #need to update alpha finally
    my_table = all_table( len(break_points_index) + 1, 2)
    locus_seq = '' 
    for i in range(len(snp_list)):
        if i in break_points_index:
            locus_seq += '|' + str(i)
            locus_seq += '_'
        else:
            locus_seq += '|' + str(i)
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
    # print (locus_set, index_set)
    # print ('mytable', my_table)
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
    # print (seq_list)
    # print (new_seq_list_set)
                
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
    os.system('tabix %s/%s.%s.rephase.vcf.gz'%(outdir, gene, index))
    exon_ref = '/home/wangmengyao/scripts/NeedleHLA/script/database/exon/hla_exon_ref.fasta'
    gene_ref = ' /%s.new.fasta'%(gene)
    for i in range(1, k+1):
        fastq2 = '/home/wangmengyao/scripts/NeedleHLA/script/bin/samtools faidx %s %s | \
        /home/wangmengyao/miniconda2/bin/bcftools consensus -H %s %s/%s.%s.rephase.vcf.gz\
         >%s/%s.%s.%s.fasta'%(exon_ref, gene, i, outdir, gene, index, outdir, gene, index,i)
        os.system(fastq2)
        # fa = '%s/%s.%s.%s.fasta'%(outdir, gene, index,i)
        # makedb = '/home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/makeblastdb -in %s -dbtype nucl\
        #  -parse_seqids -out %s'%(fa, fa)
        # blast = '/home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/blastn -query %s -out %s/%s.blast.%s.%s.out\
        #      -db %s -outfmt 7 -num_threads 4 -max_target_seqs 3500'%(fa, outdir, gene, index, i, gene_ref)
        # os.system(makedb)
        # os.system(blast)

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

def coverage_filter(bamfile, outdir):
    f = open(outdir + '/unreliable.genes.txt', 'w')
    d = open(outdir + '/mean.depth.genes.txt', 'w')
    print ('#gene mean median std', file = d)
    exon_region_list = exon_region()
    depth_file = '%s/bam.depth'%(outdir)
    depth_order = 'samtools depth %s > %s'%(bamfile, depth_file)
    os.system(depth_order)
    depth_set = []
    for gene in ['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']:
        locus_list = []
        depth_list = []
        for exon in exon_region_list:
            if exon[0] == gene:
                locus_list.append([exon[1], exon[2]])
        low_depth_window = []
        for line in open(depth_file, 'r'):
            line = line.strip()
            array = line.split()
            array[1] = float(array[1])
            array[2] = float(array[2])
            if array[0] != gene:
                continue
            
            exon_flag = False
            for locus in locus_list:
                if array[1] >= float(locus[0]) and array[1] <= float(locus[1]):
                    exon_flag = True
            if exon_flag == False:
                continue
            depth_set.append(array[2])
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
                print ('%s coverage dirty!'%(gene), low_depth_window)
                low_depth_window = []
                # break
        print (gene, round(np.mean(depth_list),1), np.median(depth_list), np.std(depth_list), file = d)
    f.close()
    d.close()
    depth_set = np.array(depth_set)
    mean_depth, std_depth = np.mean(depth_set), np.std(depth_set)
    print ("depth, mean, std, median", mean_depth, std_depth, np.median(depth_set))
    return mean_depth, std_depth
            
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
        
def choose_allele():
    f = open('/home/wangmengyao/packages/hla-polysolver/data/HLA_FREQ.txt', 'r')
    i = 1
    for line in f:
        line = line.strip()
        array = line.split()
        if i == 1:
            freq_dict = {}
            name_dict = {}
            for j in range(1, len(array)):
                freq_dict[j] = {}
                name_dict[array[j]] = j
                # print (freq_dict, name_dict)
        else:
            allele_name = array[0].upper()
            name_array = allele_name.split('_')
            allele_name =  '%s*%s:%s'%(name_array[1], name_array[2], name_array[3])            # ':'.join(allele_name.split('_'))
            for j in range(1, len(array)):
                freq_dict[j][allele_name] = float(array[j])
            # print (allele_name)
        i += 1
    # print (freq_dict)
    return freq_dict, name_dict

def read_blast():
    blast_file = '/home/wangmengyao/scripts/NeedleHLA/1000WES/noise_num/NA18943/HLA_A.blast.1.1.out'
    for line in open(blast_file, 'r'):
        if line[0] == '#':
            continue
        array = line.strip().split()
        match = re.search('(.*?)\*(\d+):(\d+)', array[1])
        allele_name = '%s*%s:%s'%(match.group(1), match.group(2), match.group(3))
        # print (line)
        # print (allele_name)
        return allele_name, float(array[2])


if __name__ == "__main__":
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})
    else: 
        # '''
        raw_bamfile,outdir,strainsNum,snp_dp,germline_flag,database_flag,indel_len,freq_bias=\
            args.bamfile,args.outdir,args.k,args.snp_dp,args.germline,\
            args.database_flag,args.indel_len,args.freq_bias           
        os.system('rm %s/*fasta'%(outdir))
        if not os.path.exists(outdir):
            os.system('mkdir '+ outdir) 
        #######realign and find break points for SV
        bamfile = raw_bamfile
        vcffile = args.vcf
        realign_and_sv_break.main_no_realign(raw_bamfile, outdir, args.rlen, 'WES')
        # fq1 = outdir +'/sample.read1.fastq.gz'
        # fq2 = outdir +'/sample.read2.fastq.gz'
        # dup_file = outdir +'/select.DRB1.seq.txt'
        sco = open(outdir + '/edge.score.txt', 'w')
        #######realign and find break points for SV
        # for gene in ['HLA_A:1504-1773']:
        # gene = args.gene
        mean_depth, std_depth = coverage_filter(bamfile, outdir)
        
        #for gene in ['HLA_A','HLA_B','HLA_C', 'HLA_DQB1','HLA_DRB1']:
        os.system('rm %s/*.rephase.vcf.gz*'%(outdir))
        os.system('rm %s/*.*.fasta'%(outdir))
        for gene in ['HLA_A','HLA_B','HLA_C','HLA_DQB1','HLA_DRB1','HLA_DQA1','HLA_DPA1','HLA_DPB1']:
        #for gene in ['HLA_A']:

            ######for SNV
            snp_list,beta_set,allele_set = read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,gene,\
                freq_bias,strainsNum,mean_depth, std_depth)     
            if germline_flag == True: 
                #For germline, the first-order genotype frequency is 0.5 for all locus, 
                #thus will not provide any info.
                args.weight = 0.0
            if len(snp_list)==0:
                #same seq for two haps, need change with >2 haps
                print ('No hete locus')
                if strainsNum==0:
                    strainsNum = 1
                no_snv_gene_phased(vcffile, outdir, gene, strainsNum)
                delta_set, seq_list, final_alpha = [], [[],[]], [0.5, 0.5] 
                print (gene, 0, 0, 0, 0, file = sco)
            else:  
                delta_set=second_beta(bamfile,snp_list)   
                r_score = edge_score(delta_set)  #just for double haps
                print (gene, r_score, len(delta_set), end= '\t', file = sco)

                for i in range(len(snp_list)-1):
                    print (snp_list[i], beta_set[i], delta_set[i])
                    # index_sort=np.argsort(delta_set[i])
                    # if index_sort[-1] + index_sort[-2] == 3:
                    #     delta_set[i][index_sort[-3]] = 0
                    # # for j in range(4):
                    # #     if delta_set[i][j] > 3:
                    # #         delta_set[i][j] = 6
                    # #     else:
                    # #         delta_set[i][j] = 0
                    # print (snp_list[i], beta_set[i], delta_set[i])

                fir_beta,sec_beta=rectify(snp_list,beta_set,delta_set,args.lambda1,args.lambda2,\
                    germline_flag,database_flag)
                for i in range(len(snp_list)-1):
                    print (snp_list[i], fir_beta[i], sec_beta[i])
                if strainsNum==0:
                    wo=Workflow(fir_beta,sec_beta,delta_set,args.weight,args.elbow,allele_set)
                    final_alpha,seq_list,loss=wo.choose_k()
                    strainsNum = len(final_alpha)
                else:
                    wo=Workflow(fir_beta,sec_beta,delta_set,args.weight,args.elbow,allele_set)
                    final_alpha,seq_list,loss = wo.given_k(strainsNum)  
                print (loss, end= '\t', file = sco)
                if len(snp_list)==0:
                    print (0, file = sco)
                else:
                    print (float(loss)/len(snp_list), file = sco)
                print ('Final alpha is %s.'%(final_alpha))
                output(outdir,final_alpha,seq_list,snp_list,gene)
            break_points_index = generate_break_points(outdir,snp_list,delta_set,strainsNum)
            print ("#", len(break_points_index), file = sco)
            if len(break_points_index) > 10:
                break_points_index = break_points_index[:10]
            all_possibilities(break_points_index, strainsNum, snp_list, seq_list, final_alpha, gene)   
            print ('%s Done!'%(gene))
            # read_blast()
            # freq_dict, name_dict = choose_allele()
        #annotation = 'perl /mnt/disk2_workspace/wangshuai/00.strain/08.NeedleHLA/wes/new_assign/select_db/anno_HLA_pop.pl %s 2 %s '%(outdir, args.popu)
        annotation = 'perl /home/wangmengyao/scripts/NeedleHLA/script/anno_HLA_pop.pl %s %s %s %s '%(args.prefix, outdir, strainsNum, args.popu)
        os.system(annotation)
        os.system('cat %s/hla.result.txt'%(outdir))
        #os.system('rm %s/*.rephase.vcf.gz*'%(outdir))
        os.system('rm %s/HLA_*.*.fasta'%(outdir))





        #         reph='perl /home/wangmengyao/scripts/NeedleHLA/script/rephase_exonV1.pl %s/%s_break_points.txt\
        #             %s %s %s/%s_break_points_phased.txt %s %s'%(outdir,gene,outdir,strainsNum,outdir,gene,\
        #             args.block_len,args.points_num)
        #         os.system(str(reph))
        #         update_seqlist=newphase(outdir,final_alpha,seq_list,snp_list,vcffile,gene) 
        # type_order = 'perl /home/wangmengyao/scripts/NeedleHLA/script/annoHLAexon.pl sample %s %s %s'%(outdir, outdir, strainsNum)
        # os.system(type_order) 
        # sco.close()
        # '''       



