#!/usr/bin/env python3
from argparse import ArgumentParser
from argparse import ArgumentTypeError
from my_imports import *
import time
import re
from itertools import combinations, permutations

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Please give right flag (True or False).')

Usage = \
"""
python3 single_tgs.py [options] 

Help information can be found by python3 single_tgs.py -h/--help, additional information can be found in \
README.MD or https://github.com/deepomicslab/HLAPro.
"""
scripts_dir=sys.path[0]+'/'
parser = ArgumentParser(description="SpecHLA.",prog='python3 single_tgs.py',usage=Usage)
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
required.add_argument("--gene",help="gene",dest='gene',metavar='', type=str)
required.add_argument("--tgs",help="PACBIO TGS fastq",dest='tgs',metavar='', type=str)
required.add_argument("-o", "--outdir",help="The output directory.",dest='outdir',metavar='')

#alternative parameter
optional.add_argument("--nanopore",help="NANOPORE TGS fastq",dest='nanopore',metavar='', type=str)
optional.add_argument("--hic_fwd",help="fwd_hic.fastq",dest='hic_fwd',metavar='', type=str)
optional.add_argument("--hic_rev",help="rev_hic.fastq",dest='hic_rev',metavar='', type=str)
optional.add_argument("--tenx",help="10X data",dest='tenx',metavar='', type=str)
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
optional.add_argument("--reads_num",help="The number of supporting reads between two adjcent loci\
     lower than this value will be regard as break points.(default is 10)",dest='reads_num',\
     metavar='',default=10, type=int)
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
    # for longshot output snps
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
        dp = 100
        if record.qual == None:
            print ('WARNING: no vcf quality value.')
            continue
        if record.chrom != chrom_name:
            continue
        if dp < snp_dp:
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
        if not (record.filter.keys()[0] == 'PASS' or record.qual > snp_qual):
            continue
        ##########after filter snp################
        snp_index_dict[record.pos] = snp_index
        snp_index += 1

        if geno != (1,1):
            if geno == (1,2):                
                snp = [record.chrom,record.pos,record.alts[0],record.alts[1],record.ref]
            else:
                snp=[record.chrom,record.pos,record.ref,record.alts[0],record.ref]

            reads_list = reads_support(samfile, snp)
            allele_dp = [len(reads_list[0]), len(reads_list[1])]
            new_dp=sum(allele_dp)
            # print (record, allele_dp)
            if new_dp == 0:
                print ('WARNING: the depth of the locus obtained by pysam is zero!',snp)
                continue

            snp_list.append(snp)
            beta_set.append(allele_dp)
            allele_set.append(2)
            record.samples[sample].phased=False
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

def freq_output(outdir, gene, final_alpha, germline_flag):
    # if germline_flag:
    #     final_alpha = [0.5, 0.5]
    ra_file=open(outdir+'/%s_freq.txt'%(gene),'w')    
    print ('# HLA\tFrequency',file=ra_file)
    for j in range(len(final_alpha)):
        print ('str-'+str(j+1),final_alpha[j],file=ra_file)
    ra_file.close()
      
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

def newphase(outdir,seq_list,snp_list,gene):
    file=outdir+'/%s_break_points_phased.txt'%(gene)
    # print (file)
    if os.path.isfile(file):
        block_dict=relate_order_other(file, snp_list)
    else:
        block_dict={gene:{}}
    seq=np.array(seq_list)
    seq=np.transpose(seq)
    snp=snp_list  
    update_seqlist=[]  
    k=2
    if gene in block_dict.keys():
        gene_dict=block_dict[gene]
    else:
        gene_dict={}
    ref_order=[]
    for orde in range(k):
        ref_order.append(orde)
    he=0
    os.system('%s/../bin/tabix -f %s/middle.vcf.gz'%(sys.path[0],outdir))
    m = VariantFile('%s/middle.vcf.gz'%(outdir))
    rephase_file = '%s/%s.rephase.vcf.gz'%(outdir,gene)
    if os.path.isfile(rephase_file):
        os.system('rm %s'%(rephase_file))
    out = VariantFile(rephase_file,'w',header=m.header)
    sample = list(m.header.samples)[0]
    for record in m.fetch():
        geno = record.samples[sample]['GT']   
        # print (record, record.samples[sample].phased) 
        if record.samples[sample].phased != True:
            if geno == (1,2):
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

def gene_phased(update_seqlist,snp_list, gene):
    gene_profile={}
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

def focus_region():
    return {'HLA_A':[1000,4503],'HLA_B':[1000,5081],'HLA_C':[1000,5304],'HLA_DPA1':[1000,10775],\
        'HLA_DPB1':[1000,12468],'HLA_DQA1':[1000,7492],'HLA_DQB1':[1000,8480],'HLA_DRB1':[1000,12229]}

class Share_reads():

    def __init__(self, deletion_region, outdir, strainsNum, gene, gene_profile, ins_seq):
        self.deletion_region = deletion_region
        self.bamfile = outdir + '/newref_insertion.bam'
        self.ins_bam = outdir + '/newins_insertion.bam'
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
                if not read.query_sequence:
                    continue
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
        # print ("normal reads support", reads_support[0], reads_support[1])
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
        samfile = pysam.AlignmentFile(self.ins_bam, "rb")
        for read in samfile.fetch('%s_%s'%(self.gene,self.deletion_region[deletion_index][0])):
            # print (read.query_name)
            for i in range(self.strainsNum):
                if read.query_name in self.reads_support[i]:
                    link_reads[i] += 1
        print ("insertion reads", deletion_index, link_reads[0], link_reads[1])
        return self.most_support(link_reads)

    def deletion_phase(self):
        for deletion_index in range(len(self.deletion_region)):
            # print (self.deletion_region[deletion_index])
            self.deletion_reads(deletion_index)

    def split_seg(self):
        normal_region, segs = self.generate_normal_region()
        id_name = {}
        gap = ''
        contain_dup_seg = ''
        for seg in segs:
            print ('seg', seg)
            if seg[2] == 'insertion':
                continue
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
                    # print ('###########insertion', seg, assign_index, self.deletion_region[deletion_index][2])
                    insertion_seg = '%s_%s'%(self.gene, seg[0])
                    if self.deletion_region[deletion_index][2] == 2:                        
                        insert_seq = self.link_diploid_insertion(insertion_seg)
                        hap_seq += insert_seq[i]                       
                    elif  i == assign_index and self.deletion_region[deletion_index][2] == 1:
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
                if len(insert_phase_result[1][-1][1]) != 1 or len(insert_phase_result[0][-1][1]) != 1:
                    continue          

        samfile = pysam.AlignmentFile(self.ins_bam, "rb")
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
        print ("-------------------------", len(insert_reads_support[0]), len(insert_reads_support[1]))
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
        print ('link sv support reads number', r00, r01)
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
        
def segment_mapping_pre(tgs, ins_seq, outdir, gene, gene_ref):
    newref=outdir+'/newref_insertion.fa'
    os.system('cp %s %s'%(gene_ref, newref))
    for seg in ins_seq.keys():
        f = open(newref, 'a')
        print ('>%s_%s\n%s'%(gene, int(seg), ins_seq[seg]), file = f)
        f.close()
        # index the ref
    alignment_order = f"""
    fq={tgs}
    ref={newref}
    outdir={outdir}
    bin={sys.path[0]}/../bin
    sample='newref_insertion'
    $bin/samtools faidx $ref
    $bin/minimap2 -a $ref $fq > $outdir/$sample.tgs.sam
    $bin/samtools view -F 2308 -b -T $ref $outdir/$sample.tgs.sam > $outdir/$sample.tgs.bam
    $bin/samtools sort $outdir/$sample.tgs.bam -o $outdir/$sample.tgs.sort.bam
    $bin/samtools index $outdir/$sample.tgs.sort.bam
    longshot --bam $outdir/$sample.tgs.sort.bam --ref $ref --out $outdir/$sample.freebayes.vcf -F
    echo longshot --bam $outdir/$sample.tgs.sort.bam --ref $ref --out $outdir/$sample.freebayes.vcf -F
    bgzip -f $outdir/$sample.freebayes.vcf 
    tabix -f $outdir/$sample.freebayes.vcf.gz
    echo alignment done....
    """

    print ('New mapping starts to link long InDels.')
    os.system(alignment_order)
    for ins in ins_seq.keys():
        ins_call = """%s/../bin/samtools faidx %s/newref_insertion.fa %s_%s |%s/../bin/bcftools consensus -H 1 %s/newref_insertion.freebayes.vcf.gz  >%s/fresh_ins.fa
        """%(sys.path[0],outdir,gene,int(ins),sys.path[0],outdir,outdir)
        print (ins_call)
        os.system(ins_call)
        ins_seq[ins] = read_fasta('%s/fresh_ins.fa'%(outdir))
    return ins_seq
  
def segment_mapping(tgs, ins_seq, outdir, gene, gene_ref):
    newref=outdir+'/newref_insertion.fa'
    os.system('cp %s %s'%(gene_ref, newref))
    for seg in ins_seq.keys():

        f = open(newref, 'a')
        print ('>%s_%s\n%s'%(gene, int(seg), ins_seq[seg]), file = f)
        f.close()
        # index the ref
    print ('New mapping starts to link long InDels.')
    map_call = f"""
    fq={tgs}
    ref={newref}
    outdir={outdir}
    bin={sys.path[0]}/../bin
    sample='newref_insertion'
    $bin/samtools faidx $ref
    $bin/minimap2 -a $ref $fq > $outdir/$sample.sam
    $bin/samtools view -F 2308 -b -T $ref $outdir/$sample.sam > $outdir/$sample.unsort.bam
    $bin/samtools sort $outdir/$sample.unsort.bam -o $outdir/$sample.bam
    $bin/samtools index $outdir/$sample.bam
    longshot --bam $outdir/$sample.bam --ref $ref --out $outdir/$sample.freebayes.vcf -F
    echo new realignment done.
    """
    os.system(map_call)
    map_call = f"""
    fq={tgs}
    ref={outdir}/fresh_ins.fa
    outdir={outdir}
    bin={sys.path[0]}/../bin
    sample='newins_insertion'
    $bin/samtools faidx $ref
    $bin/minimap2 -a $ref $fq > $outdir/$sample.sam
    $bin/samtools view -F 2308 -b -T $ref $outdir/$sample.sam > $outdir/$sample.unsort.bam
    $bin/samtools sort $outdir/$sample.unsort.bam -o $outdir/$sample.bam
    $bin/samtools index $outdir/$sample.bam
    longshot --bam $outdir/$sample.bam --ref $ref --out $outdir/$sample.ins.freebayes.vcf -F
    echo new realignment done.
    """
    os.system(map_call)

def sv_copy_number(deletion_region, sv_list):
    print (deletion_region)
    print (sv_list)
    for i in range(len(deletion_region)):
        deletion_region[i] += [0]
    for i in range(len(deletion_region)):
        for sv in sv_list:
            if deletion_region[i][0] >= int(sv[0]) and deletion_region[i][1] <= int(sv[1]):
                deletion_region[i][2] += sv[3]
                # break
        if deletion_region[i][2] > 2:
            deletion_region[i][2] = 2
        print (deletion_region[i], '### copy number', deletion_region[i][2])
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
        sv = [array[1], array[4], array[6], int(array[7])]
        if array[0] not in sv_dict:
            sv_dict[array[0]] = [sv]
        else:
            sv_dict[array[0]].append(sv)
    return sv_dict

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
            deletion_region.append([int(float(sv[0])), int(float(sv[1])), int(sv[3])])
            seg = float(sv[0])
            ins_seq[seg] = sv[2]
    print ('ordered deletion region:', deletion_region)
    return deletion_region, ins_seq

def split_vcf(gene, outdir, deletion_region):
    vcf = '%s/%s.vcf.gz'%(outdir,gene)
    os.system('%s/../bin/tabix -f %s'%(sys.path[0], vcf))
    vcf_gap = []
    start = 1001
    break_points_list = [3950]
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
    # print (vcf_gap)
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

def phase_insertion(gene, outdir, hla_ref, shdir):
    order = """
    sample=%s
    outdir=%s
    ref=%s/fresh_ins.fa
    cat $outdir/newins_insertion.ins.freebayes.vcf|grep '#'>$outdir/filter_newref_insertion.freebayes.vcf
    awk -F'\t' '{if($6>5) print $0}' $outdir/newins_insertion.ins.freebayes.vcf|grep -v '#' >>$outdir/filter_newref_insertion.freebayes.vcf
    %s/../bin/ExtractHAIRs --triallelic 1 --pacbio 1 --mbq 4 --mmq 0 --indels 1 \
    --ref $ref --bam $outdir/newref_insertion.bam --VCF $outdir/filter_newref_insertion.freebayes.vcf --out $outdir/$sample.fragment.file > spec.log 2>&1
    sort -n -k3 $outdir/$sample.fragment.file >$outdir/$sample.fragment.sorted.file
    bgzip -f $outdir/filter_newref_insertion.freebayes.vcf
    tabix -f $outdir/filter_newref_insertion.freebayes.vcf.gz
    %s/../bin/SpecHap --window_size 15000 -N --vcf $outdir/filter_newref_insertion.freebayes.vcf.gz --frag $outdir/$sample.fragment.sorted.file --out $outdir/$sample.insertion.phased.raw.vcf
    cat $outdir/$sample.insertion.phased.raw.vcf| sed -e 's/1\/1/1\|1/g'>$outdir/$sample.insertion.phased.vcf
    bgzip -f $outdir/$sample.insertion.phased.vcf
    tabix -f $outdir/$sample.insertion.phased.vcf.gz
    """%(gene, outdir, outdir, sys.path[0], shdir)
    os.system(order)
    print ('insertion phasing done.')

if __name__ == "__main__":   
    if len(sys.argv)==1:
        print (Usage%{'prog':sys.argv[0]})
    else:     
        
        bamfile,outdir,snp_dp,indel_len,freq_bias=\
            args.bamfile,args.outdir,args.snp_dp,args.indel_len,args.freq_bias           
        snp_qual,gene,vcffile = args.snp_qual,args.gene,args.vcf
        strainsNum = 2
        if not os.path.exists(outdir):
            os.system('mkdir '+ outdir) 
        sv_dict = long_InDel_breakpoints(args.sv)
        if gene in sv_dict.keys():
            sv_result = sv_dict[gene]
        else:
            sv_result = []
        
        deletion_region, ins_seq = find_deletion_region(sv_result)

        
        ######phase small variants
        snp_list,beta_set,allele_set,snp_index_dict = read_vcf(vcffile,outdir,snp_dp,bamfile,indel_len,gene,\
            freq_bias,strainsNum,deletion_region, snp_qual)   
        # hla_ref = '%s/../db/ref/hla.ref.extend.fa'%(sys.path[0])
        hla_ref = '%s/../db/ref/%s.fa'%(sys.path[0], gene)
        if len(snp_list)==0:
            print ('No heterozygous locus, no need to phase.')
            gene_profile = no_snv_gene_phased(vcffile, outdir, gene, strainsNum)
        else:
            my_new_vcf = '%s/middle.vcf.gz'%(outdir)
            os.system('%s/../bin/tabix -f %s'%(sys.path[0], my_new_vcf))
            extract_order = '%s/../bin/ExtractHAIRs --pacbio 1 --triallelic 1 --indels 1 --ref %s --bam %s --VCF %s --out %s/fragment.file'%(sys.path[0], hla_ref, bamfile, my_new_vcf, outdir)
            os.system(extract_order)  
            print (extract_order)  
            os.system('sort -n -k3 %s/fragment.file >%s/fragment.sorted.file'%(outdir, outdir))           
            spec_order='%s/../bin/SpecHap --window_size 15000 --vcf %s --frag %s/fragment.sorted.file --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene) 
            os.system(spec_order)
            seq_list = read_spechap_seq('%s/%s.specHap.phased.vcf'%(outdir, gene), snp_list)           
            convert(outdir, gene, '%s/%s.specHap.phased.vcf'%(outdir,gene), seq_list, snp_list)

            if gene == 'HLA_DRB1':
                split_vcf(gene, outdir, deletion_region)
            reph='perl %s/whole/rephaseV1.pl %s/%s_break_points_spechap.txt\
                %s %s %s/%s_break_points_phased.txt %s %s'%(sys.path[0],outdir,gene,outdir,strainsNum,outdir,\
                gene,args.block_len,args.points_num)
            os.system(str(reph))
            update_seqlist=newphase(outdir,seq_list,snp_list,gene)  
            # print (update_seqlist) 
            gene_profile=gene_phased(update_seqlist,snp_list,gene)


        ######link long indels
        if len(ins_seq) > 0:
            print ("segment_mapping_pre")
            ins_seq = segment_mapping_pre(args.tgs, ins_seq, outdir, gene, hla_ref)
            # print ("-----------------", ins_seq)
            segment_mapping(args.tgs, ins_seq, outdir, gene, hla_ref)
            print ("segment_mapping")
        else:
            os.system('cp %s/%s.bam %s/newref_insertion.bam'%(outdir, gene.split('_')[-1], outdir))
            os.system('%s/../bin/samtools index %s/newref_insertion.bam'%(sys.path[0], outdir))
            os.system('cp %s %s/newref_insertion.freebayes.vcf'%(vcffile, outdir))
        deletion_region = sv_copy_number(deletion_region, sv_result)
        if len(ins_seq) > 0:
            phase_insertion(gene, outdir, args.ref, sys.path[0])
        sh = Share_reads(deletion_region, outdir, strainsNum, gene, gene_profile, ins_seq)
        sh.split_seg()
        print (f"{gene} gene is done.")


