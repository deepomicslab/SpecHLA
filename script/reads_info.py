import pysam
import numpy as np

class Reads():
    def __init__(self,vcffile,bamfile,conn_file):
        self.vcffile=vcffile
        self.bamfile=bamfile
        self.conn_file=conn_file
    def extract_snp(self):
        #need to filter the vcf file, maybe in the calling step.
        snp_list,freq_set=[],[]
        called_num,dp_sum=0,0
        f=open(self.vcffile,'r')
        for line in f:
            if line[0] != "#":
                called_num+=1
                line=line.strip()
                array=line.split()
                # snp_tag=str(array[0])+"_"+str(array[1])
                sub_array=array[-1].split(":")
                # print (sub_array,line)
                mini_array=sub_array[1].split(",")
                alt_dp=mini_array[1]
                dp=sub_array[2]
                dp_sum+=float(dp)
                if float(dp) == 0:
                    alt_freq=0
                else:
                    alt_freq=float(alt_dp)/float(dp)            
                freq_set.append(round(alt_freq,6))
                single_snp=[]
                single_snp.append(array[0])
                single_snp.append(array[1])
                single_snp.append(array[3])
                single_snp.append(array[4])
                snp_list.append(single_snp)
        print ("mean depth is",int(float(dp_sum)/len(snp_list)))
        return freq_set,snp_list
    def extract_delta(self):
        freq_set,snp_list=self.extract_snp()
        delta_set=init_delta(len(freq_set))
        f=open(self.conn_file,'r')
        for line in f:
            line=line.strip()
            array=line.split()
            if array[0] == '1': 
                delta_index=int(array[2])-1
                geno_type=array[3]
                for i in range(len(geno_type)-1):
                # delta_set[delta_index][int(geno_type[0])][int(geno_type[1])]+=1
                    delta_set[delta_index+i][int(geno_type[i])][int(geno_type[i+1])]+=1
                #print (geno_type[0],geno_type[1])
        # print (delta_set[:10])
        for delta in delta_set:
            sum_dp=sum(delta[0])+sum(delta[1])
            # print (sum_dp)
            if sum_dp>4:
                delta[0][0]=float(delta[0][0])/sum_dp
                delta[0][1]=float(delta[0][1])/sum_dp
                delta[1][0]=float(delta[1][0])/sum_dp
                delta[1][1]=float(delta[1][1])/sum_dp
            else:
                delta[:]=[[0,0],[0,0]]
        return freq_set,delta_set

def reads_support(samfile,first):   
    reads_list=[]
    allele_num=len(first[3])+1
    for i in range(allele_num):
        reads_list.append([])
    num=0
    for read in samfile.fetch(str(first[0]),int(first[1])-1,int(first[1])):
        # if int(first[1]) == 1799:
        #     print (read.cigartuples)
        num+=1
        if int(first[1])-1 in read.get_reference_positions(full_length=True) and read.mapping_quality >10:   
                  
            reads_index=read.get_reference_positions(full_length=True).index(int(first[1])-1)
            start=reads_index
            if len(first[3][0])>1:
                end=reads_index+len(first[3][0])
                # start=reads_index-len(first[3][0])+1
                # end=reads_index+1
            elif len(first[2])>1:
                end=reads_index+len(first[2])
            else:
                end=reads_index+1
            
            # ref_allele=read.get_reference_sequence()[start:end].upper()
            read_allele=read.query_sequence[start:end].upper()
            # if int(first[1]) == 1799:
            #     print (read.cigartuples,read_allele,read.get_reference_positions(full_length=True))
            flag=True
            for i in range(len(read_allele)):
                point_flag=isin(int(first[1])-1+i,read.get_reference_positions(full_length=True))
                if point_flag == False:
                    flag=False
                
            alt=False   
            cigar_count=0
            for cig in read.cigartuples:
                if (cig[0] == 1 or cig[0] == 2) and cigar_count+cig[1] >= read.get_reference_positions().index(int(first[1])-1) and read.get_reference_positions().index(int(first[1])-1) >cigar_count:
                    alt=True
                cigar_count+=cig[1]
                # print (read_allele,read.cigartuples,alt,read.get_reference_positions().index(int(first[1])-1))

            ref_flag=False
            if read.cigartuples[0] == (0, 150):
                ref_flag=True
            # if int(first[1]) == 5659:
            #     print (alt,first[2],first[3][0],read_allele,read.cigartuples)


            if len(first[2]) >= len(first[3][0]):
                if (read_allele==first[2]  and flag == True) :#or ref==True :
                    reads_list[0].append(read.query_name)
                else:
                    for i in range(allele_num-1):
                        if read_allele==first[3][i] or len(first[2])>1:
                            reads_list[1+i].append(read.query_name)
            else:                   
                if (read_allele==first[3][0] or alt == True) and ref_flag==False:
                    reads_list[1].append(read.query_name)
                else:
                    reads_list[0].append(read.query_name)
    # if int(first[1]) == 5659:
    #     print (len(reads_list[0]),len(reads_list[1]),read_allele,num)
    return reads_list

def share_reads(samfile,left,right,new_left):
    left_num=len(left[3])+1
    right_num=len(right[3])+1
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


if __name__ == '__main__':
    bamfile='/mnt/disk2_workspace/wangshuai/00.strain/01.real_strains/sample100/gene_k2_5/gene_k2_5.ref.bam'
    samfile = pysam.AlignmentFile(bamfile, "rb")
    left=['gi|506938955|ref|NZ_AOQL01000165.1|:c1645-1',1513,'G',('C')]
    right=['gi|506938955|ref|NZ_AOQL01000165.1|:c1645-1',1521,'G',('A')]
    share_reads(samfile,left,right)
    
    # reads_support(samfile,first)