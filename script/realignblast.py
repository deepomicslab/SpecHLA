import pysam
import os
import argparse
import re 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser=argparse.ArgumentParser()
parser.add_argument("-i", required=True, help="input bamfile")
parser.add_argument("-o", required=True, help="Output bamfile")
parser.add_argument("-r", default="rematch.read.format.txt", required=True, help="rematch.read.list")




def read_rematch(readsfile):
    dicts={}
    fp = open(readsfile, 'r')
    for line in fp:
        fields = line.split("\t")
        keys = ' '.join((fields[0], fields[1]))
        value = ' '.join((fields[2], fields[3], fields[4]))
        dicts[keys] = value
    fp.close() 
    return dicts   

def remove_dup_reads(inbam, outbam, readsfile):
    dicts = read_rematch(readsfile)
    bf = pysam.AlignmentFile(inbam, 'rb')   
    outf = pysam.AlignmentFile(outbam, "wb", header=bf.header)

    #find the reads pair that contain xA tag
    dup_reads = []
    for r in bf:
        ifdup = False
        for tag in r.tags:
            if tag[0] == 'XA':
                ifdup = True
                if r.query_name not in dup_reads:
                    dup_reads.append(r.query_name)
    bf.close()


    
    nf = pysam.AlignmentFile(inbam, 'rb')
    for r in nf:
        if r.is_supplementary:
            continue
        if r.is_unmapped:
            continue
        if r.is_read1:
            tag=1
        if r.is_read2:
            tag=2
        key = ' '.join((r.query_name, str(tag)))
        if key in dicts:
            arrs = dicts[key].split(' ')
            r.reference_start = int(arrs[1]) -1
            r.cigarstring = arrs[2]
            r.reference_name = arrs[0]

        if r.reference_name != r.next_reference_name:
            continue

        #remove the pair end reads that contain XA tag
        # if r.query_name not in dup_reads:
        #     outf.write(r)

        #remove the XA tag for reads
        # new_tag = []
        # for tag in r.tags:
        #     if tag[0] != 'XA':
        #         new_tag.append(tag)
        # r.tags = new_tag
        outf.write(r)

    nf.close()
    outf.close()



if __name__ == "__main__":
    args=parser.parse_args()
    if args.i:
        inbam = args.i
    if args.o:
        outbam = args.o
    if args.r:
        readsfile = args.r
    # outbam = 'test.bam'
    remove_dup_reads(inbam, outbam, readsfile)

'''

'''

# for r in bf:
#     if r.is_supplementary:
#         continue
#     if r.is_unmapped:
#         continue
#     if r.is_read1:
#         tag=1
#     if r.is_read2:
#         tag=2
#     print (r)
'''
    key = ' '.join((r.query_name, str(tag)))
    if key in dicts:
        arrs = dicts[key].split(' ')
        r.reference_start = int(arrs[1]) -1
        r.cigarstring = arrs[2]
        r.reference_name = arrs[0]

    if r.reference_name != r.next_reference_name:
        continue



    outf.write(r)

'''
