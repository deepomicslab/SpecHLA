#!/bin/bash

###
### SpecHLA: Full resolution HLA typing with paired-end, PacBio, Nanopore,
### Hi-C, and 10X data. Supports WGS, WES, and RNASeq.
### 
###
### Usage:
###   sh SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -t <sample.pacbio.fq.gz> -p <Asian>
###
### Options:
###   -n        Sample ID. <required>
###   -1        The first fastq file. <required>
###   -2        The second fastq file. <required>
###   -o        The output folder to store the typing results.
###   -u        Choose full-length or exon typing. 0 indicates full-length, 1 means exon, 
###             default is to perform full-length typing.
###   -p        The population of the sample [Asian, Black, Caucasian, Unknown, nonuse]. 
###             Default is Unknown, meaning use mean frequency. nonuse indicates only adopting 
###             mapping score and considering alleles with frequency as zero. 
###   -t        Pacbio fastq file.
###   -e        Nanopore fastq file.
###   -c        fwd hi-c fastq file.
###   -d        rev hi-c fastq file.
###   -x        Path of folder created by 10x demultiplexing. Prefix of the filenames of FASTQs
###             should be the same as Sample ID. Please install Longranger in the system env.
###   -w        How to use linkage info from allele imbalance [0, 0.5, 1], default is 0 that means 
###             not use, 0.5 means use both reads and imbalance info, 1 means only use imbalance info.
###   -j        Number of threads [5]
###   -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
###             is 2. It should be set larger to infer novel alleles.
###   -y        The minimum different mapping score between the best- and second-best aligned gene. 
###             Discard the read if the score is lower than this value. Deault is 0.1. 
###   -v        True or False. Consider long InDels if True, else only consider short variants. 
###             Default is False. 
###   -q        Minimum variant quality. Default is 0.01. Set it larger in high quality samples.
###   -s        Minimum variant depth. Default is 5.
###   -a        Use this long InDel file if provided.
###   -r        The minimum Minor Allele Frequency (MAF), default is 0.05 for full length and
###             0.1 for exon typing.
###   -g        Whether use G-translate in annotation [1|0], default is 0.
###   -k        The mean depth in a window lower than this value will be masked by N, default is 5.
###             Set 0 to avoid masking.
###   -z        Whether only mask exon region, True or False, default is False.
###   -f        The trio infromation; child:parent_1:parent_2 [Example: NA12878:NA12891:NA12892]. 
###             Note: this parameter should be used after performing SpecHLA once.
###   -b        Whether use database for phasing [1|0], default is 1.
###   -h        Show this message.

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

while getopts ":n:1:2:p:f:m:v:q:t:a:e:x:c:d:r:y:o:j:w:u:s:g:k:z:y:f:b:" opt; do
  case $opt in
    n) sample="$OPTARG"
    ;;
    1) fq1="$OPTARG"
    ;;
    2) fq2="$OPTARG"
    ;;
    p) pop="$OPTARG"
    ;;
    m) nm="$OPTARG"
    ;;
    v) long_indel="$OPTARG"
    ;;
    q) snp_quality="$OPTARG"
    ;;
    t) tgs="$OPTARG"
    ;;
    a) sv="$OPTARG"
    ;;
    e) nanopore_data="$OPTARG"
    ;;
    x) tenx_data="$OPTARG"
    ;;
    c) hic_data_fwd="$OPTARG"
    ;;
    d) hic_data_rev="$OPTARG"
    ;;
    r) maf="$OPTARG"
    ;;
    o) given_outdir="$OPTARG"
    ;;
    j) num_threads="$OPTARG"
    ;;
    w) weight_imb="$OPTARG"
    ;;
    u) focus_exon="$OPTARG"
    ;;
    s) snp_dp="$OPTARG"
    ;;
    g) trans="$OPTARG"
    ;;
    k) mask_depth="$OPTARG"
    ;;
    z) mask_exon="$OPTARG"
    ;;
    y) mini_score="$OPTARG"
    ;;
    f) trio="$OPTARG"
    ;;
    b) use_database="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


dir=$(cd `dirname $0`; pwd)
export LD_LIBRARY_PATH=$dir/../../spechla_env/lib
bin=$dir/../../bin
db=$dir/../../db
hlaref=$db/ref/hla.ref.extend.fa

if [ ${given_outdir:-NA} == NA ]
  then
    outdir=$(pwd)/output/$sample
  else
    outdir=$given_outdir/$sample   
fi
focus_exon_flag=${focus_exon:-0} # default is perform full-length typing

echo Start profiling HLA for $sample. 
mkdir -p $outdir
# exec >$outdir/$sample.log 2>&1 #redirect log info to the outdir
group='@RG\tID:'$sample'\tSM:'$sample
echo use ${num_threads:-5} threads.

bash $dir/SpecHLA.sh -1 $fq1 -2 $fq2 -p $pop -o $given_outdir -m $nm -j $num_threads -n $sample -z $mask_exon
python3 $dir/top_allele_2_reads.py $given_outdir/$sample/
bash $dir/SpecHLA.sh -1 $fq1 -2 $fq2 -p $pop -o $given_outdir -m $nm -j $num_threads -n $sample -z $mask_exon -t $given_outdir/$sample/top_allele.fastq
