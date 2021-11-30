#!/bin/bash


###
### The extended version of SpecHLA, performs HLA assembly and HLA Typing with full-length basd on NGS and TGS data.
### This script can use PacBio, Nanopore, Hi-C, and 10X sequencing data to improve the phasing performance if provided.
###
###
### Usage:
###   sh SpecHLA_tgs.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -t <sample.pacbio.fq.gz> -p <Asian>
###
### Options:
###   -n        Sample ID. <required>
###   -1        The first fastq file. <required>
###   -2        The second fastq file. <required>
###   -t        Pacbio TGS fastq file.
###   -e        Nanopore TGS fastq file.
###   -x        Path of folder created by 10x demultiplexing. Prefix of the filenames of FASTQs
###             should be the same as Sample ID. You can regard reads as normal NGS reads and use 
###             this parameter to adopt barcode information to improve phasing.
###   -c        fwd hi-c fastq file.
###   -d        rev hi-c fastq file.
###   -p        The population of the sample: Asian, Black, or Caucasian. Use mean frequency
###             if not provided.
###   -f        True or False. The annotation database only includes the alleles with population 
###             frequency higher than zero if set True. Otherwise, it includes all alleles. Default 
###             is True.
###   -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
###             is 2. It should be set larger to infer novel alleles.
###   -v        True or False. Consider long InDels if True, else only consider short variants. 
###             Default is False.
###   -q        Minimum variant quality. Default is 0.01. Set larger in high quality samples.
###   -s        True or False. Use SpecHap to phase if True, else use PStrain. Default is True.
###             We recommend use SpecHap.
###   -a        Use this long InDel file if provided.
###   -r        The minimum Minor Allele Frequency (MAF), default is 0.05.
###   -h        Show this message.

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

while getopts ":n:1:2:p:f:m:s:v:q:t:a:e:x:c:d:r:y:" opt; do
  case $opt in
    n) sample="$OPTARG"
    ;;
    1) fq1="$OPTARG"
    ;;
    2) fq2="$OPTARG"
    ;;
    p) pop="$OPTARG"
    ;;
    s) phase="$OPTARG"
    ;;
    m) nm="$OPTARG"
    ;;
    f) annotation="$OPTARG"
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
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


dir=$(cd `dirname $0`; pwd)
bin=$dir/../../bin
db=$dir/../../db
hlaref=$db/ref/hla.ref.extend.fa
outdir=$(pwd)/output/$sample

echo Start profiling HLA for $sample. 
mkdir -p $outdir
group='@RG\tID:'$sample'\tSM:'$sample
# :<<!
$bin/python3 $dir/../uniq_read_name.py $fq1 $outdir/$sample.uniq.name.R1.gz
$bin/python3 $dir/../uniq_read_name.py $fq2 $outdir/$sample.uniq.name.R2.gz
fq1=$outdir/$sample.uniq.name.R1.gz
fq2=$outdir/$sample.uniq.name.R2.gz

echo map the reads to database to assign reads to corresponding genes.
$bin/novoalign -d $db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx -f $fq1 $fq2 -F STDFQ -o SAM -o FullNW -r All 100000 --mCPU 10 -c 10  -g 20 -x 3  | $bin/samtools view -Sb - | $bin/samtools sort -  > $outdir/$sample.novoalign.bam

$bin/samtools index $outdir/$sample.novoalign.bam
$bin/python3 $dir/../assign_reads_to_genes.py -o $outdir -b ${outdir}/${sample}.novoalign.bam -nm ${nm:-2}

$bin/python3 $dir/../check_assign.py $fq1 $fq2 $outdir

$bin/bwa mem -U 10000 -L 10000,10000 -R $group $hlaref $fq1 $fq2 | $bin/samtools view -H  >$outdir/header.sam
#hlas=(A B C)
hlas=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
for hla in ${hlas[@]}; do
        hla_ref=$db/HLA/HLA_$hla/HLA_$hla.fa
        # $bin/bwa mem -U 10000 -L 10000,10000 -O 7,7 -E 2,2 -R $group $hla_ref $outdir/$hla.R1.fq.gz $outdir/$hla.R2.fq.gz | $bin/samtools view -bS -F 0x800 -| $bin/samtools sort - >$outdir/$hla.bam
        $bin/bwa mem -U 10000 -L 10000,10000 -R $group $hla_ref $outdir/$hla.R1.fq.gz $outdir/$hla.R2.fq.gz | $bin/samtools view -bS -F 0x800 -| $bin/samtools sort - >$outdir/$hla.bam
        $bin/samtools index $outdir/$hla.bam
done
#samtools merge -f -h $outdir/header.sam $outdir/$sample.merge.bam $outdir/A.bam $outdir/B.bam $outdir/C.bam 
$bin/samtools merge -f -h $outdir/header.sam $outdir/$sample.merge.bam $outdir/A.bam $outdir/B.bam $outdir/C.bam $outdir/DPA1.bam $outdir/DPB1.bam $outdir/DQA1.bam $outdir/DQB1.bam $outdir/DRB1.bam
$bin/samtools index $outdir/$sample.merge.bam

echo start realignment.
#sh /mnt/disk2_workspace/wangmengyao/NeedleHLA/select_wgs/realign/run.assembly.realign.sh $sample $outdir/$sample.merge.bam $outdir 70
sh $dir/run.assembly.realign.sh $sample $outdir/$sample.merge.bam $outdir 70

# !
bam=$outdir/$sample.realign.sort.bam
vcf=$outdir/$sample.realign.filter.vcf


if [ ${long_indel:-False} == True ]
  then
  port=$(date +%N|cut -c5-9)
  sh $dir/../ScanIndel/run_scanindel_sample.sh $sample $bam $outdir $port
  bfile=$outdir/Scanindel/$sample.breakpoint.txt

  if [ ${tgs:-NA} != NA ]
    then
    $bin/pbmm2 align $hlaref ${tgs:-NA} $outdir/$sample.movie1.bam --sort --preset HIFI --sample $sample --rg '@RG\tID:movie1'
    $bin/pbsv discover $outdir/$sample.movie1.bam $outdir/$sample.svsig.gz
    $bin/pbsv call $hlaref $outdir/$sample.svsig.gz $outdir/$sample.var.vcf
    $bin/python3 $dir/vcf2bp.py $outdir/$sample.var.vcf $outdir/$sample.tgs.breakpoint.txt
    cat $outdir/$sample.tgs.breakpoint.txt >>$bfile
  fi
else
  bfile=nothing
fi

if [ ${sv:-NA} != NA ]
  then
  bfile=$sv
fi
echo $bfile

echo start haplotyping.

bam=$outdir/$sample.realign.sort.bam
vcf=$outdir/$sample.realign.filter.vcf

# hlas=(A)
hlas=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
for hla in ${hlas[@]}; do
hla_ref=$db/ref/HLA_$hla.fa
$bin/python3 $dir/../phase_tgs.py \
-o $outdir \
-b $bam \
-s $bfile \
-v $vcf \
--fq1 $outdir/$hla.R1.fq.gz \
--fq2 $outdir/$hla.R2.fq.gz \
--gene HLA_$hla \
-a ${phase:-True} \
--freq_bias ${maf:-0.05} \
--block_len 200 --points_num 1 --reads_num 2 --snp_qual ${snp_quality:-0.01} \
--ref $hla_ref \
--tgs ${tgs:-NA} \
--nanopore ${nanopore_data:-NA} \
--hic_fwd ${hic_data_fwd:-NA} \
--hic_rev ${hic_data_rev:-NA} \
--tenx ${tenx_data:-NA} \
--sa $sample
done

# --fq1 $fq1 \
# --fq2 $fq2 \

# hlas=(DRB1)
# for hla in ${hlas[@]}; do
# bfile=/mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20201116/Scanindel/$sample/$sample.breakpoint.txt
# #bfile=/mnt/disk2_workspace/wangmengyao/NeedleHLA/select_wgs/Scanindel/breakpoint.txt
# #bfile=$indir/Scanindel/$sample.breakpoint.txt
# hla_ref=/home/wangmengyao/scripts/PHLAT/database/ref/extend/HLA_$hla.fa
# strainsNum=2
# /home/wangshuai/softwares/Python-3.7.0/bin/python3 /mnt/disk2_workspace/wangshuai/00.strain/08.NeedleHLA/scripts/PHLA_WGS_no_realign_phase_sv8_yonghan3.py  \
# -k $strainsNum \
# -o $outdir \
# -b $bam \
# -s $bfile \
# -v $vcf \
# --fq1 $indir/$hla.R1.fq.gz \
# --fq2 $indir/$hla.R2.fq.gz \
# --gene HLA_$hla \
# -d F \
# -g T \
# --freq_bias 0.01 \
# --rlen 150 \
# --block_len 200 --points_num 1 --reads_num 2 \
# --ref $hla_ref \
# --snp_qual 5
# done


echo start annotation.
if [ ${annotation:-True} == True ]
then
  if [ ${phase:-True} == True ]
  then
    annotation_parameter=spechap
  else
    annotation_parameter=pstrain
  fi
else
  annotation_parameter=all
fi

# echo perl $dir/annoHLApop.pl $sample $outdir $outdir 2 $pop $annotation_parameter
perl $dir/annoHLApop.pl $sample $outdir $outdir 2 $pop $annotation_parameter

# sh $dir/../clear_output.sh $outdir/
cat $outdir/hla.result.txt
echo $sample is done.
