#!/bin/bash
###
### The Whole version of HLAPro, performs HLA assembly and HLA Typing with full-length.
###
### Usage:
###   sh HLAPro_whole.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -p <Asian>
###
### Options:
###   -n        Sample ID..
###   -1        The first fastq file.
###   -2        The second fastq file.
###   -p        The population of the sample: Asian, Black, or Caucasus. Use mean frequency
###             if not provided.
###   -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
###             is 2. It should be set larger to infer novel alleles.
###   -g        True or False. The ratio of two haplotype will be regarded the same if set True.
###             Otherwise, we will calculate the haplotype ratio. Default is False.
###   -h        Show this message.

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

while getopts ":n:1:2:p:" opt; do
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
    g) germline="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
#sample=$1
#fq1=$2
#fq2=$3
#pop=$4
#outdir=$5

dir=$(cd `dirname $0`; pwd)
bin=$dir/../../bin
db=$dir/../../db
hlaref=$db/ref/hla.ref.extend.fa
outdir=$(pwd)/output/$sample

echo Start profiling HLA for $sample. 
mkdir -p $outdir
group='@RG\tID:'$sample'\tSM:'$sample
#:<<!
$bin/python3 $dir/../uniq_read_name.py $fq1 $outdir/$sample.uniq.name.R1.gz

$bin/python3 $dir/../uniq_read_name.py $fq2 $outdir/$sample.uniq.name.R2.gz
fq1=$outdir/$sample.uniq.name.R1.gz
fq2=$outdir/$sample.uniq.name.R2.gz

echo map the reads to database to assign reads to corresponding genes.
$bin/novoalign -d $db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx -f $fq1 $fq2  -o SAM -o FullNW -r All 100000 --mCPU 10 -c 10  -g 20 -x 3  | $bin/samtools view -Sb - | $bin/samtools sort -  > $outdir/$sample.novoalign.bam

$bin/samtools index $outdir/$sample.novoalign.bam
#$bin/samtools view -h -F 0x800 -F 4 -F 8 $outdir/$sample.novoalign.bam | $bin/samtools sort -n -O SAM -o $outdir/$sample.novoalign.sam
$bin/python3 $dir/../assign_reads_to_genes.py -o $outdir -b ${outdir}/${sample}.novoalign.bam -nm ${nm:-2}

$bin/python3 $dir/../check_assign.py $fq1 $fq2 $outdir
#rm -rf $outdir/$sample.novoalign.sam
$bin/bwa mem -U 10000 -L 10000,10000 -R $group $hlaref $fq1 $fq2 | $bin/samtools view -H  >$outdir/header.sam
#hlas=(A B C)
hlas=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
for hla in ${hlas[@]}; do
        hla_ref=$db/HLA/HLA_$hla/HLA_$hla.fa
        $bin/bwa mem -U 10000 -L 10000,10000 -O 7,7 -E 2,2 -R $group $hla_ref $outdir/$hla.R1.fq.gz $outdir/$hla.R2.fq.gz | $bin/samtools view -bS -F 0x800 -| $bin/samtools sort - >$outdir/$hla.bam
        $bin/samtools index $outdir/$hla.bam
done
#samtools merge -f -h $outdir/header.sam $outdir/$sample.merge.bam $outdir/A.bam $outdir/B.bam $outdir/C.bam 
$bin/samtools merge -f -h $outdir/header.sam $outdir/$sample.merge.bam $outdir/A.bam $outdir/B.bam $outdir/C.bam $outdir/DPA1.bam $outdir/DPB1.bam $outdir/DQA1.bam $outdir/DQB1.bam $outdir/DRB1.bam
$bin/samtools index $outdir/$sample.merge.bam

echo start realignment.
sh $dir/run.assembly.realign.sh $sample $outdir/$sample.merge.bam $outdir 70
bam=$outdir/$sample.realign.sort.fixmate.bam
vcf=$outdir/$sample.realign.filter.vcf

port=$(date +%N|cut -c5-9)
sh $dir/../ScanIndel/run_scanindel_sample.sh $sample $bam $outdir $port

echo start assembly.
bfile=$outdir/Scanindel/$sample.breakpoint.txt
hlas=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
for hla in ${hlas[@]}; do
hla_ref=$db/HLA_$hla.fa
strainsNum=2
$bin/python3 $dir/../phase_whole.py -k $strainsNum \
-o $outdir \
-b $bam \
-s $bfile \
-v $vcf \
--fq1 $outdir/$hla.R1.fq.gz \
--fq2 $outdir/$hla.R2.fq.gz \
--gene HLA_$hla \
-d F \
-g ${germline:-False} \
--freq_bias 0.05 \
--rlen 150 \
--block_len 200 --points_num 1 --reads_num 2 --snp_qual 0.01 \
--ref $hla_ref
done

echo start assnotation.
perl $dir/annoHLApop.pl \
$sample $outdir $outdir 2 $pop
echo $sample is done.
