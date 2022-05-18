#!/bin/bash


###
### 
### SpecHLA performs haplotyping with only PacBio or Nanopore reads.
###
###
### Usage:
###   sh SpecHLA_single.sh -n <sample> -t <sample.pacbio.fq.gz> -p <Asian>
###
### Options:
###   -n        Sample ID. <required>
###   -t        Pacbio TGS fastq file.
###   -e        Nanopore TGS fastq file.
###   -p        The population of the sample: Asian, Black, or Caucasian. Use mean frequency
###             if not provided.
###   -f        True or False. The annotation database only includes the alleles with population 
###             frequency higher than zero if set True. Otherwise, it includes all alleles. Default 
###             is True.
###   -v        True or False. Consider long InDels if True, else only consider short variants. 
###             Default is True.
###   -q        Minimum variant quality. Default is 0.01. Set larger in high quality samples.
###   -s        True or False. Use SpecHap to phase if True, else use PStrain. Default is True.
###             We recommend use SpecHap.
###   -r        The minimum Minor Allele Frequency (MAF), default is 0.05.
###   -h        Show this message.

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

while getopts ":n:p:f:m:s:v:q:t:a:e:x:c:d:r:y:" opt; do
  case $opt in
    n) sample="$OPTARG"
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

python3 $dir/../general_pipeline.py $sample $tgs $outdir
$bin/minimap2 -a $hlaref $tgs | $bin/samtools view -H  >$outdir/header.sam

hlas=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
for hla in ${hlas[@]}; do
        hla_ref=$db/HLA/HLA_$hla/HLA_$hla.fa
        $bin/minimap2 -a $hla_ref $outdir/$hla.fq.gz | $bin/samtools view -bS -F 0x800 -| $bin/samtools sort - >$outdir/$hla.bam
        $bin/samtools index $outdir/$hla.bam
done
$bin/samtools merge -f -h $outdir/header.sam $outdir/$sample.merge.bam $outdir/A.bam $outdir/B.bam $outdir/C.bam $outdir/DPA1.bam $outdir/DPB1.bam $outdir/DQA1.bam $outdir/DQB1.bam $outdir/DRB1.bam
$bin/samtools index $outdir/$sample.merge.bam
longshot --bam $outdir/$sample.merge.bam --ref $db/ref/hla.ref.extend.fa --out $outdir/$sample.longshot.vcf -F
bgzip -f $outdir/$sample.longshot.vcf
$bin/tabix -f $outdir/$sample.longshot.vcf.gz

bam=$outdir/$sample.merge.bam
vcf=$outdir/$sample.longshot.vcf.gz

echo start haplotyping.

# hlas=(A)
bfile=nothing
hlas=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
for hla in ${hlas[@]}; do

if [ ${long_indel:-True} != False ]
  then
bfile=$outdir/$sample.$hla.sv.breakpoint.txt
hla_ref=$db/HLA/HLA_$hla/HLA_$hla.fa
$bin/pbmm2 align $hla_ref ${tgs:-NA} $outdir/$sample.movie1.bam --sort --preset HIFI --sample $sample --rg '@RG\tID:movie1'
$bin/pbsv discover $outdir/$sample.movie1.bam $outdir/$sample.svsig.gz
$bin/pbsv call -m 150 --max-ins-length 4000 -t DEL,INS $hla_ref $outdir/$sample.svsig.gz $outdir/$sample.var.vcf
echo $bin/pbsv call -m 150 --max-ins-length 4000 -t DEL,INS $hla_ref $outdir/$sample.svsig.gz $outdir/$sample.var.vcf
python3 $dir/vcf2bp.py $outdir/$sample.var.vcf $bfile
echo -----------------
fi

hla_ref=$db/ref/HLA_$hla.fa
python3 $dir/../single_tgs.py \
-o $outdir \
-b $bam \
-s $bfile \
-v $vcf \
--gene HLA_$hla \
--freq_bias ${maf:-0.05} \
--block_len 200 --points_num 1 --reads_num 2 --snp_qual ${snp_quality:-0.01} \
--ref $hla_ref \
--tgs $outdir/$hla.fq.gz \
--nanopore ${nanopore_data:-NA} \
--hic_fwd ${hic_data_fwd:-NA} \
--hic_rev ${hic_data_rev:-NA} \
--tenx ${tenx_data:-NA} \
--snp_qual 5 \
--sa $sample
done


echo start annotation.


rm $outdir/hla.allele.*.HLA_*.fasta.fai
perl $dir/annoHLApop.pl $sample $outdir $outdir 2 $pop spechap

# sh $dir/../clear_output.sh $outdir/
cat $outdir/hla.result.txt
echo $sample is done.
