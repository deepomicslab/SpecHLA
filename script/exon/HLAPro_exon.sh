#!/bin/bash

###
### The Exon version of HLAPro, performs HLA assembly and HLA Typing. 
###
### Usage:
###   sh HLAPro_exon.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -p <Asian>
###
### Options:
###   -n        Sample ID.
###   -1        The first fastq file.
###   -2        The second fastq file.
###   -p        The population of the sample: Asian, Black, or Caucasian. Use mean frequency 
###             if not provided.
###   -f        True or False. The annotation database only includes the alleles with population 
###             frequency higher than zero if set True. Otherwise, it includes all alleles. Default 
###             is True.
###   -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault 
###             is 3. It should be set larger to infer novel alleles.
###   -g        True or False. Split more blocks if set True; then rely more on allele database. 
###             Otherwise, split less blocks. Default is True.
###   -s        True or False. Use SpecHap to phase if True, else use PStrain. Default is True.
###             We recommend SpecHap.
###   -h        Show this message.

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

while getopts ":n:1:2:p:s:m:g:f:" opt; do
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
    g) pstrain_bk="$OPTARG"
    ;;
    f) annotation="$OPTARG"
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


mkdir -p $outdir
group='@RG\tID:'$sample'\tSM:'$sample
# :<<!
## remove the repeat read ##
$bin/python3 $dir/../uniq_read_name.py $fq1 $outdir/$sample.uniq.name.R1.gz
$bin/python3 $dir/../uniq_read_name.py $fq2 $outdir/$sample.uniq.name.R2.gz
fq1=$outdir/$sample.uniq.name.R1.gz
fq2=$outdir/$sample.uniq.name.R2.gz

## map the HLA reads to the allele database ##
$bin/novoalign -d $db/ref/hla_gen.format.filter.extend.DRB.no26789.ndx -f $fq1 $fq2  -o SAM -o FullNW -r All 100000 --mCPU 10 -c 10  -g 20 -x 3  | $bin/samtools view -Sb - | $bin/samtools sort -  > $outdir/$sample.novoalign.bam
$bin/samtools index $outdir/$sample.novoalign.bam

## assign the reads to corresponding gene ##
$bin/python3 $dir/../assign_reads_to_genes.py -o $outdir -b ${outdir}/${sample}.novoalign.bam -nm ${nm:-3}
$bin/python3 $dir/../check_assign.py $fq1 $fq2 $outdir
#rm -rf $outdir/$sample.novoalign.sam

## align the reads to corresponding gene reference ##
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

## realignment ##
sh $dir/run.assembly.realign.sh $sample $outdir/$sample.merge.bam $outdir 70
!

## phase, link blocks, calculate haplotype ratio, give typing results ##
$bin/python3 $dir/../phase_exon.py -b $outdir/$sample.realign.sort.bam -v $outdir/$sample.realign.filter.vcf \
-o $outdir/ -g ${pstrain_bk:-True} --snp_dp 0 --block_len 80 --points_num 1 --freq_bias 0.1 --reads_num 2 -a ${phase:-True} 

# \
# -e ${annotation:-True}

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

echo perl $dir/anno_HLA_pop.pl $sample $outdir 2 $pop $annotation_parameter
perl $dir/anno_HLA_pop.pl $sample $outdir 2 $pop $annotation_parameter
cat $outdir/hla.result.txt
# sh $dir/../clear_output.sh $outdir/
echo $sample is done.
