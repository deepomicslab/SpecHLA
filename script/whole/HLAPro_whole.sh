#!/bin/bash

sample=$1
fq1=$2
fq2=$3
pop=$4
#outdir=$5

dir=$(cd `dirname $0`; pwd)
bin=$dir/../../bin
db=$dir/../../db
hlaref=$db/ref/hla.ref.extend.fa
outdir=$(pwd)/output/$sample

mkdir -p $outdir
group='@RG\tID:'$sample'\tSM:'$sample
#:<<!
$bin/python3 $dir/../uniq_read_name.py $fq1 $outdir/$sample.uniq.name.R1.gz

$bin/python3 $dir/../uniq_read_name.py $fq2 $outdir/$sample.uniq.name.R2.gz
fq1=$outdir/$sample.uniq.name.R1.gz
fq2=$outdir/$sample.uniq.name.R2.gz
$bin/novoalign -d $db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx -f $fq1 $fq2  -o SAM -o FullNW -r All 100000 --mCPU 10 -c 10  -g 20 -x 3  | $bin/samtools view -Sb - | $bin/samtools sort -  > $outdir/$sample.novoalign.bam

$bin/samtools index $outdir/$sample.novoalign.bam
$bin/samtools view -h -F 0x800 -F 4 -F 8 $outdir/$sample.novoalign.bam | $bin/samtools sort -n -O SAM -o $outdir/$sample.novoalign.sam
$bin/python3 $dir/../assign_reads_to_genes.py -o $outdir -b ${outdir}/${sample}.novoalign.bam -nm 2

$bin/python3 $dir/../check_assign.py $fq1 $fq2 $outdir
rm -rf $outdir/$sample.novoalign.sam
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

sh $dir/run.assembly.realign.sh $sample $outdir/$sample.merge.bam $outdir 70
bam=$outdir/$sample.realign.sort.fixmate.bam
vcf=$outdir/$sample.realign.filter.vcf

port=$(date +%N|cut -c5-9)
sh $dir/../ScanIndel/run_scanindel_sample.sh $sample $bam $outdir $port

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
-g F \
--freq_bias 0.05 \
--rlen 150 \
--block_len 200 --points_num 1 --reads_num 2 --snp_qual 0.01 \
--ref $hla_ref
done

perl $dir/annoHLApop.pl \
$sample $outdir $outdir 2 $pop

