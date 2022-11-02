#!/bin/bash


# alignment for pacbio reads 
# param:
# 1: reads file
# 2: prefix 
# 3: reference (hg19)
function align_pacbio()
{
    minimap2 -ax map-pb $3 $1 > $prefix.sam
    samtools view -F 2308 -b -T $3 $prefix.sam > $prefix.bam
    rm $prefix.sam
}

# alignment for nanopore reads 
# param:
# 1: reads file
# 2: prefix 
# 3: reference (hg19)
function align_nanopore()
{
    minimap2 -ax map-ont $3 $1 > $prefix.sam
    samtools view -F 2308 -b -T $3 $prefix.sam > $prefix.bam
    rm $prefix.sam
}

# alignment for hic paired-end reads 
# param:
# 1: forward reads
# 2: reverse reads
# 3: prefix 
# 4: reference (hg19)
function align_hic()
{
    bwa mem -t 8 -5SP $4 $1 $2 | samtools view -F 2304 -b -T $4 | samtools sort --thread 8 -n | samtools fixmate -O bam - ${3}.fixmate.bam
    java -jar picard.jar MarkDuplicates I=./${3}.fixmate.bam O=./${3}.bam M=./${3}.markdup.txt
    rm  ./${3}.fixmate.bam
}

# alignment for ngs paired-end reads 
# param:
# 1: forward reads
# 2: reverse reads
# 3: prefix 
# 4: reference (hg19)
function align_ngs()
{
    bwa mem -t 8 -Y -M $4 $1 $2 | samtools view -F 2304 -b -T $4 | samtools sort --thread 8 -n | samtools fixmate -O bam - ${3}.fixmate.bam
    java -jar picard.jar MarkDuplicates I=./${3}.fixmate.bam O=./${3}.bam M=./${3}.markdup.txt
    rm  ./${3}.fixmate.bam
} 

# alignment for ngs paired-end reads 
# param:
# 1: fastq_dir
# 2: reference_dir
# 3: prefix 
# 4: gatk path
function align_tenx()
{
    longranger wgs \
    --id=$3 \
    --sample=$3 \
    --fastqs=$1 \
    --reference=$2 \
    --vcmode=$4 \
    --jobmode=local \
    --localcores=24 \
    --localmem=96
}
