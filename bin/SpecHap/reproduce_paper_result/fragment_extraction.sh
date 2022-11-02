#!/bin/bash

declair -a chrs=("chr1" "chr21" "chr22")
declair -a samples=("NA12878" "NA19240")

#notice that the extractHAIR_spec refers to the SpecHap refined version
mkdir NGS 10X HiC PACBIO NANOPORE
# NGS 

cd NGS
for chr in  "${chrs[@]}"
do
    extractHAIRS --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --out HG00403_${chr}.lst
    extractHAIRS_spec --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --out HG00403_${chr}_spec.lst
    sort -n -k3 HG00403_${chr}.lst > HG00403_${chr}.sorted.lst
    sort -n -k3 HG00403_${chr}_spec.lst > HG00403_${chr}_spec.sorted.lst
    mv HG00403_${chr}.sorted.lst HG00403_${chr}.lst
    mv HG00403_${chr}_spec.sorted.lst HG00403_${chr}_spec.lst
done


for sample in  "${samples[@]}"
do
    extractHAIRS --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}.lst
    extractHAIRS_spec --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}_spec.lst
    sort -n -k3 ${sample}.lst > ${sample}.sorted.lst
    sort -n -k3 ${sample}_spec.lst > ${sample}_spec.sorted.lst
    mv ${sample}.sorted.lst ${sample}.lst
    mv ${sample}_spec.sorted.lst ${sample}_spec.lst
done
cd ..


# 10x 
# bam was renamed for the longranger outputted phased_possorted_bam.bam

cd 10X

for chr in  "${chrs[@]}"
do
    extractHAIRS --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --out HG00403_${chr}.lst --10x 1
    extractHAIRS_spec --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --out HG00403_${chr}_spec.lst --10x 1
    python3 LinkFragments.py --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --fragments HG00403_${chr}.lst --out HG00403_${chr}.linked.lst
    sort -n -k6 HG00403_${chr}_spec.lst > HG00403_${chr}_spec.sorted.lst
    BarcodeExtract HG00403_${chr}.bam HG00403_${chr}.bed ${chr}
    bgzip -c HG00403_${chr}.bed > HG00403_${chr}.bed.gz
    tabix -p bed HG00403_${chr}.bed.gz
    mv HG00403_${chr}.linked.lst HG00403_${chr}.lst
    mv HG00403_${chr}_spec.sorted.lst HG00403_${chr}_spec.lst
done

for sample in  "${samples[@]}"
do
    extractHAIRS --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}.lst --10x 1
    extractHAIRS_spec --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}_spec.lst --10x 1
    python3 LinkFragments.py --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --fragments ${sample}.lst --out ${sample}.linked.lst
    sort -n -k6 ${sample}_spec.lst > ${sample}_spec.sorted.lst
    BarcodeExtract ${sample}.bam ${sample}.bed 
    bgzip -c ${sample}.bed > ${sample}.bed.gz
    tabix -p bed ${sample}.bed.gz 
    mv ${sample}.linked.lst ${sample}.lst
    mv ${sample}_spec.sorted.lst ${sample}_spec.lst
done

cd ..


# HiC
cd HiC
for sample in  "${samples[@]}"
do
    extractHAIRS --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}.lst --hic 1
    extractHAIRS_spec --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}_spec.lst --hic 1
    sort -n -k6 ${sample}.lst > ${sample}.sorted.lst
    sort -n -k6 ${sample}_spec.lst > ${sample}_spec.sorted.lst
    mv ${sample}.sorted.lst ${sample}.lst
    mv ${sample}_spec.sorted.lst ${sample}_spec.lst
done
cd ..

# PACBIO

cd PACBIO
for chr in  "${chrs[@]}"
do
    extractHAIRS --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --out HG00403_${chr}.lst --pacbio 1 --ref hg19.fa
    extractHAIRS_spec --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --out HG00403_${chr}_spec.lst --pacbio 1 --ref hg19.fa 
    sort -n -k3 HG00403_${chr}.lst > HG00403_${chr}.sorted.lst
    sort -n -k3 HG00403_${chr}_spec.lst > HG00403_${chr}_spec.sorted.lst
    mv HG00403_${chr}.sorted.lst HG00403_${chr}.lst
    mv HG00403_${chr}_spec.sorted.lst HG00403_${chr}_spec.lst
done


for sample in  "${samples[@]}"
do
    extractHAIRS --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}.lst --pacbio 1 --ref hg19.fa
    extractHAIRS_spec --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}_spec.lst --pacbio 1 --ref hg19.fa
    sort -n -k3 ${sample}.lst > ${sample}.sorted.lst
    sort -n -k3 ${sample}_spec.lst > ${sample}_spec.sorted.lst
    mv ${sample}.sorted.lst ${sample}.lst
    mv ${sample}_spec.sorted.lst ${sample}_spec.lst
done
cd ..

# ONT

cd NANOPORE
for chr in  "${chrs[@]}"
do
    extractHAIRS --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --out HG00403_${chr}.lst --ont 1 --ref hg19.fa
    extractHAIRS_spec --bam HG00403_${chr}.bam --VCF HG00403_${chr}_hetsnp_sorted.vcf --out HG00403_${chr}_spec.lst --ont 1 --ref hg19.fa 
    sort -n -k3 HG00403_${chr}.lst > HG00403_${chr}.sorted.lst
    sort -n -k3 HG00403_${chr}_spec.lst > HG00403_${chr}_spec.sorted.lst
    mv HG00403_${chr}.sorted.lst HG00403_${chr}.lst
    mv HG00403_${chr}_spec.sorted.lst HG00403_${chr}_spec.lst
done


for sample in  "${samples[@]}"
do
    extractHAIRS --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}.lst --ont 1 --ref hg19.fa
    extractHAIRS_spec --bam ${sample}.bam --VCF ${sample}_hetsnp_sorted.vcf --out ${sample}_spec.lst --ont 1 --ref hg19.fa
    sort -n -k3 ${sample}.lst > ${sample}.sorted.lst
    sort -n -k3 ${sample}_spec.lst > ${sample}_spec.sorted.lst
    mv ${sample}.sorted.lst ${sample}.lst
    mv ${sample}_spec.sorted.lst ${sample}_spec.lst
done
cd ..