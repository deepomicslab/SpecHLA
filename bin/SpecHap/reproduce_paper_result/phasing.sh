#!/bin/bash
declair -a chrs=("chr1" "chr21" "chr22")
declair -a samples=("NA12878" "NA19240")

mkdir SpecHap HapCUT2 FastHare DCHap RefHap

cd SpecHap
#SpecHap
for chr in  "${chrs[@]}"
do 
cd NGS
    SpecHap  --vcf  HG00403_${chr}_hetsnp_sorted.vcf.gz --frag HG00403_${chr}_spec.lst -o HG00403_${chr}.vcf
cd ..

cd 10X
    SpecHap  --vcf  HG00403_${chr}_hetsnp_sorted.vcf.gz --frag_stat HG00403_${chr}.bed.gz --frag HG00403_${chr}_spec.lst -o HG00403_${chr}.vcf --tenx
cd ..

cd HiC
    SpecHap  --vcf  HG00403_${chr}_hetsnp_sorted.vcf.gz --frag HG00403_${chr}_spec.lst -o HG00403_${chr}.vcf --hic
cd ..

cd PACBIO
    SpecHap  --vcf  HG00403_${chr}_hetsnp_sorted.vcf.gz --frag HG00403_${chr}_spec.lst -o HG00403_${chr}.vcf --pacbio
cd ..

cd NANOPORE
    SpecHap  --vcf  HG00403_${chr}_hetsnp_sorted.vcf.gz --frag HG00403_${chr}_spec.lst -o HG00403_${chr}.vcf --nanopore
cd ..
done 

for sample in  "${samples[@]}"
do
cd NGS
    SpecHap  --vcf  ${sample}_hetsnp_sorted.vcf.gz --frag ${sample}_spec.lst -o ${sample}.vcf
cd ..

cd 10X
    SpecHap  --vcf  ${sample}_hetsnp_sorted.vcf.gz --frag ${sample}_spec.lst -o ${sample}.vcf --frag_stat ${sample}.bed.gz --tenx
cd ..

cd HiC
    SpecHap  --vcf  ${sample}_hetsnp_sorted.vcf.gz --frag ${sample}_spec.lst -o ${sample}.vcf --hic
cd ..

cd PACBIO
    SpecHap  --vcf  ${sample}_hetsnp_sorted.vcf.gz --frag ${sample}_spec.lst -o ${sample}.vcf --pacbio
cd ..

cd NANOPORE
    SpecHap  --vcf  ${sample}_hetsnp_sorted.vcf.gz --frag ${sample}_spec.lst -o ${sample}.vcf --nanopore
cd ..
done 
cd ..


# HapCUT2 
cd HapCUT2
for chr in  "${chrs[@]}"
do 
cd NGS
    HapCUT2  --VCF  HG00403_${chr}_hetsnp_sorted.vcf.gz --f HG00403_${chr}.lst --o HG00403_${chr}.txt
cd ..

cd 10X
    HapCUT2  --VCF  HG00403_${chr}_hetsnp_sorted.vcf.gz --f HG00403_${chr}.lst --o HG00403_${chr}.txt --nf 1
cd ..

cd HiC
    HapCUT2  --VCF  HG00403_${chr}_hetsnp_sorted.vcf.gz --f HG00403_${chr}.lst --o HG00403_${chr}.txt --hic 1
cd ..

cd PACBIO
    HapCUT2  --VCF  HG00403_${chr}_hetsnp_sorted.vcf.gz --f HG00403_${chr}.lst --o HG00403_${chr}.txt
cd ..

cd NANOPORE
    HapCUT2  --VCF  HG00403_${chr}_hetsnp_sorted.vcf.gz --f HG00403_${chr}.lst --o HG00403_${chr}.txt
cd ..
done 

for sample in  "${samples[@]}"
do
cd NGS
    HapCUT2  --VCF  ${sample}_hetsnp_sorted.vcf --f ${sample}.lst --o ${sample}.txt
cd ..

cd 10X
    HapCUT2  --VCF  ${sample}_hetsnp_sorted.vcf --f ${sample}.lst --o ${sample}.txt --nf 1
cd ..

cd HiC
    HapCUT2  --VCF  ${sample}_hetsnp_sorted.vcf --f ${sample}.lst --o ${sample}.txt --hic 1
cd ..

cd PACBIO
    HapCUT2  --VCF  ${sample}_hetsnp_sorted.vcf --f ${sample}.lst --o ${sample}.txt
cd ..

cd NANOPORE
    HapCUT2  --VCF  ${sample}_hetsnp_sorted.vcf --f ${sample}.lst --o ${sample}.txt
cd ..
done 
cd ..

declare -a protocals=("NGS" "PACBIO" "NANOPORE")

# FastHare
cd FastHare
for protocal in "${protocals[@]}"
do
    cd protocal
    for chr in  "${chrs[@]}"
    do 
    cd NGS
        java -cp SingleIndividualHaplotyper/SIH.jar mpg.molgen.sih.main.SIH -a FastHare -v HG00403_${chr}.allvar -c 2 HG00403_${chr}.lst HG00403_${chr}_phased.txt
    cd ..
    done 
done 

for sample in  "${samples[@]}"
do
    for protocal in "${protocals[@]}"
    do
    cd protocal
    java -cp SingleIndividualHaplotyper/SIH.jar mpg.molgen.sih.main.SIH -a FastHare -v ${sample}.allvar -c 2 ${sample}.lst ${sample}_phased.txt
    cd ..
    done
done

cd ..


# RefHap
cd RefHap
for chr in  "${chrs[@]}"
do 
    for protocal in "${protocals[@]}"
    do
    cd protocal
    java -cp SingleIndividualHaplotyper/SIH.jar mpg.molgen.sih.main.SIH  -v HG00403_${chr}.allvar -c 2 HG00403_${chr}.lst HG00403_${chr}_phased.txt
    cd ..
    done 
done 

for sample in  "${samples[@]}"
do
    for protocal in "${protocals[@]}"
    do
    cd protocal 
    java -cp SingleIndividualHaplotyper/SIH.jar mpg.molgen.sih.main.SIH -v ${sample}.allvar -c 2 ${sample}.lst ${sample}_phased.txt
    cd ..
    done
done
cd ..


# DCHap
cd DCHap
for chr in  "${chrs[@]}"
do 
    for protocal in "${protocals[@]}"
    do
    cd protocal
    python2 Haplotype-phasing/dchap.py --reads HG00403_${chr}_dchap.lst --k 2  --parsed-reads HG00403_${chr}_parsed_fragment.txt --phase HG00403_${chr}_haplotype.txt --assignments HG00403_${chr}_assignments.txt
    python2 Haplotype-phasing/dchap-postprocess.py --parsed-reads HG00403_${chr}_parsed_fragment.txt --blocks HG00403_${chr}_haplotype.txt --assignments HG00403_${chr}_assignments.txt --corrected-blocks HG00403_${chr}_post_haplotype.txt
    
    cd ..
    done 
done 

for sample in  "${samples[@]}"
do
    for protocal in "${protocals[@]}"
    do
    cd protocal  
    python2 Haplotype-phasing/dchap.py --reads ${sample}_dchap.lst --k 2  --parsed-reads ${sample}_parsed_fragment.txt --phase ${sample}_haplotype.txt --assignments ${sample}_assignments.txt
    python2 Haplotype-phasing/dchap-postprocess.py --parsed-reads ${sample}_parsed_fragment.txt --blocks ${sample}_haplotype.txt --assignments ${sample}_assignments.txt --corrected-blocks ${sample}_post_haplotype.txt
    cd ..
    done
done
cd ..