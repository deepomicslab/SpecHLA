#!/bin/bash

# subsample from phase 3 vcf file
# para: 
# 1: phase 3 vcf.gz indexed
# 2: sample name
# 3: chromosome

function get_sample_vcf()
{
    if [[$# -eq 3]]
    then 
        bcftools view -g het -v snps -s $2 $1 $3 > ${2}_${3}.vcf
    elif [[$# -eq 2]]
    then
        bcftools view -g het -v snps -s $2 $1 > ${2}.vcf
    fi
}

# get diploid fasta with HG00403
# para:
# 1: fasta
# 2: vcf 
# 3: output prefix 
function get_diploid_fasta()
{
    bcftools consensus -H 1 -f $1 -o ${3}_p.fa $2
    bcftools consensus -H 2 -f $1 -o ${3}_m.fa $2
}

# simulation with HG00403
# wgsim with haplotype mode
# 30x, chr1, chr21, chr22 
# param:
# 1: coverage (30) 
# 2: read_len (150)
# 3: insertion_size (350)
# 4: error_rate 
# 5: prefix
# 6: paternal fasta
# 7: maternal fasta
function NGS_simulation()
{
    fasta_length=$(cat $6 | grep -v '^>' | wc -lc | awk '{print $2 - $1}')
    let "pe_reads = $fasta_length*$coverage/300"
    pe_reads=$(echo "$pe_reads 1000000" | awk '{printf "%.6f", $1/$2}')
    let "pe_length = $fasta_length*$1/2"
    wgsim -e $4 -d $3 -s 35 -N $pe_reads -1 $2 -2 $2 -r0 -R0 -X0 $6 ${5}_paternal_r1.fq ${5}_paternal_r2.fq
    wgsim -e $4 -d $3 -s 35 -N $pe_reads -1 $2 -2 $2 -r0 -R0 -X0 $7 ${5}_maternal_r1.fq ${5}_maternal_r2.fq
    cat ${5}_paternal_r1.fq ${5}_maternal_r1.fq | gzip -c > ${5}_r1.fq.gz
    cat ${5}_paternal_r2.fq ${5}_maternal_r2.fq | gzip -c > ${5}_r2.fq.gz
    rm ${5}_paternal_r1.fq ${5}_maternal_r1.fq ${5}_paternal_r2.fq ${5}_maternal_r2.fq
}

# simulation with HG00403
# 30x, 50k SLR length, chr1, chr21, chr22
# lrsim with haplotype mode
# param:
# 1: coverage (30)
# 2: read_len (150)
# 3: insertion_size (350)
# 4: error_rate 
# 5: prefix
# 6: paternal fasta
# 7: maternal fasta
# 8: molecule_len (50)
function tenx_simulation()
{
    fasta_length=$(cat $6 | grep -v '^>' | wc -lc | awk '{print $2 - $1}')
    let "pe_reads = $fasta_length*$coverage/300"
    pe_reads=$(echo "$pe_reads 1000000" | awk '{printf "%.6f", $1/$2}')
    let "pe_length = $fasta_length*$1/2"
    partition_lower=$(echo "$pe_length $8" | awk '{printf "%.6f", $1/$2/1000}')
    let "partition_lower = $(printf %.0f $partition_lower)+1"
    partition_lower=$(echo "$partition_lower 1000" | awk '{printf "%.3f", $1/$2}')

    ./simulateLinkedReads.pl -r $6 -p ${5}_paternal -d 1 -1 $fasta_length -n -i $3 -x $pe_reads -f $8 -t $partition_lower -m 1 -o -4 1 -7 1 -e $4
    ./simulateLinkedReads.pl -r $7 -p ${5}_maternal -d 1 -1 $fasta_length -n -i $3 -x $pe_reads -f $8 -t $partition_lower -m 1 -o -4 1 -7 1 -e $4

    zcat ${5}_paternal/${5}_paternal_S1_L001_R1_001.fastq.gz | awk 'NR%4==1{gsub("/1","",$0)}; {print $0}' > $${5}_paternal/${5}_paternal_L001_R1_001.fq 
    zcat ${5}_paternal/${5}_paternal_S1_L001_R2_001.fastq.gz | awk 'NR%4==1{gsub("/2","",$0)}; {print $0}' > $${5}_paternal/${5}_paternal_L001_R2_001.fq

    zcat ${5}_maternal/${5}_maternal_S1_L001_R1_001.fastq.gz | awk 'NR%4==1{gsub("/1","",$0)}; {print $0}' > $${5}_maternal/${5}_maternal_L001_R1_001.fq 
    zcat ${5}_maternal/${5}_maternal_S1_L001_R2_001.fastq.gz | awk 'NR%4==1{gsub("/2","",$0)}; {print $0}' > $${5}_maternal/${5}_maternal_L001_R2_001.fq

    cat ${5}_paternal/${5}_paternal_L001_R1_001.fq ${5}_maternal_L001_R1_001.fq | gzip -c > ${5}_L001_R1_001.fq.gz
    cat ${5}_paternal/${5}_paternal_L001_R2_001.fq ${5}_maternal_L001_R2_001.fq | gzip -c > ${5}_L001_R2_001.fq.gz
}


# simulation with HG00403
# 50x, chr1, chr21, chr22
# pbsim with pacbio clr
# param:
# 1: maternal_fasta
# 2: paternal_fasta
# 3: prefix
# 4: real pacbio clr reads for homo sapien
function pacbio_simulation()
{   
    pbsim --prefix ${3}_paternal --data-type CLR --depth 25 --sample-fastq $4 $2
    pbsim --prefix ${3}_maternal --data-type CLR --depth 25 --sample-fastq $4 $1
    # cat paternal and maternal by user
}

#simulation with HG00403 
# 50x, chr1, chr21, chr22
# oxford nanopore with pretrained model
# param:
# 1: maternal_fasta
# 2: paternal_fasta
# 3: prefix
function nanopore_simulation()
{
    ./deep_simulator.sh -i $2 -c 20 -B 2 -K 25 -o ${3}_paternal
    ./deep_simulator.sh -i $1 -c 20 -B 2 -K 25 -o ${3}_maternal
    # cat paternal and maternal by user
}