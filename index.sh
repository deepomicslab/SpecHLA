#!/bin/bash

# construct config file for ScanIndel
HLAs=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
dir=$(cd `dirname $0`; pwd)
for hla in ${HLAs[@]}; do
    config_file=$dir/db/HLA/HLA_${hla}.config.txt
    echo bwa=$dir/db/HLA/HLA_${hla}/HLA_${hla}.fa >$config_file
    echo freebayes=$dir/db/HLA/HLA_${hla}/HLA_${hla}.fa >>$config_file
    echo blat=$dir/db/HLA/HLA_${hla}/ >>$config_file
done

# index the database
$dir/bin/bowtie2/bowtie2-build $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta
$dir/bin/bowtie2/bowtie2-build $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta

license=$dir/bin/novoalign.lic
if [ -f "$license" ];then
    $dir/bin/novoindex  -k 14 -s 1 $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.ndx \
    $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta

    $dir/bin/novoindex  -k 14 -s 1 $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx \
    $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta
fi