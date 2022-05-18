#!/bin/bash

# construct config file for ScanIndel
HLAs=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
dir=$(cd `dirname $0`; pwd)
for hla in ${HLAs[@]}; do
    config_file=$dir/db/HLA/HLA_${hla}.config.txt
    echo bwa=$dir/db/HLA/HLA_${hla}/HLA_${hla}.fa >>$config_file
    echo freebayes=$dir/db/HLA/HLA_${hla}/HLA_${hla}.fa >>$config_file
    echo blat=$dir/db/HLA/HLA_${hla}/ >>$config_file
done

# index the database
./bin/novoindex  -k 14 -s 1 ./db/ref/hla_gen.format.filter.extend.DRB.no26789.ndx \
./db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta

./bin/novoindex  -k 14 -s 1 ./db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx \
./db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta
