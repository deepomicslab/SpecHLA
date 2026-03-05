#!/bin/bash
set -euo pipefail

# Ensure the dir is set relative to the actual path, links should be resolved
# so that this script can be linked somewhere else.
script_path=$(dirname $(realpath $0))
dir=$(cd $script_path; pwd)

# Source centralized path resolution
source "$dir/script/spechla_env.sh"

# construct config file for ScanIndel
HLAs=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
for hla in ${HLAs[@]}; do
    config_file=$SPECHLA_DB/HLA/HLA_${hla}.config.txt
    echo bwa=$SPECHLA_DB/HLA/HLA_${hla}/HLA_${hla}.fa >$config_file
    echo freebayes=$SPECHLA_DB/HLA/HLA_${hla}/HLA_${hla}.fa >>$config_file
    echo blat=$SPECHLA_DB/HLA/HLA_${hla}/ >>$config_file
done

# Index the database for alignment
if command -v novoalign &>/dev/null && [ -f "$(dirname $(which novoalign))/novoalign.lic" ]; then
    echo "Found novoalign on PATH with license."
    if [ -f "$SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.ndx" ] && \
       [ -f "$SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx" ]; then
        echo "Found reference index for Novoalign."
    else
        echo "Index the reference for novoalign..."
        novoindex -k 14 -s 1 $SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.ndx \
            $SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.fasta
        novoindex -k 14 -s 1 $SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx \
            $SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta
    fi
else
    echo "Novoalign not found, index the reference for bowtie2."
    bowtie2-build $SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.fasta \
        $SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.fasta
    bowtie2-build $SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta \
        $SPECHLA_DB/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta
fi

# Build SpecHap
mkdir -p $dir/bin/SpecHap/build
cd $dir/bin/SpecHap/build
cmake .. -DCMAKE_PREFIX_PATH=${CONDA_PREFIX:-/usr}
make

# Build ExtractHAIRs
mkdir -p $dir/bin/extractHairs/build
cd $dir/bin/extractHairs/build
cmake .. -DCMAKE_PREFIX_PATH=${CONDA_PREFIX:-/usr}
make

cd $dir

echo " "
echo " "
echo " "
echo "The installation is finished! Please start use SpecHLA."
echo "-------------------------------------"
