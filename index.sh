#!/bin/bash
set -euo pipefail

# construct config file for ScanIndel
HLAs=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)

# Ensure the dir is set relative to the actual path, links should be resolved
# so that this script can be linked somewhere else.
script_path=$(dirname $(realpath $0))
dir=$(cd $script_path; pwd)

# :<<!
for hla in ${HLAs[@]}; do
    config_file=$dir/db/HLA/HLA_${hla}.config.txt
    echo bwa=$dir/db/HLA/HLA_${hla}/HLA_${hla}.fa >$config_file
    echo freebayes=$dir/db/HLA/HLA_${hla}/HLA_${hla}.fa >>$config_file
    echo blat=$dir/db/HLA/HLA_${hla}/ >>$config_file
done




license=$dir/bin/novoalign.lic
if [ -f "$license" ];then
    # index the database for novoalign
    echo "Find License for Novoalign."
    if [ -f "$dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.ndx" ] && [ -f "$dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx" ];then
        echo "Find reference index for Novoalign."
    else
        echo "Index the reference for novoalign..."
        $dir/bin/novoindex  -k 14 -s 1 $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.ndx \
        $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta

        $dir/bin/novoindex  -k 14 -s 1 $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.ndx \
        $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta
    fi
else
    # index the database for bowtie2
    echo "Cannot Find License for Novoalign, index the reference for bowtie2."
    bowtie2-build $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta
    bowtie2-build $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta $dir/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta

fi

# the lib required by samtools
ln --force -s $CONDA_PREFIX/lib/libncurses.so.6 $CONDA_PREFIX/lib/libncurses.so.5
ln --force -s $CONDA_PREFIX/lib/libtinfo.so.6 $CONDA_PREFIX/lib/libtinfo.so.5
#ln -s $CONDA_PREFIX/lib/libhts.so.3 $CONDA_PREFIX/lib/libhts.so.2
# !

# install spechap
mkdir -p $dir/bin/SpecHap/build
cd $dir/bin/SpecHap/build
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
make

# install spechap
mkdir -p $dir/bin/extractHairs/build
cd $dir/bin/extractHairs/build
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
make


cd $dir

echo " "
echo " "
echo " "
echo The installation is finished! Please start use SpecHLA.
echo "-------------------------------------"
