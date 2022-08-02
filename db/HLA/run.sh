HLAs=(HLA_A HLA_B HLA_C HLA_DPA1 HLA_DPB1 HLA_DQA1 HLA_DQB1 HLA_DRB1)
sdir=/home/wangmengyao/SpecHLA/bin
dir=$( cd `dirname $0`; pwd)
for hla in ${HLAs[@]}; do
     $sdir/samtools faidx hla.ref.extend.fa $hla >$hla/$hla.fa 
     $sdir/samtools faidx $hla/$hla.fa
     $sdir/bwa index $hla/$hla.fa
     $sdir/faToTwoBit $hla/$hla.fa $hla/$hla.2bit
     $sdir/makeblastdb -in $hla/$hla.fa -dbtype nucl -parse_seqids -out $hla/$hla
#     echo "bwa=$dir/$hla/$hla.fa">./$hla.config.txt
#     echo "freebayes=$dir/$hla/$hla.fa">>./$hla.config.txt
#     echo "blat=$dir/$hla/">>./$hla.config.txt
done  

