for i in {1..50}
do
sample=child_$i
mkdir pacbio/$sample
pbsim --data-type CLR --seed 88 --accuracy-mean 0.85 --accuracy-min 0.80 --prefix pacbio/$sample/$sample --depth 10 --model_qc pacbio/model_qc_clr fasta/$sample.fasta
cat pacbio/$sample/${sample}_*fastq>pacbio/$sample/${sample}.fastq
rm pacbio/$sample/${sample}_*fastq
done
