for i in {1..50}
do
#i=19
sample=child_${i}
/mnt/e/hla_tgs/nanopore/DeepSimulator/deep_simulator.sh -i fasta/$sample.fasta -l 3000 -B 2 -c 6 -K 10 -P 0 -o ont_error/$sample
done
