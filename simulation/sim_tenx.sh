for i in {1..50}
do
sample=child_$i
cat /mnt/disk2_workspace/wangshuai/00.strain/08.NeedleHLA/other_data/batch_run/sim/fasta/$sample.fasta /home/wangshuai/softwares/LRSIM/test/ecoli.fa >$sample.ecoli.fasta 
perl /home/wangshuai/softwares/LRSIM/simulateLinkedReads.pl -r $sample.ecoli.fasta -p $sample -1 1000000 -4 1000000 -n -e 0 -E 0 -f 20 -x 1 -o -c /home/wangshuai/softwares/LRSIM/test/fragmentSizesList
done
