for i in {1..50}
do
sample=child_${i}
sim3C --dist uniform -n 10000 -l 150 -e NlaIII -m hic fasta/$sample.fasta hic/$sample.fastq --simple-reads --insert-mean 500 --anti-rate 0 --trans-rate 0 --spurious-rate 0 -r 88 
rm hic/profile.tsv
python3 /mnt/d/plot_HLA/hic/split_hic.py hic/$sample.fastq $sample hic
done
