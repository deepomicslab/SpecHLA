#!/bin/bash

# remove temporary files

outdir=$1

echo "Clean output dir."

rm -f  $outdir/HLA_*spechap.vcf.gz*
rm -f  $outdir/HLA_*.fasta
rm -f  $outdir/assembly*
rm -f  $outdir/*fa
rm -f  $outdir/sample*
rm -f  $outdir/*hap*.fasta.out
rm -f  $outdir/rematch.*
rm -f  $outdir/*uniq.name.R*.gz
rm -f  $outdir/header.sam
rm -f $outdir/newref_insertion*
rm -f $outdir/*map_database.bam*
rm -f $outdir/*merge.bam*
rm -f $outdir/*realign.bam*
rm -f $outdir/*realign.vcf.gz*
rm -f $outdir/fragment*file
rm -f  $outdir/*.R*.fq.gz
rm -f  $outdir/A.bam*
rm -f  $outdir/B.bam*
rm -f  $outdir/C.bam*
rm -f  $outdir/DPA1.bam*
rm -f  $outdir/DPB1.bam*
rm -f  $outdir/DQA1.bam*
rm -f  $outdir/DQB1.bam*
rm -f  $outdir/DRB1.bam*
rm -f  $outdir/DRB1.hla.count
rm -f  $outdir/extract.read.blast
rm -rf  $outdir/tmp/







