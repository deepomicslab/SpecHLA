while read file gene
do
  cp $file $gene
  /home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/makeblastdb -in $gene -dbtype nucl -parse_seqids -out $gene
done < exon.list
