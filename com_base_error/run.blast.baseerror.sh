dir=$1

for sample in `cat ./sample.list`
	do
            less ./select.sample.hla.list |grep $sample >./$sample.hla.list

while read sample hla1 hla2 gene
do
type1=${gene}_0
type2=${gene}_1
outdir=./$sample
hla_ref=/mnt/d/HLAPro_backup/HLAPro/db/ref/extend/hla.ref.extend.fa
if [[ ! -d $outdir ]]; then
        mkdir -p $outdir
fi

phase1=$dir/$sample/hla.allele.1.$gene.fasta
phase2=$dir/$sample/hla.allele.2.$gene.fasta
phase=$phase2
types=($type1 $type2)
hlas=($hla1 $hla2)
for type in ${types[@]}; do
	for hla in ${hlas[@]}; do
tag=`echo $type|cut -d "_" -f 3`
if [ "$tag" -eq "0" ]; then
	phase=$phase1
else
	phase=$phase2
fi
echo "$type  $phase"
samtools faidx /mnt/d/HLAPro_backup/HLAPro/hla_gen.format.filter.fasta $hla >$outdir/$hla.fasta
makeblastdb -in $outdir/$hla.fasta -dbtype nucl -parse_seqids -out $outdir/$hla
blastn -query $phase -out $outdir/$type.$hla.blast.txt -db $outdir/$hla -outfmt 7 -penalty -1 -reward 1 -gapopen 4 -gapextend 1 -strand plus 
blastn -query $phase -out $outdir/$type.$hla.blast -db $outdir/$hla -outfmt 3 -penalty -1 -reward 1 -gapopen 4 -gapextend 1 -strand plus 
#len=($(less $outdir/$hla.fasta|grep -v ">"|wc -c))
#                len=$(($len - 1))
#                java -jar /home/wangmengyao/packages/jvarkit/dist/blastn2snp.jar < $outdir/$type.$hla.blast.xml |grep -v "#"|cut -f 1,5-9|awk '{OFS="\t"}{print $1,$3,".",$5,$6}' > $outdir/tmp.vcf
#                ebase=($(less $outdir/tmp.vcf |cut -f 4|grep -v "-"|xargs|awk 'gsub(" ","")'|wc -c))
#                if [ "$ebase" -ne 0 ]; then 
#                        ebase=$(($ebase - 1))
#                fi
#                base_error=`echo "scale=5;$ebase/$len"|bc`
#
#                echo "#base_error with $hla is $base_error ($ebase:$len)" > $outdir/$type.$hla.blast.vcf
#                echo -e "#CHR\tPOS\tID\tExpected\tDetected" >> $outdir/$type.$hla.blast.vcf
#                cat $outdir/tmp.vcf >> $outdir/$type.$hla.blast.vcf
 #               rm -rf $outdir/tmp.vcf
	done
done


done < ./$sample.hla.list
done
