sample=$1
bam=$2
outdir=$3
rlen=$4
dir=$(cd `dirname $0`; pwd)
sdir=$dir/../../bin
db=$dir/../../db
hla_fa=$db/ref/hla.ref.extend.fa

rm -rf $outdir/rematch.total.read.format.txt
for pos in `cat $dir/select.region.exon.txt`
do

#pos=HLA_B:1500-1800
hla=`echo $pos|cut -d ":" -f 1`
ref=$db/HLA/$hla/$hla

samtools view -f 64 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##1 "$2,$3}' > $outdir/extract.fa
samtools view -f 128 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##2 "$2,$3}' >> $outdir/extract.fa

$sdir/fermikit/fermi.kit/fermi2.pl unitig -s3g -t4 -T 10 -2 -l $rlen -p $outdir/prefix2 $outdir/extract.fa > $outdir/prefix2.mak

make -f $outdir/prefix2.mak

if [ ! -f "$outdir/prefix2.mag.gz" ]
then
	echo "$pos"
else
        less $outdir/prefix2.mag.gz|awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}'  > $outdir/prefix2.mag.fa
        less $outdir/prefix2.mag.gz|awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' | while read L; do  echo $L|awk '{print $1"_r"}' >>$outdir/prefix2.mag.fa; read L; echo "$L" | rev | tr "ATGC" "TACG" >>$outdir/prefix2.mag.fa; done

        $sdir/blastn -query $outdir/prefix2.mag.fa -out $outdir/prefix2.mag.blast -db $ref -outfmt 6 -strand plus  -penalty -1 -reward 1 -gapopen 4 -gapextend 1 
        less $outdir/prefix2.mag.blast |cut -f 1|sort|uniq > $outdir/id.list

if [ ! -f "$outdir/id.list" ]
then
	echo "$pos"
else
        $sdir/samtools faidx $outdir/prefix2.mag.fa 
        rm -rf $outdir/assembly.fa
        for id in `cat $outdir/id.list`
        do
	     $sdir/samtools faidx $outdir/prefix2.mag.fa $id >> $outdir/assembly.fa
        done

        $sdir/blastn -query $outdir/assembly.fa -out $outdir/assembly.blast -db $ref -strand plus -penalty -1 -reward 1 -gapopen 4 -gapextend 1 -line_length 700
        perl $sdir/blast2sam.pl $outdir/assembly.blast > $outdir/assembly.blast.sam

        $sdir/bwa index $outdir/assembly.fa
        $sdir/bwa mem -L 10000,10000 -O 6,10 -E 7,7 -A 3 -B 0 $outdir/assembly.fa $outdir/extract.fa | $sdir/samtools view -Sb -F 4 - | $sdir/samtools sort - > $outdir/rematch.bam
        $sdir/samtools view -F 0x800 $outdir/rematch.bam|cut -f 1,3,4,6 > $outdir/rematch.bam.txt
        perl $dir/../rematchblast.pl $outdir/assembly.blast.sam $outdir/extract.fa $outdir/rematch.bam.txt $outdir/rematch.read.format.txt 30
        cat $outdir/rematch.read.format.txt >>$outdir/rematch.total.read.format.txt
       
fi
fi
rm -rf $outdir/prefix2* $outdir/id.list
done

python $dir/../realignblast.py -i $bam -o $outdir/$sample.realign.bam -r $outdir/rematch.total.read.format.txt
$sdir/samtools sort $outdir/$sample.realign.bam > $outdir/$sample.realign.sort.bam
java -jar $sdir/picard.jar FixMateInformation I=$outdir/$sample.realign.sort.bam O=$outdir/$sample.realign.sort.fixmate.bam
$sdir/samtools index $outdir/$sample.realign.sort.fixmate.bam

$sdir/freebayes -a -f $hla_fa -p 3 $outdir/$sample.realign.sort.fixmate.bam > $outdir/$sample.realign.vcf && \
rm -rf $outdir/$sample.realign.vcf.gz 
$sdir/bgzip $outdir/$sample.realign.vcf
$sdir/tabix $outdir/$sample.realign.vcf.gz
less $outdir/$sample.realign.vcf.gz |grep "#" > $outdir/$sample.realign.filter.vcf
$sdir/bcftools filter -R $dir/exon_extent.V2.bed $outdir/$sample.realign.vcf.gz |grep -v "#" |awk '$6>0.01' >> $outdir/$sample.realign.filter.vcf  

