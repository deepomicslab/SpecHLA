sample=$1
bam=$2
outdir=$3
rlen=$4

bedfile=$5
thread=$6
dir=$(cd `dirname $0`; pwd)
sdir=$dir/../bin
db=$dir/../db
hla_fa=$db/ref/hla.ref.extend.fa

exec >$outdir/$sample.local_assem.log 2>&1 #redirect log info to the outdir

rm -rf $outdir/rematch.total.read.format.txt
echo "readid	tag	chr	start	cigar" > $outdir/rematch.total.read.format.txt

for pos in `cat $bedfile`
do

hla=`echo $pos|cut -d ":" -f 1`
ref=$db/HLA/$hla/$hla

samtools view --threads $thread -f 64 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##1 "$2,$3}' > $outdir/extract.fa
samtools view --threads $thread -f 128 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##2 "$2,$3}' >> $outdir/extract.fa

if [ ! -f "$outdir/extract.fa" ]
then
	echo "$pos noreads"
else
$sdir/fermikit/fermi.kit/fermi2.pl unitig -s1k -t $thread  -T 10 -2 -l $rlen -p $outdir/prefix2 $outdir/extract.fa > $outdir/prefix2.mak

make -f $outdir/prefix2.mak

if [ ! -f "$outdir/prefix2.mag.gz" ]
then
	echo "$pos assemble fail"
else
        zless $outdir/prefix2.mag.gz|awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}'  > $outdir/prefix2.mag.fa
        zless $outdir/prefix2.mag.gz|awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' | while read L; do  echo $L|awk '{print $1"_r"}' >>$outdir/prefix2.mag.fa; read L; echo "$L" | rev | tr "ATGC" "TACG" >>$outdir/prefix2.mag.fa; done

        $sdir/blastn -num_threads $thread -query $outdir/prefix2.mag.fa -out $outdir/prefix2.mag.blast -db $ref -outfmt 6 -strand plus  -penalty -1 -reward 1 -gapopen 4 -gapextend 1 
        less $outdir/prefix2.mag.blast |cut -f 1|sort|uniq > $outdir/id.list

if [ ! -f "$outdir/id.list" ]
then
	echo "$pos bad contig"
else
        samtools faidx $outdir/prefix2.mag.fa 
        rm -rf $outdir/assembly.fa
        for id in `cat $outdir/id.list`
        do
	     samtools faidx $outdir/prefix2.mag.fa $id >> $outdir/assembly.fa
        done

        $sdir/blastn -num_threads $thread -query $outdir/assembly.fa -out $outdir/assembly.blast -db $ref -strand plus -penalty -1 -reward 1 -gapopen 4 -gapextend 1 -line_length 700
        perl $sdir/blast2sam.pl $outdir/assembly.blast > $outdir/assembly.blast.sam

        bwa index $outdir/assembly.fa
        bwa mem -t $thread -L 10000,10000 -O 6,10 -E 7,7 -A 3 -B 0 $outdir/assembly.fa $outdir/extract.fa | samtools view --threads $thread -Sb -F 4 - | samtools sort --threads $thread - > $outdir/rematch.bam
        samtools view --threads $thread -F 0x800 $outdir/rematch.bam|cut -f 1,3,4,6 > $outdir/rematch.bam.txt
        perl $dir/rematchblast.pl $outdir/assembly.blast.sam $outdir/extract.fa $outdir/rematch.bam.txt $outdir/rematch.read.format.txt 30
        cat $outdir/rematch.read.format.txt >>$outdir/rematch.total.read.format.txt
       
fi
fi
fi

rm -rf $outdir/prefix2* $outdir/id.list
done

$dir/../spechla_env/bin/python3 $dir/realignblast.py -i $bam -o $outdir/$sample.realign.bam -r $outdir/rematch.total.read.format.txt
samtools sort --threads $thread $outdir/$sample.realign.bam > $outdir/$sample.realign.sort.bam
#java -jar $sdir/picard.jar FixMateInformation I=$outdir/$sample.realign.sort.bam O=$outdir/$sample.realign.sort.fixmate.bam
samtools index $outdir/$sample.realign.sort.bam

