## please repace the folder to your own software folder.
bin=/home/wangmengyao/TheHumanPangenomeMHC/bin
export PATH=$PATH:/home/wangmengyao/anaconda3/envs/py37/bin/:$bin

odir=`pwd`
me=`basename $0`
dir=`dirname $0`
function usage {
	echo "TGS phased read of MHC region"
	echo "USAGE: <PATH-TO>/$me -s <sample_id> -f <mhc.read> -v <mhc.vcf> -r <refGenomedir> -t <Sequencer>"
	echo
	echo "-s     : sample name (ex: H002) [required]"
	echo
	echo "-f     : extracted MHC reads (ex: mhc.read.fastq.gz) [required]"
	echo
	echo "-v     : extracted MHC snv vcf (ex: mhc.deepvariant.vcf) [required]"
	echo
	echo "-r     : reference dir (dir/hs37d5.fa) [required]"
	echo
	echo "-t     : sequencer (ex: Pacbio/ONT/10X) [required]"
	exit 1
}


help() {
	sed -rn 's/^### ?//;T;p' "$0"
}
if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
	usage
	exit 1
fi



while getopts ":s:f:v:r:t" opt; do
	case $opt in
		s) sample="$OPTARG"
			;;
		f) read="$OPTARG"
			;;
                v) vcf="$OPTARG"
			;;
		r) ref="$OPTARG"
			;;
		t) seq="$OPTARG"
			;;
		\?) echo "Invalid option - $OPTARG" >&4
			;;
	esac
done



mkdir $odir/cram

#sample=HG002
group='@RG\tID:1\tSM:'$sample
thread=8
echo $group

#ref=/home/wangmengyao/TheHumanPangenomeMHC/ref/hs37d5.fa
#read=/home/wangmengyao/TheHumanPangenomeMHC/data/pacbio/reads/HG002.MHConly.fastq.gz
#vcf=/home/wangmengyao/TheHumanPangenomeMHC/data/pacbio/reads/HG002.mhc.het.vcf
seq="pacbio"
db=`dirname $ref`


# Generating Phased Read Sets

if [ $seq == 'pacbio' ]; then

        $bin/minimap2 -t $thread -R $group -ax asm20 $ref $read | samtools view -bS > $odir/cram/$sample.MHC.bam

fi

if [ $seq == '10X' ]; then

	$bin/minimap2 -t $thread -R $group -ax asm20 $ref $read | samtools view -bS > $odir/cram/$sample.MHC.bam

fi

if [ $seq == "ONT" ]; then

        $bin/minimap2 -t $thread -R $group -ax map-ont $ref $read | samtools view -bS - > $odir/cram/$sample.MHC.bam

fi

$bin/samtools sort -@ $thread -o $odir/cram/$sample.MHC.sorted.bam $odir/cram/$sample.MHC.bam
$bin/samtools index $odir/cram/$sample.MHC.sorted.bam



mkdir $odir/whatshap
$bin/whatshap genotype --ignore-read-groups --chromosome 6 \
	--reference $ref \
	-o $odir/whatshap/$sample.regenotype.vcf \
	$vcf $odir/cram/$sample.MHC.sorted.bam

whatshap phase --ignore-read-groups --chromosome 6 \
	        --reference $ref \
		-o $odir/whatshap/$sample.phased.vcf \
	        $odir/whatshap/$sample.regenotype.vcf $odir/cram/$sample.MHC.sorted.bam 


$bin/bgzip -f $odir/whatshap/$sample.phased.vcf
$bin/tabix -f $odir/whatshap/$sample.phased.vcf.gz



$bin/whatshap haplotag --ignore-read-groups \
	--reference $ref \
	-o $odir/whatshap/$sample.sorted.tagged-newsplit.bam \
	$odir/whatshap/$sample.phased.vcf.gz \
	$odir/cram/$sample.MHC.sorted.bam



mkdir $odir/phased_reads

$bin/samtools view $odir/whatshap/$sample.sorted.tagged-newsplit.bam | grep HP:i:1 | awk '{print $1}' | sort -u  > $odir/phased_reads/H1_reads 
$bin/samtools view $odir/whatshap/$sample.sorted.tagged-newsplit.bam | grep HP:i:2 | awk '{print $1}' | sort -u  > $odir/phased_reads/H2_reads 
$bin/samtools view $odir/whatshap/$sample.sorted.tagged-newsplit.bam | grep -v HP:i:2 | grep -v HP:i:1 | awk '{print $1}' | sort -u > $odir/phased_reads/tmp_reads

cat $odir/phased_reads/H1_reads $odir/phased_reads/H2_reads | sort | uniq > $odir/phased_reads/phased_reads
sort $odir/phased_reads/tmp_reads $odir/phased_reads/phased_reads $odir/phased_reads/phased_reads |uniq -u > $odir/phased_reads/unphased_reads


mkdir $odir/tmp
bam=$odir/whatshap/$sample.sorted.tagged-newsplit.bam
java -jar $bin/picard.jar FilterSamReads I=$bam O=$odir/tmp/$sample.filter.bam READ_LIST_FILE=$odir/phased_reads/H1_reads FILTER=includeReadList VALIDATION_STRINGENCY=LENIENT
$bin/samtools sort -n -o $odir/tmp/$sample.filter.sort.bam $odir/tmp/$sample.filter.bam
$bin/samtools fasta -0 $odir/phased_reads/H1.fa $odir/tmp/$sample.filter.sort.bam

java -jar $bin/picard.jar FilterSamReads I=$bam O=$odir/tmp/$sample.filter.bam READ_LIST_FILE=$odir/phased_reads/H2_reads FILTER=includeReadList VALIDATION_STRINGENCY=LENIENT
$bin/samtools sort -n -o $odir/tmp/$sample.filter.sort.bam $odir/tmp/$sample.filter.bam
$bin/samtools fasta -0 $odir/phased_reads/H2.fa $odir/tmp/$sample.filter.sort.bam

java -jar $bin/picard.jar FilterSamReads I=$bam O=$odir/tmp/$sample.filter.bam READ_LIST_FILE=$odir/phased_reads/unphased_reads FILTER=includeReadList VALIDATION_STRINGENCY=LENIENT
$bin/samtools sort -n -o $odir/tmp/$sample.filter.sort.bam $odir/tmp/$sample.filter.bam
$bin/samtools fasta -0 $odir/phased_reads/unphased.fa $odir/tmp/$sample.filter.sort.bam


#02_run_assembler

echo $odir/phased_reads/H1.fa > $odir/phased_reads/reads.lst
echo $odir/phased_reads/unphased.fa >> $odir/phased_reads/reads.lst
$bin/pg_run.py asm $odir/phased_reads/reads.lst 4 4 4 4 4 4 4 4 4 --shimmer-r 3 \
                --with-consensus \
                --best_n_ovlp 8 \
                --output $odir/phased_reads/asm-MHC-H1 

cat $odir/phased_reads/asm-MHC-H1/p_ctg_cns.fa > $odir/p_ctg_cns_H1.fa

echo $odir/phased_reads/H2.fa > $odir/phased_reads/reads.lst
echo $odir/phased_reads/unphased.fa >> $odir/phased_reads/reads.lst


$bin/pg_run.py asm $odir/phased_reads/reads.lst 4 4 4 4 4 4 4 4 4 --shimmer-r 3 \
                --with-consensus \
                --best_n_ovlp 8 \
                --output $odir/phased_reads/asm-MHC-H2

cat $odir/phased_reads/asm-MHC-H2/p_ctg_cns.fa > $odir/p_ctg_cns_H2.fa


#HLATyping from assmbly

$bin/minimap2 -x asm5 -c --eqx  $odir/p_ctg_cns_H1.fa  $db/hla_gen.fasta > $odir/hla_gen_H1.paf 
$bin/minimap2 -x asm5 -c --eqx  $odir/p_ctg_cns_H2.fa  $db/hla_gen.fasta > $odir/hla_gen_H2.paf 

cp $dir/HLAtyping.assmbly.py $odir

python $odir/HLAtyping.assmbly.py


