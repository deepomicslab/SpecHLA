
sample=$1
bam=$2
port=$4
dir=$3

pdir=$(cd `dirname $0`; pwd)
bin=$pdir/../../bin
db=$pdir/../../db

outdir=$dir/Scanindel
mkdir -p $outdir

export PATH=$PATH:$bin

rm -rf $outdir/$sample.breakpoint.txt $outdir/$sample.ins.fa
HLAs=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
perl $pdir/merge_breakpoint.pl $dir/$sample.realign.filter.vcf $outdir/$sample.breakpoint.txt $outdir/$sample.ins.fa          
for hla in ${HLAs[@]}; do
        hfqfile=$outdir/scanindel.fq.HLA_$hla.list
        #samtools view -h $bam HLA_$hla | samtools sort -n - |samtools fastq - -1 $outdir/$sample.$hla.1.fastq.gz -2 $outdir/$sample.$hla.2.fastq.gz
        echo "$sample.$hla $dir/$hla.R1.fq.gz $dir/$hla.R2.fq.gz" >$hfqfile
        $bin/gfServer stop localhost $port
        python2 $pdir/ScanIndel.py -F 0 -p $db/HLA/HLA_$hla.config.txt -i $outdir/scanindel.fq.HLA_$hla.list\
         -o $outdir/ --gfServer_port $port
        perl $pdir/get_breakpoint2.pl $outdir/$sample.$hla.contigs.bam $outdir/$sample.$hla.merged.indel.vcf\
         $outdir/$sample.$hla.breakpoint.txt $outdir/$sample.$hla.ins.fa
        cat $outdir/$sample.$hla.breakpoint.txt >> $outdir/$sample.breakpoint.txt
        cat $outdir/$sample.$hla.ins.fa >>$outdir/$sample.ins.fa
done

