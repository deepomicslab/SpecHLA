#!/bin/bash


###
### The Exon version of SpecHLA, performs HLA assembly and HLA Typing. 
###
### Usage:
###   sh SpecHLA_exon.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -p <Asian>
###
### Options:
###   -n        Sample ID.
###   -1        The first fastq file.
###   -2        The second fastq file.
###   -o        The output folder to store the typing results.
###   -p        The population of the sample: Asian, Black, or Caucasian. Use mean frequency 
###             if not provided.
###   -j        Number of threads [5]
###   -f        True or False. The annotation database only includes the alleles with population 
###             frequency higher than zero if set True. Otherwise, it includes all alleles. Default 
###             is True.
###   -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault 
###             is 3. It should be set larger to infer novel alleles.
###   -g        True or False. Split more blocks if set True; then rely more on allele database. 
###             Otherwise, split less blocks. Default is True.
###   -h        Show this message.

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

while getopts ":n:1:2:p:m:g:f:o:j:" opt; do
  case $opt in
    n) sample="$OPTARG"
    ;;
    1) fq1="$OPTARG"
    ;;
    2) fq2="$OPTARG"
    ;;
    p) pop="$OPTARG"
    ;;
    m) nm="$OPTARG"
    ;;
    g) pstrain_bk="$OPTARG"
    ;;
    f) annotation="$OPTARG"
    ;;
    o) given_outdir="$OPTARG"
    ;;
    j) num_threads="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
dir=$(cd `dirname $0`; pwd)
export LD_LIBRARY_PATH=$dir/../../spechla_env/lib
bin=$dir/../../bin
db=$dir/../../db
hlaref=$db/ref/hla.ref.extend.fa

if [ ${given_outdir:-NA} == NA ]
  then
    outdir=$(pwd)/output/$sample
  else
    outdir=$given_outdir/$sample   
fi

echo Start profiling HLA for $sample.
mkdir -p $outdir
# exec >$outdir/$sample.log 2>&1
group='@RG\tID:'$sample'\tSM:'$sample

# :<<!

# ################ remove the repeat read name #################
python3 $dir/../uniq_read_name.py $fq1 $outdir/$sample.uniq.name.R1.gz
python3 $dir/../uniq_read_name.py $fq2 $outdir/$sample.uniq.name.R2.gz
fq1=$outdir/$sample.uniq.name.R1.gz
fq2=$outdir/$sample.uniq.name.R2.gz
# ###############################################################

# ################### assign the reads to original gene################
## map the HLA reads to the allele database ##
license=$dir/../../bin/novoalign.lic
if [ -f "$license" ];then
    $bin/novoalign -d $db/ref/hla_gen.format.filter.extend.DRB.no26789.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
    -o FullNW -r All 100000 --mCPU ${num_threads:-5} -c 10  -g 20 -x 3  | $bin/samtools view \
    -Sb - | $bin/samtools sort -  > $outdir/$sample.map_database.bam
else
    $bin/bowtie2/bowtie2 -p ${num_threads:-5} -k 10 -x $db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta -1 $fq1 -2 $fq2|\
    $bin/samtools view -bS -| $bin/samtools sort - >$outdir/$sample.map_database.bam
fi
$bin/samtools index $outdir/$sample.map_database.bam

## assign the reads to corresponding gene ##
python3 $dir/../assign_reads_to_genes.py -1 $fq1 -2 $fq2 -n $bin -o $outdir \
-b ${outdir}/${sample}.map_database.bam -nm ${nm:-2}
# #####################################################################




# ########### align the gene-specific reads to the corresponding gene reference########
## align the reads to corresponding gene reference ##
$bin/bwa mem -U 10000 -L 10000,10000 -R $group $hlaref $fq1 $fq2 | $bin/samtools view -H  >$outdir/header.sam
#hlas=(A B C)
hlas=(A B C DPA1 DPB1 DQA1 DQB1 DRB1)
for hla in ${hlas[@]}; do
        hla_ref=$db/HLA/HLA_$hla/HLA_$hla.fa
        $bin/bwa mem -t ${num_threads:-5} -U 10000 -L 10000,10000 -O 7,7 -E 2,2 -R $group $hla_ref\
         $outdir/$hla.R1.fq.gz $outdir/$hla.R2.fq.gz | $bin/samtools view -bS -F 0x800 -| $bin/samtools sort - >$outdir/$hla.bam
        $bin/samtools index $outdir/$hla.bam
done
#samtools merge -f -h $outdir/header.sam $outdir/$sample.merge.bam $outdir/A.bam $outdir/B.bam $outdir/C.bam 
$bin/samtools merge -f -h $outdir/header.sam $outdir/$sample.merge.bam $outdir/A.bam $outdir/B.bam $outdir/C.bam $outdir/DPA1.bam $outdir/DPB1.bam $outdir/DQA1.bam $outdir/DQB1.bam $outdir/DRB1.bam
$bin/samtools index $outdir/$sample.merge.bam
# #####################################################################################



# ################################### local assembly and realignment #################################
sh $dir/../run.assembly.realign.sh $sample $outdir/$sample.merge.bam $outdir 70 $dir/select.region.exon.txt ${num_threads:-5}
echo realignment is done.
!
$bin/freebayes -a -f $hlaref -p 3 $outdir/$sample.realign.sort.bam > $outdir/$sample.realign.vcf && rm -rf $outdir/$sample.realign.vcf.gz 
bgzip -f $outdir/$sample.realign.vcf
tabix -f $outdir/$sample.realign.vcf.gz
#cp $outdir/$sample.realign.vcf $outdir/$sample.realign.filter.vcf
less $outdir/$sample.realign.vcf.gz |grep "#" > $outdir/$sample.realign.filter.vcf
$bin/bcftools filter -R $dir/exon_extent.bed $outdir/$sample.realign.vcf.gz |grep -v "#"  >> $outdir/$sample.realign.filter.vcf  
# #####################################################################################################

# !
# ###################### phase, link blocks, calculate haplotype ratio, give typing results ##############
python3 $dir/../phase_exon.py -b $outdir/$sample.realign.sort.bam -v $outdir/$sample.realign.filter.vcf \
-o $outdir/ -g ${pstrain_bk:-True} --snp_dp 0 --block_len 80 --points_num 1 --freq_bias 0.1 --reads_num 2



if [ ${annotation:-True} == True ]
then
  annotation_parameter=with_pop
else
  annotation_parameter=no_pop
fi
# #####################################################################################################


# ############################ annotation ####################################
perl $dir/anno_HLA_pop.pl $sample $outdir 2 $pop $annotation_parameter
# #############################################################################


cat $outdir/hla.result.txt
# sh $dir/../clear_output.sh $outdir/
echo $sample is done.
