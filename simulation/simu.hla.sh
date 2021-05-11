#!/bin/bash

me=`basename $0`
function usage {
	    echo "Simulation HLA capture reads for HLA Typing"
	    echo "USAGE: <PATH-TO>/$me -n <sample_size> -r <read_length> -1 <depth of haplotype1> -2 <depth of haplotype2>"
            echo
	    echo " -n        : sample size of simulation [required]"
	    echo
            echo " -r        : read length [required]"
	    echo
	    echo " -1        : depth of haplotype1 [required]"
	    echo
	    echo " -2        : depth of haplotype2 [required]"
	    echo
	    exit 1
}

help() {
	    sed -rn 's/^### ?//;T;p' "$0"
    }

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
	    usage
	    exit 1
fi

while getopts ":n:r:1:2:" opt; do
	  case $opt in
            n) nsample="$OPTARG"
	    ;;
	    r) readlen="$OPTARG"
	    ;;
	    1) depth1="$OPTARG"
	    ;;
            2) depth2="$OPTARG"
	    ;;
            \?) echo "Invalid option -$OPTARG" >&2
	    ;;
esac
done

#nsample=$1 #10
#readlen=$2 #150
#depth1=$3 #50
#depth2=$4 #50

function range_fq(){
    hla=$1
    echo $1
    shift
    new_hla_list=($@)
    
    flag=0
    dir=simu/fastq/${hla}_T_${depth1}-${depth2}
    mkdir -p $dir
    for H in ${new_hla_list[@]};do
    if [ $flag -eq 0 ] ; then
         less simu/fasta/$H.$depth1.bwa.read1.fastq.gz|gzip >> $dir/${hla}_T_${depth1}-${depth2}.read1.fastq.gz
         less simu/fasta/$H.$depth1.bwa.read2.fastq.gz|gzip >> $dir/${hla}_T_${depth1}-${depth2}.read2.fastq.gz
	 flag=1
    else
         less simu/fasta/$H.$depth2.bwa.read1.fastq.gz|gzip >> $dir/${hla}_T_${depth1}-${depth2}.read1.fastq.gz
	 less simu/fasta/$H.$depth2.bwa.read2.fastq.gz|gzip >> $dir/${hla}_T_${depth1}-${depth2}.read2.fastq.gz
         flag=0
    fi
done
}

function simu_pe_reads(){
    hla=$1
   
    fasta=simu/fasta/$hla.fasta
    echo $fasta
    ../bin/dwgsim -e 0 -E 0 -1 $readlen -2 $readlen -d 500 -C $depth1 -r 0 -o 1 $fasta simu/fasta/$hla.$depth1
    rm -rf simu/fasta/$hla.$coverage.mutations.*
    ../bin/dwgsim -e 0 -E 0 -1 $readlen -2 $readlen -d 500 -C $depth2 -r 0 -o 1 $fasta simu/fasta/$hla.$depth2
    rm -rf simu/fasta/$hla.$coverage2.mutations.*

}

for i in ${nsample};
do
    perl ./select.pl >select.hla.list
    hla_list=$(less select.hla.list| awk 'gsub("\t",";")')
    HLA_list=$(less select.hla.list|cut -f 2)
    mkdir -p simu/fasta
    HLA="HLA_"$i

    for line in ${hla_list[@]}; do
                class=`echo $line|cut -d ";" -f 1`
		hla=`echo $line|cut -d ";" -f 2`
		hla_fa="simu/fasta/hla.fasta"
		`samtools faidx ./hla_gen.format.filter.fasta $hla > $hla_fa`
		`head -1 $hla_fa > simu/fasta/$hla.fasta`
		`samtools faidx ./extend.fa HLA_$class\_1 | grep -v ">"  > simu/fasta/tmp.fa`
                `less $hla_fa|grep -v ">" >> simu/fasta/tmp.fa`
                `samtools faidx ./extend.fa HLA_$class\_2 | grep -v ">" >> simu/fasta/tmp.fa`
                `sed -i ":a;N;s/\n//g;ta" simu/fasta/tmp.fa`
                `cat simu/fasta/tmp.fa >> simu/fasta/$hla.fasta`
		rm -rf simu/fasta/tmp.fa $hla_fa
		`samtools faidx simu/fasta/$hla.fasta`
    done
     
        
        echo $HLA >>merge.hla.list
        `cat select.hla.list >>merge.hla.list`
        for hla in ${HLA_list[@]} ; do
             echo $hla 
             simu_pe_reads ${hla}
        done
        range_fq $HLA ${HLA_list[@]}

done





