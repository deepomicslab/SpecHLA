# SpecHLA: a comprehensive package for HLA haplotyping and LOH inference

## Install  
First, create the env with conda, and activate the env.
```
cd SpecHLA/
conda env create --prefix=./spechla_env -f environment.yml
conda activate ./spechla_env
```
Second, make the softwares in `bin/` executable.
```
chmod +x -R bin/*
```
Third, index the database.
```
bash index.sh
```
Note: `Novoalign` requires a license in `bin/`, if not detected, SpecHLA uses `bowtie2` automatically. 
License file of `Novoalign` should be put in `bin/` before index



## Test
Please go to the `example/` folder, and run SpecHLA with different data types. 
If we install successfully, we will find results in the `output/`.

## Basic Usage  
### Extract HLA-related reads
We should first extract HLA reads with enrichment-free data. Please map the reads to `hg19` or `hg38`, then use the `ExtractHLAread.sh` in `script/` folder to extract HLA-related reads. We use the script of Kourami with minor revision for this step. 
Extract HLA-related reads by
```
sh ExtractHLAread.sh -s <sample_id> -b <bamfile> -r <refGenome>
```

### Typing 
Full-resolution and exon HLA typing with paired-end reads. SpecHLA can use PacBio, Nanopore, Hi-C, and 10X data to improve the phasing performance if provided. 

Perform full-resolution HLA typing with `paired-end` reads by
```
sh SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/
```
Perform exon HLA typing with paired-end reads by
```
sh SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ -u 1
```
Perform full-resolution HLA typing with `paired-end` reads and `PacBio` reads by
```
sh SpecHLA.sh -n <sample> -t <sample.pacbio.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `Nanopore` reads by
```
sh SpecHLA.sh -n <sample> -e <sample.nanopore.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `Hi-C` reads by
```
sh SpecHLA.sh -n <sample> -c <sample.hic.fwd.fq.gz> -d <sample.hic.rev.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `10X` linked reads by (LongRanger should be installed in system env)
```
sh SpecHLA.sh -n <sample> -x <sample.10x.read.folder> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Consider `long Indels` and use `population information` for annotation by
```
sh SpecHLA.sh -n <sample> -v True -p <Asia> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing only with `long reads` by 
```
python3 long_read_typing.py -n <sample> -r <Pacbio or Nanopore reads> -o outdir/ 
```

### Pedigree phasing
After performing HLA typing, we can afford trio information to update the results of child by
```
sh SpecHLA.sh -n <child> -1 <child.fq.1.gz> -2 <child.fq.2.gz> -o outdir/ --trio child:parent1:parent2
```
The fresh results can be found in the trio/ folder inside the original result folder.

### Detect HLA LOH in cancer samples
If the purity of cancer sample and the ploidy of HLA region are known, we can infer HLA LOH event by
```
perl cal.hla.copy.pl -purity <Purity> -ploidy <Ploidy> -F <freq.file.list> -O <outdir>
```
Use `ls output/*/*_freq.txt >freq.list` to generate `freq.file.list`.

## Interpret output
In the denoted outdir, the results of each sample are saved in a folder named as the sample ID.  

In the directory of one specific sample, you will find the below files:
| Output | Description |
| --- | --- |
| hla.result.txt | HLA-typing results for all alleles |
| hla.result.details.txt | Detail results of mapping the allele to the database |
| hla.allele.\*.HLA_\*.fasta | Sequence of each allele |
| HLA_*_freq.txt | Haplotype frequencies of each gene |
| HLA_*.rephase.vcf.gz | Phased vcf file for each gene |

If you performed pedigree phasing, you will find below files in the output/ folder.
| Output | Description |
| --- | --- |
| trio/HLA_*.rephase.vcf.gz | Phased vcf file after pedigree phasing |
| trio/hla.allele.\*.HLA_\*.fasta | Sequence of each allele after pedigree phasing|

An example for "hla.result.txt" is as below:
```
Sample  HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2 HLA_C_1 HLA_C_2 HLA_DPA1_1      HLA_DPA1_2      HLA_DPB1_1      HLA_DPB1_2      HLA_DQA1_1      HLA_DQA1_2      HLA_DQB1_1      HLA_DQB1_2      HLA_DRB1_1      HLA_DRB1_2
NA06985 A*02:01:01G     A*03:01:01G     B*07:02:01G     B*57:01:01G     C*07:02:01G     C*06:02:01G     DPA1*01:03:01G  DPA1*01:03:01G  DPB1*04:01:01G  DPB1*04:01:01G  DQA1*01:02:01G  DQA1*01:02:01G  DQB1*06:02:01G  DQB1*06:02:01G  DRB1*15:01:01G  DRB1*15:01:01G
```
An example for "HLA_*_freq.txt" is as below:
```
# Haplotype     Frequency
hla.allele.1.HLA_A.fasta 0.532
hla.allele.2.HLA_A.fasta 0.468
# The number of heterozygote variant is 30
```


## Dependencies 
### Programming 
* python=3.8.12 or above  
* Python libraries: numpy=1.22.3, pulp=2.6.0, pysam=0.19.0, scipy=1.8.0, and biopython=1.79. 
* Perl 5 or above

### Third party packages
The softwares in the `/bin` folder can also be installed by yourself, please remind link the executable file to the `bin/` folder.
Please apply Novoalign, and put the novoalign.lic (license file) in the `bin/` folder.

* Third party packages:
    * conda 4.12.0 
    * SamTools-1.11
    * bamutil Version: 1.0.14
    * picard-tools-2.1.0 (java 11.0.1)
    * BWA-0.7.17-r1188
    * blast 2.3.0
    * BLAT v. 36x2
    * bedtools-v2.26.0
    * bcftools-Version: 1.9
    * Freebayes-v1.2.0-2-g29c4002
    * Novoalign V4.02.01
    * Tabix-Version: 1.9
    * ScanIndel v1.3
    * Fermikit-0.13
    * SpecHap v1.0.1 and ExtractHAIRs
    * Minimap 2.17-r941
    * pbmm2-1.4.0
    * pbsv Version 2.6.2
    * Bowtie 2 version 2.3.4.1
    * longshot-0.4.5


## Simulation
Please go to the `simulation/` folder, and run `simu.hla.sh` like:
```
Simulation HLA capture reads for HLA Typing
USAGE: <PATH-TO>/simu.hla.sh -n <sample_size> -r <read_length> -1 <depth of haplotype1> -2 <depth of haplotype2>

 -n        : sample size of simulation [required]

 -r        : read length [required]

 -1        : depth of haplotype1 [required]

 -2        : depth of haplotype2 [required]
```


## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).


## Testing functions

### Long-read only assembly-based typing 
To assemble the long reads, and perform HLA typing based on the contigs, please use `/assembly/MHC.TGS.mapping.sh`. 
```
USAGE: <PATH-TO>/MHC.TGS.mapping.sh -s <sample_id> -f <mhc.read> -v <mhc.vcf> -r <refGenomedir> -t <Sequencer>

-s     : sample name (ex: H002) [required]

-f     : extracted MHC reads (ex: mhc.read.fastq.gz) [required]

-v     : extracted MHC snv vcf (ex: mhc.deepvariant.vcf) [required]

-r     : reference (dir/hs37d5.fa) [required]

-t     : sequencer (ex: Pacbio/ONT/10X) [required]
```
To run MHC.TGS.mapping.sh, some dependencies and preprocessings are needed:
```
mkdir bin

conda install -c bioconda -y minimap2 samtools pbmm2 pbsv pysam whatshap bcftools bamtools
git clone https://github.com/cschin/Peregrine.git
git clone -b peregrine https://github.com/cschin/pypeFLOW.git
conda install -c conda-forge -y pypy3.6


mkdir ref
cd ref
wget http://ftp.ebi.ac.uk:/pub/databases/ipd/imgt/hla/fasta/hla_gen.fasta
less hla_gen.fasta | tr " " "_" > tmp
mv tmp hla_gen.fasta
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz


# Extract MHC reads snps

echo '{ "filters" : [ { "id" : "inHP1", "tag" : "HP:1" }, { "id" : "inHP2", "tag" : "HP:2" } ], "rule" : "!(inHP1 | inHP2)" }' > filter_script.json

samtools view -bh $bam 6:28477797-33448354 | bamtools filter -in stdin -tag HP:1 | samtools bam2fq - | bgzip -c > tmp.HP1.fastq.gz

samtools view -bh $bam 6:28477797-33448354 | bamtools filter -in stdin -tag HP:2 | samtools bam2fq - | bgzip -c > tmp.HP2.fastq.gz

samtools view -bh $bam 6:28477797-33448354 | bamtools filter -in stdin -script filter_script.json | samtools bam2fq - | bgzip -c > tmp.MHConly.UnPartitioned.fastq.gz

zcat tmp.*.fastq.gz |pigz > sample.MHConly.fastq.gz

bcftools view -g het -v snps  deepvariant/pacbio-15kb-hapsort-wgs.vcf.gz  chr6:28477797-33448354 > sample.mhc.het.vcf
```
