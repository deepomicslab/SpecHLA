# SpecHLA: a comprehensive package for HLA haplotyping and LOH inference

- [SpecHLA: a comprehensive package for HLA haplotyping and LOH inference](#spechla--a-comprehensive-package-for-hla-haplotyping-and-loh-inference)
  * [Install](#install)
  * [Run](#run)
  * [Test](#test)
  * [Dependencies](#dependencies)
    + [Programming](#programming)
    + [Third party packages](#third-party-packages)
  * [Basic Usage](#basic-usage)
    + [Extract HLA-related reads](#extract-hla-related-reads)
    + [Whole Mode](#whole-mode)
    + [Exon Mode](#exon-mode)
    + [Pedigree phasing](#pedigree-phasing)
    + [Detect HLA LOH in cancer samples](#detect-hla-loh-in-cancer-samples)
    + [TGS only alignment-based typing Testing](#tgs-only-alignment-based-typing-testing)
    + [TGS only assembly-based typing Testing](#tgs-only-assembly-based-typing-testing)
  * [Interpret output](#interpret-output)
  * [Simulation](#simulation)
  * [Getting help](#getting-help)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


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
Note: `Novoalign` requires a license in `bin/`, if it is not detected, SpecHLA will use `bowtie2`. 
## Run
Then we can perform HLA typing in exons.
```
cd SpecHLA/script/exon/
sh SpecHLA_exon.sh -h
```
Perform HLA typing in whole genes with paired-end data. Also, PacBio, Nanopore, Hi-C, or 10X could be integrated.
```
cd SpecHLA/script/whole/
sh SpecHLA_whole.sh -h
```
Adpot pedigree information to refine the typing results.
```
python3 phase_pedigree.py --trio child:mother:father --workdir output/
```
Detect LOH events
```
ls output/*/*_freq.txt >freq.list
perl SpecHLA/script/cal.hla.copy.pl  -purity <Purity> -ploidy <Ploidy> -F freq.list -O <Outdir>
```
## Test
Please go to the `example/` folder, and run `sh ./test_exon.sh` in `exon/` or `sh ./test_whole.sh` in `whole/`. 
If you install successfully, you will find results in the `output/`.


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
    * SamTools-1.3.1
    * bamutil Version: 1.0.14
    * picard-tools-2.1.0
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


## Basic Usage  
### Extract HLA-related reads
We should first extract HLA reads with enrichment-free data before profiling the HLA. Please first map the fastq reads to `hg19` or `hg38`, then use the `ExtractHLAread.sh` in `script/` folder to extract HLA-related reads. We use the script of Kourami with minor revision for this step. 
```
HLA-related reads extractor for HLA Typing
USAGE: <PATH-TO>/ExtractHLAread.sh -s <sample_id> -b <bamfile> -r <refGenome>

 -s          : desired sample name (ex: NA12878) [required]

 -b          : sorted and indexed bam (ex: NA12878.bam) [required]

 -r          : hg38 or hg19
```

### Whole Mode 
Users can integrate the PacBio, Nanopore, 10X linked-reads, Hi-C data to improve the performance.
```
Full-length HLA typing with paired-end reads. This script can use PacBio, Nanopore,
Hi-C, and 10X sequencing data to improve the phasing performance if provided.

Usage:
  sh SpecHLA_whole.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -t <sample.pacbio.fq.gz> -p <Asian>

Options:
  -n        Sample ID. <required>
  -1        The first fastq file. <required>
  -2        The second fastq file. <required>
  -o        The output folder to store the typing results.
  -t        Pacbio TGS fastq file.
  -e        Nanopore TGS fastq file.
  -x        Path of folder created by 10x demultiplexing. Prefix of the filenames of FASTQs
            should be the same as Sample ID. You can regard reads as normal NGS reads and use
            this parameter to adopt barcode information to improve phasing.
  -c        fwd hi-c fastq file.
  -d        rev hi-c fastq file.
  -p        The population of the sample: Asian, Black, or Caucasian. Use mean frequency
            if not provided.
  -f        True or False. The annotation database only includes the alleles with population
            frequency higher than zero if set True. Otherwise, it includes all alleles. Default
            is True.
  -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
            is 2. It should be set larger to infer novel alleles.
  -v        True or False. Consider long InDels if True, else only consider short variants.
            Default is False.
  -q        Minimum variant quality. Default is 0.01. Set it larger in high quality samples.
  -s        True or False. Use SpecHap to phase if True, else use PStrain. Default is True.
            We recommend use SpecHap.
  -a        Use this long InDel file if provided.
  -r        The minimum Minor Allele Frequency (MAF), default is 0.05.
  -h        Show this message.
```

### Exon Mode
```
The Exon version of SpecHLA, performs HLA assembly and HLA Typing.

Usage:
  sh SpecHLA_exon.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -p <Asian>

Options:
  -n        Sample ID.
  -1        The first fastq file.
  -2        The second fastq file.
  -o        The output folder to store the typing results.
  -p        The population of the sample: Asian, Black, or Caucasian. Use mean frequency
            if not provided.
  -f        True or False. The annotation database only includes the alleles with population
            frequency higher than zero if set True. Otherwise, it includes all alleles. Default
            is True.
  -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
            is 3. It should be set larger to infer novel alleles.
  -g        True or False. Split more blocks if set True; then rely more on allele database.
            Otherwise, split less blocks. Default is True.
  -s        True or False. Use SpecHap to phase if True, else use PStrain. Default is True.
            We recommend SpecHap.
  -h        Show this message.
```

### Pedigree phasing
After performing haplotyping, User can afford pedigree information to update the results.
```
Example: python3 phase_pedigree.py --trio NA12878:NA12891:NA12892 --workdir output/

required arguments:
  --trio           The trio infromation; give sample names in the order of
                   child:mother:father. Example: NA12878:NA12891:NA12892
  -o , --workdir   The directory consists of all samples' results. The workdir
                   running previous scipts.
```
The fresh reconstructed sequences can be found in the trio/ folder at each original sample's output folder.

### Detect HLA LOH in cancer samples
If the purity of cancer sample and the ploidy of HLA region are known, we can infer HLA LOH event.
```
Usage:
       perl cal.hla.copy.pl [Options] -purity <Purity> -ploidy <Ploidy> -F <filelist> -O <outdir>

OPTIONS:
        -purity  [f]  the purity of tumor sample. <required>
        -ploidy  [f]  the ploidy of tumor sample in HLA gene region. <required>
        -F       [s]  the filelist. <required>
                     Format, separated by table:  Sample/HLA_A_freq.txt
        -O       [s]  the output dir. <required>
```


### TGS only alignment-based typing Testing
If you only have PacBio or Nanopore reads, you can perform haplotyping with
```
sh SpecHLA_single.sh -n <sample> -t <sample.pacbio.fq.gz> -p <Asian>
```
Please install Longshot 0.4.1 and add it to the system path to call small variants with TGS data. Please use "sh SpecHLA_single.sh -h" to see detailed usage information.

### TGS only assembly-based typing Testing
To assemble the TGS reads, and perform HLA typing based on the contigs, please use `/assembly/MHC.TGS.mapping.sh`. 
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

## Interpret output
SpecHLA will automaticly build a folder named `output/` in the current path. In the `output/` folder, the results of each sample are saved in a folder named as the sample ID.  

In the directory of one specific sample, you will find the below files:
| Output | Description |
| --- | --- |
| hla.result.txt | HLA-typing results for all alleles |
| hla.allele.\*.HLA_\*.fasta | Sequence of each allele |
| HLA_*_freq.txt | Haplotype frequencies of each gene |
| HLA_*.rephase.vcf.gz | Phased vcf file for each gene |

If you performed pedigree phasing, you will find below files in the output/ folder.
| Output | Description |
| --- | --- |
| trio/HLA_*.rephase.vcf.gz | Phased vcf file after pedigree phasing |
| trio/hla.allele.\*.HLA_\*.fasta | Sequence of each allele after pedigree phasing|

An example for the formats of "hla.result.txt" is as below:
```
Sample  HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2 HLA_C_1 HLA_C_2 HLA_DPA1_1      HLA_DPA1_2      HLA_DPB1_1      HLA_DPB1_2      HLA_DQA1_1      HLA_DQA1_2      HLA_DQB1_1      HLA_DQB1_2      HLA_DRB1_1      HLA_DRB1_2
NA06985 A*02:01:01G     A*03:01:01G     B*07:02:01G     B*57:01:01G     C*07:02:01G     C*06:02:01G     DPA1*01:03:01G  DPA1*01:03:01G  DPB1*04:01:01G  DPB1*04:01:01G  DQA1*01:02:01G  DQA1*01:02:01G  DQB1*06:02:01G  DQB1*06:02:01G  DRB1*15:01:01G  DRB1*15:01:01G
Score   100     100     100             100     100     100     100     100     100     100     100     99.6376811594203        100     1       1
```
This means the typing result of each gene's haplotype. The third line indicates the match score between the assembly sequence and the typing allele in the database.   

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

