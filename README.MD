# SpecHLA: full-resolution HLA typing from sequencing data

## Install  
First, create the env with conda, and activate the env.
```
git clone https://github.com/deepomicslab/SpecHLA.git --depth 1
cd SpecHLA/
conda env create --prefix=./spechla_env -f environment.yml
conda activate ./spechla_env
```
Second, index the database and install the packages.
```
unset LD_LIBRARY_PATH && unset LIBRARY_PATH 
bash index.sh
```
Note: 1) SpecHLA requires a License of [Novoalign](https://www.novocraft.com/buy-now/) in `bin/`, if not detected, it uses `bowtie2` as replacement automatically. 
The License file of `Novoalign` should be put in the `bin/` folder `before` running `bash index.sh`. 2) SpecHLA now supports Linux and Windows WSL.

## Test
Please go to the `example/` folder, run SpecHLA with given scripts, and find results in the `output/`.

## Basic Usage  
### Extract HLA-related reads
We should first extract HLA reads with enrichment-free data. Please map the reads to `hg19` or `hg38`, then use the `ExtractHLAread.sh` in `script/` folder to extract HLA-related reads. We use the script of [Kourami](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1388-2) with minor revision for this step. 
Extract HLA-related reads by
```
USAGE: bash script/ExtractHLAread.sh -s <sample_id> -b <bamfile> -r <refGenome> -o <outdir>

 -s          : desired sample name (ex: NA12878) [required]

 -b          : sorted and indexed bam (ex: NA12878.bam) [required]

 -r          : hg38 or hg19

 -o          : folder to save extracted reads [required]
```

### HLA Typing 
Full-resolution and exon HLA typing with paired-end reads. SpecHLA can use PacBio, Nanopore, Hi-C, and 10X data to improve the phasing performance if provided. 

Perform full-resolution HLA typing with `paired-end` reads by
```
sh script/whole/SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/
```
Perform exon HLA typing with paired-end reads by
```
sh script/whole/SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ -u 1
```
Perform full-resolution HLA typing with `paired-end` reads and `PacBio` reads by
```
sh script/whole/SpecHLA.sh -n <sample> -t <sample.pacbio.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `Nanopore` reads by
```
sh script/whole/SpecHLA.sh -n <sample> -e <sample.nanopore.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `Hi-C` reads by
```
sh script/whole/SpecHLA.sh -n <sample> -c <sample.hic.fwd.fq.gz> -d <sample.hic.rev.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `10X` linked reads by (LongRanger should be installed in system env)
```
sh script/whole/SpecHLA.sh -n <sample> -x <sample.10x.read.folder> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Consider `long Indels` and use `population information` for annotation by
```
sh script/whole/SpecHLA.sh -n <sample> -v True -p <Asia> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing only with `long reads` by 
```
python3 script/long_read_typing.py -n <sample> -r <Pacbio or Nanopore reads> -o outdir/ 
```

### Full arguments
```
SpecHLA: Full-resolution HLA typing from paired-end, PacBio, Nanopore,
Hi-C, and 10X sequencing. Supports WGS, WES, and RNASeq data.


Usage:
  sh SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -t <sample.pacbio.fq.gz> -p <Asian>

Options:
  -n        Sample ID. <required>
  -1        The first fastq file. <required>
  -2        The second fastq file. <required>
  -o        The output folder to store the typing results.
  -u        Choose full-length or exon typing. 0 indicates full-length, 1 means exon,
            default is to perform full-length typing.
  -p        The population of the sample [Asian, Black, Caucasian, Unknown, nonuse].
            Default is Unknown, meaning use mean frequency. nonuse indicates only adopting
            mapping score and considering alleles with frequency as zero.
  -t        Pacbio fastq file.
  -e        Nanopore fastq file.
  -c        fwd hi-c fastq file.
  -d        rev hi-c fastq file.
  -x        Path of folder created by 10x demultiplexing. Prefix of the filenames of FASTQs
            should be the same as Sample ID. Please install Longranger in the system env.
  -w        The weight to use linkage info from genotype frequency [0-1], default is 0 that means
            not use, 1 means only use imbalance info, other values integrate use both.
  -j        Number of threads [5]
  -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
            is 2. It should be set larger to infer novel alleles.
  -y        The minimum different mapping score between the best- and second-best aligned gene.
            Discard the read if the score is lower than this value. Deault is 0.1.
  -v        True or False. Consider long InDels if True, else only consider short variants.
            Default is False.
  -q        Minimum variant quality. Default is 0.01. Set it larger in high quality samples.
  -s        Minimum variant depth. Default is 5.
  -a        Use this long InDel file if provided.
  -r        The minimum Minor Allele Frequency (MAF), default is 0.05 for full length and
            0.1 for exon typing.
  -g        Whether use G-translate in annotation [1|0], default is 0.
  -k        The mean depth in a window lower than this value will be masked by N, default is 5.
            Set 0 to avoid masking.
  -z        Whether only mask exon region, True or False, default is False.
  -f        The trio infromation; child:parent_1:parent_2 [Example: NA12878:NA12891:NA12892].
            Note: this parameter should be used after performing SpecHLA already.
  -b        Whether use database for unlinked block phasing [1|0], default is 1.
  -h        Show this message.
```

### Pedigree phasing
After performing HLA typing, we can afford trio information to update the results of child by
```
sh script/SpecHLA.sh -n <child> -1 <child.fq.1.gz> -2 <child.fq.2.gz> -o outdir/ --trio child:parent1:parent2
```
The fresh results can be found in the trio/ folder inside the original result folder.

### Detect HLA LOH in cancer samples
If the purity of cancer sample and the ploidy of HLA region are known, we can infer HLA LOH event by
```
usage:  perl script/cal.hla.copy.pl [Options] -S <samplename> -C <5> -purity <Purity> -ploidy <Ploidy> -F <filelist> -T <hla.result.txt> -O <outdir>

        OPTIONS:
        -purity  [f]  the purity of tumor sample. <required>
        -ploidy  [f]  the ploidy of tumor sample in HLA gene region. <required>
        -S       [s]  sample name <required>
        -C       [f]  the cutoff of heterogeneous snp number. <required>
        -F       [s]  the filelist. <required>     Format, separated by table: ./HLA_A_freq.txt
        -T       [s]  the hla typing result file of Spechla <required>
        -O       [s]  the output dir. <required>
```
Use `ls typing_result_dir/*_freq.txt >freq.list` to generate the file required by `-F`.

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

If you performed pedigree phasing, you will find below files.
| Output | Description |
| --- | --- |
| trio/HLA_*.trio.vcf.gz | Phased vcf file after pedigree phasing |
| hla.allele.\*.HLA_\*.fasta | Allele sequence after pedigree phasing|

An example for "hla.result.txt" is as below:
```
Sample  HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2 HLA_C_1 HLA_C_2 HLA_DPA1_1      HLA_DPA1_2      HLA_DPB1_1      HLA_DPB1_2      HLA_DQA1_1    HLA_DQA1_2      HLA_DQB1_1      HLA_DQB1_2      HLA_DRB1_1      HLA_DRB1_2
HG00118 A*24:02:01:01   A*02:01:01:01   B*07:02:01:01   B*40:01:02:02   C*03:04:01:01   C*07:02:01:03   DPA1*01:03:01:02        DPA1*01:03:01:02    DPB1*04:01:01:01        DPB1*04:01:01:01        DQA1*01:02:01:01        DQA1*01:02:01:01        DQB1*06:02:01:01        DQB1*06:02:01:01    DRB1*15:01:01:01        DRB1*15:01:01:01
```
An example for "HLA_*_freq.txt" is as below:
```
# Haplotype     Frequency
hla.allele.1.HLA_A.fasta 0.532
hla.allele.2.HLA_A.fasta 0.468
# The number of heterozygous variant is 30
```

An example for HLA LOH result is as below:
```
Sample	HLA	Allele1	Allele2	copyratio	KeptHLA	LossHLA
DNA_CRL-2314_p3_s6	A	A*03:01:01:01	A*02:05:01:01	2:0	A*03:01:01:01	A*02:05:01:01
DNA_CRL-2314_p3_s6	B	B*35:03:01:01	B*49:01:01:01	2:0	B*35:03:01:01	B*49:01:01:01
DNA_CRL-2314_p3_s6	C	C*12:03:01:01	C*07:01:01:01	2:0	C*12:03:01:01	C*07:01:01:01
DNA_CRL-2314_p3_s6	DPA1	DPA1*01:03:01:01	DPA1*01:03:01:01	2:0	DPA1*01:03:01:01	homogeneous
DNA_CRL-2314_p3_s6	DPB1	DPB1*15:01:01:01	DPB1*04:01:01:01	0:2	DPB1*04:01:01:01	DPB1*15:01:01:01
DNA_CRL-2314_p3_s6	DQA1	DQA1*01:02:01:01	DQA1*05:05:01:01	2:0	DQA1*01:02:01:01	DQA1*05:05:01:01
DNA_CRL-2314_p3_s6	DQB1	DQB1*06:02:01:01	DQB1*03:19:01:01	2:0	DQB1*06:02:01:01	DQB1*03:19:01:01
DNA_CRL-2314_p3_s6	DRB1	DRB1*11:02:01:01	DRB1*15:01:01:01	0:2	DRB1*15:01:01:01	DRB1*11:02:01:01
```


## Dependencies 

### Systematic requirement
SpecHLA requires `conda 4.12.0` and `GCC 11.3.0` for environment construction and installation.

### Programming 
* python=3.8.12 or above  
* Python libraries: numpy=1.22.3, pulp=2.6.0, pysam=0.19.0, scipy=1.8.0, and biopython=1.79. 
* Perl 5 or above

### Third party packages
SpecHLA enables automatic installation of these third party packages. 
Also, users can install the softwares by themselves, please remind link the executable file to the `bin/` folder.
To apply [Novoalign](https://www.novocraft.com/buy-now/), please put the novoalign.lic (License file) in the `bin/` folder.

* Third party packages:
    * SamTools-1.3
    * bamutil Version: 1.0.14
    * BWA-0.7.17-r1188
    * blast 2.3.0
    * BLAT v. 36x2
    * bedtools-v2.26.0
    * bcftools-Version: 1.9
    * Freebayes-v1.2.0-2-g29c4002
    * Novoalign V4.02.01
    * ScanIndel v1.3
    * Fermikit-0.13
    * SpecHap v1.0.1 and ExtractHAIRs
    * Minimap 2.17-r941
    * pbmm2-1.4.0
    * pbsv Version 2.6.2
    * Bowtie 2 version 2.3.4.1
    * longshot-0.4.5
   

## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).