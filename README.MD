# SpecHLA: full-resolution HLA typing from sequencing data

SpecHLA is a software package leveraging reads binning and local assembly to achieve accurate full-resolution HLA typing and loss-of-heterozygosity detection. 
- SpecHLA reconstructs diploid sequences of HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, and -DRB1 genes.
- SpecHLA accepts short-reads-only, long-reads-only, and short+long-reads data.
- SpecHLA supports WES, WGS, and RNA-seq data.
- SpecHLA detects HLA loss-of-heterozygosity events.


## Install  
First, create the env with conda, and activate the env.
```
git clone https://github.com/deepomicslab/SpecHLA.git --depth 1
cd SpecHLA/
conda env create --prefix=./spechla_env -f environment.yml
conda activate ./spechla_env
```
Second, make the softwares in bin/ executable.
```
chmod +x -R bin/*
```
Third, index the database and install the packages.
```
unset LD_LIBRARY_PATH && unset LIBRARY_PATH 
bash index.sh
```
Perform SpecHLA with
```
bash script/whole/SpecHLA.sh -h
```
With only long reads, run
```
python3 script/long_read_typing.py -h
```

Note:

- SpecHLA requires a License of [Novoalign](https://www.novocraft.com/buy-now/) in `bin/`. If not detected, it uses `bowtie2` as a replacement automatically. 
The License file of `Novoalign` should be put in the `bin/` folder `before` running `bash index.sh`. 
- If you want to run SpecHLA with `only` long-read data, there is no need for the `Novoalign` license and running `bash index.sh`. 
You need only construct the environment from `conda env create --prefix=./spechla_env -f long_only_env.yml`.
- SpecHLA now supports Linux and Windows WSL systems.
- SpecHLA does not accept short single-end reads.
- The conda environment should be located in the local folder.


## Test
Please go to the `example/` folder, run SpecHLA with given scripts, and find results in the `output/`.

## Basic Usage  

### Main functions
| Scripts | Description |
| --- | --- |
|script/ExtractHLAread.sh| Extract HLA reads from enrichment-free data.|
|script/whole/SpecHLA.sh| HLA typing with paired-end (PE), PE+long reads, PE+HiC, or PE+10X data.  |
|script/long_read_typing.py|HLA typing with only long-read data.|
|script/typing_from_assembly.py  |HLA Typing from diploid assemblies.|
|script/cal.hla.copy.pl|Detect HLA LOH events based on SpecHLA's typing results.|


### Extract HLA-related reads
First extract HLA reads with enrichment-free data. Otherwise, HLA typing would be slow. Map reads to `hg19` or `hg38`, then use `script/ExtractHLAread.sh` to extract HLA-related reads. We use the script of [Kourami](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1388-2) with minor revision for this step. 
Extract HLA-related reads by
```
USAGE: bash script/ExtractHLAread.sh -s <sample_id> -b <bamfile> -r <refGenome> -o <outdir>

 -s          : desired sample name (ex: NA12878) [required]

 -b          : sorted and indexed bam or cram (ex: NA12878.bam) [required]

 -r          : hg38 or hg19

 -o          : folder to save extracted reads [required]
```

### HLA Typing 
Full-resolution and exon HLA typing using SpecHLA. With Exome data like WES or RNASeq, only support exon typing. For efficient HLA typing, we strongly recommend utilizing only HLA-related reads. Specifically, for enrichment-free data, we recommend first performing the aforementioned step.

Perform full-resolution HLA typing with `paired-end` reads by
```
bash script/whole/SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/
```
Perform exon HLA typing with `paired-end` reads by
```
bash script/whole/SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ -u 1
```
Perform full-resolution HLA typing with `paired-end` reads and `PacBio` reads by
```
bash script/whole/SpecHLA.sh -n <sample> -t <sample.pacbio.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `Nanopore` reads by
```
bash script/whole/SpecHLA.sh -n <sample> -e <sample.nanopore.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `Hi-C` reads by
```
bash script/whole/SpecHLA.sh -n <sample> -c <sample.hic.fwd.fq.gz> -d <sample.hic.rev.fq.gz> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Perform full-resolution HLA typing with `paired-end` reads and `10X` linked reads by (LongRanger should be installed in system env)
```
bash script/whole/SpecHLA.sh -n <sample> -x <sample.10x.read.folder> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```
Consider `long Indels` and use `population information` for annotation by
```
bash script/whole/SpecHLA.sh -n <sample> -v True -p <Asia> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o outdir/ 
```

Full arguments can be seen in
```
SpecHLA: Full-resolution HLA typing from sequencing data.

Note:
  1) Use HLA reads only, otherwise, it would be slow. Use ExtractHLAread.sh to extract HLA reads first.
  2) WGS, WES, and RNASeq data are supported.
  3) With Exome data like WES or RNASeq, must select exon typing  (-u 0).
  4) Short single-end read data are not supported.

Usage:
  bash SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -o <outdir>

Options:
  -n        Sample ID. <required>
  -1        The first fastq file of paired-end data. <required>
  -2        The second fastq file of paired-end data. <required>
  -o        The output folder to store the typing results. Default is ./output
  -u        Choose full-length or exon typing [0|1]. 0 indicates full-length, 1 means exon,
            default is 0. With Exome or RNA data, must select 1 (i.e., exon typing).
  -p        The population of the sample [Asian, Black, Caucasian, Unknown, nonuse] for annotation.
            Default is Unknown, meaning use mean allele frequency in all populations. nonuse indicates
            only adopting mapping score and considering zero-frequency alleles.
  -j        Number of threads [5].
  -t        Pacbio fastq file.
  -e        Nanopore fastq file.
  -c        fwd hi-c fastq file.
  -d        rev hi-c fastq file.
  -x        Path of folder created by 10x demultiplexing. Prefix of the filenames of FASTQs
            should be the same as Sample ID. Please install Longranger in the system env.
  -w        The weight to use allele imbalance info for phasing [0-1]. Default is 0 that means
            not use. 1 means only use imbalance info; other values integrate reads and allele imbalance.
  -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
            is 2. It should be set larger to infer novel alleles.
  -y        The minimum different mapping score between the best and second-best aligned genes.
            Discard the read if the score is lower than this value. Deault is 0.1.
  -v        True or False. Consider long InDels if True, else only consider small variants.
            Default is False.
  -q        Minimum variant quality. Default is 0.01. Set it larger in high quality samples.
  -s        Minimum variant depth. Default is 5.
  -a        Use this long InDel file if provided.
  -r        The minimum Minor Allele Frequency (MAF), default is 0.05 for full length and
            0.1 for exon typing.
  -k        The mean depth in a window lower than this value will be masked by N, default is 5.
            Set 0 to avoid masking.
  -z        Whether only mask exon region, True or False, default is False.
  -f        The trio infromation; child:parent_1:parent_2 [Example: NA12878:NA12891:NA12892]. If provided,
            use trio info to improve typing. Note: use it after performing SpecHLA once already.
  -b        Whether use database for unlinked block phasing [0|1], default is 1 (i.e., use).
  -i        Location of the IMGT/HLA database folder, default is db.
  -l        Whether remove all tmp files [0|1], default is 1.
  -h        Show this message.
```

### HLA typing with long-read data alone
Perform HLA typing only with `long reads` by 
```
usage: python3 long_read_typing.py -h

HLA Typing with only long-read data.

Required arguments:
  -r                  Long-read fastq file. PacBio or Nanopore. (default: None)
  -n                  Sample ID (default: None)
  -o                  The output folder to store the typing results. (default: ./output)

Optional arguments:
  -p                  The population of the sample [Asian, Black, Caucasian, Unknown, nonuse] for annotation.
                        Unknown means use mean allele frequency in all populations. nonuse indicates only adopting
                        mapping score and considering zero-frequency alleles. (default: Unknown)
  -j                  Number of threads. (default: 5)
  -d                  Minimum score difference to assign a read to a gene. (default: 0.001)
  -g                  Whether use G group resolution annotation [0|1]. (default: 0)
  -m                  1 represents typing, 0 means only read assignment (default: 1)
  -k                  The mean depth in a window lower than this value will be masked by N, set 0 to avoid masking
                        (default: 5)
  -a                  Prefix of filtered fastq file. (default: long_read)
  -y                  Read type, [nanopore|pacbio]. (default: pacbio)
  --minimap_index     Whether build Minimap2 index for the reference [0|1]. Using index can reduce memory usage.
                        (default: 0)
  --db                db dir. (default: /home/wangshuai/softwares/SpecHLA/script/../db/)
  --strand_bias_pvalue_cutoff
                        Remove a variant if the allele observations are biased toward one strand (forward or reverse).
                        Recommand setting 0 to high-depth data. (default: 0.01)
  --seed              seed to generate random numbers (default: 8)
  --max_depth         maximum depth for each HLA locus. Downsample if exceed this value. (default: 2000)
  -h, --help
```

### HLA Typing from diploid assemblies
```
usage: python3 typing_from_assembly.py -h

HLA Typing from diploid assemblies.

Required arguments:
  -1        Assembly file of the first haplotype in fasta formate (default: None)
  -2        Assembly file of the second haplotype in fasta formate (default: None)
  -n        Sample ID (default: None)
  -o        The output folder to store the typing results. (default: ./output)

Optional arguments:
  -j        Number of threads. (default: 10)
  -h, --help
```
An example of running command is 
```
python SpecHLA/script/typing_from_assembly.py -j 15 -o output -n v12_HG00096 -1 v12_HG00096_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta -2 v12_HG00096_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta
```
The result is stored in `<sample>_extracted_HLA_alleles.fasta`, please see result description in the below `Interpret output` section.


### Pedigree phasing
After performing HLA typing, we can afford trio information to update the results of child by
```
bash script/SpecHLA.sh -n <child> -1 <child.fq.1.gz> -2 <child.fq.2.gz> -o outdir/ --trio child:parent1:parent2
```
The fresh results can be found in the trio/ folder inside the original result folder.

### Detect HLA LOH in tumor samples
Based on the tumor purity and ploidy, we can infer HLA LOH event by
```
usage:  perl cal.hla.copy.pl [Options] -S <samplename> -C <5> -purity <Purity> -ploidy <Ploidy> -F <filelist> -T <hla.result.txt> -O <outdir>

        OPTIONS:
        -purity  [f]  the purity of the tumor sample. <required>
        -ploidy  [f]  the ploidy of the tumor sample. <required> Note: tumor purity and ploidy can be inferred by software like ABSOLUTE and ASCAT.
        -S       [s]  sample name <required>
        -C       [f]  the cutoff of heterogeneous snp number. Default is 5.
        -F       [s]  the HLA_*_freq.txt file list. <required>  Obtained by "ls typing_result_dir/*_freq.txt >freq.list"
        -T       [s]  the hla typing result file of Spechla <required>
        -O       [s]  the output dir. <required> Need exist before running.

The result file "merge.hla.copy.txt" can be found in the output dir.
```
Use `ls typing_result_dir/*_freq.txt >freq.list` to generate the file required by `-F`. The command example is
```
perl script/cal.hla.copy.pl -purity 0.5 -ploidy 2 -S test -F freq.list -T hla.result.txt -O ./
```


## Interpret output
In the denoted outdir, the results of each sample are saved in a folder named as the sample ID.  

In the directory of one specific sample, you will find the below files:
| Output | Description |
| --- | --- |
| hla.result.txt | HLA-typing results for all alleles |
| hla.result.details.txt | All the alleles with the highest mapping score in annotation |
| hla.result.g.group.txt | G group resolution HLA type |
| hla.allele.\*.HLA_\*.fasta | Sequence of each allele (the low-depth region is masked by N) |
| HLA_*_freq.txt | Haplotype frequencies of each gene |
| HLA_*.rephase.vcf.gz | Phased vcf file for each gene |
| hlala.like.results.txt | Report all potential types like HLA*LA, only generated by long-read-only mode |

If you performed pedigree phasing, you will find below files.
| Output | Description |
| --- | --- |
| trio/HLA_*.trio.vcf.gz | Phased vcf file after pedigree phasing |
| hla.allele.\*.HLA_\*.fasta | Allele sequence after pedigree phasing|


1. **An example for `hla.result.txt` is as below:** 
```
Sample  HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2 HLA_C_1 HLA_C_2 HLA_DPA1_1      HLA_DPA1_2      HLA_DPB1_1      HLA_DPB1_2      HLA_DQA1_1    HLA_DQA1_2      HLA_DQB1_1      HLA_DQB1_2      HLA_DRB1_1      HLA_DRB1_2
HG00118 A*24:02:01:01   A*02:01:01:01   B*07:02:01:01   B*40:01:02:02   C*03:04:01:01   C*07:02:01:03   DPA1*01:03:01:02        DPA1*01:03:01:02    DPB1*04:01:01:01        DPB1*04:01:01:01        DQA1*01:02:01:01        DQA1*01:02:01:01        DQB1*06:02:01:01        DQB1*06:02:01:01    DRB1*15:01:01:01        DRB1*15:01:01:01
```
Note: 
- We handle each gene independently. Hence the result does not contain the linkage information between different genes.
- The Symbol `-` means low confidence sequence which cannot map to any allele in the database.

2. **An example for `HLA_*_freq.txt` is as below:**
```
# Haplotype     Frequency
hla.allele.1.HLA_A.fasta 0.532
hla.allele.2.HLA_A.fasta 0.468
# The number of heterozygous variant is 30
```

3. **After typing from diploid assmebly, HLA sequences are in `<sample>_extracted_HLA_alleles.fasta`, and HLA types can be got by `cat <sample>_extracted_HLA_alleles.fasta| grep >`, an example is as below:**
```
>v12_HG00096.h1.HLA-A   cluster13_contig_98:2073847-2077365     A*29:02:01:01   # version: IPD-IMGT/HLA 3.51.0
>v12_HG00096.h1.HLA-B   cluster13_contig_98:657830-661911       B*44:03:01:01   # version: IPD-IMGT/HLA 3.51.0
...
>v12_HG00096.h2.HLA-A   cluster13_scaffold_119:2516943-2520446  A*01:01:01:01   # version: IPD-IMGT/HLA 3.51.0
...
```
Interpret each column in the annotation line as
|Column| Description |
| --- | --- |
|first| gene and haplotype, h1 means the first haplotype|
|  second |  the region to extract the sequence from the assembly |
| third  | the best matched IMGT allele, i.e., the typing result  |
|  four | IMGT version  |

4. **An example for HLA LOH result file: `merge.hla.copy.txt` is as below:**
```
Sample  HLA     Allele1                 Allele2              copyratio       KeptHLA                 LossHLA                 Freq1   Freq2   Purity  Het_num LOH
test    A       A*24:02:01:01           A*02:01:01:01           1:1          A*24:02:01:01           A*02:01:01:01           0.507   0.493   0.5     100     N
test    B       B*07:02:01:01           B*40:01:02:01           1:1          B*07:02:01:01           B*40:01:02:01           0.502   0.498   0.5     32      N
test    C       C*03:04:01:01           C*07:02:01:03           1:1          C*07:02:01:03           C*03:04:01:01           0.433   0.567   0.5     99      N
test    DPA1    DPA1*01:03:01:02        DPA1*01:03:01:02        2:0          DPA1*01:03:01:02        homogeneous             1       0       0.5     0       N
test    DPB1    DPB1*04:01:01:01        DPB1*04:01:01:01        2:0          DPB1*04:01:01:01        DPB1*04:01:01:01        0.909   0.091   0.5     1       N
test    DQA1    DQA1*01:02:01:01        DQA1*01:02:01:01        1:1          DQA1*01:02:01:01        DQA1*01:02:01:01        0.548   0.452   0.5     1       N
test    DQB1    DQB1*06:02:01:01        DQB1*06:02:01:01        2:0          DQB1*06:02:01:01        homogeneous             1       0       0.5     0       N
test    DRB1    DRB1*15:01:01:01        DRB1*15:01:01:01        2:0          DRB1*15:01:01:01        DRB1*15:01:01:01        0.65    0.35    0.5     2       N
```
Interpret each column as
|Header| Description |
| --- | --- |
|Sample| Sample Name|
|HLA| HLA locus|
|Allele1| HLA type of the first allele|
|Allele2| HLA type of the second allele|
|copyratio| Copy numbers of two alleles|
|KeptHLA| The remaining allele|
|LossHLA| The possible lost allele|
|Freq1| Frequency of the first allele|
|Freq2|Frequency of the second allele|
|Purity| Tumor purity|
|Het_num|Number of heterozygous variant|
|LOH| Whether LOH exists (Y or N)|

## Update to latest IMGT/HLA database
To nomenclature the HLA alleles using the latest IMGT/HLA database, please execute the following command:
```
cd SpecHLA/bin
perl renew_HLA_annotation_db.pl
```
This script will download the latest IMGT/HLA database from Github, and preprocess it for SpecHLA to use.

## Dependencies 

### Systematic requirement
SpecHLA requires `conda 4.12.0+`, `cmake 3.16.3+`, and `GCC 9.4.0+` for environment construction and software installation.

### Programming 
* python=3.8.12 or above  
* Python libraries: numpy=1.22.3, pulp=2.6.0, pysam=0.19.0, scipy=1.8.0, and biopython=1.79. 
* Perl 5 or above

### Third party packages
SpecHLA enables automatic installation of these third party packages using `conda`. 
Also, we can install the softwares by ourselves, and add their locations to system path (Exception: put bcftools, fermikit, and novoalign to `bin/` folder).
After applying for the authority of [Novoalign](https://www.novocraft.com/buy-now/), please put the novoalign.lic (License file) in the `bin/` folder.

* Third party packages:
    * SamTools-1.3.1
    * bamutil Version: 1.0.14
    * BWA-0.7.17-r1188
    * blast 2.12.0
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
   
## Citation
[Wang, Shuai, Mengyao Wang, Lingxi Chen, Guangze Pan, Yanfei Wang, and Shuai Cheng Li. "SpecHLA enables full-resolution HLA typing from sequencing data." Cell Reports Methods (2023).](https://doi.org/10.1016/j.crmeth.2023.100589)

## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).
