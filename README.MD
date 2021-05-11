# HLAPro: a hybrid framework for novel HLA allele discovery and allele-specific HLA LOH inference

## Install  
```
* Exon Mode
cd HLAPro/script/exon/
sh HLAPro_exon.sh -h

* Whole Mode
cd HLAPro/script/whole/
sh HLAPro_whole.sh -h

* Detect LOH events
ls output/*/*_freq.txt >freq.list
perl HLAPro/script/cal.hla.copy.pl  -purity <Purity> -ploidy <Ploidy> -F freq.list -O <Outdir>

```
## Dependencies 
* Python 3.7.0 or above  
* Python libraries:  
    * Numpy  
    * Scipy  
    * Pandas  
    * Pulp
    * Pysam  
* Perl 5 or above

HLAPro needs a HLA allele database. We have uploaded the database in the db/ folder to the link https://drive.google.com/drive/folders/1M5GSF-oDTJ-tWgHt30TskOTWZvGBrD0Y?usp=sharing, you can directly download the db/ folder, and put it in the same path with the script/ folder. HLAPro will use them by default. 

HLAPro relies on these softwares, we have uploaded the binary file of them to the folder bin/ in the link https://drive.google.com/drive/folders/1M5GSF-oDTJ-tWgHt30TskOTWZvGBrD0Y?usp=sharing, you can directly download the bin/ folder, and put it in the same path with the script/ folder. HLAPro could use them by default. Also, you can install them by yourself, and build a folder named bin/ in the same folder of script/, and link the binery files to the bin/ folder. HLAPro will use them by default.   
* Third party pipelines:
    * SamTools-1.3.1
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

Then you can use HLAPro's codes in the script/ folder. 


## Basic Usage  
### Extract HLA-related reads
Before profiling the HLA, we should first extract HLA reads. Then we use the extracted HLA reads to do the following steps. Please first map the fastq reads to hg19 or hg38, then use the ExtractHLAread.sh in script/ folder to extract HLA-related reads.
```
USAGE: 
    sh ExtractHLAread.sh <sample_id> <bamfile> <refGenome>

    sample_id        : desired sample name (ex: NA12878) [required]

    bamfile          : sorted and indexed bam (ex: NA12878.bam) [required]

    ref              : hg38 or hg19

```
### Exon Mode (profile HLA in exons)
```
Usage:
  sh HLAPro_exon.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -p <Asian>

Options:
  -n        Sample ID.
  -1        The first fastq file.
  -2        The second fastq file.
  -p        The population of the sample: Asian, Black, or Caucasus. Use mean frequency
            if not provided.
  -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
            is 2. It should be set larger to infer novel alleles.
  -g        True or False. The ratio of two haplotype will be regarded the same if set True.
            Otherwise, we will calculate the haplotype ratio. Default is False.
  -h        Show this message.
```
### Whole Mode (profile HLA with full length)
The Whole version of HLAPro, performs HLA assembly and HLA Typing with full-length.
```
Usage:
  sh HLAPro_whole.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -p <Asian>

Options:
  -n        Sample ID..
  -1        The first fastq file.
  -2        The second fastq file.
  -p        The population of the sample: Asian, Black, or Caucasus. Use mean frequency
            if not provided.
  -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
            is 2. It should be set larger to infer novel alleles.
  -g        True or False. The ratio of two haplotype will be regarded the same if set True.
            Otherwise, we will calculate the haplotype ratio. Default is False.
  -h        Show this message.
```

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

## Interpret output
HLAPro will automaticly build a folder named output/ in the current path. In the output/ folder, the results of each sample are saved in a folder named as the sample ID.  

In the directory of one specific sample, you will find the below files:
```
"hla.result.txt" contains the Typing results of each allele.
"result.HLA_*.fasta" contains the Assembly results of each allele.
"HLA_*_freq.txt" contains the Haplotype Ratios of each gene.

```
An example for the formats of "hla.result.txt" is as below:
```
Sample  HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2 HLA_C_1 HLA_C_2 HLA_DPA1_1      HLA_DPA1_2      HLA_DPB1_1      HLA_DPB1_2      HLA_DQA1_1      HLA_DQA1_2      HLA_DQB1_1      HLA_DQB1_2      HLA_DRB1_1      HLA_DRB1_2
NA06985 A*02:01:01G     A*03:01:01G     B*07:02:01G     B*57:01:01G     C*07:02:01G     C*06:02:01G     DPA1*01:03:01G  DPA1*01:03:01G  DPB1*04:01:01G  DPB1*04:01:01G  DQA1*01:02:01G  DQA1*01:02:01G  DQB1*06:02:01G  DQB1*06:02:01G  DRB1*15:01:01G  DRB1*15:01:01G
Score   100     100     100             100     100     100     100     100     100     100     100     99.6376811594203        100     1       1
```
This means the typing result of each gene's haplotype. The third line indicates the match score between the assembly sequence and the typing allele in the database.   

## Test
Please go to the example/ folder, and run the test code in exon/ or whole/. 
```
sh ./test_exon.sh
or
sh ./test_whole.sh
```
If you install successfully, you will find results in the example/exon/output/ or example/whole/output/.

## Simulation
Please go to the simulation/ folder, and run simu.hla.sh like:
```
Simulation HLA capture reads for HLA Typing
USAGE: <PATH-TO>/simu.hla.sh -n <sample_size> -r <read_length> -1 <depth of haplotype1> -2 <depth of haplotype2>

 -n        : sample size of simulation [required]

 -r        : read length [required]

 -1        : depth of haplotype1 [required]

 -2        : depth of haplotype2 [required]
```

## Getting help
Should you have any queries, please feel free to contact me, I will reply as soon as possible (swang66-c@my.cityu.edu.hk).

