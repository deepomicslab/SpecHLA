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
### Programming 
* Python 3.7.0 or above  
* Python libraries:  
    * Numpy  
    * Scipy  
    * Pandas  
    * Pulp
    * Pysam  
* Perl 5 or above

### Reference database
HLAPro needs a HLA allele database. Please download db/ folder from https://drive.google.com/drive/folders/1M5GSF-oDTJ-tWgHt30TskOTWZvGBrD0Y?usp=sharing, and put it in the same path with the script/ folder. HLAPro will use it by default. 

### Third party packages
Please download the folder bin/ from https://drive.google.com/drive/folders/1M5GSF-oDTJ-tWgHt30TskOTWZvGBrD0Y?usp=sharing, and put it in the same path with the script/ folder. 

After that, please install Novoalign and Tabix by yourself. And link the binary files of novoalign, novoalign.lic, Tabix, python2, and python3 to the db/ folder. HLAPro will use them by default.

* Third party packages:
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
    * SpecHap v1.0.1

Then you can test HLAPro in the example/ folder. 

Also, you can install them by yourself, and link the binary files to the bin/ folder. 

## Test
Please go to the example/ folder, and run the test code in exon/ or whole/. 
```
sh ./test_exon.sh
or
sh ./test_whole.sh
```
If you install successfully, you will find results in the example/exon/output/ or example/whole/output/.

## Basic Usage  
### Extract HLA-related reads
Before profiling the HLA, we should first extract HLA reads. Then we use the extracted HLA reads to do the following steps. Please first map the fastq reads to hg19 or hg38, then use the ExtractHLAread.sh in script/ folder to extract HLA-related reads. We use the script of Kourami with minor revision for this step. 
```
USAGE: 
    sh ExtractHLAread.sh <sample_id> <bamfile> <refGenome>

    sample_id        : desired sample name (ex: NA12878) [required]

    bamfile          : sorted and indexed bam (ex: NA12878.bam) [required]

    ref              : hg38 or hg19

```
### Exon Mode (profile HLA in exons)
```
The Exon version of HLAPro, performs HLA assembly and HLA Typing.

Usage:
  sh HLAPro_exon.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -p <Asian>

Options:
  -n        Sample ID.
  -1        The first fastq file.
  -2        The second fastq file.
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
### Whole Mode (profile HLA with full length)
The Whole version of HLAPro, performs HLA assembly and HLA Typing with full-length.
```
Usage:
  sh HLAPro_whole.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -p <Asian>

Options:
  -n        Sample ID.
  -1        The first fastq file.
  -2        The second fastq file.
  -p        The population of the sample: Asian, Black, or Caucasian. Use mean frequency
            if not provided.
  -f        True or False. The annotation database only includes the alleles with population
            frequency higher than zero if set True. Otherwise, it includes all alleles. Default
            is True.
  -m        The maximum mismatch number tolerated in assigning gene-specific reads. Deault
            is 2. It should be set larger to infer novel alleles.
  -s        True or False. Use SpecHap to phase if True, else use PStrain. Default is True.
            We recommend use SpecHap.
  -h        Show this message.
```

### Extended Mode (integrating TGS, 10X, Hi-C reads)
Users can integrate the PacBio, Nanopore, 10X linked-reads, Hi-C data to improve performance in this mode.
```

The extended version of SpecHLA, performs HLA assembly and HLA Typing with full-length basd on NGS and TGS data.
This script can use PacBio, Nanopore, Hi-C, and 10X sequencing data to improve the phasing performance if provided.


Usage:
  sh SpecHLA_tgs.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -t <sample.pacbio.fq.gz> -p <Asian>

Options:
  -n        Sample ID. <required>
  -1        The first fastq file. <required>
  -2        The second fastq file. <required>
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
  -q        Minimum variant quality. Default is 0.01. Set larger in high quality samples.
  -s        True or False. Use SpecHap to phase if True, else use PStrain. Default is True.
            We recommend use SpecHap.
  -a        Use this long InDel file if provided.
  -r        The minimum Minor Allele Frequency (MAF), default is 0.05.
  -h        Show this message.
```
To use PacBio reads, you should install Minimap 2.17-r941, pbmm2-1.4.0, and pbsv Version 2.6.2 in the bin/ folder.
To use Nanopore ONT reads, you need to install Minimap 2.17-r941.
You need LongRanger2.1.2 to adopt 10X reads.


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

## Interpret output
HLAPro will automaticly build a folder named output/ in the current path. In the output/ folder, the results of each sample are saved in a folder named as the sample ID.  

In the directory of one specific sample, you will find the below files:
```
"hla.result.txt" contains the Typing results of each allele.
"hla.allele.*.HLA_*.fasta" contains the Haplotyping sequence results of each allele.
"HLA_*_freq.txt" contains the Haplotype Ratios in each gene.

```
An example for the formats of "hla.result.txt" is as below:
```
Sample  HLA_A_1 HLA_A_2 HLA_B_1 HLA_B_2 HLA_C_1 HLA_C_2 HLA_DPA1_1      HLA_DPA1_2      HLA_DPB1_1      HLA_DPB1_2      HLA_DQA1_1      HLA_DQA1_2      HLA_DQB1_1      HLA_DQB1_2      HLA_DRB1_1      HLA_DRB1_2
NA06985 A*02:01:01G     A*03:01:01G     B*07:02:01G     B*57:01:01G     C*07:02:01G     C*06:02:01G     DPA1*01:03:01G  DPA1*01:03:01G  DPB1*04:01:01G  DPB1*04:01:01G  DQA1*01:02:01G  DQA1*01:02:01G  DQB1*06:02:01G  DQB1*06:02:01G  DRB1*15:01:01G  DRB1*15:01:01G
Score   100     100     100             100     100     100     100     100     100     100     100     99.6376811594203        100     1       1
```
This means the typing result of each gene's haplotype. The third line indicates the match score between the assembly sequence and the typing allele in the database.   

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
Should you have any queries, please feel free to contact us, we will reply as soon as possible (swang66-c@my.cityu.edu.hk).

