# SpecHap

SpecHap is an ultra fast phasing algorithm based on spectral graph analysis. SpecHap currently support general WGS sequencing, Hi-C, 10X linked-reads, PacBio SMRT and Oxford Nanopore .


## Getting Started

To build SpecHap, run the following command:

```
cd /path/to/SpecHap/
mkdir build
cd build
cmake ..
make && make install
```

To see the option SpecHap support, run the following command:

```
SpecHap --help
```

A modified utility software ExtractHair, originally from [HAPCUT2](https://github.com/vibansal/HapCUT2), is needed for fragment processing. To install, run
```
cd /path/to/SpecHap/
git submodule init && git submodule update
cd submodules/htslib
git checkout 26229a3
cd ../samtools
git checkout 255f97d
cd ../../hair-src
make
```
### Prerequisites

SpecHap relies on ARPACK for Eigen-calculation. To gain stable utilization, we recommend [arapack-ng](https://github.com/opencollab/arpack-ng).

Arpack-ng can be easily compile with cmake, ensure you have BLAS and LAPACK installed before compiling.

To build [arapack-ng](https://github.com/opencollab/arpack-ng), try run
```
cd /path/to/arpack-ng/
sh bootstrap
./configure --enable-icb
make
make && make install
```

[Htslib](https://github.com/samtools/htslib) is also required. To install [htslib](https://github.com/samtools/htslib), simply run
```
cd /path/to/htslib
autoheader  #required if htslib is cloned from github
autoconf    #required if htslib is cloned from github
./configure
make && make install
```



### Using SpecHap

#### Data preprocessing
SpecHap requires at least a fragment file and a bgziped and indexed VCF to perform phasing. Ensure your VCF is sorted by position. The fragment file is also required to be sorted by the same order.

To generate the fragment file, run the following command
```
extractHAIRS --bam /your/bam/file --VCF /your/vcf/file --out fragment_file
``` 
With Hi-C sequenced file, try 
```
extractHAIRS --bam /your/bam/file --VCF /your/vcf/file --out fragment_file --hic 1
```
With 10X linked reads, try 
```
extractHAIRS --bam /your/bam/file --VCF /your/vcf/file --out fragment_file --10x 1
```
You also need a bed file indicating each barcode's inferred spanning range. You can use the BarcodeExtract to do your job
```
BarcodeExtract /you/bam/file barcode_spnanning.bed
bgzip -c barcode_spanning.bed > barcode_spanning.bed.gz
tabix -p bed barcode_spanning.bed.gz
```
With PacBio SMRT:
```
extractHAIRS --pacbio 1 --bam /your/bam/file --VCF /your/vcf/file --out fragment_file --ref /reference/file
```
Similarly with Nanopore, change the --pacbio into --ont  

With a fragment file, you can sort it with following command, if your are phasing with 10X, Hi-C or the new-format is specified 
```
sort -n -k6 in.frag > sorted.frag
```  
If paired-ended NGS, PacBio SMRT or Oxford Nanopore is used with default format, use
```
sort -n -k3 in.frag > sorted.frag
```

#### Run SpecHap
The detailed usage can be found by
```
SpecHap --help
```

For instance, to phase PacBio SMRT reads:
```
SpecHap --vcf /your/gzvcf/file --frag /your/fragment/file --out /phased/vcf --pacbio
```

### Script for VCF handling.
You may find bunch of scripts that we use to benchmark the accuracy and completeness of assembled haplotype under folder ./scripts. 

### Author
SpecHap is developed by DeepOmics lab under the supervision of Dr. Li Shuaicheng, City University of Hong Kong, Hong Kong, China.

To contact us, send email to [yonghanyu2@cityu.edu.hk](yonghanyu2@cityu.edu.hk)

## Built With

* [htslib](https://github.com/samtools/htslib)
* [arapack-ng](https://github.com/opencollab/arpack-ng)
* [Eigen3](https://eigen.tuxfamily.org/dox/)
* [Lean Mean C++ Option Parser](https://eigen.tuxfamily.org/dox/)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


