## Getting started

You can copy-paste the following command lines in your terminal from a Linux
machine (or replace `k8-linux` with `k8-darwin` on Mac):
```sh
# download the k8 executable
wget -O- http://sourceforge.net/projects/biobin/files/devtools/k8-0.2.1.tar.bz2/download \
	| bzip2 -dc | tar xf -
# download test VCF files
wget ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/vcf-flt/CHM1.mem.hc.flt.vcf.gz
wget ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/vcf-flt/NA12878.mem.hc.flt.vcf.gz
# acquire the evaluation script
git clone https://github.com/lh3/hapdip.git
# evaluate (it is fast, 0.5 minute)
./k8-linux hapdip/hapdip.js eval CHM1.mem.hc.flt.vcf.gz NA12878.mem.hc.flt.vcf.gz
```

The output is
```
hapdip  SNP     FP      39167
hapdip  SNP     TP      2095440
hapdip  INDEL   FP      46043
hapdip  INDEL   TP      460382
hapdip  INDEL   FP      21319   INDEL-1bp
hapdip  INDEL   TP      196706  INDEL-1bp
```

Or to evaluate excluding variants overlapping low-complexity regions (LCRs;
replace `-B` with `-b` to evaluate in LCRs only):
```sh
wget ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/LCR-hs37d5.bed.gz
./k8-linux hapdip/hapdip.js eval -B LCR-hs37d5.bed.gz \
	CHM1.mem.hc.flt.vcf.gz NA12878.mem.hc.flt.vcf.gz
```
The output is (note the much lower FP for INDELs)
```
hapdip  SNP     FP      29882
hapdip  SNP     TP      2028817
hapdip  INDEL   FP      5003
hapdip  INDEL   TP      189720
hapdip  INDEL   FP      1767    INDEL-1bp
hapdip  INDEL   TP      98315   INDEL-1bp
```

If you want to try a new variant caller, you need to download the reference
genome and the BWA-MEM alignments with
```sh
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
wget ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/CHM1.mem.bam
wget ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/NA12878.mem.bam
```
The alignments, 240GB in total, contain all the raw reads. The CHM1 alignment
has been processed with [MarkDuplicates][dedup].

## Methodology

HapDip is a benchmark suite for evaluating the accuracy of small variant
calling. The core evaluation command `hapdip.js eval` takes two [VCF files][vcf]
as input, one called from the CHM1 haploid genome and the other from the NA12878
diploid genome ([link to data][ftp]). If we assume CHM1 is truly haploid,
heterozygotes (hets) called from CHM1 can be considered as false positives (FP).
Because NA12878 has similar ancestry and is of similar coverage and quality to
CHM1, we would expect the variant caller to make similar number of errors on the
NA12878 data. Then the number of false NA12878 het calls should equal the number
of CHM1 het calls, and true positives (TP) of hets should equal the rest of
NA12878 hets. This way, we get an estimate of TP and FP without comparing
variant calls. For details and more discussions, please see [the paper][varcmp]
or [the preprint][arxiv].

This method not only evaluates variant callers, but also evaluates the reference
genome. A good genome should yield fewer het calls from CHM1 and more het calls
from NA12878.

## Usage

This repository contains a single script `hapdip.js` consisting of several
subcommands. It is modified from [the companion scripts][script] with the paper.
The script is written in a dialect of Javascript and requires the `k8`
javascript shell to run. The k8 executables for Linux and Mac are available at
[the same ftp site][ftp].

For evaluation, we only need to run the `eval` subcommand. It can be simply
invoked as
```sh
k8 hapdip.js eval CHM1-calls.vcf NA12878-calls.vcf
```
Other subcommands filter VCFs or have miscellaneous functionalities.
In particular, the filtered VCFs at the FTP site were all generated with:
```sh
# de-overlap, re-estimate GT and annotate
k8 hapdip.js deovlp raw.vcf | k8 hapdip.js upd1gt | k8 hapdip.js anno > anno.vcf
# apply preset hard filters
k8 hapdip.js filter -A anno.vcf > flt.vcf
```
Note that `anno` only works with specific versions of GATK, SAMtools, FreeBayes
and Platypus as it needs to extract information from caller-specific tags. If
it does not work for your VCF, please filter with your own program.

[varcmp]: http://bioinformatics.oxfordjournals.org/content/early/2014/07/03/bioinformatics.btu356.abstract
[vcf]: http://vcftools.sourceforge.net/specs.html
[ftp]: ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/
[arxiv]: http://arxiv.org/abs/1404.0929
[script]: https://github.com/lh3/varcmp/tree/master/scripts
[dedup]: http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates
[giab]: ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/
