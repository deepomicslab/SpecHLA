## Introduction

This project started as a set of [companion scripts][script] used in [my
review][varcmp] (or [preprint][arxiv]). After the publication, I polished
scripts and separated them out in this independent repo. With time, the project
grows with more functionalities particularly focusing on VCF processing. This
README gives an overview of key routines. [README-hapdip.md][hapdip-md]
explains the original hapdip benchmark in more details.

As most of my programs, this `hapdip.js` script consists of multiple
subcommands. It is very fast and seamlessly reads both plain and gzip'd text
files. In addition, for VCF processing, the script properly breaks MNPs and
complex events down to atomic SNPs and INDELs. It is aware of the CIGAR INFO
field produced by [FreeBayes][freebayes].

This script requires the k8 javascript shell to run, which can be found
[here][k8].

## Functionalities

### Filtering single-sample VCF files

VCF produced from whole-genome deep resequencing can be filtered with:
```sh
k8 hapdip.js deovlp raw.vcf.gz | k8 hapdip.js anno | gzip -1 > anno.vcf.gz
k8 hapdip.js filter anno.vcf.gz | bgzip > filtered.vcf.gz
```
Here `deovlp` chooses the best VCF record among overlapping records, `anno`
adds generic depth-related INFO fields and `filter` uses the fields to filter
variants with hard thresholds in [my review][varcmp]. In my experience,
for single-sample calling, these hard filters are better than GATK's VQSR (see
the [FermiKit manuscript][fermikit-manu]).

### Converting VCF to unary BED

```sh
k8 hapdip.js vcf2bed var.vcf.gz > var.bed
```
In the output, each alternate allele takes one BED line, which consists of chr,
start, end, allele length (0 for SNPs, positive for insertions and negative for
deletions), ref allele, alt allele, count of alt allele, total called
haplotypes and filters in VCF. Complex events are broken down to atomic SNPs
and INDELs.

### Subsetting VCF and reordering/replacing sample names

```sh
echo -e "old2\tnew1\nold4\tnew2" > list.txt
k8 hapdip.js vcfsub list.txt var.vcf > var-sub.vcf
```
The command chooses a subset of samples and reorder the samples according to
the input. If the input `list.txt` has a second column, the command replace
the old sample name in the first column to the new one on the second column.

### Evaluating calls from the CHM1-NA12878 pair

```sh
k8 hapdip.js eval CHM1.vcf.gz NA12878.vcf.gz
```
Please see [README-hapdip.md][hapdip-md] for details.

### Evaluating against a truth set

The `distEval` command evaluates a call set against a GIAB-like truth dataset,
which consists of one BED file giving high-confidence regions and one
VCF file giving variants.
```sh
# download GIAB-2.18 truth dataset (or from http://bit.ly/giab218)
wget ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/GIAB/N+P.auto.bed.gz
wget ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/GIAB/P.vcf.gz
# download an evaluation call set
wget ftp://hengli-data:lh3data@ftp.broadinstitute.org/hapdip/vcf-flt/NA12878.mem.hc.flt.vcf.gz
# evaluate
k8 hapdip.js distEval -d 10 -b N+P.auto.bed.gz P.vcf.gz NA12878.mem.hc.flt.vcf.gz
```
The output is:
```
distEval        SNP     N+P     2083729399
distEval        SNP     TP      2634728
distEval        SNP     FN      27555
distEval        SNP     FP      263
distEval        INDEL   N+P     2083729399
distEval        INDEL   TP      164730
distEval        INDEL   FN      2074
distEval        INDEL   FP      1762
```
where `N+P` gives the total length of the non-overlapping regions in
`N+P.auto.bed.gz`, `TP` is the number of true SNPs (INDELs) within *d*-bp from
a called variant (including both SNPs and INDELs), `FN` the number of true SNPs
(INDELs) *not* within *d*-bp from any called variants and `FP` the number of
called SNPs (INDELs) *not* within *d*-bp from any true variants. Distance-based
evaluation is robust to the multi-representation of variants.



[varcmp]: http://bioinformatics.oxfordjournals.org/content/early/2014/07/03/bioinformatics.btu356.abstract
[arxiv]: http://arxiv.org/abs/1404.0929
[script]: https://github.com/lh3/varcmp/tree/master/scripts
[hapdip-md]: https://github.com/lh3/hapdip/blob/master/README-hapdip.md
[fermikit-manu]: https://github.com/lh3/fermikit/tree/master/tex
[giab218]: http://bit.ly/giab218
[freebayes]: https://github.com/ekg/freebayes
[k8]: https://sourceforge.net/projects/biobin/files/devtools/
