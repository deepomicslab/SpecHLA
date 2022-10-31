## Introduction

HTSbox is a fork of early [HTSlib][htslib]. It is a collection of small
*experimental* tools manipulating HTS-related files. While some of these tools
are already part of the [official SAMtools][samtools] package, others are for
niche use cases in my own work, so are maintained by myself. Please keep in mind
that HTSbox is **NOT** the official repository. It lags far behind HTSlib in
terms of features, activity, clarity and robustness. If you are looking for a
high-quality library and related tools, the [SAMtools organization
repositories][github] are the right place.

## Usage by examples

1. Summary pileup:
   ```sh
   htsbox pileup -f ref.fa sorted1.bam sorted2.bam
   ```
   This gives the number of observed alleles, including InDels, in each input
   BAM. Additional options can be applied for filtering:
   ```sh
   htsbox pileup -f ref.fa -Q20 -q30 -l90 sorted.bam
   ```
   This filters alignments shorter than 90bp and with mapping quality lower than
   20 and filters bases with quality lower than 20.

2. Naive variant calling:
   ```sh
   htsbox pileup -f ref.fa -Q20 -q30 -cs3 sorted.bam
   ```
   The output in VCF gives positions with ALT alleles appearing 3 or more times.
   Note that variant calling in this way is very crude. It is usually not
   recommended to use this for whole-genome and exon variant calling from short
   reads. Nonetheless, `pileup' is a proper tool to call variants from contigs,
   in particular unitigs produced by fermi:
   ```sh
   htsbox pileup -cuf ref.fa unitig.bam
   ```
   Option `-u` enables multiple settings.

3. Generate consensus FASTA:
   ```sh
   htsbox pileup -f ref.fa -Q20 -q30 -Fs3 sorted.bam
   ```

4. Pairwise alignment summary:
   ```sh
   htsbox samview -p aln.bam
   ```
   In the output, each line gives QName, QLen, QStart, QEnd, Strand, RName,
   RLen, RStart, REnd, PerBaseDivergence, MapQ and semicolon-delimited misc
   information.

5. Count alignment break points (mainly used to evaluate misassemblies):
   ```sh
   htsbox abreak name-srt.sam.gz
   ```
   or call structural variantions:
   ```sh
   htsbox abreak -u name-srt.sam.gz
   ```

6. Reduce quality resolution with Illumina binning:
   ```sh
   htsbox qualbin -t2 -bm7 in.bam
   ```
   This is the only command so far that explores multi-threading.

7. Evaluate empirical base quality, similar to MAQ's mapcheck:
   ```sh
   htsbox mapchk sorted.bam ref.fa
   ```
   The output is a bit complicated. Let me know if you are interested (I will
   see if it is worth documenting the output).

8. Generate per-base read depth:
   ```sh
   htsbox depth -p0 sorted.bam
   ```
   The last two columns give the total read depth and depth of reads mapped with
   high mapQ (threshold defaults to 20). When `-p` takes a value between 0 and
   1, it will use the [CODOC][codoc] strategy to find windows with relatively
   stable read depth. The output can be much smaller as it loses some
   information.

[htslib]: https://github.com/samtools/htslib
[samtools]: https://github.com/samtools/samtools
[github]: https://github.com/samtools/
[codoc]: http://www.ncbi.nlm.nih.gov/pubmed/24872424
