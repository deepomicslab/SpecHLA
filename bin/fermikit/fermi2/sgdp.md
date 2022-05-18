This page describes how to assemble unitigs from high-coverage human samples
sequenced with Illumina 100bp reads, and to construct the FM-index for a
collection of human samples. The following procedure is used to assemble and
index ~300 samples from [Simons Genome Diversity Project][sgdp] (SGDP).

1. Download and unpack precompiled binaries:

        wget -O- http://sf.net/projects/biobin/files/unclassified/sgdp-1.1_x64-linux.tar.bz2/download \
            | bzip2 -dc | tar xvf -

   This will create a directory `sgdp-1.0_x64-linux`. For simplicity, we assume
   you have put the directory on the `$PATH`.

2. Generate Makefiles:

        fermi2.pl unitig -t8 -p foo "trimadap foo.fq.gz" > foo.mak

   This requests 8 threads at the maximum and to write the results with prefix
   `foo`. `trimadap` trims typical adapters used in Illumina resequencing. It
   is part of [seqtk][seqtk]. If the input is BAM, you may replace with:

        fermi2.pl unitig -t8 -p foo "htscmd bam2fq foo.bam | trimadap" > foo.mak

   If you are running this pipeline the first time, it is recommended to have a
   look at `foo.mak`. It is simple.

3. Assemble unitigs:

        make -f foo.mak

   For 100bp human data at ~30X coverage, it usually take 2.5-3.5 days in real
   time with peak memory around 85GB. File `foo.mag.gz` contains the final
   unitigs. You can delete the rest of files, though it is usually a good idea
   to keep `*.log` files for timing and trouble shooting.

4. Append `*.mag.gz` in the working directory to an existing index `bar.fmr`:

        fermi2.pl mag2fmr -i bar.fmr *.mag.gz > all.mak
        make -f all.mak
        k8 sgdp-1.0_x64-linux/fermi2.js log2tbl *.fmr.log > all.tbl

   This appends each set of unitigs to the previous FM-index. It takes about 30
   minutes to add one sample. The last `.fmr` file is the final FM-index you
   should keep. If you do not have an existing index, simply drop the `-i
   bar.fmr` option. `all.tbl` keeps the sample information.

5. Variant calling:

        bwa index ref.fa; samtools faidx ref.fa
        bwa mem -B9 -O16 -L10 -t8 ref.fa foo.mag.gz | gzip -1 > foo.sam.gz
		samtools view -uS foo.sam.gz | samtools sort - foo
		htsbox pileup -cq40 -Q5 -f ref.fa foo.bam > foo.vcf

   where `htsbox` is available [via github][htsbox].

[sgdp]: http://www.simonsfoundation.org/life-sciences/simons-genome-diversity-project/
[seqtk]: https://github.com/lh3/seqtk
[htsbox]: https://github.com/lh3/htsbox
