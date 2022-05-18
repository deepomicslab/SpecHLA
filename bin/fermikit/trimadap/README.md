Trimadap is a small tool to trim adapter sequences from Illumina data. It
performs SSE2-SW between each read and each adapter sequence and identifies
adapter sequences with a few heuristic rules which can be found in the
`ta_trim1()` function in `trimadap-mt.c`. The default adapters it uses are
included in `illumina.txt`. These are typical Illumina adapters from
paired-end sequencing.

Trimadap is designed as an on-the-fly stream filter. It is very fast. In the
multi-threading mode, it is as fast as reading through a gzip-compressed FASTQ
file. On the other hand, trimadap is very conservative. It is not good in terms
of accuracy as of now. I will probably fine tune the heuristic rules in
future. This should not be hard in principle, but it takes development time.
