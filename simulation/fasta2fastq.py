from Bio import SeqIO
import sys


for r in SeqIO.parse(sys.argv[1], "fasta"):
    r.letter_annotations["solexa_quality"] = [40] * len(r)
    print(r.format("fastq"), end='')
