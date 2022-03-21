
def get_allele_to_hits(fn):
    allele_to_hits = {}
    with open(fn)  as f:
        for row in f:
            if row[0] in ["#","@"]:
                continue
            row = row.strip().split("\t")
            HLA_allele = row[0].split("_")[1].split("*")[0]

            data = [c.split(":") for c in row[13:]]
            data = dict( ((c[0],c[2]) for c in data) )

            allele_to_hits.setdefault(HLA_allele, [])
            allele_to_hits[HLA_allele].append( (float(data["de"]), row[0], data["cg"], row[5], int(row[7]), row[4] ) ) 
    return allele_to_hits

Alleles = {}
finalAlleles1 = {}
finalAlleles2 = {}
with open("HLA-ALN-H1.tsv","w") as f:
    for allele, hits in get_allele_to_hits("hla_gen_H1.paf").items():
        Alleles[allele] = 1
        hits.sort(key=lambda x: x[0])
        for h in hits[:5]:
            print(allele, *[str(c) for c in h], sep="\t", file=f)
            if allele in finalAlleles1:
                #print(allele)
            else:
                finalAlleles1[allele] = h[1]
                #print(allele, h[1])


with open("HLA-ALN-H2.tsv","w") as f:
    for allele, hits in get_allele_to_hits("hla_gen_H2.paf").items():
        Alleles[allele] = 1
        hits.sort(key=lambda x: x[0])
        for h in hits[:5]:
            print(allele, *[str(c) for c in h], sep="\t", file=f)
            if allele in finalAlleles2:
                #print(allele)
            else:
                finalAlleles2[allele] = h[1]
                #print(allele, h[1])

ff = open("HLAtyping.result.txt", "w")
for allele in Alleles:
    A1 = "NA"
    A2 = "NA"
    if allele in finalAlleles1:
        A1 = finalAlleles1[allele].split("_")[1]
    if allele in finalAlleles2:
        A2 = finalAlleles2[allele].split("_")[1]
    out = "\t".join((allele, A1, A2))
    ff.write(out + "\n")
ff.close()

