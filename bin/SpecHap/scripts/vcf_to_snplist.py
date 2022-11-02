import vcf
import argparse

def main():
    parser = argparse.ArgumentParser(prog='vcf_to_snplist.py')
    parser.add_argument('-v', '--vcf', help='VCF file', required=True)
    parser.add_argument('-o', '--out', help='SNP list', required=True)
    options = parser.parse_args()
    
    inf = vcf.Reader(filename=options.vcf)
    out = open(file=options.out, mode='w')

    record : vcf.model._Record
    for record in inf:
        chromo = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = record.ALT
        stype = "SNP"
        out.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (chromo, pos, pos, ref, alt[0], stype))

    out.close()   
    return


if __name__ == '__main__':
    main()