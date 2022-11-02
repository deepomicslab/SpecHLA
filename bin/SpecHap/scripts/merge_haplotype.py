import argparse
import vcf
from record import Record, PhaseSet, ChromosomoHaplotype


def merge_chromosome_haplotype(chromo_info1: ChromosomoHaplotype, chromo_info2: ChromosomoHaplotype):
    chromo_info1.construct_connection_graph(chromo_info2)
    chromo_info1.merge_chromo_haplotype(chromo_info2)
    chromo_info1.finalize_chromosome_haplotype()


def write_chromosome(in_vcf: vcf.Reader, out_vcf: vcf.Writer, chromo_haplotype: ChromosomoHaplotype, contig: str):
    rec:vcf.model._Record
    for rec in in_vcf.fetch(contig):
        het = rec.samples[0].gt_type
        if het != 1:        # not het loci
            out_vcf.write_record(rec)
        else:
            record = chromo_haplotype.chromo_record[rec.POS]
            record.finalize_record(rec)
            out_vcf.write_record(rec)

        
def merge_haplotype(in_vcfs: list, out_vcf: vcf.Writer):
    contigs = in_vcfs[0].contigs.keys()
    for contig in contigs:
        chromo_infos = list()
        try:
            in_vcfs[0].fetch(contig)
        except:
            continue
        for in_vcf in in_vcfs:
            chromo_infos.append(ChromosomoHaplotype(in_vcf, str(contig)))
        for i in range(1, len(chromo_infos)):
            merge_chromosome_haplotype(chromo_infos[0], chromo_infos[i])
        write_chromosome(in_vcfs[0], out_vcf, chromo_infos[0], str(contig))



def main():
    parser = argparse.ArgumentParser("merge_haplotype.py")
    parser.add_argument('-v', '--vcf', nargs='+', help='VCF file(s) in gz format, indexed, seperated by space', required=True)
    parser.add_argument('-o', '--out', help='Haplotype merged vcf file', required=True)
    options = parser.parse_args()
    
    print (options.vcf)

    in_vcfs = list()
    for in_name in options.vcf:
        in_vcfs.append(vcf.Reader(filename=in_name))

    if len(in_vcfs) == 1:
        return

    out_vcf = vcf.Writer(open(options.out, 'w'), in_vcfs[0])
    merge_haplotype(in_vcfs, out_vcf)
    out_vcf.close()


if __name__ == '__main__':
    main()