import pysam
import np
import sys
def allele_imba(beta_set, snp_index_dict, snp_list, outdir):
    """
    Utilizing allelic imbalance information to phase
    get the linkage info from allele frequencies at each variant locus
    return the linkage matrix with a SpecHap acceptable format
    """
    f = open(outdir + '/fragment.imbalance.file', 'w')
    locus_num = len(beta_set)
    for i in range(locus_num - 1):
        j = i + 1
        same = max([ beta_set[i][0] *  beta_set[j][0], beta_set[i][1] *  beta_set[j][1] ])
        reverse = max([ beta_set[i][0] *  beta_set[j][1], beta_set[i][1] *  beta_set[j][0] ])

        edge_same = max(np.log(same/reverse), 0)
        edge_reverse = max(np.log(reverse/same), 0)
        first_locus = snp_index_dict[snp_list[i][1]] - 1 # 0 index
        second_locus = snp_index_dict[snp_list[j][1]] - 1
        print(second_locus, first_locus, edge_same, edge_reverse, edge_reverse, edge_same, file = f)
    f.close()


def parse_ad_info_tag(vcf_file_path, outdir):
    vcf_file = pysam.VariantFile(vcf_file_path)
    snp_index_dict = {}
    snp_list = []
    beta_set = []
    snp_idx = 0
    ad_values = []
    for variant in vcf_file:
        snp_list.append([variant.chrom, variant.pos, variant.id])
        snp_index_dict[variant.pos] = snp_idx
        snp_idx += 1
        for sample in variant.samples.values():
            if 'AD' in sample:
                ad_sum = sum(sample['AD'])
                beta_set.append([sample['AD'][0]/ad_sum, sample['AD'][1]]/ad_sum)

    allele_imba(beta_set, snp_index_dict, snp_list, outdir)

if __name__ == '__main__':
    # 1: input vcf, 2: output dir
    parse_ad_info_tag(sys.argv[1], sys.argv[2])
