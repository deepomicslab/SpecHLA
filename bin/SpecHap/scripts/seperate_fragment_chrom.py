import argparse

def parse_fragment(inf, prefix, chromo_snp_count:dict, chromos:list):
    line_new = ""
    chromo_idx = 0
    prev_snp_count = 0
    curr_chromo = chromos[0]
    outf = open(prefix + "_"+curr_chromo+".lst", "w")
    curr_snp_count = chromo_snp_count[curr_chromo]
    for line in inf:
        line_new = ""
        buff = line.split()
        nb = int(buff[0])
        name = buff[1]
        line_new = buff[0] + " " + name 
        for i in range(0, nb):
            idx = int(buff[2 + 2*i])
            if (idx > curr_snp_count):
                outf.close()
                chromo_idx += 1
                curr_chromo = chromos[chromo_idx]
                prev_snp_count = curr_snp_count
                curr_snp_count = chromo_snp_count[curr_chromo]
                outf = open(prefix + "_"+curr_chromo+".lst", "w")
            line_new = line_new + " " + str(idx-prev_snp_count) + " " + buff[2 + 2*i + 1]
        line_new = line_new + " " +  buff[-1] + "\n"
        outf.write(line_new)
    outf.close()

def read_snplist(inf, prefix):
    curr_chromo = ""
    snp_count = 0
    chromo_snp_count = dict()
    chromos = list()
    outf = None
    for line in inf:
        buff = line.split()
        chromo = buff[0]
        if chromo != curr_chromo:
            if outf is not None:
                outf.close()
                chromo_snp_count[curr_chromo] = snp_count 
                chromos.append(curr_chromo)  
            outf = open(prefix + "_" + chromo + ".allvar", "w")
            curr_chromo = chromo
        outf.write(line)
        snp_count += 1
    outf.close()
    chromos.append(curr_chromo)
    chromo_snp_count[curr_chromo] = snp_count 
    return chromo_snp_count, chromos

def main():
    parser = argparse.ArgumentParser('phaseset_to_vcf.py')
    parser.add_argument('-i', '--inf',  help='in fragment', required=True)
    parser.add_argument('-s', '--snplst', help='input snplist', required=True)
    parser.add_argument('-o', '--out', help='output prefix', required=True)
    options = parser.parse_args()

    inf = open(options.inf, "r")
    snplist = open(options.snplst, "r") 

    chromo_snp_count, chromos = read_snplist(snplist, options.out)
    parse_fragment(inf, options.out, chromo_snp_count, chromos)
    inf.close()
    snplist.close()
    return

if __name__ == '__main__':
    main()