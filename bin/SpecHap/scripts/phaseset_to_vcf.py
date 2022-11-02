import vcf
import argparse
from record import Record

def read_phased_record(inf):
    ps = 0
    update_ps = False
    phased_record = dict()
    for line in inf:
        if line[0] == 'B':
            update_ps = True
        elif line[0] == '*':
            continue
        else:
            buff = line.split()
            pos = int(buff[0])
            hap0 = buff[1]
            hap1 = buff[2]
            if hap0 == '-' or hap1 == '-':
                continue
            if update_ps:
                ps = pos
                update_ps = False
            phased_record[pos] = Record(pos, int(hap0), int(hap1), ps)
    return phased_record

def reset_phasing_info(record: vcf.model._Record):
    gt_str = record.samples[0]['GT']
    new_gt_str = gt_str[0] + '/' + gt_str[2]
    record.samples[0].gt_nums = new_gt_str
    if 'PS' in record.FORMAT.split(':'):
        record.samples[0].data = record.samples[0].data._replace(PS=0)

    record.samples[0].data = record.samples[0].data._replace(GT=new_gt_str)
 
def update_phasing_info_dchap(in_vcf: vcf.Reader, out_vcf: vcf.Writer, phased_record:dict):
    rec: vcf.model._Record
    count = 0
    for rec in in_vcf:
        count += 1
        pos = rec.POS
        if count in phased_record.keys(): #phased
            record = phased_record[count]
            record.finalize_record(rec)
        else:                           #not phased
            reset_phasing_info(rec)

        out_vcf.write_record(rec)

def update_phasing_info(in_vcf: vcf.Reader, out_vcf: vcf.Writer, phased_record:dict):
    rec: vcf.model._Record
    for rec in in_vcf:
        pos = rec.POS
        if pos in phased_record.keys(): #phased
            record = phased_record[pos]
            record.finalize_record(rec)
        else:                           #not phased
            reset_phasing_info(rec)

        out_vcf.write_record(rec)


def main():
    parser = argparse.ArgumentParser('phaseset_to_vcf.py')
    parser.add_argument('-i', '--inf',  help='in phaseset, refhap format', required=True)
    parser.add_argument('-v', '--vcf', help='input vcf', required=True)
    parser.add_argument('-o', '--out', help='output vcf', required=True)
    parser.add_argument('--dchap', help='dchap format', action="store_true")
    options = parser.parse_args()

    inf = open(file=options.inf, mode="r")
    in_vcf = vcf.Reader(filename=options.vcf)
    outf = vcf.Writer(open(options.out, "w"), template=in_vcf)
    

    phased_record = read_phased_record(inf)
    if options.dchap:
        update_phasing_info_dchap(in_vcf, outf, phased_record)
    else:
        update_phasing_info(in_vcf, outf, phased_record)

    inf.close()
    outf.close()

    return

if __name__ == '__main__':
    main()