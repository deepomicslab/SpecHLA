import vcf
import sys, getopt

def vcf_read_write(inn, out):
    inf = vcf.Reader(open(inn, 'r'))
    header = inf.formats
    
    
    return



def main(argv):
    inputfile = ''
    outputfile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:o:",["invcf=","outvcf="])
    except getopt.GetoptError:
        print ('SVaba_vcf_concordant.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--invcf"):
            inputfile = arg
        elif opt in ("-o", "--outvcf"):
            outputfile = arg
    vcf_read_write(inputfile, outputfile)
    
    return

if __name__ == '__main__':
    main(sys.argv[1:])