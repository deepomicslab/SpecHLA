"""
python /mnt/d/HLAPro_backup/HLAPro/evaluation/visualize_trio/visualize.py --trio HG002:HG003:HG004 --workdir ./ --outfile test.file --outinf test.info
"""

from argparse import ArgumentParser
from argparse import ArgumentTypeError
import os
import sys

parser = ArgumentParser(description="preprocess for making figures to visualize trio.",prog='vis',usage="-h")
optional=parser._action_groups.pop()
required=parser.add_argument_group('required arguments')
flag_parser = parser.add_mutually_exclusive_group(required=False)
flag_data = parser.add_mutually_exclusive_group(required=False)
#necessary parameter
required.add_argument("--trio",help="The trio infromation; give sample names in the order of child:mother:father.\
 Example: NA12878:NA12891:NA12892. The order of mother and father can be ambiguous.",dest='trio',metavar='', type=str)
required.add_argument("--workdir",help="The directory consists of all samples' results.\
 The workdir running previous scipts.",dest='workdir',metavar='')
required.add_argument("--outfile",help="The outfile for making figure.",dest='outfile',metavar='')
required.add_argument("--outinf",help="output sampleinfo.",dest='outinf',metavar='')
parser._action_groups.append(optional)
args = parser.parse_args()

def main():
    sample_list = args.trio.split(":")
    vcf_list = args.workdir + "/trio.vcf.list"
    i = 0
    for sample in sample_list:
        if i == 0:
            command = "ls %s/%s/HLA_*.rephase.vcf.gz >%s"%(args.workdir, sample, vcf_list)
        else:
            command = "ls %s/%s/HLA_*.rephase.vcf.gz >>%s"%(args.workdir, sample, vcf_list)
        os.system(command)
        i += 1
    command = "perl %s/convert.tro.pl %s %s %s"%(sys.path[0], vcf_list, args.outfile, args.outinf)
    # print (command)
    os.system(command)
    os.system("rm %s"%(vcf_list))

if __name__ == "__main__":
    main()