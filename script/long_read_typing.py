"""
Function 1 : assign long reads to the gene
Function 2 : typing with only long reads
"""

import sys
import os
import numpy
import pysam
import gzip
import argparse

interval_dict={}
gene_list=[]

class Read_Obj():
    def __init__(self, read):
        mis_NM = 0
        for ta in read.get_tags():
            if ta[0] == 'NM':
                mis_NM = ta[1]  

        self.read_length = 0
        self.match_num = 0  
        for ci in read.cigar:
            self.read_length += ci[1]
            if ci[0] == 0:
                self.match_num += ci[1]
        #    elif ci[0] == 4 and ci[1] > 50:
        #        self.read_length -= ci[1]
        #    elif ci[0] == 5:
        #        self.read_length -= ci[1]
        # self.read_length = len(read.query_sequence) 
        #print (read.query_name, read.reference_name, self.read_length, self.match_num, mis_NM)
        self.read_name = read.query_name
        self.allele_name = read.reference_name
        self.mismatch_num = mis_NM
        #self.mismatch_rate = round(float(mis_NM)/self.read_length, 6)
        self.mismatch_rate = round(float(mis_NM)/self.match_num, 6)
        self.match_rate = round(self.match_num / self.read_length, 6)
        if self.match_num < 650 and self.allele_name[0:3] == "KIR":
            self.match_rate = 0
        #elif self.allele_name[0:6] == "CYP8A1":
        #    print(mis_NM)
        #elif self.allele_name[0:6] == "CYP4B1":
        #    print(mis_NM)
        elif self.match_num < 800 and self.allele_name[0:3] == "CYP":
            self.match_rate = 0
        elif self.match_num < 600:
            self.match_rate = 0
        #self.match_rate = round(float(self.match_num - 2 * mis_NM)/self.read_length, 6)
        #self.match_rate = round(self.match_num /self.read_length , 6)
        self.loci_name = self.allele_name.split("*")[0]
        if self.loci_name == "KIR2DL5A" or self.loci_name == "KIR2DL5B":
            self.loci_name = "KIR2DL5"

class Score_Obj():
    def __init__(self):
        self.loci_score = {}
        self.loci_mismatch_score = {}
        self.reads_len_dict = {}
        self.read_loci = {}
    
    def add_read(self, read_obj):
        self.reads_len_dict[read_obj.read_name] = read_obj.read_length
        if gene_class == "HLA":
            score = round(numpy.exp(read_obj.match_rate) * (1 - read_obj.mismatch_rate), 6)
        elif gene_class == "KIR":
            score = round(numpy.exp(read_obj.match_rate) * (1 - read_obj.mismatch_rate), 6)
        else:
            score = round(1 - read_obj.mismatch_rate, 6)
        
        if read_obj.read_name not in self.loci_score:
            self.loci_score[read_obj.read_name] = {}
            self.loci_score[read_obj.read_name][read_obj.loci_name] = score
        elif read_obj.loci_name not in self.loci_score[read_obj.read_name]:
            self.loci_score[read_obj.read_name][read_obj.loci_name] = score
        else:
            if score > self.loci_score[read_obj.read_name][read_obj.loci_name]:
                self.loci_score[read_obj.read_name][read_obj.loci_name] = score
    
    def assign(self, assign_file):
        f = open(assign_file, 'w')
        # print (len(self.loci_score))
        for read_name in self.loci_score:
            gene_score = sorted(self.loci_score[read_name].items(), key=lambda item: item[1], reverse = True)
            if gene_score[0][1] < Min_score:
                print (read_name, 'MinScore', gene_score)
                continue
            if len(gene_score) == 1:
                assigned_locus = gene_score[0][0]
                print (read_name, 'OK', gene_score)
            else:
                if gene_score[0][1] - gene_score[1][1] >= Min_diff:
                    assigned_locus = gene_score[0][0]
                    print (read_name, 'OK', gene_score)
                else:
                    print (read_name, 'MinDiff', gene_score)
                    continue
            print (read_name, assigned_locus, file = f)
            self.read_loci[read_name] = assigned_locus
        f.close()
        return self.read_loci


class Pacbio_Binning():

    def __init__(self):
        self.infq = f"{parameter.raw_fq}"
        if gene_class == "HLA":
            self.db = f"{sys.path[0]}/../db/HLA/ref/HLA.select.fasta"
        if gene_class == "KIR":
            self.infq = f"{parameter.outdir}/{parameter.sample}.split.fastq" 
            self.db = f"{sys.path[0]}/../db/KIR/ref/KIR.extend.select.fasta"
            #self.db = f"{sys.path[0]}/../db/KIR/ref/KIR.extend.fasta"
        if gene_class == "CYP":
            self.db = f"{sys.path[0]}/../db/CYP/ref/CYP.select.fasta"
        self.map2db()
        self.sam = f"{parameter.outdir}/{parameter.sample}.db.sam"
        self.bamfile = pysam.AlignmentFile(self.sam, 'r')   
        self.assign_file = f"{parameter.outdir}/{parameter.sample}.assign.txt"

    def map2db(self):
        if gene_class == "KIR":
            align_split = f"""
            ref={sys.path[0]}/../db/KIR/ref/KIR.ref.fasta
            {sys.path[0]}/../bin/minimap2 -t {parameter.threads} -p 0.5 -H -O 50,60 -N 100000 -a $ref {parameter.raw_fq} > {parameter.outdir}/{parameter.sample}.split.sam
            perl {sys.path[0]}/split.read.pl {parameter.raw_fq} {parameter.outdir}/{parameter.sample}.split.sam {parameter.outdir}/{parameter.sample}.split.fastq
            echo split fastq done
            """
            os.system(align_split)
        
        alignDB_order = f"""
        ref={self.db}
        fq={self.infq}
        outdir={parameter.outdir}
        bin={sys.path[0]}/../bin
        sample={parameter.sample}
        $bin/minimap2 -t {parameter.threads} -p 0.5 -H -O 50,60 -N 100000 -a $ref $fq > $outdir/$sample.db.sam
        echo alignment done.
        """
        #print(alignDB_order)
        os.system(alignDB_order)

    def read_bam(self):
        scor = Score_Obj()
        for read in self.bamfile:
            if read.is_unmapped:
                continue
            # print (read)
            read_obj = Read_Obj(read)
            genename = read_obj.allele_name.split('*')[0]
            if read_obj.match_rate == 0:
                continue
            if gene_class == "KIR" and read_obj.mismatch_rate > 0.02:
                continue
            if gene_class == "HLA" and read_obj.mismatch_rate > 0.03:
                continue
            #if genename == "CYP4F2" and read_obj.mismatch_rate >0.01:
            #    continue
            if gene_class == "CYP" and read_obj.mismatch_rate > 0.02 and genename != "CYP3A5":
                continue
            if read_obj.allele_name.split('*')[0] == "CYP3A5" and read_obj.mismatch_rate > 0.1:
                continue

            print (read_obj.read_name, read_obj.match_rate,read_obj.mismatch_rate, read_obj.allele_name, read_obj.match_num, read_obj.mismatch_num, read_obj.read_length)
            scor.add_read(read_obj)
        read_loci = scor.assign(self.assign_file)
        for gene in gene_list:
            self.filter_fq(gene, read_loci)
        print ("reads-binning done.")

    def filter_fq(self, gene, dict):
        i = 0
        self.infq = f"{parameter.raw_fq}"
        if gene_class == "KIR":
            self.infq = f"{parameter.outdir}/{parameter.sample}.split.fastq"
        outfile = parameter.outdir + '/%s.%s.fq'%(gene, args["a"])
        out = open(outfile, 'w')
        flag = False
        if self.infq.split(".")[-1] == "gz":
            f = gzip.open(self.infq,'rt')
        else:
            f = open(self.infq)
        for line in f:
            line = line.strip()
            if i % 4 == 0:
                read_name = line.split()[0][1:]
                if read_name in dict.keys() and dict[read_name] == gene:
                    flag = True
                    num = 1
                    print (line, file = out)
            elif flag:
                print (line, file = out)
                num += 1
                if num == 4:
                    flag = False
            i += 1
        f.close()
        out.close()
        os.system('gzip -f %s'%(outfile))

class Parameters():

    def __init__(self):

        self.sample = args["n"]
        self.raw_fq = args["r"]
        outdir = args["o"]
        self.population = args["p"]
        self.threads = args["j"]
        self.db = "%s/../db/"%(sys.path[0])
        self.bin = "%s/../bin/"%(sys.path[0])      
        self.outdir = "%s/%s/"%(outdir, self.sample)
        self.whole_dir = "%s/"%(sys.path[0])

        if not os.path.exists(args["o"]):
            os.system("mkdir %s"%(args["o"]))
        if not os.path.exists(self.outdir):
            os.system("mkdir %s"%(self.outdir))


class Fasta():

    def vcf2fasta(self, gene):
        for index in range(2):
            #parameters="-p 0.5 -B 10 -ax asm20"
            #if gene == "CYP8A1" or gene == "CYP4B1":
            #    parameters="-a"
            
            order = """
            sample=%s
            bin=%s
            db=%s
            outdir=%s
            gene=%s
            i=%s
            j=%s
            gene_ref=$db/split_ref/$gene.fasta
            $bin/minimap2 -t %s -p 0.5 -B 10 -ax asm20 $gene_ref $outdir/$gene.%s.fq.gz | $bin/samtools view -bS -F 0x800 -| $bin/samtools sort - >$outdir/$gene.bam
            $bin/samtools index $outdir/$gene.bam

            longshot -F -S -c 2 -q 10 -Q 2 -a 5 -y 20 -e 2 -E 0.1 --hom_snv_rate 0.01 --het_snv_rate 0.01 --bam $outdir/$gene.bam --ref $gene_ref --out $outdir/$sample.$gene.longshot.vcf 
            bgzip -f $outdir/$sample.$gene.longshot.vcf
            tabix -f $outdir/$sample.$gene.longshot.vcf.gz

            if [ %s == 1 ];then 
                $bin/bcftools filter -R $db/$gene_class/$gene_class.ref.exon.extend.txt $outdir/$sample.$gene.longshot.vcf.gz >$outdir/$sample.$gene.phased.vcf
            else 
                zcat $outdir/$sample.$gene.longshot.vcf.gz >$outdir/$sample.$gene.phased.vcf
            fi            
            bgzip -f $outdir/$sample.$gene.phased.vcf
            tabix -f $outdir/$sample.$gene.phased.vcf.gz

            $bin/samtools faidx $gene_ref %s |$bin/bcftools consensus -H $i $outdir/$sample.$gene.phased.vcf.gz >$outdir/allele.$i.$gene.raw.fasta
            echo ">%s_$j" >$outdir/allele.$i.$gene.fasta
            cat $outdir/allele.$i.$gene.raw.fasta|grep -v ">" >>$outdir/allele.$i.$gene.fasta
            
            $bin/samtools faidx $outdir/allele.$i.$gene.fasta        
            """%(parameter.sample, parameter.bin, parameter.db, parameter.outdir, gene, index+1, index,  parameter.threads, args["a"], args["u"], interval_dict[gene], gene)
            os.system(order)
            # -S -A -Q 10 -E 0.3 -e 5


    def get_fasta(self):
        for gene in gene_list:
            self.vcf2fasta(gene)
        # self.annotation()

    def annotation(self):
        anno = """
        perl %s/annoHLA.pl -s %s -i %s -p %s -r tgs -g %s -c %s
        cat %s/allele.phase.result.txt
        """%(parameter.whole_dir, parameter.sample, parameter.outdir, parameter.population, args["g"], args["i"], parameter.outdir)
        # print (anno)
        os.system(anno)



    # def phase(self):
    #     hairs = f"$bin/ExtractHAIRs --triallelic 1 --pacbio 1 --indels 1 --ref $ref\
    #      --bam $outdir/$sample.tgs.sort.bam --VCF %s --out $outdir/fragment.tgs.file"
    # order = '%s/../bin/SpecHap -P --window_size 15000 --vcf %s --frag %s/fragment.sorted.file\
    #  --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene)

if __name__ == "__main__":   

    parser = argparse.ArgumentParser(description="HLA Typing with long-reads", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-r", type=str, help="fastq file for the long-reads sample", metavar="\b")
    required.add_argument("-n", type=str, help="Sample ID", metavar="\b")
    required.add_argument("-o", type=str, help="outdir", metavar="\b", default="./output")
    required.add_argument("-i", type=str, help="HLA,KIR,CYP",metavar="\b", default="HLA")
    optional.add_argument("-p", type=str, help="Population information", metavar="\b", default="Unknown")
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=5)
    optional.add_argument("-d", type=float, help="Minimum score difference to assign a read to a gene.", metavar="\b", default=0.001)
    optional.add_argument("-m", type=int, help="1 represents typing, 0 means only read assignment", metavar="\b", default=1)
    optional.add_argument("-a", type=str, help="prefix of filtered fastq file.", metavar="\b", default="long_read")
    optional.add_argument("-g", type=int, help="Whether use G-translate in annotation [1|0], default is 0.", metavar="\b", default=1)
    optional.add_argument("-u", type=str, help="Choose full-length or exon typing. 0 indicates full-length, 1 means exon.", metavar="\b", default="0")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    parameter = Parameters()
    # Min_score = 0.1  #the read is too long, so the score can be very low.
    Min_score = 0  #the read is too long, so the score can be very low.
    Min_diff = args["d"]  #0.001
    gene_class = args["i"]
    if gene_class == "HLA":
        gene_list = [ 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-K', 'HLA-L', 'HLA-P', 'HLA-V', 'HLA-DQA2', 'HLA-DPA2', 'HLA-N', 'HLA-S', 'HLA-T', 'HLA-U', 'HLA-W', 'MICA', 'MICB', 'TAP1', 'TAP2', 'HFE' ]
        interval_dict = {"HFE": "HFE:301-8261","HLA-A": "HLA-A:301-3802","HLA-B": "HLA-B:301-4381","HLA-C": "HLA-C:301-4618","HLA-DMA": "HLA-DMA:301-5310","HLA-DMB": "HLA-DMB:301-7040","HLA-DOA": "HLA-DOA:301-3953","HLA-DOB": "HLA-DOB:301-5086","HLA-DPA1": "HLA-DPA1:301-10075","HLA-DPA2": "HLA-DPA2:301-7043","HLA-DPB1": "HLA-DPB1:301-11826","HLA-DPB2": "HLA-DPB2:301-18134","HLA-DQA1": "HLA-DQA1:301-6784","HLA-DQA2": "HLA-DQA2:301-6152","HLA-DQB1": "HLA-DQB1:301-7402","HLA-DRA": "HLA-DRA:301-6005","HLA-DRB1": "HLA-DRB1:301-11380","HLA-DRB3": "HLA-DRB3:301-13888","HLA-DRB4": "HLA-DRB4:301-15764","HLA-DRB5": "HLA-DRB5:301-13745","HLA-E": "HLA-E:301-4122","HLA-F": "HLA-F:301-3848","HLA-G": "HLA-G:301-3438","HLA-H": "HLA-H:301-3810","HLA-J": "HLA-J:301-3844","HLA-K": "HLA-K:301-3852","HLA-L": "HLA-L:301-4070","HLA-N": "HLA-N:301-935","HLA-P": "HLA-P:301-3231","HLA-S": "HLA-S:301-1174","HLA-T": "HLA-T:301-2787","HLA-U": "HLA-U:301-1030","HLA-V": "HLA-V:301-2203","HLA-W": "HLA-W:301-3272","MICA": "MICA:301-13027","MICB": "MICB:301-12616","TAP1": "TAP1:301-9570","TAP2": "TAP2:301-10907"}
    if gene_class == "CYP":
        gene_list = [ 'CYP19A1', 'CYP1A1', 'CYP1B1', 'CYP26A1', 'CYP2A13', 'CYP2A6', 'CYP2B6', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP2D6', 'CYP2F1', 'CYP2J2', 'CYP2R1', 'CYP2S1', 'CYP2W1', 'CYP3A4', 'CYP3A43', 'CYP4A22', 'CYP4B1', 'CYP4F2', 'CYP8A1', 'CYP3A5', 'CYP3A7' ]
        interval_dict = {"CYP4B1": "CYP4B1:1-20353","CYP3A4": "CYP3A4:1-34205","CYP3A7": "CYP3A7:1-37162","CYP2F1": "CYP2F1:1-20929","CYP2A13": "CYP2A13:1-14732","CYP2C9": "CYP2C9:1-58934","CYP26A1": "CYP26A1:1-11410","CYP3A5": "CYP3A5:1-38805","CYP1A1": "CYP1A1:1-7878","CYP3A43": "CYP3A43:1-45538","CYP2A6": "CYP2A6:1-13910","CYP19A1": "CYP19A1:1-47775","CYP2S1": "CYP2S1:1-21330","CYP1B1": "CYP1B1:1-12177","CYP2B6": "CYP2B6:1-34098","CYP8A1": "CYP8A1:1-64264","CYP2C19": "CYP2C19:1-99871","CYP4F2": "CYP4F2:1-27051","CYP2D6": "CYP2D6:1-11312","CYP4A22": "CYP4A22:1-12568","CYP2R1": "CYP2R1:1-21197","CYP2C8": "CYP2C8:1-39726","CYP2W1": "CYP2W1:1-13442","CYP2J2": "CYP2J2:1-40444"}
    if gene_class == "KIR":
        gene_list = ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5', 'KIR2DP1', 'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3', 'KIR3DP1', 'KIR3DS1']
        interval_dict = {"KIR2DL1": "KIR2DL1:301-15041","KIR2DL2": "KIR2DL2:301-15082","KIR2DL3": "KIR2DL3:301-15068","KIR2DL4": "KIR2DL4:301-11434","KIR2DL5": "KIR2DL5:301-10072","KIR2DP1": "KIR2DP1:301-13426","KIR2DS1": "KIR2DS1:301-15020","KIR2DS2": "KIR2DS2:301-14878","KIR2DS3": "KIR2DS3:301-15403","KIR2DS4": "KIR2DS4:301-16392","KIR2DS5": "KIR2DS5:301-15548","KIR3DL1": "KIR3DL1:301-14846","KIR3DL2": "KIR3DL2:301-17301","KIR3DL3": "KIR3DL3:301-12699","KIR3DP1": "KIR3DP1:301-4540","KIR3DS1": "KIR3DS1:301-15232"}


    ###assign reads
    if args["m"] == 10086:
        print ("skip assignment, just for testing")
    else:
        pbin = Pacbio_Binning()
        pbin.read_bam()        

    if args["m"] != 0:
        fa = Fasta()
        fa.get_fasta()
        fa.annotation()





