import sys
import os
import pysam
import gzip

Min_score = 0.1  #the read is too long, so the score can be very low.
Min_diff = 0.001

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
        # print (read.query_name, read.reference_name, read_length, mis_NM)   
        # self.read_length = len(read.query_sequence) 

        self.read_name = read.query_name
        self.allele_name = read.reference_name
        self.mismatch_num = mis_NM
        self.mismatch_rate = round(float(mis_NM)/self.read_length, 6)
        self.match_rate = round(float(self.match_num)/self.read_length, 6)
        self.loci_name = self.allele_name.split("*")[0]

class Score_Obj():
    def __init__(self):
        self.loci_score = {}
        self.loci_mismatch_score = {}
        self.reads_len_dict = {}
        self.read_loci = {}
    
    def add_read(self, read_obj):
        self.reads_len_dict[read_obj.read_name] = read_obj.read_length
        score = round(read_obj.match_rate * (1 - read_obj.mismatch_rate), 6)
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
        for read_name in self.loci_score:
            gene_score = sorted(self.loci_score[read_name].items(), key=lambda item: item[1], reverse = True)
            # if gene_score[0][0] == "B" or gene_score[0][0] == "C":
            #     print (read_name, gene_score)
            if gene_score[0][1] < Min_score:
                continue
            if len(gene_score) == 1:
                assigned_locus = gene_score[0][0]
            else:
                if gene_score[0][1] - gene_score[1][1] >= Min_diff:
                    assigned_locus = gene_score[0][0]
                else:
                    continue
            print (read_name, assigned_locus, file = f)
            self.read_loci[read_name] = assigned_locus
        f.close()
        return self.read_loci


class Pacbio_Binning():

    def __init__(self, raw_fq, outdir, ID):
        
         
        self.fastq = raw_fq
        self.outdir = outdir
        self.ID = ID
        self.db = f"{sys.path[0]}/../db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"
   
        self.map2db()
        self.sam = f"{self.outdir}/{self.ID}.db.sam"
        self.bamfile = pysam.AlignmentFile(self.sam, 'r')   
        self.assign_file = f"{self.outdir}/{self.ID}.assign.txt"


    def map2db(self):
        alignDB_order = f"""
        fq={self.fastq}
        ref={self.db}
        outdir={self.outdir}
        bin={sys.path[0]}/../bin
        sample={self.ID}
        $bin/minimap2 -t 12 -p 0.1 -N 100000 -a $ref $fq > $outdir/$sample.db.sam
        echo alignment done.
        """
        os.system(alignDB_order)

    def read_bam(self):
        scor = Score_Obj()
        for read in self.bamfile:
            if read.is_unmapped:
                continue
            # print (read)
            read_obj = Read_Obj(read)
            scor.add_read(read_obj)
            # print (read_obj.read_name, read_obj.mismatch_rate, read_obj.allele_name )
        read_loci = scor.assign(self.assign_file)
        for gene in ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']:
            self.filter_fq(gene, read_loci)
        print ("reads-binning done.")

    def filter_fq(self, gene, dict):
        i = 0
        #gene = 'A'
        outfile = self.outdir + '/%s.fq'%(gene)
        out = open(outfile, 'w')
        flag = False
        if self.fastq.split(".")[-1] == "gz":
            f = gzip.open(self.fastq,'rt')
        else:
            f = open(self.fastq)
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



class PacBio():

    def __init__(self, gene, outdir, ID, db):
        self.fq = self.outdir + '/%s.fq.gz'%(gene)
        self.outdir = outdir
        self.ID = ID
        self.ref = f"{sys.path[0]}/../db/HLA/HLA_{gene}/HLA_{gene}.fa"
        self.db = db

    def alignment(self):
        alignment_order = f"""
        fq={self.fq}
        ref={self.ref}
        outdir={self.outdir}
        bin={sys.path[0]}/../bin
        sample={self.ID}
        $bin/minimap2 -a $ref $fq > $outdir/$sample.tgs.sam
        $bin/samtools view -F 2308 -b -T $ref $outdir/$sample.tgs.sam > $outdir/$sample.tgs.bam
        $bin/samtools sort $outdir/$sample.tgs.bam -o $outdir/$sample.tgs.sort.bam
        $bin/samtools index $outdir/$sample.tgs.sort.bam
        echo alignment done.
        """
        return alignment_order

    def run(self):
        # alignment_order = self.alignment()
        # os.system(alignment_order)
        alignDB_order = self.map2db()
        os.system(alignDB_order)

    # def smallVariants(self):


    # def phase(self):
    #     hairs = f"$bin/ExtractHAIRs --triallelic 1 --pacbio 1 --indels 1 --ref $ref\
    #      --bam $outdir/$sample.tgs.sort.bam --VCF %s --out $outdir/fragment.tgs.file"
    # order = '%s/../bin/SpecHap -P --window_size 15000 --vcf %s --frag %s/fragment.sorted.file\
    #  --out %s/%s.specHap.phased.vcf'%(sys.path[0],my_new_vcf, outdir, outdir,gene)

if __name__ == "__main__":   
# fq = "/mnt/d/HLAPro_backup/insert/high_pacbio/child_1/child_1.fastq"
    # ID = "child_1"
    # raw_fq = "/mnt/d/HLAPro_backup/insert/single_pacbio/fq/child_1_A.fastq"
    # outdir = "/mnt/d/HLAPro_backup/insert/single_pacbio"
    ID = sys.argv[1]
    raw_fq = sys.argv[2]
    outdir = sys.argv[3]
    
    # ref = "/mnt/d/HLAPro_backup/HLAPro/db/ref/HLA_A.fa"


    ###assign reads
    pbin = Pacbio_Binning(raw_fq, outdir, ID)
    pbin.read_bam()

    # pac = PacBio(gene, outdir, ID, ref, db)
    # pac.run()
