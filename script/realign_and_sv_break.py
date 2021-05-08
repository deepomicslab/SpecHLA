import os

def realign_B(raw_bamfile, outdir, version_num, the_para, rlen):
    realign_dir = outdir + '/realign'
    if not os.path.exists(realign_dir):
        os.system('mkdir '+ realign_dir)  
    # test_order = """ %s """
    final_order = "sh /mnt/disk2_workspace/wangmengyao/NeedleHLA/WES/assembly.realign.sh sample %s %s %s"%(rlen, raw_bamfile, realign_dir, version_num, the_para)
    '''
    head_order = "bam=%s\noutdir=%s\nsample=sample\n"%(raw_bamfile, realign_dir)   
    fixed_order_1 = r"""
        if [[ ! -d $outdir ]]; then
                mkdir -p $outdir
        fi

        hla_fa=/home/wangmengyao/scripts/NeedleHLA/ref/hla.ref.extend.fa

        rm -rf $outdir/rematch.total.read.format.txt
        for pos in `cat /home/wangmengyao/scripts/NeedleHLA/script/select.region.txt`
        do
        hla=`echo $pos|cut -d ":" -f 1`
        ref=/home/wangmengyao/scripts/NeedleHLA/ref/$hla
        samtools view -f 64 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##1 "$2,$3}' > $outdir/extract.fa
        samtools view -f 128 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##2 "$2,$3}' >> $outdir/extract.fa
        """
    fixed_order_2="""
        /home/wangmengyao/packages/fermikit/fermi.kit/fermi2.pl unitig -s3g -t4 -l %s -p $outdir/prefix $outdir/extract.fa > $outdir/prefix.mak
        """%(rlen)
    fixed_order_3=r"""
        make -f $outdir/prefix.mak

        if [ ! -f "$outdir/prefix.mag.gz" ]
        then
                echo "$pos"
        else
                less $outdir/prefix.mag.gz|awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}'  > $outdir/prefix.mag.fa
                less $outdir/prefix.mag.gz|awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' | while read L; do  echo $L|awk '{print $1"_r"}' >>$outdir/prefix.mag.fa; read L; echo "$L" | rev | tr "ATGC" "TACG" >>$outdir/prefix.mag.fa; done

                /home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/blastn -query $outdir/prefix.mag.fa -out $outdir/prefix.mag.blast -db $ref -outfmt 6 -strand plus  -penalty -1 -reward 1 -gapopen 4 -gapextend 1
                less $outdir/prefix.mag.blast |cut -f 1|sort|uniq > $outdir/id.list

        if [ ! -f "$outdir/id.list" ]
        then
                echo "$pos"
        else
                samtools faidx $outdir/prefix.mag.fa
                rm -rf $outdir/assembly.fa
                for id in `cat $outdir/id.list`
                do
                    samtools faidx $outdir/prefix.mag.fa $id >> $outdir/assembly.fa
                done

                /home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/blastn -query $outdir/assembly.fa -out $outdir/assembly.blast -db $ref -strand plus -penalty -1 -reward 1 -gapopen 4 -gapextend 1 -line_length 700
                perl /home/BIOINFO_TOOLS/mutation_tools/SamTools/SamTools-1.9/misc/blast2sam.pl $outdir/assembly.blast > $outdir/assembly.blast.sam

                bwa index $outdir/assembly.fa
                bwa mem $outdir/assembly.fa $outdir/extract.fa | samtools view -Sb -F 4 - | samtools sort - > $outdir/rematch.bam
                samtools view -F 0x800 $outdir/rematch.bam|cut -f 1,3,4,6 > $outdir/rematch.bam.txt
        """
    fixed_order = fixed_order_1 + fixed_order_2 + fixed_order_3 
    diff_order = """
                perl /mnt/disk2_workspace/wangmengyao/NeedleHLA/realign/rematchblast.V2.%s.pl $outdir/assembly.blast.sam $outdir/extract.fa $outdir/rematch.bam.txt $outdir/rematch.read.format.txt %s
        """%(version_num, the_para)
    final_order = r"""
                cat $outdir/rematch.read.format.txt >>$outdir/rematch.total.read.format.txt

        fi
        fi
        rm -rf $outdir/prefix* $outdir/id.list
        done
        echo hi
        echo $outdir
        echo $outdir/$sample.realign.bam
        python /mnt/disk2_workspace/wangmengyao/NeedleHLA/realign/realignblast.py -i $bam -o $outdir/$sample.realign.bam -r $outdir/rematch.total.read.format.txt
        samtools sort $outdir/$sample.realign.bam > $outdir/$sample.realign.sort.bam
        java -jar /home/BIOINFO_TOOLS/alignment_tools/Picard/picard-tools-2.1.0/picard.jar FixMateInformation I=$outdir/$sample.realign.sort.bam O=$outdir/$sample.realign.sort.fixmate.bam
        samtools index $outdir/$sample.realign.sort.fixmate.bam
        less /home/wangmengyao/scripts/NeedleHLA/script/vcfhead.txt > $outdir/$sample.realign.vcf
        /home/wangmengyao/packages/freebayes/bin/freebayes -f $hla_fa -p 3 $outdir/$sample.realign.sort.fixmate.bam | awk '$6>10'| grep -v "##" >> $outdir/$sample.realign.vcf && \
        rm -rf $outdir/$sample.realign.vcf.gz
        bgzip $outdir/$sample.realign.vcf
        echo end
        /home/wangmengyao/miniconda2/bin/bcftools filter -t HLA_A:1000-4503,HLA_B:1000-5081,HLA_C:1000-5304,HLA_DPA1:1000-10775,HLA_DPB1:1000-12468,HLA_DQA1:1000-7492,HLA_DQB1:1000-8480,HLA_DRB1:1000-12229 $outdir/$sample.realign.vcf.gz -o $outdir/$sample.realign.filter.vcf"""
    final_order = head_order + fixed_order + diff_order + final_order
    # f = open('test_fermikit.sh','w')
    # print (final_order, file = f)
    '''
    os.system(final_order)

def realign(raw_bamfile, outdir, rlen, version):
    realign_dir = outdir + '/realign'
    if not os.path.exists(realign_dir):
        os.system('mkdir '+ realign_dir)  
    if version == 'WES':
        final_order = "sh /home/wangmengyao/scripts/NeedleHLA/script/run.assembly.realign.wes.sh sample %s %s %s"%( raw_bamfile, realign_dir, rlen)
    elif version == 'WGS':
        final_order = "sh /home/wangmengyao/scripts/NeedleHLA/script/run.assembly.realign.sh sample %s %s %s"%( raw_bamfile, realign_dir, rlen)
    os.system(final_order)

def bam2fq(raw_bamfile, outdir):
    order = "samtools sort -n %s -o %s/sample.for.fq.bam\nsamtools fastq  %s/sample.for.fq.bam -1 %s/sample.read1.fastq -2 %s/sample.read2.fastq\ngzip -f %s/*fastq\nrm %s/sample.for.fq.bam"\
        %(raw_bamfile, outdir, outdir, outdir, outdir, outdir, outdir)
    os.system (order)
    f = open(outdir + '/fq.list','w')
    print ('sample\t%s/sample.read1.fastq.gz\t%s/sample.read1.fastq.gz'%(outdir, outdir), file = f)
    f.close()

def break_point_sv(outdir):
    order = """export PATH=$PATH:/home/wangmengyao/packages/ScanIndel/tools/
        samtools=/home/BIOINFO_TOOLS/mutation_tools/SamTools/SamTools-0.1.18/samtools
        /home/wangmengyao/packages/ScanIndel/tools/gfServer stop localhost 50000
        python /home/wangmengyao/packages/ScanIndel/ScanIndel.py -F 0 -p \
        /home/wangmengyao/packages/ScanIndel/config.txt \
        -i %s/fq.list -o %s
        bam=%s/sample.contigs.bam
        sample=`basename $bam .contigs.bam`
        dir=`dirname $bam`
        group='@RG\tID:'$sample'\tSM:'$sample
        outdir=%s
        perl /mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20200318/dup/ScanIndel/get_breakpoint2.pl $bam $dir/$sample.merged.indel.vcf $outdir/$sample.breakpoint.txt $outdir/$sample.ins.fa
        """%(outdir, outdir, outdir, outdir)
    os.system (order)

def main(raw_bamfile, outdir, rlen, version):
    if os.path.isfile('%s/realign/sample.realign.vcf.gz'%(outdir)) \
        and os.path.isfile('%s/realign/sample.realign.sort.fixmate.bam'%(outdir)):
        print ('realign files exists already!') 
    else:
        realign(raw_bamfile, outdir, rlen, version)
    if os.path.isfile('%s/sample.read1.fastq.gz'%(outdir)):
        print ('fq exists already!')
    else:
        bam2fq(raw_bamfile, outdir)
    # if os.path.isfile('%s/sample.breakpoint.txt'%(outdir)):
    #     print ('break points file for SV exists already!')
    # elif version == 'WGS':
    #     break_point_sv(outdir)

def main_no_realign(raw_bamfile, outdir, rlen, version):
    if os.path.isfile('%s/sample.read1.fastq.gz'%(outdir)):
        print ('fq exists already!')
    else:
        bam2fq(raw_bamfile, outdir)

def main_WES(raw_bamfile, outdir, rlen):
    if os.path.isfile('%s/realign/sample.realign.vcf.gz'%(outdir)):
        print ('realign files exists already!') 
    else:
        realign(raw_bamfile, outdir, rlen, 'WES')
    if os.path.isfile('%s/sample.read1.fastq.gz'%(outdir)):
        print ('fq exists already!')
    else:
        bam2fq(raw_bamfile, outdir)

def typing(outdir, strainsNum):
    order = """perl /home/wangmengyao/scripts/NeedleHLA/script/annoHLA.pl \
        sample %s %s %s"""%(outdir, outdir, strainsNum)
    os.system(order)

def dup_region_type(outdir, strainsNum):
    order = r"""
        bam=%s/realign/sample.realign.sort.fixmate.bam
        outdir=%s
        k=%s
        pos=HLA_DRB1:3898-4400        
        ref=/mnt/disk2_workspace/wangmengyao/NeedleHLA/GA_rich/DRB1/bwa/DRB1_dup_extract_ref.fasta
        samtools view -f 64 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##1 "$2,$3}' > $outdir/extract.fa
        samtools view -f 128 $bam $pos| cut -f 1,6,10|sort|uniq |awk '{OFS="\n"}{print ">"$1"##2 "$2,$3}' >> $outdir/extract.fa
        #/home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/makeblastdb -in DRB1_dup_extract_ref.fasta -dbtype nucl -parse_seqids -out DRB1_dup_extract_ref.fasta
        /home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/blastn -query $outdir/extract.fa -out $outdir/extract.read.blast -db $ref -outfmt 6 -strand plus  -penalty -1 -reward 1 -gapopen 4 -gapextend 1
        perl /mnt/disk2_workspace/wangmengyao/NeedleHLA/GA_rich/DRB1/bwa/count.read.pl $outdir
        less $outdir/DRB1.hla.count| sort -k3,3nr -k4,4nr | head -n $k |awk '$3>0.8'|awk '$4>5' >$outdir/select.DRB1.seq.txt
        """%(outdir, outdir, strainsNum)
    os.system(order)

if __name__ == "__main__":
    sample = 'HLA_2_T_50-50'
    raw_bamfile = '/mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20200318/simu/bwa/%s/%s.bam'%(sample, sample)
    # outdir = 'test_pipeline'
    outdir = 'test_3_18'
    # realign(raw_bamfile, outdir)
    # bam2fq(raw_bamfile, outdir)
    # break_point_sv(outdir)
    # main(raw_bamfile, outdir)
    # typing(outdir, 2)
    dup_region_type(outdir, 2)
