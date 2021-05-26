#!/usr/bin/perl
use FindBin qw($Bin);
#perl annoHLA.pl HLA_10_T_50-50 /mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20200318/assembly/HLA_10_T_50-50/phase/HLA_10_T_50-50 ./ 2
my ($sample,$fadir,$workdir,$k,$pop, $version) = @ARGV;

my $db="$Bin/../../db/HLA";
my $bin="$Bin/../../bin";

my @hlas = ("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1");
my (%hash,%hashp);
open FIN, "$db/HLA_FREQ_HLA_I_II.txt" or die "$!\n";
while(<FIN>){
    chomp;
    my ($gene,$c,$b,$a) = (split);
    if(!$pop){$hashp{$gene} = ($a+$b+$c)/3}
    if($pop eq "Asian"){$hashp{$gene} = $a}
    elsif($pop eq "Black"){$hashp{$gene} = $b}
    elsif($pop eq "Caucasian"){$hashp{$gene} = $c}
    else{$hashp{$gene} = ($a+$b+$c)/3}
}
close FIN;
my ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000);
if($version eq "pstrain"){
     ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (30,500,200,100,100,100,100,70);
}
elsif($version eq "spechap"){
     ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (500,10000,10000,10000,10000,10000,10000,70);
}

my %hashm;
open INL, "$db/hla_nom_g.txt" or die "$!\n";
while(<INL>){
        chomp;
        next if(/^#/);
        my ($tag, $line, $hh) = (split /;/, $_)[0,1,2];
        if($hh){
                my $value = "$tag"."$hh";
                my @arrs = (split /\//,$line);
                foreach my $kk(@arrs){
                        my $hla = "$tag"."$kk";
                        $hashm{$hla} = $value;
                }
        }else{
                my $hla = "$tag"."$line";
                $hashm{$hla} = $hla;
        }

}
close INL;

foreach my $class(@hlas){
	#my $ref="/home/wangmengyao/scripts/NeedleHLA/script/db/HLA_$class.fasta";
	my $ref="$db/whole/HLA_$class";
	if($class eq "DRB1"){$ref="$db/whole/HLA_DRB1.exon";}
	#if($class eq "B"){$ref="/home/wangmengyao/scripts/PHLAT/database/ref/HLA_B";}
	for(my $i=1;$i<=$k;$i++){

               my $fa="$fadir/hla.allele.$i.HLA_$class.fasta";
	       if($class eq "DQB1"){
		       my $j = $i -1;
		       `$bin/samtools faidx $fa HLA_$class\_$j:500-2400 HLA_$class\_$j:5200-7300 >$fadir/HLA_$class.temp.fasta`;
		       $fa = "$fadir/HLA_$class.temp.fasta";
	       }
	       if($class eq "DQA1"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa HLA_$class\_$j:500-900 HLA_$class\_$j:4600-6100 >$fadir/HLA_$class.temp.fasta`;
                       $fa = "$fadir/HLA_$class.temp.fasta";
               }
               if($class eq "DPA1"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa HLA_$class\_$j:400-700 HLA_$class\_$j:4200-5500 >$fadir/HLA_$class.temp.fasta`;
                       $fa = "$fadir/HLA_$class.temp.fasta";
               }
               if($class eq "DPB1"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa HLA_$class\_$j:300-600 HLA_$class\_$j:5000-5300 HLA_$class\_$j:9200-10600 >$fadir/HLA_$class.temp.fasta`;
                       $fa = "$fadir/HLA_$class.temp.fasta";
               }

	       if($class eq "DRB1"){
		       my $j = $i -1;
		       #`samtools faidx $fa HLA_$class\_$j:5950-6260 >$fadir/HLA_$class.temp.fasta`;
		       `$bin/samtools faidx $fa HLA_$class\_$j:100-11000 >$fadir/HLA_$class.temp.fasta`;
                       $fa = "$fadir/HLA_$class.temp.fasta";
	       }
	       if($class eq "A"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa HLA_$class\_$j:100-3300 >$fadir/HLA_$class.temp.fasta`;
                       $fa = "$fadir/HLA_$class.temp.fasta";
               }

	       if($class eq "B"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa HLA_$class\_$j:150-4000 >$fadir/HLA_$class.temp.fasta`;
                       $fa = "$fadir/HLA_$class.temp.fasta";
               }
                 if($class eq "C"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa HLA_$class\_$j:400-3500 >$fadir/HLA_$class.temp.fasta`;
                       $fa = "$fadir/HLA_$class.temp.fasta";
               }
               my $tag = "$class"."_"."$i";
	       `$bin/makeblastdb -in $fa -dbtype nucl -parse_seqids -out $fa`;
	       `$bin/blastn -query $fa -out $workdir/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 1000 `;
	       `$bin/blastn -query $ref.fasta -out $workdir/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 1000 `;

               my (%hash11,%hash12, %hash21, %hash22,$gene,$score);
	      
               open IN1, "$workdir/$tag.blast.out1" or die "$!\n";
               while(<IN1>){
                       chomp;
                       next if(/^#/);
                       my ($hla, $t, $m,$d) = (split)[1,3,4,5];
		       #next if($hla =~ /[N|Q]$/);
		       next if($t <250);
                       $hash11{$hla} += $t;
                       $hash12{$hla} += $m + $d;
              }
              close IN1;

              open IN2, "$workdir/$tag.blast.out2" or die "$!\n";
              while(<IN2>){
	              chomp;
                      next if(/^#/);
         	      my ($hla, $t, $m,$d) = (split)[0,3,4,5];
		      #next if($hla =~ /[N|Q]$/);
		      next if($t <250);
	              $hash21{$hla} += $t;
	              $hash22{$hla} += $m + $d;
              }
              close IN2;
	      #`rm -rf $workdir/blast.out1 $workdir/blast.out2`;
              $score=90;
	      my $ff=0;
              foreach my $key(sort keys %hash11){
		      next if(!exists $hash21{$key});
		      my @tt = (split /:/, $key);
		      my $kid = "$tt[0]".":"."$tt[1]";
		      my $fre=0;
		      if(exists $hashp{$kid}){$fre=$hashp{$kid}}
		      #next if($fre == 0);
	              my $s1 = 100 * (1 - $hash12{$key}/$hash11{$key});
	              my $s2 = 100 * (1 - $hash22{$key}/$hash21{$key});
		      my $s  = $s1;
		      #my $s = ($s1 + $s2)/2;
	              next if($hash11{$key} < 250);
	              next if($hash21{$key} < 250);
		      my ($ts,$tf) = ($s,$fre);
		      my $ttf = $tf;
		      next if(($tf == 0) && ($ts < 99.5));
		      if(($tf > 0) && ($ts==100)){$ttf = 1}
                      if($tf == 0){$ttf = -1}
                      my $scorel=0;
                      if($class eq "B"){$scorel = $ts * (2**($ttf/$C_b));}
		      elsif($class eq "DRB1"){$scorel = $ts * (2**($ttf/$C_drb))}
                      elsif($class eq "A"){
                                  if($tf==0){$ttf=-10}
                                  $scorel = $ts * (2**($ttf/$C_a))
                      }
                      elsif($class eq "C"){
                                 if($tf==0){$ttf=-10}
                                 $scorel = $ts * (2**($ttf/$C_c))
                      }else{$scorel = $ts * (2 ** ($tff/$C_dpa))}

	              if($scorel > $score){
		             $score = $scorel;
		             $gene = $key;
			     $ff=$fre;
	              }
		      #print "$tag\t$key\t$score\t$scorel\t$ts\t$ttf\n";
             }
	     #print "$tag\t$gene\t$score\n";
	     # $hash{$tag} = $hashm{$gene};
	     $hash{$tag} = $gene;
	     `rm -rf $fadir/HLA_$class.temp*`;
     }
}
	     

open OUT, ">$workdir/hla.result.txt";
my $line1 = "Sample";
my $line2 = $sample;
foreach my $kk(sort keys %hash){
	$line1 .= "\t$kk";
        my $hla = $hash{$kk};
	if($hla =~ /DRB1\*14:01/){
		`$bin/samtools  mpileup -r HLA_DRB1:9519-9519 -t DP -t SP -uvf $db/hla.ref.extend.fa $fadir/$sample.merge.bam --output $fadir/snp.vcf`;
		open TE, "$fadir/snp.vcf" or die "$!\n";
		while(<TE>){
			chomp;
		        next if(/^#/);													                                  my $alt = (split)[4];
	                if($alt =~ /T/){}
			else{$hla = "DRB1*14:54"}
		}
		close TE;
	}
	$line2 .= "\t$hla";
}
print OUT "$line1\n$line2\n";

