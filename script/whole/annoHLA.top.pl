#!/usr/bin/perl -w
use FindBin qw($Bin);
use Getopt::Long;

my ($sample, $dir, $pop, $mode, $wxs, $help);

GetOptions(
           "s=s"     =>      \$sample,
           "i=s"     =>      \$dir,
           "p=s"     =>      \$pop,
           "m=s"     =>      \$mode,
           "r=s"     =>      \$wxs,
           "h"       =>      \$help
);
my $usage = <<USE;
Usage:
description: HLAtyping annotation of SpecHLA phased sequence
usage: perl $0 [options]
        Common options:
        -s       <tr>    sample name
        -i       <tr>    the directory of phased sequence
        -p       <tr>    population information "Asian|Black|Caucasian|Unknown|nonuse"
        -m       <tr>    Spechla mode "pstrain|spechap"
        -r       <tr>    focus region "exon|whole" ("exon" is suitable for WES or RNAseq; "whole" is suitable for WGS )
        -help|?           print help information
e.g.:
        perl $0 -s samplename -i indir -p Unknown -m spechap -r exon
USE
die $usage unless ($sample && $dir && $pop && $mode && $wxs) ;

print "parameter:\tsample:$sample\tdir:$dir\tpop:$pop\tmode:$mode\twxs:$wxs\n";

my $k = 2;
my (%hashp, %hashpp, %hashg, %hashc, %hash,%hashdd);
my $db="$Bin/../../db/HLA";
my $bin="$Bin/../../bin";
my @hlas = ("HLA_A","HLA_B","HLA_C","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1","HLA_DRB1");
my $fadir=$dir;
my $workdir = "$dir/tmp";
`mkdir  -p $workdir`;
#population frequency
open FIN, "$db/HLA_FREQ_HLA_I_II.txt" or die "$!\n";
<FIN>;
while(<FIN>){
    chomp;
    my ($gene,$c,$b,$a) = (split);
    $a = sprintf "%.3f",$a;
    $b = sprintf "%.3f",$b;
    $c = sprintf "%.3f",$c;
    $hashpp{$gene} = "$c;$b;$a";
    if($pop eq "Unknown"){$hashp{$gene} = ($a+$b+$c)/3}
    if($pop eq "Asian"){$hashp{$gene} = $a}
    if($pop eq "Black"){$hashp{$gene} = $b}
    if($pop eq "Caucasian"){$hashp{$gene} = $c}
    if($pop eq "nonuse"){$hashp{$gene} = 0}
}
close FIN;

my ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000);
if($wxs eq "whole" && $mode eq "pstrain"){
     ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (30,500,200,100,100,100,100,70);
}
elsif($wxs eq "whole" && $mode eq "spechap"){
     ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (500,10000,10000,10000,10000,10000,10000,70);
}
elsif($wxs eq "exon" && $mode eq "pstrain"){
     ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (20,8,200,100,100,100,100,50);
}
elsif($wxs eq "exon" && $mode eq "spechap"){
     ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (20,8,200,100,100,100,100,50);
}
#Convert HLA nomenclature
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
                        $hashg{$hla} = $value;
                }
        }else{
                my $hla = "$tag"."$line";
                $hashg{$hla} = $hla;
        }

}
close INL;
#convertion of HLA id
open IN, "$db/Allelelist.txt" or die "$!\n";
while(<IN>){
        chomp;
        next if(/^#/);
        my ($id, $name) = (split /,/, $_)[0,1];
        my $key = "HLA:"."$id";
        my $hla=$name;
   #     if(exists $hashg{$name}){$hla=$hashg{$name}}
        $hashc{$key} = $hla;
}
close IN;

sub exon_blast{
    foreach my $class(@hlas){
        my $ref="$db/exon/$class.fasta";
        my (%hashs, %idks);
        for(my $i=1;$i<=$k;$i++){
                my $fa = "$dir/hla.allele.$i.$class.fasta";
                `echo ">$class.$i" >$workdir/$class.$i.temp.fasta`;
                `less $fa |grep -v ">" >>$workdir/$class.$i.temp.fasta`;
                my $seq = `less $fa|grep -v ">"`; chomp $seq; $seq=~s/\s+//g;
                `$bin/blastn -query $workdir/$class.$i.temp.fasta -out $workdir/$class.$i.blast.out1 -db $ref -outfmt 6 -num_threads 4 -max_target_seqs 1000 `;
                open BIN, "$workdir/$class.$i.blast.out1" or die "$!\n";
                my (%hash1, %hash2,%hash3, %hashas, %hashhk, %hash4);
                while(<BIN>){
                     chomp;
                     my ($id, $hla, $t, $m, $i, $si) = (split)[0,1,3,4,5,11];
                     next if($hla eq "HLA:HLA10778");
                     my $nhla = $hashc{$hla};
                     $hash1{$nhla} += $m + $i;
                     $hash2{$nhla} += $t - $i;
                     $hash4{$nhla} += $si;
                 }
                 my $tag = "$class"."_"."$i";
                 my (%hash_max, %hash_score);
                 my ($hh,$mscorel) = ("",0);
                 foreach my $hla (sort keys %hash1){
                         my $mis = $hash1{$hla};
                         my $len = $hash2{$hla};
                         my $score = 100 * (1 - $mis/$len);
                         if($class =~ /DRB1/){$score = $hash4{$hla} / 499}
                    
                         my @tt = (split /:/, $hla);
                         my $kid = "$tt[0]".":"."$tt[1]";
                         my $fre=0;
                         if(exists $hashp{$kid}){$fre=$hashp{$kid}}
                         $hash3{$hla} = $score;
                         my ($ts,$tf) = ($score,$fre);
                         next if(($tf == 0) && ($ts < 99.5) && $pop ne "nonuse" && (!$class =~ /DRB1/));
                         my $ttf = $tf;
                         if(($tf > 0) && ($ts==100)){$ttf = 1}
                         if($tf == 0 && $pop ne "nonuse"){$ttf = -1}
                         if($pop eq "nonuse"){$ttf = 0}
                         my $scorel=0;
                         if($class eq "HLA_B"){$scorel = $ts * (2**($ttf/$C_b));}
                         elsif($class eq "HLA_DRB1"){
                                 #if($tf >0){$scorel = $ts ;}else{$scorel=70}
                                 $scorel = $ts * (2 ** ($ttf/$C_drb))
                         }
                         elsif($class eq "HLA_A"){
                                  if($tf==0){$ttf=-10}
                                  $scorel = $ts * (2**($ttf/$C_a))
                         }
                         elsif($class eq "HLA_C"){
                                 if($tf==0){$ttf=-10}
                                 $scorel = $ts * (2**($ttf/$C_c))
                         }else{$scorel = $ts * (2**($ttf/$C_dpb));}
                         $hash_score{$scorel} .= "$hh;$ts\t";      
                         if($scorel>=$mscorel){$mscorel = $scorel;$hh=$hla;$hash_max{$mscorel} .= "$hla;$ts\t"}
                 }
                 #exclude region HLA_B:670-730
                 if(($hh =~ /^B/) && ($hash3{$hh} <99.26)){
                         my $nseq1 = substr($seq,0,182);
                         my $nseq2 = substr($seq,248,);
                         my $nseq = "$nseq1"."$nseq2";
                         open OTT, ">$workdir/B.tmp.fa";
                         print OTT ">$hh\n$nseq\n";
                         close OTT;
                         `$bin/blastn -query $workdir/B.tmp.fa -out $workdir/B.blast.tmp -db $ref -outfmt 6 -num_threads 4 -max_target_seqs 500`;
                         open BT, "$workdir/B.blast.tmp" or die "$!\n";
                         my (%hashbt,%hashbf);
                         while(<BT>){
                                 chomp;
                                 my ($b1,$b2) = (split)[1,11];
                                 my $nhla = $hashc{$b1};
                                 my @tt = (split /:/, $nhla);
                                 my $kid = "$tt[0]".":"."$tt[1]";
                                 my $fre=0;
                                 if(exists $hashp{$kid}){$fre=$hashp{$kid}}
                                 next if($fre<=0 && $pop ne "nonuse");
                                 $hashbt{$b1} += $b2;
                                 $hashbf{$b1} = $fre;
                         }
                         close BT;
                         my ($shla,$btc) = ("",0);
                         foreach my $b3(sort keys %hashbt){
                                 my $ss = ($hashbt{$b3} /10 )* (2 ** ($hashbf{$b3}/7));
                              #          print "$hashc{$b3}\t$ss\t$hashbt{$b3}\t$hashbf{$b3}\n";
                                 if($ss >= $btc){$btc=$ss;$shla=$b3}
                         }
                         $hh = $hashc{$shla};
                         $hash_max{$mscorel} = "$hh;$mscorel";
                         $hash_score{$mscorel} = "$hh;$mscorel";
                 }
                 if($hh =~ /DRB1\*14:01/){
                          ` $bin/samtools  mpileup -r HLA_DRB1:9519-9519 -t DP -t SP -uvf $db/hla.ref.extend.fa $dir/$sample.merge.bam --output $workdir/snp.vcf`;
                          open TE, "$workdir/snp.vcf" or die "$!\n";
                          while(<TE>){
                                 chomp;
                                 next if(/^#/);
                                 my $alt = (split)[4];
                                 if($alt =~ /T/){} else{$hh = "DRB1*14:54";
                                          $hash_max{$mscorel} = "$hh;$mscorel";
                                          $hash_score{$mscorel} = "$hh;$mscorel";
                                 }      
                          }
                          close TE;            
                }
                my $n=10; my $tout;
                foreach my $cs(sort {$b <=> $a} keys %hash_score){
                   if($n>0){$tout .= "$hash_score{$cs}\t";$n -= 1;}
                }
                #$hash{$tag} = $hash_max{$mscorel};
                $hash{$tag} = $tout;

         }
#        `rm -rf $workdir/*`;
          
     }
}


sub whole_blast{
    foreach my $class(@hlas){
        my $ref="$db/whole/$class";
        if($class eq "DRB1"){$ref="$db/whole/HLA_DRB1.exon";}
        for(my $i=1;$i<=$k;$i++){
               my $fa="$fadir/hla.allele.$i.$class.fasta";
               if($class eq "HLA_DQB1"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa $class\_$j:500-2400 $class\_$j:5200-7300 >$fadir/$class.temp.fasta`;
                       $fa = "$fadir/$class.temp.fasta";
               }
               if($class eq "HLA_DQA1"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa $class\_$j:500-900 $class\_$j:4600-6100 >$fadir/$class.temp.fasta`;
                       $fa = "$fadir/$class.temp.fasta";
               }
               if($class eq "HLA_DPA1"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa $class\_$j:400-700 $class\_$j:4200-5500 >$fadir/$class.temp.fasta`;
                       $fa = "$fadir/$class.temp.fasta";
               }
               if($class eq "HLA_DPB1"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa $class\_$j:300-600 $class\_$j:5000-5300 $class\_$j:9200-10600 >$fadir/$class.temp.fasta`;
                       $fa = "$fadir/$class.temp.fasta";
               }
               if($class eq "HLA_DRB1"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa $class\_$j:100-11000 >$fadir/$class.temp.fasta`;
                       $fa = "$fadir/$class.temp.fasta";
               }
               if($class eq "HLA_A"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa $class\_$j:100-3300 >$fadir/$class.temp.fasta`;
                       $fa = "$fadir/$class.temp.fasta";
               }

               if($class eq "HLA_B"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa $class\_$j:150-4000 >$fadir/$class.temp.fasta`;
                       $fa = "$fadir/$class.temp.fasta";
               }
                 if($class eq "HLA_C"){
                       my $j = $i -1;
                       `$bin/samtools faidx $fa $class\_$j:400-3500 >$fadir/$class.temp.fasta`;
                       $fa = "$fadir/$class.temp.fasta";
               }
               my $tag = "$class"."_"."$i";
               `$bin/makeblastdb -in $fa -dbtype nucl -parse_seqids -out $fa`;
               `$bin/blastn -query $fa -out $workdir/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 1000 `;
               `$bin/blastn -query $ref.fasta -out $workdir/$tag.blast.out2 -db $fa -outfmt 7 -num_threads 4 -max_target_seqs 1000 `;
               my (%hash_max,%hash11,%hash12, %hash21, %hash22,$gene,$score);
              
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
              $score=90;
              my $ff=0;
              foreach my $key(sort keys %hash11){
                      next if(!exists $hash21{$key});
                      my @tt = (split /:/, $key);
                      my $kid = "$tt[0]".":"."$tt[1]";
                      my $fre=0;
                      if(exists $hashp{$kid}){$fre=$hashp{$kid}}
                      my $s1 = 100 * (1 - $hash12{$key}/$hash11{$key});
                      my $s2 = 100 * (1 - $hash22{$key}/$hash21{$key});
                      my $s  = $s1;
                      #my $s = ($s1 + $s2)/2;
                      next if($hash11{$key} < 250);
                      next if($hash21{$key} < 250);
                      my ($ts,$tf) = ($s,$fre);
                      my $ttf = $tf;
                      next if(($tf == 0) && ($ts < 95));
                      if(($tf > 0) && ($ts==100)){$ttf = 1}
                      if($tf == 0){$ttf = -1}
                      my $scorel=0;
                      if($class eq "HLA_B"){$scorel = $ts * (2**($ttf/$C_b));}
                      elsif($class eq "HLA_DRB1"){$scorel = $ts * (2**($ttf/$C_drb))}
                      elsif($class eq "HLA_A"){
                                  if($tf==0){$ttf=-10}
                                  $scorel = $ts * (2**($ttf/$C_a))
                      }
                      elsif($class eq "HLA_C"){
                                 if($tf==0){$ttf=-10}
                                 $scorel = $ts * (2**($ttf/$C_c))
                      }else{$scorel = $ts * (2 ** ($ttf/$C_dpa))}

                      if($scorel >= $score){
                             $score = $scorel;
                             $gene = $key;
                             $ff=$fre;
                             $hash_max{$scorel} .= "$gene;$ts\t";
                      }
                      #print "$tag\t$key\t$score\t$scorel\t$ts\t$ttf\n";
             }
             $hash{$tag} = $hash_max{$score};
             `rm -rf $fadir/$class.temp*`;
        }
    }
}

open COUT, ">$dir/hla.result.txt";
print COUT "Sample\tHLA_A_1\tHLA_A_2\tHLA_B_1\tHLA_B_2\tHLA_C_1\tHLA_C_2\tHLA_DPA1_1\tHLA_DPA1_2\tHLA_DPB1_1\tHLA_DPB1_2\tHLA_DQA1_1\tHLA_DQA1_2\tHLA_DQB1_1\tHLA_DQB1_2\tHLA_DRB1_1\tHLA_DRB1_2\n";

open OUT, ">$dir/hla.result.details.txt";
print OUT "Gene\tG_best\tallele\tdetails:allele;Score;Caucasian;Black;Asian\n";
if($wxs eq "exon"){
       &exon_blast; 
}
if($wxs eq "whole"){
       &whole_blast;
}
sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
my $hout = $sample;
foreach my $hla(@hlas){
       for(my $i=1;$i<=$k;$i++){
             my $id = "$hla"."_"."$i";
             my @arrs = (split /\t/,$hash{$id});
             my ($line1,$line2,$line3,%ggs);
             foreach my $oo(@arrs){
                 if($oo =~ /:/){
                 my ($allele,$score) = (split /;/,$oo)[0,1];
                 $score = sprintf "%.3f", $score;
                 my @tes = (split /:/,$allele);
                 $alle = "$tes[0]".":"."$tes[1]";
                 $oo = "$alle".";"."$score";
                 my @tt = (split /:/, $allele);
                 my $kid = "$tt[0]".":"."$tt[1]";
                 $line2 .= "$allele;";
                 if(exists $hashg{$allele}){$ggs{$hashg{$allele}} += 1}
                 else{$ggs{$allele} += 1}
                 if(exists $hashpp{$kid}){$line3 .= "$oo;$hashpp{$kid}\t"}
                 else{$line3 .= "$oo;-;-;-\t";} 
             }}
             my $max = 0; my $out="-";
             foreach my $gg(sort {$ggs{$b} <=> $ggs{$a}} keys %ggs ){
                  if($ggs{$gg} >= $max){
                       $max = $ggs{$gg};$line1 .= "$gg;";
                        if($out eq "-"){$out = $gg;}
                   }
             }
             $hout .= "\t$out";
             my @lines3 = (split /\t/,$line3); @lines3 = uniq(@lines3);
             $line3 = join("\t",@lines3);
             print OUT "$id\t$line1\t$line2\t$line3\n";
       }    
}
close OUT;
print COUT "$hout\n";
close COUT;
