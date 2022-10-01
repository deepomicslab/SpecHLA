#!/usr/bin/perl -w
use FindBin qw($Bin);
use Getopt::Long;

my ($sample, $dir, $pop, $wxs, $g_nom, $help);
$g_nom=0;
GetOptions(
           "s=s"     =>      \$sample,
           "i=s"     =>      \$dir,
           "p=s"     =>      \$pop,
           "r=s"     =>      \$wxs,
           "g=s"     =>      \$g_nom,
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
        -r       <tr>    focus region "exon|whole" ("exon" is suitable for WES or RNAseq; "whole" is suitable for WGS )
        -g       <tr>     G-translate 1|0"
        -help|?           print help information
e.g.:
        perl $0 -s samplename -i indir -p Unknown -r exon -g 1
USE
die $usage unless ($sample && $dir && $pop && $wxs );
print "parameter:\tsample:$sample\tdir:$dir\tpop:$pop\twxs:$wxs\tG_nom:$g_nom\n";

my $k = 2;
my (%hashp,%hashp2, %hashpp, %hashg, %hashc, %hash,%hashdd);
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
    if($wxs eq "exon"){
       $a = sprintf "%.3f",$a;
       $b = sprintf "%.3f",$b;
       $c = sprintf "%.3f",$c;
    }
    if($wxs eq "whole"){
       $a = sprintf "%.8f",$a;
       $b = sprintf "%.8f",$b;
       $c = sprintf "%.8f",$c;
    }
    $hashpp{$gene} = "$c;$b;$a";
    $hashp2{$gene} = ($a+$b+$c)/3;
    if($pop eq "Unknown"){$hashp{$gene} = ($a+$b+$c)/3}
    if($pop eq "Asian"){$hashp{$gene} = $a}
    if($pop eq "Black"){$hashp{$gene} = $b}
    if($pop eq "Caucasian"){$hashp{$gene} = $c}
    if($pop eq "nonuse"){$hashp{$gene} = 0}
}
close FIN;

#Convert HLA nomenclature
open INL, "$db/hla_nom_g.txt" or die "$!\n";
while(<INL>){
        chomp;
        next if(/^#/);
        my ($tag, $line, $hh) = (split /;/, $_)[0,1,2];
        if($hh && $g_nom){
                my $value = "$tag"."$hh";
                my @arrs = (split /\//,$line);
                foreach my $kk(@arrs){
                        my $hla = "$tag"."$kk";
                        # next if($hla =~ /DQB1\*02:02/);
                        # next if($hla =~ /DQB1\*03/);
                        # next if($hla =~ /DQA1\*05/);
                        # next if($hla =~ /DQA1\*03/);
                        # next if($hla =~ /C\*07/);
                        # next if($hla =~ /DRB1\*14/);
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
                system("$bin/blastn -query $workdir/$class.$i.temp.fasta -out $workdir/$class.$i.blast.out1 -db $ref -outfmt 6 -num_threads 4 -max_target_seqs 6000") ;           
     ## Read blast result
                open BIN, "$workdir/$class.$i.blast.out1" or die "$!\n";
                my (%hash1, %hash2,%hash3, %hashas, %hashhk, %hash4);
                my $blastcount =0;
                while(<BIN>){
                     chomp;
                     my ($id, $hla, $t, $m, $i, $si) = (split)[0,1,3,4,5,11];
                     next if($hla eq "HLA:HLA17631");
                     next if($hla eq "HLA:HLA08917");
                     next if($t<50);
                     my $nhla = $hashc{$hla};
                     $hash1{$nhla} += $m + $i;
                     $hash2{$nhla} += $t - $i;
                     $hash4{$nhla} += $si;
                     $blastcount += 1;
                 }
                 close BIN;
                 my $tag = "$class"."_"."$i";
                 #print "$tag\t$blastcount\n";
                 my %hash_max;
                 my ($hh,$mscorel) = ("",0);
                 foreach my $hla (sort keys %hash1){
                         my $mis = $hash1{$hla};
                         my $len = $hash2{$hla};
                         my $l = 1; 
                         if($class =~ /DQB/){$l = int($hash4{$hla}/500 + 0.5) + 1;}
                         my $score = 100 * (1 - $mis/$len) * $l;
                         if($class =~ /DRB1/){$score = $hash4{$hla} }
                #print "$hla\t$score\t$mis\t$len\t$hash4{$hla}\n";        
                         my @tt = (split /:/, $hla);
                         my $kid = "$tt[0]".":"."$tt[1]";
                         my ($fre,$fre2)=(0,0);
                         if(exists $hashp{$kid}){$fre=$hashp{$kid}; $fre2=$hashp2{$kid}} #fre:population frequency of 4 digit HLA allele
                         $hash3{$hla} = $score;
                         next if($pop ne "nonuse" && $fre2 == 0);
                         if($score>=$mscorel){$mscorel = $score;$hh=$hla;$hash_max{$mscorel} .= "$hla;$score\t"}
                 }

                $hash{$tag} = $hash_max{$mscorel};
         }
#        `rm -rf $workdir/*`;
          
     }
}


sub whole_blast{
    foreach my $class(@hlas){
        my $ref="$db/whole/$class";
        if($class =~ "DRB1" && $pop ne "nonuse"){$ref="$db/whole/HLA_DRB1.exon";}
        for(my $i=1;$i<=$k;$i++){
               my $fa="$fadir/hla.allele.$i.$class.fasta";
               #extract the diversity region for annotation
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
               system("$bin/blastn -query $fa -out $workdir/$tag.blast.out1 -db $ref -outfmt 7 -num_threads 4 -max_target_seqs 6000 ");
               my (%hash_max,%hash11,%hash12, %hash21, %hash22,$gene,$score);
               ## read blast score
               open IN1, "$workdir/$tag.blast.out1" or die "$!\n";
               my $blastcount = 0;
               while(<IN1>){
                       chomp;
                       next if(/^#/);
                       my ($hla, $t, $m,$d,$s) = (split)[1,3,4,5,8];
                       $hash11{$hla} += $t;
                       $hash12{$hla} += $m + $d;
                       $blastcount += 1;
              }
              close IN1;
              next if($blastcount == 0); ###no blast result
              $score=50;
              my $ff=0;
              foreach my $key(sort keys %hash11){
                      #next if(!exists $hash21{$key});
                      my @tt = (split /:/, $key);
                      my $kid = "$tt[0]".":"."$tt[1]";
                      my ($fre,$fre2)=(0,0);
                      if(exists $hashp{$kid}){$fre=$hashp{$kid};$fre2=$hashp2{$kid}} ## population frequency of 4 digit hla allele
                      my $s = 100 * (1 - $hash12{$key}/$hash11{$key}); #blast score
                      next if($pop ne "nonuse" && $fre2 == 0);
                      my $scorel=$s;

                      if($scorel >= $score){
                             $score = $scorel;
                             $gene = $key;
                             $ff=$fre;
                             $hash_max{$scorel} .= "$gene;$s\t";
                      }
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

#print the alleles annotation of best score for each HLA gene
my $hout = $sample;
foreach my $hla(@hlas){
       for(my $i=1;$i<=$k;$i++){
             my $id = "$hla"."_"."$i";
             my ($line1,$line2,$line3,$pout,$out) = ("","","","","-");
             my (%ggs, %hashmm);
             if(!exists $hash{$id}){ $hout .= "\t-"; print OUT "$id\t-\t-\t-\n"}
             next if(!exists $hash{$id});
             my @arrs = (split /\t/,$hash{$id});
             my ($max,$pfre) = (0,0);
             foreach my $oo(@arrs){
                 my ($allele,$score) = (split /;/,$oo)[0,1];
                 $score = sprintf "%.3f", $score;
                 #DRB1*14:01 and DRB1*14:54 differ in HLA_DRB1:9519
                 if($allele =~ /DRB1\*14:01/){
                          system("$bin/samtools  mpileup -r HLA_DRB1:9519-9519 -t DP -t SP -uvf $db/hla.ref.extend.fa $dir/$sample.merge.bam --output $workdir/snp.vcf");
                          open TE, "$workdir/snp.vcf" or die "$!\n";
                          while(<TE>){
                                 chomp;
                                 next if(/^#/);
                                 my $alt = (split)[4];
                                 if($alt =~ /T/){print "$allele\n"} else{$allele = "DRB1*14:54";}
                          }
                          close TE;
                 }
                 #C*07:01 and C*07:18 differ in HLA_C:4061
                 if($allele =~ /C\*07:01/ && $wxs eq "exon"){
                          system("$bin/samtools  mpileup -r HLA_C:4061-4061 -t DP -t SP -uvf $db/hla.ref.extend.fa $dir/$sample.merge.bam --output $workdir/snp.vcf");
                          open TE, "$workdir/snp.vcf" or die "$!\n";
                          while(<TE>){
                                 chomp;
                                 next if(/^#/);
                                 my $alt = (split)[4];
                                 if($alt =~ /T/){$allele = "C*07:18:01:01";}
                          }
                          close TE;
                 }

                 my @tt = (split /:/, $allele);
                 my $kid = "$tt[0]".":"."$tt[1]";
                 #$hashmm{$kid} .= "$allele\t";
                 $line2 .= "$allele;";
                 $oo = "$kid".";"."$score";
                 if(exists $hashg{$allele}){
                      my @ttt = (split /:/,$hashg{$allele});
                      my $tid = "$ttt[0]".":"."$ttt[1]";
                      $ggs{$tid} += 1;
                      $hashmm{$tid} .= "$allele\t";
                 }
                 else{$ggs{$kid} += 1; $hashmm{$kid} .= "$allele\t"}
                 if(exists $hashpp{$kid}){$line3 .= "$oo;$hashpp{$kid}\t"; 
                     if($pfre < $hashp{$kid}){$pfre = $hashp{$kid};$pout = $allele;
                     }
                 }
                 else{$line3 .= "$oo;-;-;-\t";} 
             }
             my @lines3 = (split /\t/,$line3); @lines3 = uniq(@lines3);
             my @lines2 = (split /;/,$line2); $line3 = join("\t",@lines3);
             if($pop eq "nonuse" || !$pout){
               foreach my $gg(sort {$ggs{$b} <=> $ggs{$a}} keys %ggs ){
                   my $agg = (split /\t/, $hashmm{$gg})[0];
                   if($ggs{$gg} >= $max){
                       $max = $ggs{$gg};$line1 .= "$agg;";
                        if($out eq "-"){$out = $agg;}
                   }
                   if($g_nom == 1 && exists $hashg{$out}){$out = $hashg{$out}}
                }
             }else{ 
                   if(exists $hashg{$pout}){$pout = $hashg{$pout};}
             #      else{$pout = $lines2[0]}
                   $out = $pout;   $line1 = $pout; 
             }
             if(!$out){$out="-"}
             $hout .= "\t$out";
             print OUT "$id\t$line1\t$line2\t$line3\n";
       }    
}
close OUT;
print COUT "$hout\n";
close COUT;
