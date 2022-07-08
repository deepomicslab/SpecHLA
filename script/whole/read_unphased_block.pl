#!/usr/bin/perl -w
use FindBin qw($Bin);
#perl  read_unphased_block.pl  break_points.txt ./ 2 prephase_breakpoints.txt 300 5
my ($bfile, $vdir, $k, $outfile, $wxs) = @ARGV;

my $db="$Bin/../../db/HLA";
my $bin="$Bin/../../bin";
my $hla_ref="$db/hla.ref.extend.fa";
my $outdir="$vdir/tmp";
`rm -rf $outdir`;
`mkdir $outdir`;
my %hash;
my $count =0;
open IN, "$bfile" or die "$bfile\t$!\n";
<IN>;
while(<IN>){
      chomp;
      my $gen = (split)[0];
      $hash{$gen} .= "$_\n";
      $count += 1;
}
close IN;
#ouput block region
my %hashlen=('HLA_A'=>'3503', 'HLA_B'=>'4081', 'HLA_C'=>'4304', 'HLA_DPA1'=>'9775','HLA_DPB1'=>'11468','HLA_DQA1'=>'6492','HLA_DQB1'=>'7480','HLA_DRB1'=>'11229');
my %hashe=('HLA_A'=>'1504-1773;2015-2290;2870-3145', 'HLA_B'=>'1486-1755;2001-2276;2851-3126', 'HLA_C'=>'1699-1968;2215-2490;3078-3353', 'HLA_DPA1'=>'5208-5453;5794-6075','HLA_DPB1'=>'6002-6265;10217-10498','HLA_DQA1'=>'5600-5848;6262-6543','HLA_DQB1'=>'3073-3342;6232-6513','HLA_DRB1'=>'6972-7241');
my (%hasha,%hashb);
foreach my $key(sort keys %hashe){
     my @rrs = (split /;/,$hashe{$key});
     foreach my $rr(@rrs){
         my ($s,$e) = (split /-/,$rr)[0,1];
         my $o1 = "$key\t$s";
         my $o2 = "$key\t$e";
         $hasha{$o1} = $key;
         $hashb{$o2} = $key;
     }
}
my %hashr;
open OUT, ">$outfile";
my ($ref,$region1,$region2,$break1,$break2,$vcf,$start1,$start2,$end1,$end2,$gene,$n,$score1,$score2);
foreach my $ge(sort keys %hash){
        $gene = $ge;
        my @lines = (split /\n/, $hash{$gene});
        $vcf = "$vdir/"."$gene".".spechap.vcf.gz";  # bug
        $ref = "$db/whole/$gene";
        ($start1,$n) = (1001,0);
        while($n<=$#lines){               
                my @oarrs = (split /\s/, $lines[$n]);
                $break1 = $oarrs[2];
                my $out = join("\t", @oarrs);
                $n += 1;
                if($n>$#lines){ $break2=$hashlen{$gene} + 1000; $n+=1; }
                else{($break2) = (split /\s/, $lines[$n])[2];}
                $end1 = $break1;
                $start2 = $break1 + 1;
                $end2 = $break2;
                #skip the long indel enriched region
                if(($gene eq "HLA_DQB1") && ($break2 >4230) && ($break2 <4600)){$end2 = 4230}
                elsif(($gene eq "HLA_DQB1") && ($break1 < 4230) && ($break2>4230) && ($break2 <4600)){$end2 = 4230}
                elsif(($gene eq "HLA_DQB1") && ($break1 >=4230) && ($break1<4600) && ($break2>4600)){$end1 = 4230; $start2 =4600}
                elsif(($gene eq "HLA_DQB1") && ($break1 < 4230) && ($break2 >4600)){$end1 = 4230; $start2=4600 }
                if(($gene eq "HLA_DRB1") && ($break2 >3950) && ($break2 <4300)){$end2 = 3950}
                elsif(($gene eq "HLA_DRB1") && ($break1 < 3950) && ($break2>3950) && ($break2 <4300)){$end2 = 3950}
                elsif(($gene eq "HLA_DRB1") && ($break1 >=3950) && ($break1<4300) && ($break2>4300)){$end1=3950; $start2 =4300}
                elsif(($gene eq "HLA_DRB1") && ($break1 < 3950) && ($break2>4300)){$end1 = 3950; $start2=4300 }
                $region1 = "$gene".":"."$start1"."-"."$end1";
                $region2 = "$gene".":"."$start2"."-"."$end2";
                $hashr{$region1} = $start1;
                $hashr{$region2} = $start2;
                $start1 = $start2;
       }
}
my %hashrr;
#add exon boundary for wes data
if($wxs eq "wes"){
       $ref = "$db/exon/$gene.fasta";
       my @arrs = (split /;/, $hashe{$gene});
       my $start = (split /-/,$arrs[0])[0];
       my $end = (split /-/,$arrs[-1])[1];
       foreach my $region(sort {$hashr{$a} <=> $hashr{$b}} keys %hashr){
                      #print "$region\t$hashr{$region}\n";
                      my ($g,$rr) = (split /:/,$region)[0,1];
                      my ($s,$e) = (split /-/,$rr)[0,1];
                      next if($e <= $start || $s >= $end);
                      my ($re,$ore) = ($region,"");
                      foreach my $arr (@arrs){
                           my ($es,$ee) = (split /-/,$arr)[0,1];
                           next if($ee <= $s);
                           if($s < $es && $ee <= $e){$ore .= "$g".":"."$arr\t";}
                           
                           if(($e >= $es && $e <= $ee) || ($s>=$es && $s <= $ee)){
                                   $re = "$gene".":"."$s"."-"."$e";
                                   my $ss = $s; my $se = $ee;
                                   if($s<$es && $e >$es){$ss = $es}
                                   if($e > $es && $e < $ee){$se = $e}
                                   $ore .= "$g".":"."$ss"."-"."$se\t";    
                                  
                           }
                            
                      }
                      print "$re\t$ore\n";
                      $hashrr{$re} = $ore;
        }

}else{%hashrr = %hashr;}
$n=0;
#return combination of each two hap
foreach my $region1(sort keys %hashrr){ 
        foreach my $region2(sort keys %hashrr){
                next if($region1 eq $region2);
                my ($re1,$re2); 
                if($wxs eq "wes"){$re1 = $hashrr{$region1};$re2 = $hashrr{$region2}}
                else{$re1 = $region1; $re2 = $region2}
                for(my $j=1; $j<=$k; $j++){
		     my $gene = (split /:/,$re1)[0];
                     `echo ">allele$j.break1" > $outdir/allele$j.break1.fa`;
                     `echo ">allele$j.break2" > $outdir/allele$j.break2.fa`;
                     `$bin/samtools faidx $hla_ref $re1 | $bin/bcftools consensus -H $j $vcf.gz | grep -v ">" >> $outdir/allele$j.break1.fa`;
                     `$bin/samtools faidx $hla_ref $re2 | $bin/bcftools consensus -H $j $vcf.gz | grep -v ">" >> $outdir/allele$j.break2.fa`;
               }
                `cat $outdir/allele1.break1.fa $outdir/allele1.break2.fa $outdir/allele2.break1.fa $outdir/allele2.break2.fa > $outdir/allele.break.merge.$n.fa`;
                `$bin/blastn -query $outdir/allele.break.merge.$n.fa -out $outdir/allele.break.merge.blast.$n -db $ref -outfmt 6 -max_target_seqs 10000 -num_threads 4 -strand plus`; 
                
                 open TE, "$outdir/allele.break.merge.blast.$n" or die "blast\t$!\n";
                 my (%hash1,%hash2,%hash11,%hash12,%hash21,%hash22); my ($score1,$score2) = (0,0);
		 my (%hashm11,%hashm12,%hashm21,%hashm22);
                 while(<TE>){
                      chomp;
                      next if(/^#/);
                      my ($id,$hla,$t,$m,$i) = (split)[0,1,3,4,5];
		      next if($t < 250 && $gene eq "HLA_DRB1" && $wxs ne "wes");
                      if($id eq "allele1.break1"){$hash11{$hla} += $t; $hashm11{$hla} += $m+$i}
                      if($id eq "allele1.break2"){$hash12{$hla} += $t; $hashm12{$hla} += $m+$i}
                      if($id eq "allele2.break1"){$hash21{$hla} += $t; $hashm21{$hla} += $m+$i}
                      if($id eq "allele2.break2"){$hash22{$hla} += $t; $hashm22{$hla} += $m+$i}
                 }
                 close TE;
                 ##return the score of hap 00/01;
                 my $hhla; my $max = 0;
                 foreach my $hh(sort keys %hash11){
                      my ($s11,$s12,$s21,$s22) = (0,0,0,0);
                      if(exists $hash11{$hh}){$s11 = 100*(1- $hashm11{$hh}/$hash11{$hh});}
                      if(exists $hash12{$hh}){$s12 = 100*(1- $hashm12{$hh}/$hash12{$hh});}
                      if(exists $hash21{$hh}){$s21 = 100*(1- $hashm21{$hh}/$hash21{$hh});}
                      if(exists $hash22{$hh}){$s22 = 100*(1- $hashm22{$hh}/$hash22{$hh});}
                      my $ss1 = ($s11 + $s12)/2;
                      my $ss2 = ($s21 + $s22)/2;
                      my $rss1 = ($s11 + $s22)/2;
                      my $rss2 = ($s12 + $s21)/2;
                      my ($S1,$S2);
                      if($ss1 >= $ss2){$S1 = $ss1}else{$S1 = $ss2}
                      if($rss1 >= $rss2){$S2 = $rss1}else{$S2 = $rss2}
                      if($S1 >= $max && $S1 >= $S2){
                           $hhla=$hh; $score1 = $S1; $score2 = $S2;
                           $max = $S1; $hash1{$max} .= "\t$score1;$score2;$hhla";
                      }
                      if( $S2 >= $max && $S2 >= $S1){
                           $hhla=$hh; $score2 = $S2; $score1 = $S1;
                           $max = $S2; $hash2{$max} .= "\t$score1;$score2;$hhla";
                      }
                           
       #if($region1 eq "HLA_DRB1:2433-2601" && $region2 eq "HLA_DRB1:2602-3808" ){print "$region1\t$region2\t$hh\t$ss1\t$ss2\t$rss1\t$rss2\t$S1\t$S2\t$score1\t$score2\t$max\n"}
                 }
                 $n += 1;
                 my $out = "";
                 if(exists $hash1{$max}){$out .= $hash1{$max}}
                 if(exists $hash2{$max}){$out .= $hash2{$max}}
                 print OUT "$region1\t$region2\t$out\n";
        }
}                
close OUT;

