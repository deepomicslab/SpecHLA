#!/usr/bin/perl -w
#perl rephase.pl break_points.txt ./ 2 prephase_breakpoints.txt 300 5
my ($bfile, $vdir, $k, $outfile,$len, $point) = @ARGV;

my $db="../../db/HLA";
my $bin="../../bin";

my $hla_ref="$db/hla.ref.extend.fa";
my $outdir="$vdir/tmp";
# ($len, $point) = (200,3);
`rm -rf $outdir`;
`mkdir $outdir`;
my $count=0;
`ls $vdir/HLA_DRB1_part_*.vcf >$outdir/vcf.list`;
open COUT, ">$outfile";
print COUT "gene\tlocus\t00\t01\t10\t11\tphase0\tphase1\n";

my %hashf;
open FI, "$outdir/vcf.list" or die "$!\n";
while(<FI>){
	chomp;
        my $file = $_;
	s/\.vcf//;
        my $line = (split /\//,$_)[-1];
	my ($c,$s,$e) = (split /_/,$line)[3,4,5];
	#next if(($e==12229) && ($s>10000));
	$hashf{$c} = "$s\t$e\t$file";
	$count += 1;
}
close FI;
print "Count:$count\n";
my %hash;
my $ccount=0;
open IN, "$bfile" or die "$!\n";
<IN>;
while(<IN>){
	chomp;
	$ccount += 1;
	my ($gen,$pos) = (split)[0,1];
	next if($gen ne "HLA_DRB1");
	foreach my $key(sort keys %hashf){
             my ($s,$e,$file) = (split /\t/, $hashf{$key})[0,1,2];
             if(($s<=$pos) && ($e>$pos)){
		     $hash{$key} .= "$_\n";
	     }
	}
}
close IN;

if($ccount > 0){
open OUT, ">$outdir/DRB1.temp.prephase_breakpoints.txt";
my $pp;
for(my $i=0; $i<$k;$i++){
	$pp .= "\t$i";
}
my ($idss,$ref, $vcf,$tag1,$tag2, $break1,$break2,$start1,$start2,$end1,$end2,$n,$gene,$num1,$num2);
($gene,$tag1,$tag2) = ("HLA_DRB1","y","y");
foreach my $contig(sort keys %hashf){
    my ($s,$e,$v) = (split /\t/,$hashf{$contig})[0,1,2];
    my $cc = $contig -1;
    if(($s>3900) && ($s<4100)){$s=4300}
    my ($s2,$e2,$v2);
    if(exists $hashf{$cc}){($s2,$e2,$v2) = (split /\t/, $hashf{$cc})[0,1,2];}
    # if(!exists $hash{$contig}){print OUT "$contig\tHLA_DRB1\t$s2\t$e2\t$s\t$e"."$pp\n"}
    #else{
      
    if(exists $hash{$contig}){
        my @lines = (split /\n/, $hash{$contig});
        $vcf= $v;
	$ref = "$db/whole/HLA_DRB1";
        $start1 = $s;
	`$bin/bcftools sort $vcf -O z -o $vcf.gz`;
	`$bin/tabix $vcf.gz`;
	$n=0;
	while($n<=$#lines){
	        ($break1,$num1) = (split /\s/, $lines[$n])[1,6];
		my @oarrs = (split /\s/, $lines[$n])[0,1,2,3,4,5];
                my $out = join("\t", @oarrs);
		$n += 1;
		next if($break1 - $start1 < $len);
		#print "$gene\t$start1\t$break1\tregion1\n";
                &rephase_blast1;
		next if($tag1 =~ /noid/);
                my $tag = 1; $num2 = 0;
		while($tag && $n <= $#lines+1){
	            if($n>$#lines){$num2=$point + 1; $break2=$e; $n+=1; $tag=0;}else{
		    my $m;
	            ($break2,$m) = (split /\s/, $lines[$n])[1,6];
		    $num2 += $m;
		    if (($break2 - $break1 < $len) || ($num2 <$point)){$n += 1; $tag=1; next;}
		    else{$tag = 0;}
	           }
		   #   print "$gene\t$break1\t$break2\tregion2\n";
		   next if($break2-$break1 < $len - 50) ;
		    &rephase_blast2;
		    next if(!$idss);
		    next if($idss =~ /10/);
		    next if($tag2 eq "mix");	       
               }
	        next if($idss =~ /10/);
		next if($start1 >= $start2);
		#print OUT "$contig\t$out\t$idss\n";
		
                print OUT "$contig\t$gene\t$start1\t$end1\t$start2\t$end2\t$idss\n";
                print "$gene\t$start1\t$end1\t$gene\t$start2\t$end2\t$num2\t$idss\n";
                $start1 = $start2;
		$end1 = $end2;
        }
   }else{print OUT "$contig\t$gene\t$start2\t$end2\t$s\t$e"."$pp\n";$start2=$s;$end2=$e}
}
close OUT;

my (%hashl,%hashtt, %hashpt, %hashtg);
open CIN, "$outdir/DRB1.temp.prephase_breakpoints.txt" or die "$!\n";
while(<CIN>){
    chomp;
    my @arrs = (split /\t/, $_);
    my $contig = shift @arrs;
    my $out = join("\t",@arrs);
    $hashpt{$contig} .= "$_\n";
    my $gene = shift @arrs;
    my $vcf = (split /\t/,$hashf{$contig})[2];
    `$bin/bcftools sort $vcf -O z -o $vcf.gz`;
    `$bin/tabix -f $vcf.gz`;

    $hashtt{$contig} += 1;
    my ($start1,$end1,$start2,$end2,$p1,$p2) = @arrs[0,1,2,3,4,5];
    my $name1 = "$gene"."-"."$contig"."-"."1";
    my $name2 = "$gene"."-"."$contig"."-"."2";
    if($hashtt{$contig} == 1){
	    my $ta1 = `$bin/samtools faidx $hla_ref $gene:$start1-$end1 | $bin/bcftools consensus -H 1 $vcf.gz |grep -v ">"`;
            chomp $ta1;
	    $hashl{$name1} = $ta1;
	    my $ta2 = `$bin/samtools faidx $hla_ref $gene:$start1-$end1 | $bin/bcftools consensus -H 2 $vcf.gz |grep -v ">"`;
	    chomp $ta2;
	    $hashl{$name2} = $ta2;
	   
	    # if($ta1 eq $ta2){print "$gene:$start1-$end1\t$vcf\tsameta\n"}else{print "$gene:$start1-$end1\t$vcf\tdiffta\n"}
    }

    my $fa1=`$bin/samtools faidx $hla_ref $gene:$start2-$end2 | $bin/bcftools consensus -H 1 $vcf.gz |grep -v ">"`;
    my $fa2=`$bin/samtools faidx $hla_ref $gene:$start2-$end2 | $bin/bcftools consensus -H 2 $vcf.gz |grep -v ">"`;
    chomp $fa2;
    chomp $fa1;
    $hashtg{$contig} .= "$p1";
    #if($fa1 eq $fa2){print "$gene:$start2-$end2\t$vcf\tsamefa\n"}
    #else{print "$gene:$start2-$end2\t$vcf\tdifffa\n"}
    if($p1==0){
           $hashl{$name1} .= $fa1;
	   $hashl{$name2} .= $fa2;
    }else{
	   $hashl{$name1} .= $fa2;
	   $hashl{$name2} .= $fa1;
	   my $tmp = $hashl{$name1};
	   $hashl{$name1} = $hashl{$name2};
	   $hashl{$name2} = $tmp;
    }
    
}
close CIN;


my @blocks;
&generate_block;
my %hashb;
foreach my $block(@blocks){
     my $blockr;
     my $len = length($block);
     my ($fa1,$fa2);
     for(my $i=0;$i<$len;$i++){
	     my $str = substr($block,$i,1);
	     my $strr = abs($str - 1);
	     $blockr .= $strr;
	     my $ptag = 0;
	     if(exists $hashtg{$i}){$ptag= ($hashtg{$i} =~ tr/1/1/);}
             my $pp = $str + 1; 
	     my $pr = $strr + 1;
	     my $id1 = "$gene"."-"."$i"."-"."$pp";
	     my $id2 = "$gene"."-"."$i"."-"."$pr";
	     next if(!exists $hashl{$id1});
	     if($ptag %2 ==0){
	          $fa1 .= $hashl{$id1};
	          $fa2 .= $hashl{$id2};
	     }else{
	          $fa1 .= $hashl{$id2};
		  $fa2 .= $hashl{$id1};
	     }
     } 
     $fa1 =~ s/\s//g;
     $fa2 =~ s/\s//g;
     #print TT ">DRBA_$block\n$fa1\n";
     `echo '>DRB1_$block\n$fa1' > $outdir/DRB1_$block.fa`;
     `echo '>DRB1_$blockr\n$fa2' > $outdir/DRB1_$blockr.fa`;
     `$bin/blastn -query $outdir/DRB1_$block.fa -out $outdir/DRB1_$block.blast -db $ref -outfmt 7 -max_target_seqs 1000 -num_threads 4 -strand plus`;
     `$bin/blastn -query $outdir/DRB1_$blockr.fa -out $outdir/DRB1_$blockr.blast -db $ref -outfmt 7 -max_target_seqs 1000 -num_threads 4 -strand plus`;
     my $max = 0;
     my $okey;
     my @bbs = ($block,$blockr);
     foreach my $bb(@bbs){
	   my $max1 = 80;
           my (%hasht1,%hasht2,%hasht3);
           open MIN, "$outdir/DRB1_$bb.blast" or die "$!\n";
           while(<MIN>){
	     chomp;
	     next if(/^#/);
	     my ($r1,$r2,$r3,$r4,$r5) = (split /\t/,$_)[0,1,3,4,11];
	     next if($r3 < $len);
	     my $k1 = "$r1\t$r2";
             $hasht1{$k1} += $r5;
             $hasht2{$k1} += $r3;
             $hasht3{$k1} += $r4;
          }
          close MIN;
          foreach my $key1(sort keys %hasht1){
		  #	  next if($hasht1{$key1} < 200);
		  my $score = 100 * (1 - ($hasht3{$key1}/$hasht2{$key1}));
		  if($score > $max1){$max1 = $score;$okey=$key1}
	  }
	  $max += $max1;
	  print "$bb\t$max1\t";
    }
    $hashb{$block} = $max;
    print "$block\t$max\n";
}

my $select_block;
my $max=90;
foreach my $key(sort keys %hashb){
     my $score = $hashb{$key};
     #   print "$score\t$key\n";
     if($score>=$max){$max = $score;$select_block=$key}
}
print "$select_block\t$max\n";
#open CIN, "$outdir/DRB1.temp.prephase_breakpoints.txt" or die "$!\n";
my %thash;
for(my $i=0;$i<$count;$i++){
      my $contig = $i;
      my ($s,$e,$v) = (split /\t/,$hashf{$contig})[0,1,2];
      my $tag = (substr($select_block,0,$contig+1) =~ tr/1/1/);
      my $ptag = 0;
      $ptag = substr($select_block,$contig,1);
      if(!exists $hashpt{$contig}){my $gt = abs(1-$ptag);print COUT "$gene\t$s\t+\t+\t+\t+\t$ptag\t$gt\n"}
      next if(!exists $hashpt{$contig});
      my @arrs = (split /\n/, $hashpt{$contig});
      my $line1 =  $arrs[0];
     
      my $s1 = (split /\t/,$line1)[3];
      my $j=0;
      my $pp;
      foreach my $line(@arrs){
              my ($start1,$end1,$start2,$end2,$p1,$p2) = (split /\t/,$line)[2,3,4,5,6,7];
              $j += 1;
              $pp .= "$p1";
	      my $GT;
	      if(($i>0) && ($#arrs>0) && ($j==1)){
                      if($ptag==0){$GT=0}else{$GT=1}
                      my $gt = abs($GT-1);
                      print COUT "$gene\t$start1\t+\t+\t+\t+\t$GT\t$gt\n";
              }

              if($#arrs==0){
                     if($ptag==0){$GT=0}
                     if($ptag==1){$GT=1}
		     my $gt = abs($GT-1);
		     print COUT "$gene\t$start2\t+\t+\t+\t+\t$GT\t$gt\n";
              }
	      else{
                     my $np = ($pp =~ tr/1/1/);
                     if(($ptag ==0) && ($np%2==0)){$GT = 0}
                     if(($ptag ==0) && ($np%2==1)){$GT = 1}
                     if(($ptag ==1) && ($np%2==0)){$GT = 1}
                     if(($ptag ==1) && ($np%2==1)){$GT = 0}
                     my $pos = $end1 + 2;
	             my $gt = abs($GT-1);
	             print COUT "$gene\t$pos\t-\t-\t-\t-\t$GT\t$gt\n";
              }
              
      }

}
close COUT;


sub generate_block{
     my @arrs = ("0");
     for(my $i=1;$i<$count;$i++){
           my @arrst;
           foreach my $arr(@arrs){
                   my $arr0 = "$arr"."0";
                   my $arr1 = "$arr"."1";
                   push @arrst, $arr0;
                   push @arrst, $arr1;
           }
           @arrs = @arrst;
     }
     @blocks = @arrs;
}

sub rephase_blast1{
                $end1 = $break1 -2;
		$tag1 = "y";
		next if ($start1>=$end1);
		#print "region1\t$gene\t$start1\t$end1\n";
		for(my $j=1; $j<=$k; $j++){
			`$bin/samtools faidx $hla_ref $gene:$start1-$end1 | $bin/bcftools consensus -H $j $vcf.gz > $outdir/allele$j.break1.fa`;
			`$bin/blastn -query $outdir/allele$j.break1.fa -out $outdir/allele$j.break1.blast -db $ref -outfmt 7 -max_target_seqs 1000 -num_threads 4 -strand plus`;	
			my $max = 200; my @fix; my %ha;
			open TMP, "$outdir/allele$j.break1.blast" or die "$!\n";
			while(<TMP>){
				chomp;
				next if(/^#/);
				my ($id, $pre, $rlen, $mis, $score) = (split)[1,2,3,4,11];
				next if($rlen < ($end1 - $start1)/3);
				#next if($pre < 95);
				#next if($mis > 55);
				if($score > $max){$max = $score}
				if(($pre == 100) && ($rlen > 100)){ push @fix, $id;}
				$ha{$score} .= "$id\n";
			}
                        close TMP;
			`rm -rf $outdir/allele$j.hla.fasta`;
                        my @hlas;
		        if(exists $ha{$max}){@hlas = (split /\n/, $ha{$max});}
			push @hlas, @fix;
			my %has; my @uhlas=grep{++$has{$_}<2}@hlas;
			if(!@uhlas){$tag1 .= "noid";}else{
				$tag1 .= "y";
				#	print "$gene\t$start1\t$end1\t@uhlas\t$tag1\n";
			        foreach my $hla(@uhlas){
				    chomp $hla;
				    `$bin/samtools faidx $ref.fasta $hla >> $outdir/allele$j.hla.fasta`;
			    }  
			}
			#print "$#uhlas\t$tag1\n";
		}
	}
sub rephase_blast2{
         $start2 = $break1 -1;
         $end2 = $break2 -2;
	 my $len2 = $end2 - $start2;
	 #next if($len2 < 100);
	 next if($end2<=$start2);
         my $region="$gene:$start2-$end2";
	 # print "region2\t$gene\t$start2\t$end2\n";
	 `rm -rf $outdir/break2.fa`;
	 my ($oid,$oscore)=(9,0);
         my (@ids, @scores,@sabs);
         for(my $j=1; $j<=$k; $j++){	 
		         `$bin/samtools faidx $hla_ref $region | $bin/bcftools consensus -H $j $vcf.gz | sed  "s/$region\$/allele$j.break2/" >> $outdir/break2.fa`;}
	for($j=1; $j<=$k; $j++){ 
	                `$bin/makeblastdb -in $outdir/allele$j.hla.fasta -dbtype nucl -parse_seqids -out $outdir/allele$j.hla -logfile $outdir/log.txt`;
			`$bin/blastn -query $outdir/break2.fa -out $outdir/allele$j.break2.blast -db $outdir/allele$j.hla -outfmt 7 -max_target_seqs 1000 -num_threads 4 -strand plus `;
                        my %hasht; my %hasht2; my %hasht3;
			open TE, "$outdir/allele$j.break2.blast" or die "$!\n";
			while(<TE>){
				chomp;
				next if(/^#/);
				my ($r1,$r2,$r3,$r4,$r5) = (split /\t/,$_)[0,1,3,4,11];
                                my $k1 = "$r1\t$r2";
				#next if($r3 < 200);
				$hasht{$k1} += $r5;
				$hasht2{$k1} += $r3;
                                $hasht3{$k1} += $r4;
			}
			close TE;
			my ($id,$iscore)= (11,80);
			my $tagi = 50;
			my ($a1,$a2)=(0,0);
			foreach my $key1(sort keys %hasht){
				# if($hasht{$key1} > $tagi){
					 my $tags = $hasht{$key1};
					 next if ($tags < 50);
					 $key1 =~ /^allele(\d)/;
					 my $score = 100 * (1 - ($hasht3{$key1}/$hasht2{$key1}));
                                         if($score > $iscore){$id = $1; $iscore = $score}
					 if(($1 == 1) && ($score >= $iscore )){$a1=$score}
					 if(($1 == 2) && ($score >= $iscore )){$a2=$score}
					 #if(($score == $iscore) && ($tags>$tagi)){$id = $1; $tagi = $tags}
				         
					 # }
			 }
			if(($a1 == $iscore) && ($a2 == $iscore)){$id=11}
			print "$a1\t$a2\t$id\n";
			if(!$id){$id=11}
			$id = $id -1;
			$ids[$j-1] = $id;
			$scores[$j-1] = $iscore;
			$sabs[$j-1] = abs($a1-$a2);
			$oid = $id;
			$oscore = $iscore;

		}
		my $cc = 0; my %hl; my $tt; my @nn;
		print "$region\t@ids\t@scores\tnew\n";
		for(my $a=0;$a<$k;$a++){
			if($ids[$a] == 10){$cc += 1}
			else{$hl{$ids[$a]} += 1}
		}
		if($cc >1){$tag2 = "mix"}else{$tag2 = "y"}
		for($a=0;$a<$k;$a++){
			if(! exists $hl{$a}){$tt = $a}
			if(exists $hl{$ids[$a]}){
		         	if($hl{$ids[$a]} > 1){push @nn, $a}
			}
		}
		my (@idx, @sx);
		for($a=0;$a<$k;$a++){
			if(($cc==1) && ($ids[$a] == 10)){$ids[$a] = $tt}
			if($#nn ==1){
			    if(($a==$nn[0]) || ($a==$nn[1])){push @idx, $a; push @sx, $scores[$a]}
		        }
		}
                if($#nn == 1){
		      if($sx[1]>=$sx[0]){$ids[$idx[0]] = $tt}
		      else{$ids[$idx[1]] = $tt}}
		$idss = join("\t",@ids);
}

}
