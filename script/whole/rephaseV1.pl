#!/usr/bin/perl -w
#perl rephase.pl break_points.txt ./ 2 prephase_breakpoints.txt 300 5
my ($bfile, $vdir, $k, $outfile,$len, $point) = @ARGV;

my $db="../../db/HLA";
my $bin="../../bin";

my $hla_ref="$db/hla.ref.extend.fa";
my $outdir="$vdir/tmp";
#my ($len, $point) = (300,5);
`rm -rf $outdir`;
`mkdir $outdir`;
my %hash;
my $count =0;
open IN, "$bfile" or die "$!\n";
<IN>;
my $head="#gene\tlocus\t00\t01\t10\t11";
print "input_rephase_para:$len\t$point\n";
while(<IN>){
	chomp;
	my $gen = (split)[0];
	$hash{$gen} .= "$_\n";
	$count += 1;
}
close IN;

my %hashlen=('HLA_A'=>'3503', 'HLA_B'=>'4081', 'HLA_C'=>'4304', 'HLA_DPA1'=>'9775','HLA_DPB1'=>'11468','HLA_DQA1'=>'6492','HLA_DQB1'=>7480,'HLA_DRB1'=>'11229');

open OUT, ">$outfile";
print OUT "$head";

for(my $i=0; $i<$k;$i++){
	print OUT "\tphase$i";
}
print OUT "\n";
my ($idss,$ref, $tag1,$tag2, $break1,$break2,$vcf, $start1,$start2,$end1,$end2,$n,$gene,$num1,$num2);
($tag1,$tag2) = ("y","y");
foreach my $ge(sort keys %hash){
	$gene = $ge;
        my @lines = (split /\n/, $hash{$gene});
	$vcf = "$vdir/"."$gene".".vcf";
	$ref = "$db/whole/$gene";

	if($gene eq "HLA_A"){$len=100}
	#if($gene eq "HLA_DQB1"){$len=400}
	#if($gene eq "HLA_DRB1"){$len=1000}
	if($gene eq "HLA_DPB1"){$point=0;$len=150}
	#`/home/wangmengyao/miniconda2/bin/bcftools sort $vcf -O z -o $vcf.gz`;
	`tabix -f $vcf.gz`;
	($start1, $n) = (1001, 0);
	if($count == 1){$start1 = 1501}
	while($n<=$#lines){
		
	        ($break1,$num1) = (split /\s/, $lines[$n])[1,6];
		my @oarrs = (split /\s/, $lines[$n])[0,1,2,3,4,5];
                my $out = join("\t", @oarrs);
		if($count==1){$break1 += 10}	
		$n += 1;
		print "$break1\n";
		if(($gene eq "HLA_DQB1") && ($break1 >4230) && ($break1 < 4600)){$break1 =4230}
		#if(($gene eq "HLA_DQB1") && ($start1 <4599) && ($break1 > 4600)){$break1 =4230}
		if(($gene eq "HLA_DQB1") && ($start1 <4230) && ($break1 > 4600)){$break1 =4230}
		next if($break1 - $start1 < $len);
		print "$gene\t$start1\t$break1\tregion1\n";
                &rephase_blast1;
		next if($tag1 =~ /noid/);
                my $tag = 1; $num2 = 0;
		while($tag && $n <= $#lines+1){
	            if($n>$#lines){$num2=$point + 1; $break2=$hashlen{$gene} + 1000; $n+=1; $tag=0;}else{
                    my $m;
	            ($break2,$m) = (split /\s/, $lines[$n])[1,6];
		    $num2 += $m;
		     print "$break1\t$break2\n";
		    if(($gene eq "HLA_DQB1") && ($break2 >4230) && ($break2 <4600)){$break2 = 4230}
		    elsif(($gene eq "HLA_DQB1") &&  ($break1 <4230) && ($break2>4230) && ($break2 <4600)){$break2 = 4230}
		    elsif(($gene eq "HLA_DQB1") && ($break1 >=4230) && ($break1<4600) && ($break2>4600)){$break1 =4600}
		    elsif(($gene eq "HLA_DQB1") && ($break1 < 4230) && ($break2 >4600)){$break1 = 4600 }
		    if (($break2 - $break1 < $len) || ($num2 <=$point)){$n += 1; $tag=1; next;}
		    else{$tag = 0;}
	           }
		   # if($n>$#lines){$num2=$point + 1; $break2=`less $vcf|grep "#"|grep $gene  |cut -d "=" -f 4|awk 'gsub(">","")'` + 1000;}
		    print "$gene\t$break1\t$break2\tregion2\n";
                    &rephase_blast2;
		    next if(!$idss);
		    next if($idss =~ /10/);
		    next if($tag2 eq "mix");
	       
               }
	        next if($idss =~ /10/);
		print OUT "$out\t$idss\n";

                print "$gene\t$start1\t$end1\t$gene\t$start2\t$end2\t$num2\t$idss\n";
                $start1 = $start2;
		$end1 = $end2;

        }

}
close OUT;


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
			#	if( !exists $ha{$max} ){ $tag1 = "noid"; print "noid\t$gene\t$start1\t$end1\t$max\n";}else{$tag1 = "y";}
			#next if( !exists $ha{$max});
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
		      if($sx[1]>$sx[0]){$ids[$idx[0]] = $tt}
		      else{$ids[$idx[1]] = $tt}}
		$idss = join("\t",@ids);
}
