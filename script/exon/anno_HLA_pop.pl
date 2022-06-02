#!/usr/bin/perl
use FindBin qw($Bin);

my ($sample,$dir,$k,$pop,$version) = @ARGV;
my ($line1,$line2,$line3);

my $db="$Bin/../../db/HLA";
my $bin="$Bin/../../bin";
my @hlas = ("HLA_A","HLA_B","HLA_C","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1","HLA_DRB1");

my %hash;
open FIN, "$db/HLA_FREQ_HLA_I_II.txt" or die "$!\n";
while(<FIN>){
    chomp;
    my ($gene,$c,$b,$a) = (split);
    if(!$pop){$hash{$gene} = ($a+$b+$c)/3}
    if($pop eq "Asian"){$hash{$gene} = $a}
    elsif($pop eq "Black"){$hash{$gene} = $b}
    elsif($pop eq "Caucasian"){$hash{$gene} = $c}
    else{$hash{$gene} = ($a+$b+$c)/3}
}
close FIN;
my ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000);
if($version eq "with_pop"){
     ($C_a,$C_b,$C_c,$C_dpa,$C_dpb,$C_dqa,$C_dqb,$C_drb) = (20,8,200,100,100,100,100,50);
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
my %hashg;
open IN, "$db/Allelelist.txt" or die "$!\n";
while(<IN>){
        chomp;
        next if(/^#/);
        my ($id, $name) = (split /,/, $_)[0,1];
        my $key = "HLA:"."$id";
        my $hla=$name;
        if(exists $hashm{$name}){$hla=$hashm{$name}}
        $hashg{$key} = $hla;
}
close IN;
my %hash_hla;
foreach my $class(@hlas){
	my $ref="$db/exon/$class.fasta";
        my %hashs;
	my %idks;
	`ls $dir/$class.*.fasta >$dir/tfile.list`;
	open TE, "$dir/tfile.list" or die "$!\n";
	open OUT, ">$dir/$class.total.fasta";
	open SOUT, ">$dir/result.$class.fasta";
	while(my $file=<TE>){
		chomp $file;
		my $id = (split /\//,$file)[-1];
		$id =~ s/\.fasta//;
                print OUT ">$id\n";
		my @tes = (split /\./,$id); pop @tes;my $idk = join(".",@tes);$idks{$idk}=$idk;
		my $seq;
		open TT, "$file" or die "$!\n";
		while(my $line=<TT>){
			chomp $line;
			next if($line=~/^>/);
			$seq .= $line;
		}
		close TT;
		print OUT "$seq\n";
		$hashs{$id} = $seq;
	}
	close TE;
	close OUT;
	`$bin/blastn -query $dir/$class.total.fasta -out $dir/$class.blast.out1 -db $ref -outfmt 6 -num_threads 4 -max_target_seqs 10 `;
        open BIN, "$dir/$class.blast.out1" or die "$!\n";
	my (%hash1, %hash2,%hash3, %hashas, %hashhk, %hash4);
	while(<BIN>){
		chomp;
		my ($id, $hla, $t, $m, $i, $si) = (split)[0,1,3,4,5,11];
		next if($hla eq "HLA:HLA10778");
		#next if($hla eq "HLA:HLA03388");
                $hash1{$id}{$hla} += $m + $i;
		$hash2{$id}{$hla} += $t - $i;
		$hash4{$id}{$hla} += $si;
	}
	close BIN;
	my %hashh;
        foreach my $id (sort keys %hash1){
		my (%hashs,$mscore);
		foreach my $hla(sort keys %{$hash1{$id}}){
				my $mis = $hash1{$id}{$hla};
				my $len = $hash2{$id}{$hla};
				my $score = 100 * (1 - $mis/$len);
				if($class =~ /DRB1/){$score = $hash4{$id}{$hla} / 499}
				my $nhla = $hashg{$hla};
				my @tt = (split /:/, $nhla);
				my $kid = "$tt[0]".":"."$tt[1]";
				my $fre=0;
				if(exists $hash{$kid}){$fre=$hash{$kid}}
				#my $scorel = $score * (2 ** ($fre/20));
				my $tk = "$id\t$nhla";
				$hashhk{$tk} = "$score\t$fre";
				$hash3{$id}{$nhla} = $score;
                                if(($score >= $mscore) && ($fre > 0)){$mscore = $score}
                                $hashs{$mscore} .= "\t$hla";
				$hashas{$id} .= "$nhla\t";
		}
                $hashh{$id} = "$mscore\t$hashs{$mscore}";
	}
	my %hashss;
	my $mscore = 0;
	foreach my $idk (sort keys %idks){
		 my $ts = 0;
                 for(my $i=1;$i<=$k;$i++){
			 my $idd = "$idk"."."."$i";
                         my $s = (split /\t/, $hashh{$idd})[0];
			 $ts += $s;
		 }
		 $ts = $ts / $k;
		  print "$idk\t$ts\n";
		 if($ts >= $mscore){$mscore = $ts}
		 $hashss{$ts} .= "$idk\t";
	 }
         my @tids = (split /\t/, $hashss{$mscore});
	 my (%fhash,$ffre,%thash);
	 
         foreach my $tid(@tids){
		 my $tfre;
		 for(my $i=1;$i<=$k;$i++){
			 my $hkey = "$class"."_"."$i";
			 my $idd = "$tid"."."."$i";
			 my @arrs = (split /\t+/, $hashh{$idd});	
			 my $aa = shift @arrs;
			 my (%nhash,$mfre);
			 foreach $od(@arrs){
				 my $nd = $hashg{$od};
				 my @tt = (split /:/, $nd);
				 my $id = "$tt[0]".":"."$tt[1]";
				 my $fre=0;
				 if(exists $hash{$id}){$fre=$hash{$id}}
				 $nhash{$fre} .= "$nd\t";
				 #	 	 print "$idd\t$od\t$nd\t$fre\n";
				 if($fre >= $mfre){$mfre = $fre}
			 }
			 $tfre += $mfre;
			 if($mfre == 0){$tfre = $tfre/2;}
			 $thash{$idd} = $nhash{$mfre};
			 #print "$idd\t$tfre\t$aa\n";

		 }
		 if($tfre >=$ffre){$ffre = $tfre}
                 $fhash{$tfre} = $tid;
  	 }
         my $select_id = $fhash{$ffre};
	 print "Select_id\t$select_id\n";
         my ($hout,$scout,$nout);
	 for($i=1;$i<=$k;$i++){
		 $nout .= "\t$class"."_"."$i";
		 my $idd = "$select_id"."."."$i";
		 my $seq = $hashs{$idd};
		 print SOUT ">$idd\n$seq\n";
		 my ($hh,$mscorel) = ("",0);
		 $hh = (split /\t/,$thash{$idd})[0];

		 #}else{
		     my @hlas = (split /\t/, $hashas{$idd});
		     foreach my $hla (sort @hlas){
			 my $hk = "$idd\t$hla";
			 my ($ts,$tf) = (split /\t/,$hashhk{$hk});
			 next if(($tf == 0) && ($ts < 99.5) && (!$class =~ /DRB1/));
			 my $ttf = $tf;
			 if(($tf > 0) && ($ts==100)){$ttf = 1}
			 if($tf == 0){$ttf = -1}
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
			  print "$idd\t$hla\t$scorel\t$ts\t$tf\n";
			 if($scorel>=$mscorel){$mscorel = $scorel;$hh=$hla}
		     }
		     # }   
		     #exclude region HLA_B:670-730
		 if(($hh =~ /^B/) && ($hash3{$idd}{$hh} <99.26)){
			 my $nseq1 = substr($seq,0,182);
			 my $nseq2 = substr($seq,248,);
			 my $nseq = "$nseq1"."$nseq2";
			 open OTT, ">$dir/B.tmp.fa";
			 print OTT ">$idd\n$nseq\n";
			 close OTT;
			 `$bin/blastn -query $dir/B.tmp.fa -out $dir/B.blast.tmp -db $ref -outfmt 6 -num_threads 4 -max_target_seqs 500`;
			 open BT, "$dir/B.blast.tmp" or die "$!\n";
			 my (%hashbt,%hashbf);
			 while(<BT>){
				 chomp;
				 my ($b1,$b2) = (split)[1,11];
				 my $nhla = $hashg{$b1};
                                 my @tt = (split /:/, $nhla);
                                 my $kid = "$tt[0]".":"."$tt[1]";
                                 my $fre=0;
                                 if(exists $hash{$kid}){$fre=$hash{$kid}}
				 next if($fre<=0);
				 $hashbt{$b1} += $b2;
				 $hashbf{$b1} = $fre;
			 }
			 close BT;
			 my ($shla,$btc);
			 foreach my $b3(sort keys %hashbt){
				 my $ss = ($hashbt{$b3} /10 )* (2 ** ($hashbf{$b3}/7));
				 	print "$hashg{$b3}\t$ss\t$hashbt{$b3}\t$hashbf{$b3}\n";
				 if($ss >= $btc){$btc=$ss;$shla=$b3}
			 }
			 $hh = $hashg{$shla};
			 
		 }

		 #$hout .= "\t$hh";
		 $scout .= "\t$hash3{$idd}{$hh}";
	 
	 #$hash_hla{$class} = "$hout".";"."$scout";
                 if($hh =~ /DRB1\*14:01/){
	                  ` $bin/samtools  mpileup -r HLA_DRB1:9519-9519 -t DP -t SP -uvf $db/hla.ref.extend.fa $dir/$sample.merge.bam --output $dir/snp.vcf`;
                          open TE, "$dir/snp.vcf" or die "$!\n";
                          while(<TE>){
		                 chomp;
	                         next if(/^#/);
	                         my $alt = (split)[4];
	                         if($alt =~ /T/){} else{$hh = "DRB1*14:54"}      
	                  }
                          close TE;	       
	        }
		$hout .= "\t$hh";
         }
         $line1 .= "$nout";
	 $line2 .= "$hout";
	 $line3 .= "$scout";
	 `rm -rf $dir/$class.total.fasta`;
	 close SOUT;
 
}
open OUT, ">$dir/hla.result.txt";
print OUT "Sample"."$line1\n$sample"."$line2\nScore"."$line3\n";
close OUT;


