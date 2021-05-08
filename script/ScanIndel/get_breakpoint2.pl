#!/usr/bin/perl -w

my ($insam, $vcf, $outfile, $outfa) = @ARGV;
my %hash;
my $region;
open OUT2, ">$outfa";
open OUT, ">$outfile";
open IN, "samtools view $insam|" or die "$!\tno samfile\n";
while(<IN>){
	chomp;
	next if(/^#/);
	my ($start,$end);
	my $insert = ".";
        my($hla, $pos, $qual, $cigar, $seq) = (split)[2,3,4,5,9];
	next if($pos < 1000);
	next if($qual < 5);
	next if(($hla eq "HLA_DRB1") && ($pos <3000) && ($pos>2000));
	if($cigar =~ /^(\d+)S/){
		my $len = $1;
                s/SA:Z:(HLA_\S+)//;
                my ($nhla,$npos) = (split /,/,$1)[0,1];
		next if(($nhla eq "HLA_DRB1") && ($npos >3800) && ($npos <4400));
                next if(($nhla =~ /HLA/) && ($nhla ne $hla));
		if($len >150){
			my $ins = substr($seq, 0, $len);
			$region = "$hla\t$pos\t$pos";
			$hash{$region} = "$hla\t$pos\t"."+"."\t$hla\t$pos\t"."+"."\t$ins\t$seq\n";}
	}
	elsif($cigar =~ /(\d+)S$/){
		my $len = $1;
		s/SA:Z:(HLA_\S+)//;
		my ($nhla,$npos) = (split /,/,$1)[0,1];
	        next if(($nhla eq "HLA_DRB1") && ($npos >3800) && ($npos <4400));	
		next if(($nhla =~ /HLA/) && ($nhla ne $hla));
		if($len>150){
			my $ins = substr($seq,-$len);
			my $a = length($seq) - $len;
			$pos += $a;
			$region = "$hla\t$pos\t$pos";
			$hash{$region} = "$hla\t$pos\t"."+"."\t$hla\t$pos\t"."+"."\t$ins\t$seq\n";}
	}
	
	elsif($cigar =~ /^(\S+M)(\d\d\d+)([I|D])\d+M/){
		my $len = $2;
		my $type = $3;
		my $line = $1;
		my $k;
		if($len > 150){
                      while($line){
			      $line =~ s/^(\d+)(\w)//;
			      my($n, $a) = ($1,$2);
			      $k .= $a x $n;
		      }
		      $k =~ s/D//g;
		      my $kl = length($k);
		      $start = $pos + $kl ;
		      $end = $start + $len;
		      $insert = substr($seq,$kl-1,$len);
		      if($type eq "D"){
		            $insert = ".";
		      	      $region = "$hla\t$start\t$end";
			      # print "$region\t$len\n";
		      	      $hash{$region} = "$hla\t$start\t"."+"."\t$hla\t$end\t"."+"."\t$insert\n";
		       }
		       else{
		             
		      $region = "$hla\t$start\t$start";
		      #print "$region\t$len\n";
                      $hash{$region} = "$hla\t$start\t"."+"."\t$hla\t$start\t"."+"."\t$insert\t$seq\n";
		    
			       }

		}
	}
}
close IN;

open IN2, "$vcf" or die "$!\tnovcf\t$vcf\n";
while(<IN2>){
	chomp;
	next if(/^#/);
        my ($hla,$pos,$qual,$line) = (split /\t/,$_)[0,1,5,7];
	my ($type,$cigar,$len) = (split /;/, $line)[32,33,36];
	next if($qual eq ".");
        next if($pos >1010 && $qual <100);
	$type =~ s/TYPE=//;
	$cigar =~ s/CIGAR=//;
	$len =~ s/LEN=//;
       	if($type =~ /del/){
		#print "$hla\t$pos\t$qual\t$type\t$cigar\t$len\n";
	if($cigar =~ /,/){
		my @cigars = split /,/, $cigar;
		my @lens = split /,/, $len;
		for(my $i=0; $i<=$#lens; $i++){
			next if($lens[$i] <50);
			$cigars[$i] =~ /^(\d+)M/;
			my $start = $pos + $1;
			my $end = $start + $lens[$i];
			$region = "$hla\t$start\t$end";
                        $hash{$region} = "$hla\t$start\t"."+"."\t$hla\t$end\t"."+"."\t".".\n";
			#	print "$region\tvcf\n";
			$i = 10;
		}

	}
	else{
	     next if($len < 50);
	     $cigar =~ /^(\d+)M/;
             my $start = $pos + $1;
	     my $end = $start + $len;
	     $region = "$hla\t$start\t$end";
	     $hash{$region} = "$hla\t$start\t"."+"."\t$hla\t$end\t"."+"."\t".".\n";
     }
     # print "$region\t$len\n";
     }

}
foreach my $k1(sort keys %hash){
	print "$k1\n";
}
foreach my $reg(sort keys %hash){
	my ($hla,$start,$end) = (split /\t/,$reg)[0,1,2];
	next if(! $hash{$reg} );
	for(my $i=0; $i<=5; $i++){
		for(my $j=0; $j<=5; $j++){
		       next if(($i==0)&($j==0));
		       my $ns = $start + $i;
		       my $ns2 = $start -$i;
		       my $ne = $end + $j;
		       my $ne2 = $end - $j;
		       my $nreg1 = "$hla\t$ns\t$ne";
		       my $nreg2 = "$hla\t$ns\t$ne2";
		       my $nreg3 = "$hla\t$ns2\t$ne";
		       my $nreg4 = "$hla\t$ns2\t$ne2";
		       if(exists $hash{$nreg1}){$hash{$nreg1} = "";}
		       if(exists $hash{$nreg2}){$hash{$nreg2} = "";}
		       if(exists $hash{$nreg3}){$hash{$nreg3} = "";}
		       if(exists $hash{$nreg4}){$hash{$nreg4} = "";}
	        }
	}
}
foreach my $key(sort keys %hash){
       if( $hash{$key}){
	       my ($hla,$s,$e) = (split /\t/, $key)[0,1,2];
	       my $seq = (split /\t/, $hash{$key})[6];
	       my $len = length($seq);
	       next if(($hla eq "HLA_DRB1") && ($s >10000) && ($len>1000));
	       next if(($hla eq "HLA_DRB1") && ($s > 3800) && ($s < 4400));
	       next if(($hla eq "HLA_DRB1") && ($e > 3800) && ($e < 4400));
	       next if(($hla eq "HLA_DRB1") && ($s == 12144));
	       next if(($s == $e) && ($e<1100));
	       next if(($hla eq "HLA_DQA1") && ($s == $e) && ($e =~ /^74/));
	       next if(($hla eq "HLA_DQB1") && ($s == $e) && ($e =~ /^84/));
               next if(($hla eq "HLA_DPB1") && ($s == $e) && ($e =~ /^124/));
	       next if(($hla eq "HLA_A") && ($s == $e) && ($e =~ /^45/));
	       next if(($hla eq "HLA_B") && ($s == $e) && ($e =~ /^50/));
	       next if(($hla eq "HLA_C") && ($s == $e) && ($e =~ /^53/));
               next if(($hla eq "HLA_DPA1") && ($s == $e) && ($e =~ /^107/));
               next if(($hla eq "HLA_DRB1") && ($s == $e) && ($e =~ /^122/));
               next if(($hla eq "HLA_A") && ($e>4510));
               next if(($hla eq "HLA_B") && ($e>5100));
	       next if(($hla eq "HLA_C") && ($e>5310));
	       next if(($hla eq "HLA_DPA1") && ($e>10780));
	       next if(($hla eq "HLA_DPB1") && ($e>12470));
	       next if(($hla eq "HLA_DQA1") && ($e>7500));
	       next if(($hla eq "HLA_DQB1") && ($e>8480));
	       next if(($hla eq "HLA_DRB1") && ($e>12230));
	        print OUT "$hash{$key}";
	       if($seq =~ m/[A|T|C|G]/ ){print OUT2 ">$hla"."_"."$s"."_"."$e\n$seq\n";}
       }

}
close OUT;
#`less ./tmp.txt|sort|uniq|sort -k1,1 -k2n,2 >$outfile`;
#`rm -rf ./tmp.txt`;
