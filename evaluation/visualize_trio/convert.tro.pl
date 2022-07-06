#!/usr/bin/perl
my ($vcflist, $output, $sampleinfo) = @ARGV;
my (%hash,%hashr,%hasha,%hashaa,%hashp);
sub uniq {
	my %seen;
        grep !$seen{$_}++, @_;
}


open OUT, ">$output";
print OUT "#Sample\talt_type\tchr\tpos\tref_allele\talt_allele\tgene\tphase\tCN\tM\tm\n";
open LI, "$vcflist" or die "$!\n";
while(my $file=<LI>){
	chomp $file;
	my $sa = (split /\//, $file)[-2];
        open (IN, "gunzip -c $file|") or die "$!\n";
	while(<IN>){
		chomp;
		next if(/^#/);
		my ($gene,$pos, $ref, $alt, $le) = (split /\t/, $_)[0,1,3,4,9];
		my $phase = (split /:/, $le)[0];
		my $type = "SNP";
		my @refs = split //, $ref;
		if($phase =~ /2/){
			my @aalt = split /,/, $alt;
			my @alts1 = split //, $aalt[0];
			my @alts2 = split //, $aalt[1];
			my $len = 0;
			if(($#refs <= $#alts1) && ($#refs <= $#alts2)){$len = $#refs; }
			if(($#alts1 <= $#refs) && ($#alts1 <= $#alts2)){$len = $#alts1 }
			if(($#alts2 <= $#refs) && ($#alts2 <= $#alts1)){$len = $#alts2 }
			for(my $i=0;$i<=$len;$i++){
				my ($p, $r, $a);
				my $ph = $phase;
				$p = $pos + $i;
				$r = $refs[$i];
				next if(($refs[$i] eq $alts1[$i]) && ($alts1[$i] eq $alts2[$i]));

				if($alts1[$i] eq $alts2[$i]){
					$ph = "1|1";
					$a = $alts1[$i];
					
				}elsif($refs[$i] eq $alts2[$i]){
					$ph =~ s/2/0/;
					$a = $alts1[$i];
			           
				}elsif($refs[$i] eq $alts1[$i]){
					$ph =~ s/1/0/;
					$ph =~ s/2/1/;
					$a = $alts2[$i];
				}else{
					$a="$alts1[$i]".","."$alts2[$i]";
					
				}
				$hashr{"$gene\t$p"} = $r;
				$hasha{"$gene\t$p"} .= "$a,";
				my $kk = "$sa\t$gene\t$p";
				$hashp{$kk} = "$r\t$a\t$ph";
				#print OUT "$sa\t$type\t$gene\t$p\t$r\t$a\t$gene\t$ph\t2\t1\t1\n";
				my $key = "$sa\t$gene";
				$hash{$key} += 1;
			}
		}else{
			my @alts = split //, $alt;
                        my $len = $#alts;
			if($#refs <= $#alts){$len = $#refs}
			for(my $i=0;$i<=$len;$i++){
				next if($refs[$i] eq $alts[$i]);
				my $p = $pos + $i;
				$hashr{"$gene\t$p"} = $refs[$i];
				$hasha{"$gene\t$p"} .= "$alts[$i],";
				my $kk = "$sa\t$gene\t$p";
				$hashp{$kk} = "$refs[$i]\t$alts[$i]\t$phase";
				#print OUT "$sa\t$type\t$gene\t$p\t$refs[$i]\t$alts[$i]\t$gene\t$phase\t2\t1\t1\n";
				my $key = "$sa\t$gene";
				$hash{$key} += 1;
			}
			
		}
		#print OUT "$sa\t$type\t$gene\t$pos\t$ref\t$alt\t$gene\t$phase\t2\t1\t1\n";
		#my $key = "$sa\t$gene";
		#$hash{$key} += 1;
	}
	close IN;
}
close LI;

foreach my $k(sort keys %hasha){
	my @alts = split /,/, $hasha{$k};
        my @aas = uniq(@alts);
	my $aa = join(",", @aas);
	$hashaa{$k} = $aa;
}

foreach my $kk(sort keys %hashp){
        my ($sa,$gene,$pos) = (split /\t/, $kk)[0,1,2];
        my ($r,$a,$ph) = (split /\t/, $hashp{$kk})[0,1,2];
        my @alts = split /,/, $hashaa{"$gene\t$pos"};
	if($ph =~ /2/){
		my ($p1,$p2) = (split /\|/,$ph)[0,1];
		my ($a1,$a2) = (split /,/,$a)[0,1];
		my ($k1,$k2);
		for(my $i=0;$i<=$#alts;$i++){
			if($alts[$i] eq $a1){$k1 = $i + 1}
			if($alts[$i] eq $a2){$k2 = $i + 1}
		}
		if($k2<$k1){$a = "$a2,$a1";$ph="$k1"."|"."$k2";}
		
	}
	else{
		my $k = 1;
		for(my $i=0;$i<=$#alts;$i++){
			if($alts[$i] eq $a){$k = $i + 1;}
		}
		if($k == 2){
		      $ph =~ s/1/$k/g;
		      $a = "$alts[0],$alts[1]";
	        }
		
	}
	print OUT "$sa\tSNP\t$gene\t$pos\t$r\t$a\t$gene\t$ph\t2\t1\t1\n";
      

}

open OUT2, ">$sampleinfo";
print OUT2 "gene\tsample\tphased\tmut_num\tmut_density\n";
foreach my $key(sort keys %hash){
	my ($sa,$gene) = (split /\t/,$key)[0,1];
	my $d = 0;
	my $m = $hash{$key};
	if($gene eq "HLA_A"){$d = $m /3503}
	if($gene eq "HLA_B"){$d = $m /4081}
	if($gene eq "HLA_C"){$d = $m /4304}
	if($gene eq "HLA_DPA1"){$d = $m /9775}
	if($gene eq "HLA_DPB1"){$d = $m /11468}
	if($gene eq "HLA_DQA1"){$d = $m /6492}
	if($gene eq "HLA_DQB1"){$d = $m /7480}
	if($gene eq "HLA_DRB1"){$d = $m /11229}
	print OUT2 "$gene\t$sa\tTRUE\t$hash{$key}\t$d\n";
}
close OUT2;


