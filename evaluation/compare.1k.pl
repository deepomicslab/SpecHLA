#!/usr/bin/perl -w

my ($software,$infile, $outfile) = @ARGV;
my $k=3;

my @hlas = ("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1");
my (%hash1,%hash2,%hashs);
open OUT, ">$outfile";
print OUT "Gene\tRightCount\tTotal\t$software\n";
open IN, "../20181129_HLA_types_full_1000_Genomes_Project_panel.txt" or die "$!\n";
<IN>;
while(<IN>){
	chomp;
	s/\*//g;
	my @arrs = (split /\t/, $_);
	shift @arrs; shift @arrs;
	my $sa = shift @arrs;
	my @tts = ("A","B","C","DQB1","DRB1"); 
	for(my $i=0;$i<5;$i++){
		my $j = 2 * $i;
		my $gene = $tts[$i];
		my $hla1 = $arrs[$j];
		my $hla2 = $arrs[$j+1];
		my $key = "$sa\t$gene";
		next if(!$hla1 || !$hla2);
		my $b1 = (split /:/,$hla1)[0];
		my $b2 = (split /:/,$hla2)[0];
		if($hla1 =~ /\//){my @hlas = (split /\//,(split /:/,$hla1)[1]);$hla1=""; foreach my $tt(@hlas){$hla1 .="$b1".":"."$tt\t"; }}
		if($hla2 =~ /\//){my @hlas = (split /\//,(split /:/,$hla2)[1]);$hla2=""; foreach my $tt(@hlas){$hla2 .="$b2".":"."$tt\t"; }}
		$hla1 =~ s/\t$//; $hla2 =~ s/\t$//;
		$hash1{$key} = $hla1;
		$hash2{$key} = $hla2;
	}
}
close IN;

my @ids = ("A_1", "A_2 ", "B_1", "B_2", "C_1","C_2", "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2", "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2", "DRB1_1", "DRB1_2");
my (%hashc,%hasha,%hashb,%hasht);
open CI, "$infile" or die "$!\n";
while(<CI>){
	chomp;
	next if(/Sample/ && /A_1/);
	my @arrs = (split /\t/,$_);
        my $sa = shift @arrs;
	$hashs{$sa} = 1;
        for(my $i=0;$i<=$#arrs;$i++){
                my $id = $ids[$i];
                my $hla = $arrs[$i];
		next if(!$hla || $hla eq "-");
                my ($gene,$k) = (split /_/,$id)[0,1];
                $hla =~ s/[;|,|\/]/\t/g;
		my $nhla = "";
		my @hlas = split /\t/,$hla;
		foreach my $hh(@hlas){
			my @hhs = (split /:/, (split /\*/,$hh)[-1]);
			my $thla = "$hhs[0]".":"."$hhs[1]";
			$nhla .= "$thla\t";
		}
		my $key = "$sa\t$gene";
                if($id =~ /_1/){$hasha{$key} = $nhla}
                else{$hashb{$key} = $nhla}
        }
}
close CI;


foreach my $sa(sort keys %hashs){

	foreach my $gene(@hlas){
		my $key = "$sa\t$gene";
		my $m = 0;
		next if(!exists $hash1{$key});
		$hasht{$gene} += 2;
		next if(!exists $hasha{$key} && !exists $hashb{$key});
		my ($hla11,$hla12,$hla21,$hla22) = ($hash1{$key}, $hash2{$key}, $hasha{$key}, $hashb{$key});
		#print "$key\ttruth\t$hla11\t$hla12\ttyping\t$hla21\t$hla22\n";	
		#if(($hla11 eq $hla12) && ($hla21 ne $hla22)){
		#        if(($hla11 eq $hla21) || ($hla11 eq $hla22)){$m+=1}
		#}else{
                #        if(($hla11 eq $hla21) || ($hla11 eq $hla22)){$m+=1;}
		#        if(($hla12 eq $hla21) || ($hla12 eq $hla22)){$m+=1;}
		#}
		my @hhs11 = (split /\t/,$hla11);
		my @hhs12 = (split /\t/,$hla12);
		$hla21 =~ s/[;|,|\/]/\t/g;
		$hla22 =~ s/[;|,|\/]/\t/g;
		my @hhs21 = split /\t+/,$hla21;
		my @hhs22 = split /\t+/,$hla22;
		if($#hhs21 >= $k){splice (@hhs21,$k)}
		if($#hhs22 >= $k){splice (@hhs22,$k)}
		my (@isect1,@isect2,@isect3,@isect4);
		my %hhs11 = map{$_=>1} @hhs11; 
		my %hhs12 = map{$_=>1} @hhs12;
		my %hhs21 = map{$_=>1} @hhs21;
		my %hhs22 = map{$_=>1} @hhs22;
		@isect1 = grep($hhs11{$_},@hhs21);
		@isect2 = grep($hhs12{$_},@hhs22);
		@isect3 = grep($hhs11{$_},@hhs22);
		@isect4 = grep($hhs12{$_},@hhs21);
		#print "11\t$#isect1\t$#isect2\t10\t$#isect3\t$#isect4\t";
		if($#isect1 >=0 && $#isect2>=0){$m=2}
		elsif($#isect3 >=0 && $#isect4>=0){$m=2}
		elsif($#isect1 >=0 || $#isect2>=0 || $#isect3>=0 || $#isect4>=0){$m=1}
		else{$m=0}
		#if($m != 2 && $gene eq "A"){print "$key\ttruth\t$hla11\t$hla12\ttyping\t$hla21\t$hla22\n";}	
		$hashc{$gene} += $m;
	}
}


foreach my $gene(sort keys %hashc){
	my $count = $hashc{$gene};
	my $total = $hasht{$gene};
	my $fre = sprintf("%.3f", ($count / $total));
	print OUT "$gene\t$hashc{$gene}\t$total\t$fre\n"
}
close OUT;
