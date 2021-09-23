#!/usr/bin/perl -w
my (%hash1, %hash2);
open IN, "merge.base.error.txt1" or die "$!\n";
my $head = <IN>;
chomp $head;
my @hlas = split /\t/, $head;
shift @hlas;
while(<IN>){
	chomp;
	my $sa = (split)[0];
	$hash1{$sa} = $_;
}
close IN;


open IN1, "merge.base.error.txt2" or die "$!\n";
<IN1>;
while(<IN1>){
	chomp;
	my $sa = (split)[0];
	$hash2{$sa} = $_;
}
close IN1;

my %hash = ('A'=>3503,'B'=>4081,'C'=>4304,'DPA1'=>9775,'DPB1'=>11468,'DQA1'=>6492,'DQB1'=>7480,'DRB1'=>11229);

open OUT, ">merge.base.error.fre.txt";
print OUT "$head\n";
foreach my $sa(sort keys %hash1){
	my @arr1 = split /\t/, $hash1{$sa};
	my @arr2 = split /\t/, $hash2{$sa};
	shift @arr1;
	shift @arr2;
	print OUT "$sa";
	for(my $i=0;$i<8;$i++){
		my $e1 = $arr1[$i];
		my $e2 = $arr2[$i];
		my $hla = $hlas[$i];
		my $len = $hash{$hla};
		my $fre;
		if($e1 =~ /-/ || $e2 =~ /-/){$fre = "-"}
		else{$fre = ($e1+$e2)/(2 * $len);}
	        print OUT "\t$fre";
	}
	print OUT "\n";
}
close OUT;

