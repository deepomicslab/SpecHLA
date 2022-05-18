#!/usr/bin/perl
my $infile = shift;
my $outfile = `basename $infile`;
my %hash;
open IN, "/home/wangmengyao/scripts/NeedleHLA/script/HLA_FREQ_HLA_I_II.txt" or die "$!\n";
<IN>;
while(<IN>){
	chomp;
	my ($hla,$a,$b,$c) = split;
	my $t = $a+$b+$c;
	next if($t==0);
	$hash{$hla} = $t;
}
close IN;
my %hasha;
my $key;
open CIN, "$infile" or die "$!\n";
while(<CIN>){
	chomp;
        if(/^>/){my ($id1,$id2) = (split /\s/,$_)[0,1];$key = "$id1\t$id2";}
	else{$hasha{$key} .= $_;}
}
close CIN;

open OUT, ">$outfile";
foreach my $key(sort keys %hasha){
	my ($id1,$id2) = (split /\t/,$key)[0,1];
	my @arr = (split /:/,$id2);
	my $k = "$arr[0]".":"."$arr[1]";
	next if (! exists $hash{$k});
	print OUT "$key\n$hasha{$key}\n";
}
close OUT;

