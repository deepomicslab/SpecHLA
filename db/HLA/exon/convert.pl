#!/usr/bin/perl
my $infile = "/mnt/disk2_workspace/wangmengyao/NeedleHLA/IMGTHLA/xml/HLA-DRB1.exon.fa";
my $outfile = "HLA_DRB1.fasta";
my %hash;
open IN, "/home/wangmengyao/scripts/NeedleHLA/script/HLA_FREQ_HLA_I_II.txt" or die "$!\n";
<IN>;
while(<IN>){
	chomp;
	my ($hla,$a,$b,$c) = split;
	my $t = $a+$b+$c;
	if($hla =~ /DRB1\*01:23/){$hash{$hla} = $t}
	next if($t==0);
	$hash{$hla} = $t;
}
close IN;
my %hashk;
open LI, "/mnt/disk2_workspace/wangmengyao/NeedleHLA/IMGTHLA/Allelelist.txt" or die "$!\n";
while(<LI>){
	chomp;
	if(/^HLA/){
		my ($id,$vv) = (split /,/,$_)[0,1];
		$hashk{$vv} = $id;
	}
}
close LI;

my %hasha;
my $key;
open CIN, "$infile" or die "$!\n";
while(<CIN>){
	chomp;
        if(/^>/){s/^>HLA-//g;$key = $_;}
	else{$hasha{$key} .= $_;}
}
close CIN;

open OUT, ">$outfile";
foreach my $key(sort keys %hasha){
	#my ($id1,$id2) = (split /\t/,$key)[0,1];
	my @arr = (split /:/,$key);
	my $k = "$arr[0]".":"."$arr[1]";
	my $hh = $hashk{$key};
	print "$hh\t$k\t$key\n";
	next if (! exists $hash{$k});
	print OUT ">HLA:$hh\t$key\n$hasha{$key}\n";
}
close OUT;
`/home/wangmengyao/packages/ncbi-blast-2.9.0+/bin/makeblastdb -in HLA_DRB1.fasta -dbtype nucl -parse_seqids -out HLA_DRB1.fasta`;

