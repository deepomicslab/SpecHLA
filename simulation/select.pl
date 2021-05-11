#!/use/bin/perl
my %hash;

open IN, "./hla_gen.format.filter.fasta.fai" or die "$!\n";
while(<IN>){
      chomp;
      my $hla = (split)[0];
      my $class = (split /_/,$hla)[0];
      $hash{$class} .= "$hla\t";
      
}
close IN;

foreach my $key(sort keys %hash){
	my @arr = split /\t/,$hash{$key};
        pop @arr;
        my $s = $#arr+1;
        my $i = int(rand($s));
        my $j = int(rand($s));
        my ($hla1,$hla2) = ($arr[$i],$arr[$j]);
	print "$key\t$hla1\n$key\t$hla2\n";
}
