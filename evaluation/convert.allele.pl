#!/usr/bin/perl -w

# map the HLA types to the most recent IMGT version

my ($infile,$outfile) = @ARGV;
my %hash;
open IN, "../Allelelist_history.txt" or die "$!\n"; # https://github.com/ANHIG/IMGTHLA/Allelelist_history.txt
while(<IN>){
        chomp;
        next if(/^#/);
        next if(/^HLA_ID/);
        my @arrs = (split /,/,$_);
        my $id = shift @arrs;
        my $hla = $arrs[0];
        next if($hla eq "NA");
	foreach my $aa(@arrs){
                $hash{$aa} = $hla;
        }
}
close IN;

open CI, "$infile" or die "$!\n";
open OUT, ">$outfile";
while(<CI>){
        chomp;
        if(/^Sample/){print OUT "$_\n"}
        else{
                my @arrs = (split);
                my $sa = shift @arrs;
                my $out = $sa;
                foreach my $aa(@arrs){
                        my $nn = $aa;
                        if(exists $hash{$aa}){$nn = $hash{$aa};} # if($aa ne $nn){print "$nn\t$aa\n"}}
                        else{print "$aa\n"}
                        $out .= "\t$nn";
                }
                print OUT "$out\n";
        }
}
close CI;
close OUT;


