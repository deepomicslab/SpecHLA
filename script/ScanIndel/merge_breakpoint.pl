#!/usr/bin/perl -w

my ($fvcf, $outfile, $outfa) = @ARGV;
my %hash;
my $region;
open OUT2, ">$outfa";
open OUT, ">$outfile";

open IN3, "$fvcf" or die "$!\trealign.vcf\n";
while(<IN3>){
        chomp;
        next if(/^#/);
        my ($hla,$pos,$ref,$alt,$qual,$line) = (split /\t/,$_)[0,1,3,4,5,7];
        my ($type,$cigar,$len) = (split /;/, $line)[40,6,13];
        next if($qual <100);
        $type =~ s/TYPE=//;
        $cigar =~ s/CIGAR=//;
        $len =~ s/LEN=//;
	my $alen = length($ref);
        if($cigar =~ /,/){
                my @cigars = split /,/, $cigar;
                my @lens = split /,/, $len;
		my @alts = split /,/, $alt;
                for(my $i=0; $i<=$#lens; $i++){
                        next if($lens[$i] <50);
			my $leni = length($alts[$i]);
			next if ($leni < 30);
			my $p1 = substr($alts[$i],0,$alen);
			my $p2 = substr($alts[$i],-$alen);
			#	print "$p1\t$p2\t$ref\t$_\n";
			if($p1 eq $ref){
				my $pose = $pos + $alen;
				my $seq = substr($alts[$i], $alen-1);
				print OUT "$hla\t$pose\t"."+"."\t$hla\t$pose\t"."+"."\t$seq\t$seq\n";
				print OUT2 ">$hla"."_"."$pose"."_"."$pose\n$seq\n";
			}
			if($p2 eq $ref){
                                my $pose = $pos;
				my $seq = substr($alts[$i], 0, $leni - $alen);
				print OUT "$hla\t$pose\t"."+"."\t$hla\t$pose\t"."+"."\t$seq\t$seq\n";
				print OUT2 ">$hla"."_"."$pose"."_"."$pose\n$seq\n";
			}
		}
        }
}

close IN3;
close OUT;
