#!/usr/bin/perl -w
use strict;
use Cwd qw/abs_path/;
use Getopt::Long;


my ($purity, $ploidy, $filelist, $outdir);
my $Help;
GetOptions(
	   "purity=f"  =>      \$purity,
           "ploidy=f"  =>      \$ploidy,
	   "F=s"       =>      \$filelist,
	   "O=s"       =>      \$outdir,
	   "h"         =>      \$Help
);


#my ($purity, $ploidy) = (0.4, 1);
# my $filelist:  Sample\tHLA_A\tHLA_A_freq.txt\n
#
my @hlas = ("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1");

open OUT, ">$outdir/merge.hla.copy.txt";
#my $head = join("\t",@hlas);
print OUT "Sample\tHLA\tcopy_number\n";

open IN, "$filelist" or die "$!\tfilelist\n";
while(<IN>){
	chomp;
	my $file = $_;
	my ($sample, $hla) = (split /\//,$file)[-2,-1];
	$hla =~ s/_freq\.txt$//;
	open LI, "$file" or die "$!\tfile\n";
	<LI>;
	my $line1 = <LI>; my $line2 = <LI>;
	chomp $line1; chomp $line2;
	my $fre1 = (split /\s+/, $line1)[1];
	my $fre2 = (split /\s+/, $line2)[1];
        print OUT "$sample\t$hla";
        my $tfre;
	if($fre1 ==0 || $fre2 ==0){$ploidy = sprintf("%.0f",$ploidy); print OUT "\t$ploidy:0";}
	else{
		     if($fre1 > $fre2){$tfre = $fre1 / $fre2}
		     else{$tfre = $fre2 / $fre1}
                
		     my ($a,$b);
                     $a = ($tfre * $ploidy + ($tfre -1) * (1/$purity - 1)) / ( 1 + $tfre );
		     #$a = sprintf("%.0f", $a);
	             $b = $ploidy - $a;
		     if($b <0) { $b=0; $a = $ploidy}
		     #if($b <0){ 
        	     #        $b = ($ploidy + (1/$purity - 1) * (1 - $tfre)) / ( 1 + $tfre);
		     #        $b = sprintf("%.0f", $b);
		     #        $a = $ploidy - $b;
		     #}
		     $a = sprintf("%.0f",$a);
		     $b = sprintf("%.0f",$b);
		     print OUT "\t$a:$b";
	}
	
	print OUT "\n";
}
close IN;
close OUT;

sub para_alert{
	my $alert;
	if($Help){
		$alert = 1;
	}
	elsif(!$purity || $ploidy){
		$alert = 1;
	}
	elsif(!$filelist || !-e $filelist){
		$alert = 1;
	}
	elsif(!$outdir || !-e $outdir){															  $alert = 1;																    }
        die "
	USAGE:
	perl $0 [Options] -purity <Purity> -ploidy <Ploidy> -F <filelist> -O <outdir>
	
	OPTIONS:									     
     	-purity  [f]  the purity of tumor sample. <required>
        -ploidy  [f]  the ploidy of tumor sample in HLA gene region. <required>	
	-F       [s]  the filelist. <required>															     Format, separated by table:  Sample/HLA_A_freq.txt				
	-O       [s]  the output dir. <required>
       	\n" if($alert);
}

