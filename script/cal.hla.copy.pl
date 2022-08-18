#!/usr/bin/perl -w
use strict;
use Cwd qw/abs_path/;
use Getopt::Long;


my ($sample, $het_cutoff, $purity, $ploidy, $filelist, $tfile, $outdir);
my $Help;
GetOptions(
	   "purity=f"  =>      \$purity,
           "ploidy=f"  =>      \$ploidy,
	   "F=s"       =>      \$filelist,
           "T=s"       =>      \$tfile,  ## typing result
	   "O=s"       =>      \$outdir,
	   "S=s"       =>      \$sample,
	   "C=f"       =>      \$het_cutoff,
	   "h"         =>      \$Help
);


my @hlas = ("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1");
my %hasha;
open RIN, "$tfile" or die "$!\n";
my $line1 = <RIN>; my $line2 = <RIN>;
chomp $line1; chomp $line2;
my @array1 = split /\t/,$line1;
my @array2 = split /\t/,$line2;
for(my $i=1;$i<=$#array1;$i++){
     my $id = $array1[$i];
     my $k = (split /_/,$id)[1];
     my $alle = $array2[$i];
     $hasha{$k} .= "$alle\t";
}

open OUT, ">$outdir/merge.hla.copy.txt";
print OUT "Sample\tHLA\tAllele1\tAllele2\tcopyratio\tKeptHLA\tLossHLA\tFreq1\tFreq2\tPurity\tHet_num\tLOH\n";

open IN, "$filelist" or die "$!\tfilelist\n";
while(<IN>){
	chomp;
	my $file = $_;
	my $hla = (split /\//,$file)[-1];
	$hla =~ s/_freq\.txt$//;
        $hla =~ s/HLA_//;
	open LI, "$file" or die "$!\tfile\n";
	<LI>;
	my $line1 = <LI>; my $line2 = <LI>; my $line3 = <LI>;
	chomp $line1; chomp $line2; chomp $line3; close LI;
	my $fre1 = (split /\s+/, $line1)[1];
	my $fre2 = (split /\s+/, $line2)[1];
	my $het = (split /\s/,$line3)[-1];
        print OUT "$sample\t$hla\t$hasha{$hla}";
        my $tfre; my $tag=0; my $A=0; my $B=0;
        if($fre1 ==0){$ploidy = sprintf("%.0f",$ploidy); print OUT "0:$ploidy"; $B= $ploidy;}
	if($fre2 ==0){$ploidy = sprintf("%.0f",$ploidy); print OUT "$ploidy:0"; $A = $ploidy;}
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
		      $A = $a; $B = $b;
		     if($fre1 > $fre2){print OUT "$a:$b"; $tag =2}
                     else{print OUT "$b:$a"; $tag=1}
	}
     	my ($hla1,$hla2) = (split /\t/,$hasha{$hla})[0,1];
	my $loh = "Y";
	if($A == $B || $het < $het_cutoff ){$loh = "N"}
        if($tag ==0 ){print OUT "\t$hla1\thomogeneous\t$fre1\t$fre2\t$purity\t$het\tN\n"}
	if($tag ==1 ){print OUT "\t$hla2\t$hla1\t$fre1\t$fre2\t$purity\t$het\t$loh\n"}
        if($tag ==2 ){print OUT "\t$hla1\t$hla2\t$fre1\t$fre2\t$purity\t$het\t$loh\n"}
}
close IN;
close OUT;

sub para_alert{
	my $alert;
	if($Help){
		$alert = 1;
	}
	elsif(!$purity || !$ploidy || !$sample || !$het_cutoff){
		$alert = 1;
	}
	elsif(!$filelist || !-e $filelist){
		$alert = 1;
	}
        elsif(!$tfile || !-e $tfile){
                $alert = 1;
        }
	elsif(!$outdir || !-e $outdir){
                $alert = 1;
        }
        die "
	USAGE:
	perl $0 [Options] -S <samplename> -C <5> -purity <Purity> -ploidy <Ploidy> -F <filelist> -T <hla.result.txt> -O <outdir>
	
	OPTIONS:									     
     	-purity  [f]  the purity of tumor sample. <required>
        -ploidy  [f]  the ploidy of tumor sample in HLA gene region. <required>	
	-S       [s]  sample name <required>
	-C       [f]  the cutoff of heterogeneous snp number. <required>
	-F       [s]  the filelist. <required>     Format, separated by table: ./HLA_A_freq.txt
	-T       [s]  the hla typing result file of Spechla <required>
        -O       [s]  the output dir. <required>
       	\n" if($alert);
}

