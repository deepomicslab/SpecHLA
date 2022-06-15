#!/usr/bin/perl -w
use FindBin qw($Bin);
use Getopt::Long;

my ($sample, $gene, $dir, $pop, $help);

GetOptions(
           "s=s"     =>      \$sample,
           "g=s"     =>      \$gene,
           "i=s"     =>      \$dir,
           "p=s"     =>      \$pop,
           "h"       =>      \$help
);
my $usage = <<USE;
Usage:
description: combination selection
usage: perl $0 [options]
        Common options:
        -s       <tr>    sample name
        -g       <tr>    gene name "HLA_A|HLA_B|HLA_C|HLA_DPA1|HLA_DPB1|HLA_DQA1|HLA_DQB1|HLA_DRB1"
        -i       <tr>    the directory of phased sequence
        -p       <tr>    population information "Asian|Black|Caucasian|Unknown|nonuse"
        -help|?           print help information
e.g.:
        perl $0 -s samplename -g HLA_A -i indir -p Unknown 
USE
die $usage unless ($sample && $dir && $pop && $gene) ;

print "parameter:\tsample:$sample\tdir:$dir\tpop:$pop\tgene:$gene\n";

my $k = 2;
my (%hashp, %hashpp, %hashg, %hashc, %hash,%hashdd);
my $db="$Bin/../../db/HLA";
my $bin="$Bin/../../bin";
my $fadir=$dir;
my $workdir = "$dir/tmp";
`mkdir  -p $workdir`;
#population frequency
open FIN, "$db/HLA_FREQ_HLA_I_II.txt" or die "$!\n";
<FIN>;
while(<FIN>){
    chomp;
    my ($gene,$c,$b,$a) = (split);
    $a = sprintf "%.3f",$a;
    $b = sprintf "%.3f",$b;
    $c = sprintf "%.3f",$c;
    $hashpp{$gene} = "$c;$b;$a";
    if($pop eq "Unknown"){$hashp{$gene} = ($a+$b+$c)/3}
    if($pop eq "Asian"){$hashp{$gene} = $a}
    if($pop eq "Black"){$hashp{$gene} = $b}
    if($pop eq "Caucasian"){$hashp{$gene} = $c}
    if($pop eq "nonuse"){$hashp{$gene} = 0}
}
close FIN;

#Convert HLA nomenclature
open INL, "$db/hla_nom_g.txt" or die "$!\n";
while(<INL>){
        chomp;
        next if(/^#/);
        my ($tag, $line, $hh) = (split /;/, $_)[0,1,2];
        if($hh){
                my $value = "$tag"."$hh";
                my @arrs = (split /\//,$line);
                foreach my $kk(@arrs){
                        my $hla = "$tag"."$kk";
                        $hashg{$hla} = $value;
                }
        }else{
                my $hla = "$tag"."$line";
                $hashg{$hla} = $hla;
        }

}
close INL;
#convertion of HLA id
open IN, "$db/Allelelist.txt" or die "$!\n";
while(<IN>){
        chomp;
        next if(/^#/);
        my ($id, $name) = (split /,/, $_)[0,1];
        my $key = "HLA:"."$id";
        my $hla=$name;
        if(exists $hashg{$name}){$hla=$hashg{$name}}
        $hashc{$key} = $hla;
}
close IN;

my $class = $gene;
my $ref="$db/exon/$class.fasta";
my (%hashs,%hash_max);
my %idks;
`ls $dir/$class.*.fasta >$dir/tfile.list`;
open TE, "$dir/tfile.list" or die "$!\n";
open OUT, ">$dir/$class.total.fasta";
open SOUT, ">$dir/result.$class.fasta";
while(my $file=<TE>){
                chomp $file;
                my $id = (split /\//,$file)[-1];
                $id =~ s/\.fasta//;
                print OUT ">$id\n";
                my @tes = (split /\./,$id); pop @tes;my $idk = join(".",@tes);$idks{$idk}=$idk;
                my $seq;
                open TT, "$file" or die "$!\n";
                while(my $line=<TT>){
                        chomp $line;
                        next if($line=~/^>/);
                        $seq .= $line;
                }
                close TT;
                print OUT "$seq\n";
                $hashs{$id} = $seq;
        }
        close TE;
        close OUT;
        `$bin/blastn -query $dir/$class.total.fasta -out $dir/$class.blast.out1 -db $ref -outfmt 6 -num_threads 4 -max_target_seqs 10 `;
        open BIN, "$dir/$class.blast.out1" or die "$!\n";
        my (%hash1, %hash2,%hash3, %hashas, %hashhk, %hash4);
        while(<BIN>){
                chomp;
                my ($id, $hla, $t, $m, $i, $si) = (split)[0,1,3,4,5,11];
                next if($hla eq "HLA:HLA10778");
                $hash1{$id}{$hla} += $m + $i;
                $hash2{$id}{$hla} += $t - $i;
                $hash4{$id}{$hla} += $si;
        }
        foreach my $id (sort keys %hash1){
                my %hashs; my $mscore=0;
                foreach my $hla(sort keys %{$hash1{$id}}){
                                my $mis = $hash1{$id}{$hla};
                                my $len = $hash2{$id}{$hla};
                                my $score = 100 * (1 - $mis/$len);
                                if($class =~ /DRB1/){$score = $hash4{$id}{$hla} / 499}
                                my $nhla = $hashc{$hla};
                                my @tt = (split /:/, $nhla);
                                my $kid = "$tt[0]".":"."$tt[1]";
                                my $fre=0;
                                if(exists $hashp{$kid}){$fre=$hashp{$kid}}
                                #my $scorel = $score * (2 ** ($fre/20));
                                my $tk = "$id\t$nhla";
                                $hashhk{$tk} = "$score\t$fre";
                                $hash3{$id}{$nhla} = $score;
                                if($pop eq "nonuse" && $score >= $mscore){$mscore = $score}
                                if(($score >= $mscore) && ($fre > 0) && $pop ne "nonuse"){$mscore = $score}
                                $hashs{$mscore} .= "\t$hla";
                                $hashas{$id} .= "$nhla\t";
                }
                $hashh{$id} = "$mscore\t$hashs{$mscore}";
        }
        my %hashss;
        my $mscore = 0;
        foreach my $idk (sort keys %idks){
                 my $ts = 0;
                 for(my $i=1;$i<=$k;$i++){
                         my $idd = "$idk"."."."$i";
                         my $s = (split /\t/, $hashh{$idd})[0];
                         $ts += $s;
                 }
                 $ts = $ts / $k;
                 print "$idk\t$ts\n";
                 if($ts >= $mscore){$mscore = $ts}
                 $hashss{$ts} .= "$idk\t";
         }
         my @tids = (split /\t/, $hashss{$mscore});
         my (%fhash,%thash); my $ffre=0;
         foreach my $tid(@tids){
                 my $tfre;
                 for(my $i=1;$i<=$k;$i++){
                         my $hkey = "$class"."_"."$i";
                         my $idd = "$tid"."."."$i";
                         my @arrs = (split /\t+/, $hashh{$idd});        
                         my $aa = shift @arrs;
                         my %nhash; my $mfre = 0;
                         foreach $od(@arrs){
                                 my $nd = $hashc{$od};
                                 my @tt = (split /:/, $nd);
                                 my $id = "$tt[0]".":"."$tt[1]";
                                 my $fre=0;
                                 if(exists $hashp{$id}){$fre=$hashp{$id}}
                                 $nhash{$fre} .= "$nd\t";
                                 #               print "$idd\t$od\t$nd\t$fre\n";
                                 if($fre >= $mfre){$mfre = $fre}
                         }
                         $tfre += $mfre;
                         if($mfre == 0 && $pop ne "nonuse"){$tfre = $tfre/2;}
                         $thash{$idd} = $nhash{$mfre};
                         #print "$idd\t$tfre\t$aa\n";

                 }
                 if($tfre >=$ffre){$ffre = $tfre}
                 $fhash{$tfre} = $tid;
         }
         my $select_id = $fhash{$ffre};
         print "Select_id\t$select_id\n";
         for($i=1;$i<=$k;$i++){
                 my $tag = "$class"."_"."$i";
                 my $idd = "$select_id"."."."$i";
                 my $seq = $hashs{$idd};
                 print SOUT ">$idd\n$seq\n";

         }
`rm -rf $dir/$class.total.fasta`;
 close SOUT;
     


