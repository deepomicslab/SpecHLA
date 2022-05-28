#!/usr/bin/perl -w
use FindBin qw($Bin);

my ($dir) = @ARGV;
my $db="$Bin/../db/ref";
my %hash;
open LI, "$db/DRB1_dup_extract_ref.fasta.fai" or die "$!\tfile1\n";
#open LI, "/mnt/disk2_workspace/wangmengyao/NeedleHLA/GA_rich/DRB1/bwa/DRB1_dup_extract_ref.fasta.fai" or die "$!\n";
while(<LI>){
	chomp;
	my ($hla, $len) = (split /\t/, $_)[0,1];
        $hash{$hla} = $len;
}
close LI;
my %hash1;
my %hashfa;
open CI, "$db/DRB1_dup_extract_ref.fasta" or die "$!\tfile2\n";
#open CI, "/mnt/disk2_workspace/wangmengyao/NeedleHLA/GA_rich/DRB1/DRB1_dup_extract.fasta" or die "$!\n";
while(<CI>){
	chomp;
	if(/^>HLA/){
             s/^>//;
	     my @temp = (split /\t/, $_);
	     my @arrs;
	     foreach my $te(@temp){
		     my $tt = substr($te,0,12);
		     push @arrs, $tt;
	     }
	     my $k = $arrs[0];
	     my $line = join(";",@arrs);
	     $hash1{$k} = $line;
	     my $seq = <CI>;
	     chomp $seq;
	     $hashfa{$k} = $seq;
     }    
}
close CI;
my %hashl;
my %hashr;
my (%hasha,%hashb);
open IN, "$dir/extract.read.blast" or die "$!\tfile3\n";
while(<IN>){
	chomp;
	my ($id,$hla,$score) = (split)[0,1,2];
	if(exists $hasha{$id}){
		if($score > $hasha{$id}){$hashb{$id} = $_;$hasha{$id} = $score}
	}else{$hasha{$id} = $score; $hashb{$id} = $_}
}
close IN;
foreach my $kk(sort keys %hashb){
	my @temp = split /\t/,$hashb{$kk};
	#print "$hashb{$kk}\n";
	next if($temp[2] < 95);
	#next if($temp[3] < 100);
	next if($temp[4] >4);
	next if($temp[5] >2);
        my $reg = "$temp[8]"."_"."$temp[9]";
	$hashl{$temp[1]} .= "\t$reg";
	$hashr{$temp[1]} .= "$temp[0]".";";
}

open OUT, ">$dir/DRB1.hla.count";
foreach my $hla(sort keys %hash){
	if(!exists $hashl{$hla}){
		#print "$hla\t$hash{$hla}\t0\n";
	}
	else{ 
              my @arrs = (split /\t/, $hashl{$hla});
	      shift @arrs;
	      my $count=0;
	      my $len = $hash{$hla};
	      my $a = 0;
	      for(my $i=1; $i<=$len;$i++){
		      my $tag = 0;
		      foreach my $reg(@arrs){
			      my ($s,$e) = (split /_/, $reg)[0,1];
			      if(($s<=$i) && ($i<=$e)){$tag = 1; $a+=1;}
		      }
		      if($tag == 1){$count += 1;}
	      }
	      my $fre = $count / $len;
	      my $depth = $a / $len;
	      $fre = sprintf "%.2f", $fre;
	      $depth = sprintf "%.2f", $depth;
	     
	      print OUT "$hla\t$len\t$fre\t$depth\t$hash1{$hla}\t$hashfa{$hla}\t$hashr{$hla}\n";

	}
}
close OUT;
