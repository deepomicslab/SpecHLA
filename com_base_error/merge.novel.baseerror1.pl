#!/usr/bin/perl -w

my %hash;
my @hlas = ("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1");
#open IN, "/mnt/disk2_workspace/wangmengyao/NeedleHLA/simu_data/simu_20201116/sample.fq.list" or die "$!\n";
open IN, "./sample.list" or die "$!\n";
open OUT, ">merge.base.error.txt1";
my $head = join("\t",@hlas);
print OUT "Sample\t$head\n";
while(<IN>){
	chomp;
	my $sample = (split)[0];
        print OUT "$sample";
        foreach my $hla(@hlas){
	     my $vcf1 = `ls $sample/HLA_$hla\_0.*.blast.txt |head -1`;
	     my $vcf2 = `ls $sample/HLA_$hla\_0.*.blast.txt |tail -1`;
	     chomp $vcf1;
	     chomp $vcf2;
	     my $e;
	     
	     if(-e $vcf1){
		     #my $ee1 = `less $vcf1 |grep base_error|cut -d " " -f 6`;
		     #my $ee2 = `less $vcf2 |grep base_error|cut -d " " -f 6`;
		     #$ee1 =~ m/(\d+):/;
		     #my $e1 = $1;
		     #$ee2 =~ m/(\d+):/;
		     my ($ee1,$ee2);
		     open V1, "$vcf1";
		     while(<V1>){
			     chomp;
			     next if(/^#/);
			     my ($l,$t,$i) = (split)[3,4,5];
			     next if($l<400);
			     $ee1 += $t+$i;
		     }
		     close V1;
		     open V2, "$vcf2";
		     while(<V2>){
			     chomp;
			     next if(/^#/);
			     my ($l,$t,$i) = (split)[3,4,5];
			     next if($l < 400);
			     $ee2 += $t+$i;
		     }
		     close V2;
	             
	             if($ee1 <= $ee2){$e = $ee1}else{$e = $ee2}
	             if(($ee1 ==0 ) && ($ee2 ==0)){print OUT "\t-";}
	             else{print OUT "\t$e";}
              }else{print OUT "\t-"}
      }
       print OUT "\n";
}
close IN;
close OUT;

