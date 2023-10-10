#!/usr/bin/perl -w
##step1 download the latest IMGT-HLA database
`git clone https://github.com/ANHIG/IMGTHLA.git --depth 1`;

my $infile = "IMGTHLA/hla_gen.fasta";
my $infile2 = "IMGTHLA/hla_nuc.fasta";
my $db_dir = "../db";
my $bin = "./";

`cp IMGTHLA/wmda/hla_nom_g.txt $db_dir/HLA`;
`cp IMGTHLA/Allelelist.txt $db_dir/HLA`;

my (%hash,%hashe,$id);
open IN, "$infile" or die "$!\tplease input IMGTHLA/hla_gen.fasta\n";
while(<IN>){
	chomp;
	if(/^>/){$id = (split /\s+/,$_)[1];}
	else{$hash{$id} .= $_;}
}
close IN;

open IN1, "$infile2" or die "$!\tplease input IMGTHLA/hla_nuc.fasta\n";
while(<IN1>){
        chomp;
        if(/^>/){$id = $_}
        else{$hashe{$id} .= $_;}
}
close IN1;

my @hlas = ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1");
foreach my $hla(@hlas){
	open OUT, ">$db_dir/HLA/whole/HLA_$hla.fasta";
	foreach my $key (sort keys %hash){
		my $hh = (split /\*/,$key)[0];
		if($hh eq $hla){print OUT ">$key\n$hash{$key}\n"}
	}
	close OUT;
	`samtools faidx $db_dir/HLA/whole/HLA_$hla.fasta`;
	`makeblastdb -in $db_dir/HLA/whole/HLA_$hla.fasta -dbtype nucl -parse_seqids -out $db_dir/HLA/whole/HLA_$hla`;

}

foreach my $hla(@hlas){
        open OUT, ">$db_dir/HLA/exon/HLA_$hla.fasta";
        foreach my $key (sort keys %hashe){
                my $hh = (split /\*/, (split /\s+/,$key)[1])[0];
                if($hh eq $hla){print OUT "$key\n$hashe{$key}\n"}
        }
        close OUT;
        `samtools faidx $db_dir/HLA/exon/HLA_$hla.fasta`;
        `makeblastdb -in $db_dir/HLA/exon/HLA_$hla.fasta -dbtype nucl -parse_seqids -out $db_dir/HLA/exon/HLA_$hla.fasta`;

}
`cp $db_dir/HLA/exon/HLA_DRB1.fasta $db_dir/HLA/whole/HLA_DRB1.exon.fasta`;
`samtools faidx $db_dir/HLA/whole/HLA_DRB1.exon.fasta`;
`makeblastdb -in $db_dir/HLA/whole/HLA_DRB1.exon.fasta -dbtype nucl -parse_seqids -out $db_dir/HLA/whole/HLA_DRB1.exon`;



