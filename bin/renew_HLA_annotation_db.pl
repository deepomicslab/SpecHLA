#!/usr/bin/perl -w
##step1 download the latest IMGT-HLA database
`git clone https://github.com/ANHIG/IMGTHLA.git --depth 1`;

my $infile = "IMGTHLA/hla_gen.fasta";
my $infile_zip = "IMGTHLA/hla_gen.fasta.zip";
my $infile2 = "IMGTHLA/hla_nuc.fasta";
my $db_dir = "../db";
my $bin = "./";

### check if $infile exists, if not, check if $infile_zip exists, and unzip $infile_zip
if(!-e $infile){
	if(-e $infile_zip){
		`unzip $infile_zip -d IMGTHLA/`;
	}else{
		die "$!\tno IMGTHLA/hla_gen.fasta found in IMGTHLA/\n";
	}
}

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

## In whole mode, DRB1 is annotated against exon 2 only (whole/HLA_DRB1.exon), named by allele.
## Do not use the full CDS here: the exon-2 alignment of a full-length record extends a few
## bases into exon 1 and collects mismatches, so partial records would outscore the true allele.
## Exon boundaries are only given in hla.dat.
my $dat = "IMGTHLA/hla.dat";
my $dat_zip = "IMGTHLA/hla.dat.zip";
if(-e $dat){open DAT, "$dat" or die "$!\tcannot open $dat\n";}
elsif(-e $dat_zip){open DAT, "unzip -p $dat_zip |" or die "$!\tcannot unzip $dat_zip\n";}
else{die "no IMGTHLA/hla.dat or IMGTHLA/hla.dat.zip found in IMGTHLA/\n";}

open OUT, ">$db_dir/HLA/whole/HLA_DRB1.exon.fasta" or die "$!\n";
my ($name, $seq, $inseq, $es, $ee, $start, $end, $count) = ("", "", 0, 0, 0, 0, 0, 0);
while(<DAT>){
	if(/^DE\s+HLA-(DRB1\*[^,\s]+)/){$name = $1;}
	next if(!$name); # not a DRB1 record
	if(/^FT   (\S+)\s+(\S+)/){ # feature line; remember the location if it is an exon
		($es, $ee) = (0, 0);
		if($1 eq "exon" && $2 =~ /^<?(\d+)\.\.>?(\d+)$/){($es, $ee) = ($1, $2);}
	}
	elsif(/^FT\s+\/number="2"/ && $es){($start, $end) = ($es, $ee);}
	elsif(/^SQ/){$inseq = 1;}
	elsif(/^\/\//){
		if($start){
			print OUT ">$name\n", uc(substr($seq, $start - 1, $end - $start + 1)), "\n";
			$count += 1;
		}elsif($seq){print STDERR "warning: no exon 2 found for $name, skipped\n";} # deleted alleles have no sequence
		($name, $seq, $inseq, $es, $ee, $start, $end) = ("", "", 0, 0, 0, 0, 0);
	}
	elsif($inseq){s/[\s\d]//g; $seq .= $_;}
}
close DAT;
close OUT;
die "no DRB1 exon 2 sequence extracted from hla.dat\n" if($count == 0);

`samtools faidx $db_dir/HLA/whole/HLA_DRB1.exon.fasta`;
`makeblastdb -in $db_dir/HLA/whole/HLA_DRB1.exon.fasta -dbtype nucl -parse_seqids -out $db_dir/HLA/whole/HLA_DRB1.exon`;



