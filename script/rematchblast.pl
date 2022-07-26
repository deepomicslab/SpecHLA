#!/usr/bin/perl -w
#perl rematch.pl assembly.blast.sam extract.fa rematch.bam.txt rematch.read.format.txt 30
my ($blastsam, $fasta,$rematchfile, $outfile, $m) = @ARGV;
if(!$m){$m=10}
my %hash2;
my %hashe;
open SAM, "$blastsam" or die "$!\n";
while(<SAM>){
        chomp;
	my $id =(split)[0];
        $hashe{$id} .= "$_\n";
}
close SAM;
foreach my $key(sort keys %hashe){
	my @lines = split /\n/, $hashe{$key};
	#pop @lines;
	my $line = shift @lines;
	my $ipos = (split /\t/, $line)[3];
	if($#lines >=1){
	     foreach my $tmp (@lines){
                 my $tpos = (split /\t/, $tmp)[3];
		 if($tpos < $ipos){$line = $tmp; $ipos = $tpos}

	    }
        }
	if($#lines == 0){
		my $line1 = $line;
		my $line2 = shift @lines;
		
		my ($id1,$hla1,$pos1,$cigar1) = (split /\t/, $line1)[0,2,3,5];
		my ($id2,$hla2,$pos2,$cigar2) = (split /\t/, $line2)[0,2,3,5];
		next if($hla1 ne $hla2);
		my ($a,$b) = ($pos1,$cigar1);
		if($pos1 > $pos2){ $pos1 = $pos2; $pos2=$a; $cigar1=$cigar2; $cigar2=$b; }
                my ($cigarstring1, $cigarstring2);
		next if(!$cigar2 =~ /^(\d+)H/);
		$cigar1 =~ s/\d+S$//;
		if($cigar1 =~ /^(\d+)H/){my $base = $1; $cigar1 =~ s/H/M/; $pos1=$pos1-$base ;}
		my $h = 0;
		if($cigar2 =~ s/^(\d+)H//){$h =$1}
		
		my $k = $pos2 - $pos1;
		
		while($cigar1){$cigar1 =~ s/^(\d+)(\w)//; $cigarstring1 .= $2 x $1;}
		while($cigar2){$cigar2 =~ s/^(\d+)(\w)//; $cigarstring2 .= $2 x $1;}
		my $d1 = $cigarstring1 =~ tr/D/D/;
		my $d12 = $cigarstring1 =~ tr/I/I/; 
		my $d2 = $cigarstring2 =~ tr/D/D/;
		
		my $inv = length($cigarstring1) - $d1 - $h - 1;
		my $break = $pos2 + $inv;
		my $end1 = $pos1 + length($cigarstring1) -$d12 -1;
                if($end1 > $pos2){
			my $k1 = $k;
		        my $lx= substr($cigarstring1,0, $k1);
		        my $dx = $lx =~ tr/MI/MI/;
		        if($dx != $k){
				$k1 += $k1 - $dx;
				$lx = substr($cigarstring1, 0, $k1);
				$dx = $lx =~ tr/MI/MI/;
			}
		        my $dk = substr($cigarstring1,$k1) =~ tr/D/D/;
                        my $di = substr($cigarstring1,$k1) =~ tr/I/I/;
		        my $dd =  $pos1 - $pos2 - $d12 + $d1 + $h + 1 + $d1-$dk ; 
			next if(($end1-$pos2) > length($cigarstring2));
			$cigarstring2 = substr($cigarstring2,$end1-$pos2);
		        next if($dd < 50);	
			my $ins = "I" x $dd;
			my $newstring = "$cigarstring1"."$ins"."$cigarstring2";
                        $hash2{$id1} = "$hla1\t$pos1\t$newstring";
		        print "$hla1\t$dd\t$dk\t$di\t$d2\t$id1\t$d1\t$d12\t$pos1\t$break\t$end1\t$pos2\t$newstring\n";
		}
		if($end1 < $pos2){
			my $dd = $pos2 - $end1;
			next if($dd > 100);
			next if($dd < 50);
			my $del = "D" x $dd;
			my $newstring = "$cigarstring1"."$del"."$cigarstring2";
			$hash2{$id1} = "$hla1\t$pos1\t$newstring";
			print "$hla1\t$dd\t$pos1\t$end1\t$pos2\n";
		}

	}
	else{
        my ($id, $hla, $pos, $cigar) = (split /\t/, $line)[0,2,3,5];
	my $tcigar =  $cigar;
	my $base = 0;
	if($tcigar =~ /^(\d+)H/){ 
		$base = $1; 
		$tcigar =~ s/H/M/;
	}
	my $cigarstring;
        while($tcigar){
		$tcigar =~ s/^(\d+)(\w)//;
		my($n, $a) = ($1,$2);
                $cigarstring .= $a x $n;
	}
	$pos = $pos - $base;
	next if($base >$m);
	next if(($cigarstring =~ tr/M/M/) <=100 );
	next if((exists $hash2{$id}));
	
        $hash2{$id} = "$hla\t$pos\t$cigarstring";
       }
	
}


my %hash;
open BI, "$fasta" or die "$!\n";
while(<BI>){
	chomp;
	next if(!/^>/);
	s/^>//;
        my ($id, $cigar) = (split);
	$hash{$id} = $cigar;
}
close BI;

open IN, "$rematchfile" or die "$!\n";
open OUT, ">$outfile";
#print OUT "readid\ttag\tchr\tstart\tcigar\n";
while(<IN>){
	chomp;
	my ($start,$cigar);
	my ($read, $id, $nstart, $ncigar) = split;
        my ($readid, $tag) = (split /##/, $read)[0,1]; # tag=1 read1; tag=1 read2;
	#print "$read\n";	
	next if(!exists $hash2{$id});
	next if($ncigar =~ /I|D/);
	my $oldcigar = $hash{$read};
	next if($oldcigar =~ /\*/);
	my $oldcigarstring;
	while($oldcigar){
                $oldcigar =~ s/^(\d+)(\w)//;
                my($n, $a) = ($1,$2);
                $oldcigarstring .= $a x $n;
        }

	my ($hla, $rstart, $cigarstring) = (split /\t/, $hash2{$id})[0,1,2];
	$ncigar =~ /(\d+)M/;
	my $len = $1;
	next if($nstart + 50 >= ($cigarstring =~ tr/MI/MI/));
	my $lo = substr($cigarstring, 0, $nstart -1) =~ tr/D/D/;

	my $lsr = substr($cigarstring, $nstart-1, $lo) =~ tr/D/D/;
	
	while($lsr){
		my $llsr = $lsr;
		$lsr = substr($cigarstring, $nstart+$lo-1,$lsr) =~ tr/D/D/;
		$lo += $llsr;

	}
	my $el=0;
	if(substr($cigarstring, 0, $nstart -1 ) =~ /(D+)$/){$el = length($1);}
        my $li = substr($cigarstring, 0, $nstart + $lo -1) =~ tr/I/I/;
	my $string;
	
	if(($len+$nstart+$lo-1) > length($cigarstring)){
		$string = substr($cigarstring, $nstart-1+$lo);
		my $b = $len -length($string); 
		$string .= "S" x $b;}
	else{$string = substr($cigarstring, $nstart-1+$lo, $len);}

	$start=$rstart+$nstart-1+$lo - $li;
	#print "$read\t$start\t$lo\t$li\t$string\n";	
	if($string =~ /^(D+)/){ 
		$string = substr($cigarstring, $nstart-1+$lo+length($1), $len); 
		$start=$rstart+$nstart-1+$lo + length($1)-$li;
	}

	my $aa = $string =~ tr/D/D/;
	my $len2 = $len;
        while(((length($string) - $aa) < $len) & (length($cigarstring) > ($len2+$nstart+$lo))){
		$len2 += 1;
		$string = substr($cigarstring, $nstart-1+$lo, $len2);
		$aa = $string =~ tr/D/D/;
		
		#if(($len2+$nstart+$lo) == length($cigarstring)){my $b = $len + $aa -length($string); $string .= "S" x $b;}
	}
	
	if(($len2+$nstart+$lo) >= length($cigarstring)){my $b = $len + $aa -length($string); $string .= "S" x $b;}
	#print "$read\t$start\t$string\n";
	if($ncigar =~ /^(\d+)S/){
		my $a0 = $1;
		my $tmp = substr($oldcigarstring, 0, $a0) ; # prefix $1S need to fill the cigar info of normal mapped reads
		my $a1 = $tmp =~ tr/D/D/;
		if($tmp =~ /S/){my $t = "S" x $a0; $string = "$t"."$string";}
		else{
			my $ntmp = substr($oldcigarstring, 0, $a0+$a1);
		        my $a2 = $ntmp =~ tr/D/D/;
		        while($a2 != $a1){
			      $a1 = $a2;
			      $ntmp = substr($oldcigarstring, 0, $a0+$a1);
                              $a2 = $ntmp =~ tr/D/D/;
		       }
                       $string = "$ntmp"."$string";
		       my $ni = $ntmp =~ tr/I/I/;
		       $start = $start - $a0 - $a1 + $ni;
	       }
	}
	if($ncigar =~ /(\d+)S$/){
                my $a0 = $1;
		my $tmp = substr($oldcigarstring, -$a0);
		my $a1 = $tmp =~ tr/D/D/;
		if($tmp =~ /S/){$string .= "S" x $a0 ;}
		else{
			my $ntmp = substr($oldcigarstring, -($a0+$a1));
		        my $a2 = $ntmp =~ tr/D/D/;
		        while($a2 != $a1){
			      $a1 = $a2;
			      $ntmp = substr($oldcigarstring, -($a0+$a1));
			      $a2 = $ntmp =~ tr/D/D/;
	                } 

		        $string .= $ntmp ;
	       }

	}
        if($string =~ /[M|I|D](S+)[M|I|D]/){$string =~ s/S/M/g;}
	my $tcigar = "1".substr($string,0,1);
        for my $i (1..length($string)-1){
                my $char = substr($string, $i, 1);
                $tcigar =~ /(\d+)(\w)$/;
                if($char eq $2){
                        my $n = $1 + 1;
                        $tcigar =~ s/(\d+)(\w)$//;
                        $tcigar = "$tcigar"."$n"."$2";
                }else{
                        $tcigar .= "1"."$char";
                }
        }
        $cigar= $tcigar;

	if($cigar =~ /^(\d+)D/){
		#	$start = $start + $1 + $el;
		$cigar =~ s/^\d+D//;
	}


	print OUT "$readid\t$tag\t$hla\t$start\t$cigar\n";
}
close IN;
close OUT;
