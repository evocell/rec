#print "type input filename and press ENTER: ";
$infile = $ARGV[0];
#print "type output filename and press ENTER: ";
$outfile = $ARGV[1];

open(iF, $ARGV[2]); ### read in sequences
while($line = <iF>){
    @line = split("",$line);
    if($line[0] eq ">"){
	$line =~ s/ //g;
	$line =~ s/\n//g;
	push @names, $line;
    }
    if(@line[0] ne ">"){
	$line =~ s/ //g;
	$line =~ s/\n//g;
	$RSS .= $line;
	$RSS .= ">";
    }
}
close(iF);

@RSS = split(">",$RSS); 
$num = @RSS; ### number of RSS
$num--; ### last index

@temp = split("",@RSS[0]);
$length = @temp; ### length of the RSS
$length--; ### index for last base

$Compression[0] = "6 7 18";
$Compression[1] = "15 16 17";
$Compression[2] = "19 20 21";
$Compression[3] = "12 13 22";
$Compression[4] = "2 14 24";
$Compression[5] = "3 4";
$Compression[6] = "8 25";
$Compression[7] = "10 26";
$Compression[8] = "9 11";
$Compression[9] = "5 27";
@Compression0 = split(" ",$Compression[0]);
@Compression1 = split(" ",$Compression[1]);
@Compression2 = split(" ",$Compression[2]);
@Compression3 = split(" ",$Compression[3]);
@Compression4 = split(" ",$Compression[4]);
@Compression5 = split(" ",$Compression[5]);
@Compression6 = split(" ",$Compression[6]);
@Compression7 = split(" ",$Compression[7]);
@Compression8 = split(" ",$Compression[8]);
@Compression9 = split(" ",$Compression[9]);

$Nuc[0] = "a";
$Nuc[1] = "g";
$Nuc[2] = "c";
$Nuc[3] = "t";

## Maximize the Mean

$Score01 = 0;
$Score11 = 0;

%MarginalCounts = (); # 4 nucs by iPos
@MarginalTotals = ();

%PairCounts = (); # 4 nucs by 4 nucs by iPos by jPos; all iPos and jPos except @Compression
@PairTotals = (); # iPos and jPos

%TripletCounts = (); # 4 nucs by 4 nucs by 4 nucs by iPos by jPos by kPos except @Compression
@TripletTotals = (); # iPos, jPos and kPos

############ COUNTS ##############

foreach $seq (@RSS){ ### for each seq in temp set
    @seq = split("",$seq);
    for $iPos (2..$length){ ### marginal counts 
	for $iNuc (0..3){
	    if($Nuc[$iNuc] eq $seq[$iPos]){
		$MarginalCounts{$Nuc[$iNuc]}[$iPos]++;
	    }
	}
    } # $iPos MC
    
    for $iNuc (0..3){
	if($Nuc[$iNuc] eq $seq[0]){
	    $MarginalCounts{$Nuc[$iNuc]}[0]++;
	}
	if($Nuc[$iNuc] eq $seq[1]){
	    $MarginalCounts{$Nuc[$iNuc]}[1]++;
	}
    }
    
    $iCompression[5] = $Compression5[0]; # previous pair counts
    $jCompression[5] = $Compression5[1];
    $iCompression[6] = $Compression6[0]; 
    $jCompression[6] = $Compression6[1];
    $iCompression[7] = $Compression7[0]; 
    $jCompression[7] = $Compression7[1];
    $iCompression[8] = $Compression8[0]; 
    $jCompression[8] = $Compression8[1];
    $iCompression[9] = $Compression9[0]; 
    $jCompression[9] = $Compression9[1];
    for $k (5..9){
	for $iNuc (0..3){
	    if($Nuc[$iNuc] eq $seq[$iCompression[$k]]){ 
		for $jNuc (0..3){
		    if($Nuc[$jNuc] eq $seq[$jCompression[$k]]){ 
			$PairCounts{$Nuc[$iNuc]}{$Nuc[$jNuc]}[$iCompression[$k]][$jCompression[$k]]++;
		    } # if jNuc
		} # for jNuc
	    } # if iNuc
	} # for iNuc
    } # $k
    
    $iCompression[0] = $Compression0[0]; # previous triplet counts
    $jCompression[0] = $Compression0[1];
    $kCompression[0] = $Compression0[2];
    $iCompression[1] = $Compression1[0]; 
    $jCompression[1] = $Compression1[1];
    $kCompression[1] = $Compression1[2];
    $iCompression[2] = $Compression2[0]; 
    $jCompression[2] = $Compression2[1];
    $kCompression[2] = $Compression2[2];
    $iCompression[3] = $Compression3[0]; 
    $jCompression[3] = $Compression3[1];
    $kCompression[3] = $Compression3[2];
    $iCompression[4] = $Compression4[0]; 
    $jCompression[4] = $Compression4[1];
    $kCompression[4] = $Compression4[2];
    for $k (0..4){
	for $iNuc (0..3){
	    if(@seq[$iCompression[$k]] eq $Nuc[$iNuc]){
		for $jNuc (0..3){
		    if(@seq[$jCompression[$k]] eq $Nuc[$jNuc]){
			for $kNuc (0..3){
			    if(@seq[$kCompression[$k]] eq $Nuc[$kNuc]){
				$TripletCounts{$Nuc[$iNuc]}{$Nuc[$jNuc]}{$Nuc[$kNuc]}[$iCompression[$k]][$jCompression[$k]][$kCompression[$k]]++;
			    } # if kNuc
			} # for kNuc
		    } # if jNuc
		} # for jNuc
	    } # if iNuc
	} # for iNuc
    } # $k
    
} # $seq in @LOO 

####################### Totals and Frequencies ####################

for $iPos (2..$length){ ### compute marginal totals and frequencies
    for $iNuc (0..3){
	$MarginalTotals[$iPos] = $MarginalTotals[$iPos] + $MarginalCounts{$Nuc[$iNuc]}[$iPos];
    } # iNuc
    for $iNuc (0..3){
	$MarginalFrequency{$Nuc[$iNuc]}[$iPos] = ($MarginalCounts{$Nuc[$iNuc]}[$iPos] + .5)/($MarginalTotals[$iPos] + 2);
    } # iNuc
} # $iPos

print ($MarginalCounts{$Nuc[$iNuc]}[0]);
for $iNuc (0..3){
    $MarginalTotals[0] = $MarginalTotals[0] + $MarginalCounts{$Nuc[$iNuc]}[0];
    $MarginalTotals[1] = $MarginalTotals[1] + $MarginalCounts{$Nuc[$iNuc]}[1];
} # iNuc
for $iNuc (0..3){
    $MarginalFrequency{$Nuc[$iNuc]}[0] = ($MarginalCounts{$Nuc[$iNuc]}[0])/($MarginalTotals[0] + 2);
    $MarginalFrequency{$Nuc[$iNuc]}[1] = ($MarginalCounts{$Nuc[$iNuc]}[1])/($MarginalTotals[1] + 2);
} # iNuc

$iCompression[5] = @Compression5[0]; # previous pair totals and freqs
$jCompression[5] = @Compression5[1];
$iCompression[6] = @Compression6[0];
$jCompression[6] = @Compression6[1];
$iCompression[7] = @Compression7[0];
$jCompression[7] = @Compression7[1];
$iCompression[8] = @Compression8[0];
$jCompression[8] = @Compression8[1];
$iCompression[9] = @Compression9[0];
$jCompression[9] = @Compression9[1];
for $k (5..9){
    for $iNuc (0..3){ 
	for $jNuc (0..3){
	    $PairTotals[$iCompression[$k]][$jCompression[$k]] = $PairTotals[$iCompression[$k]][$jCompression[$k]] + $PairCounts{$Nuc[$iNuc]}{$Nuc[$jNuc]}[$iCompression[$k]][$jCompression[$k]];
	} # for jNuc
    } # for iNuc
    for $iNuc (0..3){
	for $jNuc (0..3){
	    $PairFrequency{$Nuc[$iNuc]}{$Nuc[$jNuc]}[$iCompression[$k]][$jCompression[$k]] = ($PairCounts{$Nuc[$iNuc]}{$Nuc[$jNuc]}[$iCompression[$k]][$jCompression[$k]] + .125)/($PairTotals[$iCompression[$k]][$jCompression[$k]] + 2);
	} # for jNuc
    } # for iNuc
} # $k

$iCompression[0] = $Compression0[0]; # previous triplet totals and frequencies
$jCompression[0] = $Compression0[1];
$kCompression[0] = $Compression0[2];
$iCompression[1] = $Compression1[0]; 
$jCompression[1] = $Compression1[1];
$kCompression[1] = $Compression1[2];
$iCompression[2] = $Compression2[0]; 
$jCompression[2] = $Compression2[1];
$kCompression[2] = $Compression2[2];
$iCompression[3] = $Compression3[0]; 
$jCompression[3] = $Compression3[1];
$kCompression[3] = $Compression3[2];
$iCompression[4] = $Compression4[0]; 
$jCompression[4] = $Compression4[1];
$kCompression[4] = $Compression4[2];
for $k (0..4){
    for $iNuc (0..3){
	for $jNuc (0..3){
	    for $kNuc (0..3){
		$TripletTotals[$iCompression[$k]][$jCompression[$k]][$kCompression[$k]] = $TripletTotals[$iCompression[$k]][$jCompression[$k]][$kCompression[$k]] + $TripletCounts{$Nuc[$iNuc]}{$Nuc[$jNuc]}{$Nuc[$kNuc]}[$iCompression[$k]][$jCompression[$k]][$kCompression[$k]];
	    } # for kNuc
	} # for jNuc
    } # for iNuc
    for $iNuc (0..3){
	for $jNuc (0..3){
	    for $kNuc (0..3){
		$TripletFrequency{$Nuc[$iNuc]}{$Nuc[$jNuc]}{$Nuc[$kNuc]}[$iCompression[$k]][$jCompression[$k]][$kCompression[$k]] = ($TripletCounts{$Nuc[$iNuc]}{$Nuc[$jNuc]}{$Nuc[$kNuc]}[$iCompression[$k]][$jCompression[$k]][$kCompression[$k]] + .03125)/($TripletTotals[$iCompression[$k]][$jCompression[$k]][$kCompression[$k]] + 2);
	    } # for kNuc
	} # for jNuc
    } # for iNuc
} # $k

################ SCORE TO TEST ADDITION OF A PAIR #############
@names = ();
open(iF, "$infile"); ### read in sequences
while($line = <iF>){
    $line =~ s/ //g;
    $line =~ s/\n//g;
    @line = split("",$line);
    if(@line[0] eq ">"){
	push @names, $line;
	$whole .= ">";
    }
    if(@line[0] ne ">"){
	$whole .= $line;
    }
}
close(iF);

@set = split(">",$whole);
$size = @set;
$size--;
print "$size sequences read in.\n";

for $i (1..$size){
    $target = $set[$i];
    $target =~ tr/A-Z/a-z/;
    
    @target = split("",$target);
    $Score01 = 0;
    
    $iCompression[0] = $Compression0[0]; # start score with previous triplets 
    $jCompression[0] = $Compression0[1];
    $kCompression[0] = $Compression0[2];
    $iCompression[1] = $Compression1[0]; 
    $jCompression[1] = $Compression1[1];
    $kCompression[1] = $Compression1[2];
    $iCompression[2] = $Compression2[0]; 
    $jCompression[2] = $Compression2[1];
    $kCompression[2] = $Compression2[2];
    $iCompression[3] = $Compression3[0]; 
    $jCompression[3] = $Compression3[1];
    $kCompression[3] = $Compression3[2];
    $iCompression[4] = $Compression4[0]; 
    $jCompression[4] = $Compression4[1];
    $kCompression[4] = $Compression4[2];
    for $k (0..4){
	for $iNuc (0..3){
	    if(@target[$iCompression[$k]] eq $Nuc[$iNuc]){
		for $jNuc (0..3){
		    if(@target[$jCompression[$k]] eq $Nuc[$jNuc]){
			for $kNuc (0..3){
			    if(@target[$kCompression[$k]] eq $Nuc[$kNuc]){
				$Score01 = $Score01 + log($TripletFrequency{$Nuc[$iNuc]}{$Nuc[$jNuc]}{$Nuc[$kNuc]}[$iCompression[$k]][$jCompression[$k]][$kCompression[$k]]);
			    } # if kNuc
			} # for kNuc
		    } # if jNuc
		} # for jNuc
	    } # if iNuc
	} # for iNuc
    } #k
    
    $iCompression[5] = $Compression5[0]; # add previous pair 
    $jCompression[5] = $Compression5[1];
    $iCompression[6] = $Compression6[0]; 
    $jCompression[6] = $Compression6[1];
    $iCompression[7] = $Compression7[0]; 
    $jCompression[7] = $Compression7[1];
    $iCompression[8] = $Compression8[0]; 
    $jCompression[8] = $Compression8[1];
    $iCompression[9] = $Compression9[0]; 
    $jCompression[9] = $Compression9[1];
    for $k (5..9){
	for $iNuc (0..3){
	    if(@target[$iCompression[$k]] eq $Nuc[$iNuc]){
		for $jNuc (0..3){
		    if(@target[$jCompression[$k]] eq $Nuc[$jNuc]){
			$Score01 = $Score01 + log($PairFrequency{$Nuc[$iNuc]}{$Nuc[$jNuc]}[$iCompression[$k]][$jCompression[$k]]);
		    }
		}
	    }
	}
    } # $k
    $Score11 = $Score01;
    for $iMarginal (2..$length){ ### add all the marginals
	$flag2 = "no";
	for $k (0..2){
	    if($iMarginal == $Compression0[$k]){
		$flag2 = "yes";
	    }
	    if($iMarginal == $Compression1[$k]){
		$flag2 = "yes";
	    }
	    if($iMarginal == $Compression2[$k]){
		$flag2 = "yes";
	    }
	    if($iMarginal == $Compression3[$k]){
		$flag2 = "yes";
	    }
	    if($iMarginal == $Compression4[$k]){
		$flag2 = "yes";
	    }
	}
	for $k (0..1){
	    if($iMarginal == $Compression5[$k]){
		$flag2 = "yes";
	    }
	    if($iMarginal == $Compression6[$k]){
		$flag2 = "yes";
	    }
	    if($iMarginal == $Compression7[$k]){
		$flag2 = "yes";
	    }
	    if($iMarginal == $Compression8[$k]){
		$flag2 = "yes";
	    }
	    if($iMarginal == $Compression9[$k]){
		$flag2 = "yes";
	    }
	}
	if($flag2 eq "no"){
	    for $iNuc (0..3){
		if(@target[$iMarginal] eq $Nuc[$iNuc]){
		    $Score11 = $Score11 + log($MarginalFrequency{$Nuc[$iNuc]}[$iMarginal]);
		}
	    }
	} # $flag2
    } # $iMarginal
    if(@target[0] eq "c"){
	$Score11 = $Score11 + log($MarginalFrequency{$Nuc[2]}[0]);
    }
    if(@target[1] eq "a"){
	$Score11 = $Score11 + log($MarginalFrequency{$Nuc[0]}[1]);
    }

    if(@target[0] ne "c" or @target[1] ne "a"){
	$Score11 = - 1000;
    }

    $namesindex = $i-1;

    $temp_name = $names[$namesindex];
    $temp_name =~ s/\r//g;
    $temp_seq = $set[$i];
    $temp_seq =~ s/\r//g;
    open(oF, ">>$outfile");
    print oF "$temp_name\t$temp_seq\t$Score11\n";
    close(oF);
}
