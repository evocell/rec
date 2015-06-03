#print "type input filename and press ENTER: ";
$infile = $ARGV[0];
#print "type output filename and press ENTER: ";
$outfile = $ARGV[1];

$icomp[0] = 7;
$jcomp[0] = 8;
$kcomp[0] = 20;
$icomp[1] = 6;
$jcomp[1] = 23;
$kcomp[1] = 24; 

$icomp[2] = 27;
$jcomp[2] = 28;
$kcomp[2] = 20; # place holder
$icomp[3] = 16;
$jcomp[3] = 17;
$kcomp[3] = 20; # place holder
$icomp[4] = 12;
$jcomp[4] = 21;
$kcomp[4] = 20; # place holder
$icomp[5] = 19;
$jcomp[5] = 25;
$kcomp[5] = 20; # place holder
$icomp[6] = 10;
$jcomp[6] = 11;
$kcomp[6] = 20; # place holder
$icomp[7] = 14;
$jcomp[7] = 22;
$kcomp[7] = 20; # place holder
$icomp[8] = 34;
$jcomp[8] = 35;
$kcomp[8] = 20; # place holder
$icomp[9] = 9;
$jcomp[9] = 15;
$kcomp[9] = 20; # place holder
$icomp[10] = 3;
$jcomp[10] = 13;
$kcomp[10] = 20; # place holder
$icomp[11] = 4;
$jcomp[11] = 38;
$kcomp[11] = 20; # place holder
$icomp[12] = 33;
$jcomp[12] = 37;
$kcomp[12] = 20; # place holder

$compF = 18;
$compG = 26;
$compH = 29;
$compI = 30;
$compJ = 31;
$compK = 32;
$compL = 36;

open(iF, $ARGV[2]);
while($line = <iF>){
    $line =~ s/ //g;
    $line =~ s/\n//g;
    @line = split("",$line);
    if($line[0] eq ">"){
	push @names, $line;
    }
    if($line[0] ne ">"){
	push @RSS, $line;
    }
}
$number1 = @RSS;
$number1--;

@temp = split("",$RSS[1]);
$length = @temp;
$length--;

$Nuc[0] = "a";
$Nuc[1] = "g";
$Nuc[2] = "c";
$Nuc[3] = "t";

@MC = ();
@MT = ();
@MF = ();

@PC = ();
@PT = ();
@PF = ();

@TC = ();
@TT = ();
@TF = ();

@SeC = ();
@SeT = ();
@SeF = ();

$score = 0;
foreach $seq (@RSS){
    @LOOseq = split("",$seq);
    
    for $ilength (0..$length){ # marginals
	for $inuc (0..3){
	    if($LOOseq[$ilength] eq $Nuc[$inuc]){
		$MC[$inuc][$ilength]++;
	    }
	}
    }
    
    for $distance (1..$length){ # test pairs previous pairs
	$end = $length - $distance;
	for $ilength (0..$end){
	    $jlength = $ilength + $distance;
	    for $inuc (0..3){
		if($LOOseq[$ilength] eq $Nuc[$inuc]){
		    for $jnuc (0..3){
			if($LOOseq[$jlength] eq $Nuc[$jnuc]){
			    $PC[$inuc][$jnuc][$ilength][$jlength]++;
			}
		    }
		}
	    }
	}
    }
    
    for $cn (0..1){ # previous triplet
	for $inuc (0..3){
	    if($LOOseq[$icomp[$cn]] eq $Nuc[$inuc]){
		for $jnuc (0..3){
		    if($LOOseq[$jcomp[$cn]] eq $Nuc[$jnuc]){
			for $knuc (0..3){
			    if($LOOseq[$kcomp[$cn]] eq $Nuc[$knuc]){
				$TC[$inuc][$jnuc][$knuc][$icomp[$cn]][$jcomp[$cn]][$kcomp[$cn]]++;
			    }
			}
		    }
		}
	    }
	}
    }
    
    for $fnuc (0..3){ # previous septuplet
	if($LOOseq[$compF] eq $Nuc[$fnuc]){
	    for $gnuc (0..3){
		if($LOOseq[$compG] eq $Nuc[$gnuc]){
		    for $hnuc (0..3){
			if($LOOseq[$compH] eq $Nuc[$hnuc]){
			    for $inuc (0..3){
				if($LOOseq[$compI] eq $Nuc[$inuc]){
				    for $jnuc (0..3){
					if($LOOseq[$compJ] eq $Nuc[$jnuc]){
					    for $knuc (0..3){
						if($LOOseq[$compK] eq $Nuc[$knuc]){
						    for $lnuc (0..3){
							if($LOOseq[$compL] eq $Nuc[$lnuc]){
							    $SeC[$fnuc][$gnuc][$hnuc][$inuc][$jnuc][$knuc][$lnuc][$compF][$compG][$compH][$compI][$compJ][$compK][$compL]++;
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

for $ilength (2..$length) { # marginals
    for $inuc (0..3){
	$MT[$ilength] = $MT[$ilength] + $MC[$inuc][$ilength];
    }
    for $inuc (0..3){
	$MF[$inuc][$ilength] = ($MC[$inuc][$ilength] + .5)/($MT[$ilength] + 2);
    }
}

for $inuc (0..3){
    $MT[0] = $MT[0] + $MC[$inuc][0];
    $MT[1] = $MT[1] + $MC[$inuc][1];
}
for $inuc (0..3){
    $MF[$inuc][0] = ($MC[$inuc][0])/($MT[0] + 2);
     $MF[$inuc][1] = ($MC[$inuc][1])/($MT[1] + 2);
 }

 for $distance (1..$length){ # test pairs and previous pairs
     $end = $length - $distance;
     for $ilength (0..$end){
	 $jlength = $ilength + $distance;
	 for $inuc (0..3){
	     for $jnuc (0..3){
		 $PT[$ilength][$jlength] = $PT[$ilength][$jlength] + $PC[$inuc][$jnuc][$ilength][$jlength];
	     }
	 }
	 for $inuc (0..3){
	     for $jnuc (0..3){
		 $PF[$inuc][$jnuc][$ilength][$jlength] = ($PC[$inuc][$jnuc][$ilength][$jlength] + .125)/($PT[$ilength][$jlength] + 2);
	     }
	 }
     }
 }

 for $cn (0..1){ # previous triplet
     for $inuc (0..3){
	 for $jnuc (0..3){
	     for $knuc (0..3){
		 $TT[$kcomp[$cn]] = $TT[$kcomp[$cn]] + $TC[$inuc][$jnuc][$knuc][$icomp[$cn]][$jcomp[$cn]][$kcomp[$cn]];
	     }
	 }
     }
     for $inuc (0..3){
	 for $jnuc (0..3){
	     for $knuc (0..3){
		 $TF[$inuc][$jnuc][$knuc][$icomp[$cn]][$jcomp[$cn]][$kcomp[$cn]] = ($TC[$inuc][$jnuc][$knuc][$icomp[$cn]][$jcomp[$cn]][$kcomp[$cn]] + .03125)/($TT[$kcomp[$cn]] + 2);
	     }
	 }
     }
 }

 for $fnuc (0..3){
     for $gnuc (0..3){
	 for $hnuc (0..3){
	     for $inuc (0..3){
		 for $jnuc (0..3){
		     for $knuc (0..3){
			 for $lnuc (0..3){
			     $SeT[$compF][$compG][$compH][$compI][$compJ][$compK][$compL] = $SeT[$compF][$compG][$compH][$compI][$compJ][$compK][$compL] + $SeC[$fnuc][$gnuc][$hnuc][$inuc][$jnuc][$knuc][$lnuc][$compF][$compG][$compH][$compI][$compJ][$compK][$compL];
			 }
		     }
		 }
	     }
	 }
     }
 }
 for $fnuc (0..3){
     for $gnuc (0..3){
	 for $hnuc (0..3){
	     for $inuc (0..3){
		 for $jnuc (0..3){
		     for $knuc (0..3){
			 for $lnuc (0..3){
			     $SeF[$fnuc][$gnuc][$hnuc][$inuc][$jnuc][$knuc][$lnuc][$compF][$compG][$compH][$compI][$compJ][$compK][$compL] = ($SeC[$fnuc][$gnuc][$hnuc][$inuc][$jnuc][$knuc][$lnuc][$compF][$compG][$compH][$compI][$compJ][$compK][$compL] + .00012207031)/($SeT[$compF][$compG][$compH][$compI][$compJ][$compK][$compL] + 2);
			 }
		     }
		 }
	     }
	 }
     }
 }

 ####################### SCORE to add pair

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
    $sequence = $set[$i];
    $sequence =~ tr/A-Z/a-z/;
    @target = split("",$sequence);
    
    $score = 0;
    
    for $fnuc (0..3){ #  sept
	if($target[$compF] eq $Nuc[$fnuc]){
	    for $gnuc (0..3){ 
		if($target[$compG] eq $Nuc[$gnuc]){
		    for $hnuc (0..3){
			if($target[$compH] eq $Nuc[$hnuc]){
			    for $inuc (0..3){ 
				if($target[$compI] eq $Nuc[$inuc]){
				    for $jnuc (0..3){
					if($target[$compJ] eq $Nuc[$jnuc]){
					    for $knuc (0..3){
						if($target[$compK] eq $Nuc[$knuc]){
						    for $lnuc (0..3){
							if($target[$compL] eq $Nuc[$lnuc]){
							    $score = $score + log($SeF[$fnuc][$gnuc][$hnuc][$inuc][$jnuc][$knuc][$lnuc][$compF][$compG][$compH][$compI][$compJ][$compK][$compL]);
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    
    for $cn (2..12){ 
	for $inuc (0..3){ # add previous pair
	    if($target[$icomp[$cn]] eq $Nuc[$inuc]){
		for $jnuc (0..3){
		    if($target[$jcomp[$cn]] eq $Nuc[$jnuc]){
			$score = $score + log($PF[$inuc][$jnuc][$icomp[$cn]][$jcomp[$cn]]);
		    }
		}
	    }
	}
    }
    
    for $cn (0..1){ # add previous trip
	for $inuc (0..3){
	    if($target[$icomp[$cn]] eq $Nuc[$inuc]){
		for $jnuc (0..3){
		    if($target[$jcomp[$cn]] eq $Nuc[$jnuc]){
			for $knuc (0..3){
			    if($target[$kcomp[$cn]] eq $Nuc[$knuc]){
				$score = $score + log($TF[$inuc][$jnuc][$knuc][$icomp[$cn]][$jcomp[$cn]][$kcomp[$cn]]);
			    }
			}
		    }
		}
	    }
	}
    }
    
    for $iM (2..$length){
	$flag = "no";
	for $cn (0..12){
	    if($iM == $icomp[$cn]){
		$flag = "yes";
	    }
	    if($iM == $jcomp[$cn]){
		$flag = "yes";
	    }
	    if($iM == $kcomp[$cn]){
		$flag = "yes";
	    }
	}
	if($iM == $compF){
	    $flag = "yes";
	}
	if($iM == $compG){
	    $flag = "yes";
	}
	if($iM == $compH){
	    $flag = "yes";
	}
	if($iM == $compI){
	    $flag = "yes";
	}
	if($iM == $compJ){
	    $flag = "yes";
	}
	if($iM == $compK){
	    $flag = "yes";
	}
	if($iM == $compL){
	    $flag = "yes";
	}
	if($flag eq "no"){
	    for $inuc (0..3){
		if($target[$iM] eq $Nuc[$inuc]){
		    $score = $score + log($MF[$inuc][$iM]);
		}
	    }
	}
    }
    
    if(@target[0] eq "c"){
	$score = $score + log($MF[2][0]);
    }
    if(@target[1] eq "a"){
	$score = $score + log($MF[0][1]);
    }
    
    if(@target[0] ne "c" or @target[1] ne "a"){
	$score = - 1000;
    }
    
    $namesindex = $i - 1;
    
    open(oF, ">>$outfile");
    print oF "$names[$namesindex]\t$sequence\t$score\n";
    close(oF);
    
} # $iTarget









