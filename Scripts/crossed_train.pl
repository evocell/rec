#!/usr/bin/perl -w

#CRoSSeD v1.1

use warnings;
use strict;
use POSIX qw(floor ceil);

#Configure parameters for training
my$trainfile = 'trainset.txt';
my$modelfile = 'crossedmodel';
our$template = 'crossedtemplate'.$$;
my$correction = 1;
my$optimization = 1;
my$cut = 0;
my$debug = 0;
my$threshold = 0.9;
my$svar = 2;
for (my$i = 0;$i <= $#ARGV;$i++){
	if(($ARGV[$i] eq '-tf') and (defined($ARGV[$i+1]))){
		#trainfile with nucleotide sequence + any number of structural proprties
		$trainfile = $ARGV[$i+1];
		$i++;
	} elsif(($ARGV[$i] eq '-mf') and (defined($ARGV[$i+1]))){
		#output of this script: a crossed model file and crf file (named $modelfile.crf)
		$modelfile = $ARGV[$i+1];
		$i++;
	} elsif(($ARGV[$i] eq '-tpl') and (defined($ARGV[$i+1]))){
		#define template name for crf algorithm
		$template = $ARGV[$i+1];
		$i++;
	} elsif(($ARGV[$i] eq '-corr') and ($ARGV[$i+1] eq '0')){
		#correction extension off
		$correction = 0;
		$i++;
	} elsif(($ARGV[$i] eq '-opt') and ($ARGV[$i+1] eq '0')){
		#optimization extension off
		$optimization = 0;
		$i++;
	} elsif(($ARGV[$i] eq '-cut') and (defined($ARGV[$i+1]))){
		#cut from the sides of the motif to create a shorter data set: useful if no structure could be predicted for the edges (default is 0)
		$cut = $ARGV[$i+1];
		$i++;
	} elsif($ARGV[$i] eq '-d') {
		#debug mode, prints almost everything
		$debug = 1;
	} elsif($ARGV[$i] eq '-v') {
		#verbose mode is exactly the same as debug mode
		$debug = 1;
	} elsif(($ARGV[$i] eq '-thr') and (defined($ARGV[$i+1]))){
		#threshold for the optimization extension when improvement is significant (default is 10%)
		$threshold = $ARGV[$i+1];
		$i++;
	} elsif(($ARGV[$i] eq '-svar') and (defined($ARGV[$i+1]))){
		#maximum deviation considered for the correction extension (default is 2)
		$svar = $ARGV[$i+1];
		$i++;
	}
}
my$keepfiles = 0;

#Parse training data
open (INPUT,$trainfile) or die "Cannot open file with training data: ".$trainfile;
my$numoffeatures;
my$numofpossamples = 0;
my$numofnegsamples = 0;
my(@labels,@posdata,@negdata);
my$temprange = 0;
my$range;
my$input1 = 0;
my$input2 = 0;
while(<INPUT>){
	chomp $_;
	if ((not($input1))and(not($input2))){
		if ($_ =~ m/^[A-Z]\t[A-Z]/){
			$input2 = 1;
		} elsif ($_ =~ m/^[0-9]/){
			$input1 = 1;
		}
	}
	if($input1){
		if(($_ =~ m/^[A-Z]\t/) or($_ =~ m/^[A-Z]$/)){
			my@rowdata = split("\t",$_);
			unless (defined($numoffeatures)){
				$numoffeatures = $#rowdata;
			}
			if ($labels[$#labels] == 1){
				for my$feat (0 .. $numoffeatures){
					push @{$posdata[$numofpossamples][$feat]}, $rowdata[$feat];
				}
			} else {
				for my$feat (0 .. $numoffeatures){
					push @{$negdata[$numofnegsamples][$feat]}, $rowdata[$feat];
				}
			}
			unless(defined($range)){
			$temprange++;
			}
		} elsif ($_ = m/([0-9])/){
			push @labels, $1;
		} else {
			if ($labels[$#labels] == 1){
				$numofpossamples++;
			}else{
				$numofnegsamples++;
			}
			unless(defined($range)){
				if($temprange > 0){
					$range = ($temprange-1)/2;
					if((($temprange-1) % 2) > 0){
						die "Motif length must be odd\n";
					}
				}
			}
		}
	}
	if($input2){
		if ($_ =~ m/^[A-Z]\t[A-Z]/){
			my@rowdata = split("\t",$_);
			unless(defined($range)){
				my$motiflength = -1;
				while($rowdata[$motiflength+1] =~ m/[A-Z]/){
					$motiflength++;
				}
				$range = ($motiflength)/2;
				if((($motiflength ) % 2) > 0){
					die "Motif length (".$motiflength.") must be odd\n";
				}
			}
			unless (defined($numoffeatures)){
				$numoffeatures = $#rowdata/($range*2 + 1) - 1; #Nucleotide data is always feature 0, so we count one too many
				if(($#rowdata % ($range*2 + 1)) > 0){
					die "Invalid motif length in one (or more) of the features\nTotal length (excl label):".$#rowdata."\nNumber of features?:".floor($numoffeatures + 1)."\n";
				}
			}
			push @labels, $rowdata[$#rowdata];
			for my$feat (0 .. $numoffeatures){
				if ($labels[$#labels] == 1){
					push @{$posdata[$numofpossamples]},[@rowdata[$feat*($range*2 + 1) .. ($feat+1)*($range*2 + 1)-1]];
				} elsif ($labels[$#labels] == 0) {
					push @{$negdata[$numofnegsamples]},[@rowdata[$feat*($range*2 + 1) .. ($feat+1)*($range*2 + 1)-1]];
				} else {
					die "Unknown label ".$labels[$#labels]."\n";
				}
			}
			if ($labels[$#labels] == 1){
			$numofpossamples++;
			} elsif ($labels[$#labels] == 0) {
			$numofnegsamples++;
			}
		}
	}
}

if ((not($input1))and(not($input2))){
	die "Invalid training file\n"
}

#Output some data to assure the user that everything is running smoothly
print "Training file: ".$trainfile."\n";
print "Model file: ".$modelfile."\n";
if($correction){
	print "Correction extension ON\n";
} else {
	print "Correction extension OFF\n";
}
if($optimization){
	print "Optimization extension ON\n";
} else {
	print "Optimization extension OFF\n";
}
print "Number of positive samples found: ".$numofpossamples."\n";
print "Number of negative samples found: ".$numofnegsamples."\n";
print "Number of properties found (excl nucl): ".$numoffeatures."\n";
print "Motif length: ".($range*2+1)."\n";


#Write Template
open(TEMPLATE, '>'.$template);
my$count = 0;
for my$i (0 .. (($range*2+1)*($numoffeatures+1))-1){
	print TEMPLATE "U".$count.":%x[0,".$i."]\n";
	$count++;
	unless ($i == (($range*2+1)*($numoffeatures+1)-1)){
		print TEMPLATE "U".$count.":%x[0,".$i."]%x[0,".($i+1)."]\n";
		$count++;
	}
	if (($i - floor($i/($range*2+1))*($range*2+1)) < $range){
		print TEMPLATE "U".$count.":%x[0,".$i."]%x[0,".((floor($i/($range*2+1))+1)*($range*2+1)-($i - floor($i/($range*2+1))*($range*2+1))-1)."]\n";
		$count++;
	}

}
close TEMPLATE;


#Run first crossvalidation
my($prevscore,$hitsonvarref) = cv3fold([expand($svar,@posdata)],[@negdata],$numoffeatures,$range,$svar);
#This array will keep tabs on the highest scoring positions relative to the binding site

my@hitsonvar;
for my$pospos (0 .. $#{$hitsonvarref}){
	for my $var (0 .. $svar*2){
		if (defined(${$hitsonvarref}[$pospos][$var])){
			$hitsonvar[$pospos][$var] = 1;
		} else {
			$hitsonvar[$pospos][$var] = 0;
		}
	}
}


if($debug){print "Initial score: ".$prevscore."\n";}

#This next part is the optimisation extension
my(@smoptima,@binsoptima);
if($optimization){
if ($debug){print "Starting optimization extension...\n"}
for my$feat (1 .. $numoffeatures){
	unless ($prevscore == 0){
	
	#Optimilization of Number of bins and Smoothing value
	my$binloc = 5;
	my$smloc = 2;
	my$optimum = 0;
	
	my@binsums = (0, 1, -1, 2, -2, 3);
	my@smsums = (0, 1, -1, 2, -2);
	
	my%smposaoa;
	my%smnegaoa;
	my(%smmax,%smmin);
	my@scoremat;
	my(%binaoa,%negbinaoa);
	my@bintemp;
	
	while (not($optimum)){
		my$lowestscore;
		my$lowestnumofbins;
		my$lowestsmvalue;
		if ($debug){print "Starting at ".$smloc."-".$binloc." for feature ".$feat."\n";}
		foreach my$binsum (@binsums){
		foreach my$smsum (@smsums){
		my@bins;
		unless (abs($binsum*$smsum) > 1){
		my$numofbins = $binloc + $binsum;
		my$smvalue = $smloc + $smsum;
		
		unless ( ($numofbins < 2) or ($smvalue < 0)){
			if($debug){print $smvalue."-".$numofbins." : ";}
			unless (defined($scoremat[$smvalue][$numofbins])){
			
			my($max,$min);
			if (defined ($smmax{$smvalue})){
				$max = $smmax{$smvalue};
				$min = $smmin{$smvalue};
			} else {
			$max = -10000;
			$min = 10000;
			}
			
			#Do smoothing unless it's already been done before
			my@smposdata;
			if (defined @{$smposaoa{$smvalue}}){
				@smposdata = @{$smposaoa{$smvalue}};
			} else {
				for my$x (0 .. $#posdata){
					#if($debug){print "LOWESS value: ".$smvalue."\n"};
					my@smarray;
					unless ($smvalue == 0){
						@smarray = lowess($smvalue,@{$posdata[$x][$feat]});
					} else {
						@smarray =@{$posdata[$x][$feat]};
					}
					push @smposdata, [@smarray];
					foreach my$ele (@smarray[$cut .. $#smarray-$cut]){
						unless ($ele eq '+' ){
						if ($ele > $max){

							$max = $ele;
						} elsif ($ele < $min){

							$min = $ele;
						}
						}

					}	
				}
				$smposaoa{$smvalue} = [@smposdata];
			}
			my@smnegdata;
			if (defined @{$smnegaoa{$smvalue}}){
				@smnegdata = @{$smnegaoa{$smvalue}};
			} else {
				for my$x (0 .. $#negdata){
					my@smarray;
					unless ($smvalue == 0){
						@smarray = lowess($smvalue,@{$negdata[$x][$feat]});
					} else {
						@smarray =@{$negdata[$x][$feat]};
					}
					push @smnegdata, [@smarray];
					foreach my$ele (@smarray){
						unless ($ele eq '+' ){
						if ($ele > $max){
							$max = $ele;
						} elsif ($ele < $min){
							$min = $ele;
						}
						}
					}	
				}
			$smnegaoa{$smvalue} = [@smnegdata];
			$smmax{$smvalue} = $max;
			$smmin{$smvalue} = $min;
			}
			
			#Calculate bin limits
			for my$b (0 .. $numofbins-1){
				$bins[$b] = $min + (($max - $min)/($numofbins))*($b+1);
			}

			#With use of bin limits discretize dataset
			my@disposdata;
			for my$x (0 .. $#smposdata){
				for my$y (0 .. $#{$smposdata[$x]}){
					$disposdata[$x][$y] = where_in_bin($smposdata[$x][$y],@bins);
				}
			}
			
			my@disnegdata;
			for my$x (0 .. $#smnegdata){
				for my$y (0 .. $#{$smnegdata[$x]}){
					$disnegdata[$x][$y] = where_in_bin($smnegdata[$x][$y],@bins);
				}
			}
			
			#Run crossvalidation
			my($score,$hitsonvartemp) = cv3fold([expand($svar,@posdata)],[@negdata],$numoffeatures,$range,$svar,$feat,[expand_single($svar,@disposdata)],[@disnegdata]);
			if($debug){print $score."\t";}

			for my$pospos (0 .. $#{$hitsonvartemp}){
				for my $var (0 .. $svar*2){
					if (defined(${$hitsonvartemp}[$pospos][$var])){
						$hitsonvar[$pospos][$var]++;
					} 
				}
			}
			
			$scoremat[$smvalue][$numofbins] = $score;
			$binaoa{$smvalue}{$numofbins} = [@disposdata];
			$negbinaoa{$smvalue}{$numofbins} = [@disnegdata];
			} else { 
				if($debug){print $scoremat[$smvalue][$numofbins]."\t";}
			}
			if (not(defined($lowestscore))){
				$lowestscore = $scoremat[$smvalue][$numofbins];
				$lowestnumofbins = $numofbins;
				$lowestsmvalue = $smvalue;
				if (defined($bins[0])){
				undef(@bintemp);
				@bintemp = @bins;
				}
				
			} elsif ($scoremat[$smvalue][$numofbins] < $lowestscore){
				$lowestscore = $scoremat[$smvalue][$numofbins];
				$lowestnumofbins = $numofbins;
				$lowestsmvalue = $smvalue;
				if (defined($bins[0])){
				undef(@bintemp);
				@bintemp = @bins;
				}
				
			}
		}
		}
		#Skip untill here if the sm or bin value is invalid
		}
		}
		print "\n";
		#End of foreach Binsum/Smsum loop
		
		#Check if a new optimum was found, if so repeat the Binsum/smsum loop
		#Else save the refound optimal arrays (if it improves more than 10%) and start on the next feature
		if ( ($binloc == $lowestnumofbins) and ($smloc == $lowestsmvalue) ){
			$optimum = 1;
			if($debug){print "Optimum found at ".$smloc."-".$binloc." !\n"}
			if ($lowestscore < $prevscore*$threshold){
				if($debug){print "Score lowered! ".$prevscore."->".$lowestscore."\n"}
				push @smoptima, $smloc;
				push @binsoptima, [@bintemp];
				foreach my$seqnr (0 .. $#{$binaoa{$smloc}{$binloc}}){
					$posdata[$seqnr][$feat] = ${$binaoa{$smloc}{$binloc}}[$seqnr];
				}
				foreach my$seqnr (0 .. $#{$negbinaoa{$smloc}{$binloc}}){
					$negdata[$seqnr][$feat] = ${$negbinaoa{$smloc}{$binloc}}[$seqnr];
				}
				$prevscore = $lowestscore;
			} else {
				if($debug){print "No Improvement... ".$prevscore."->".$lowestscore."\n"}
				push @smoptima, "0";
				push @binsoptima, "0";
			}
		} else {
			$binloc = $lowestnumofbins;
			$smloc = $lowestsmvalue;
		}
	}
	} else {
		if($debug){print "No further improvement can be found!\n"}
			push @smoptima, "0";
			push @binsoptima, "0";
	}
}
#End of optimization algorithm
} else {
	for (0 .. $numoffeatures){
		push @smoptima, "0";
		push @binsoptima, "0";
	}
}


#Start training  model with the found features
if($debug){print "Training complete model\n"}
my$crftrain = 'crftrain'.$$;
open(OPTTRAIN, '>'. $crftrain);
@posdata = expand($svar,@posdata);
for my$x (0 .. (($#posdata+1)/($svar*2+1)-1)){
	my$highestvar;
	if($correction){
		$highestvar = findhighestpos(@{$hitsonvar[$x]});
		unless(defined($highestvar)){die "No optimizing steps could be completed\n"}
	} else {
		$highestvar = $svar;
	}
	if($debug){print " ".($x*($svar*2+1)+$highestvar)." ".$hitsonvar[$x][$highestvar]."\t";}
	for my$feat (0 .. $numoffeatures){
		print OPTTRAIN join("|\t",@{$posdata[$x*($svar*2+1)+$highestvar][$feat]}),"|\t"; #add | for easier parsing of the crf model
	}
	print OPTTRAIN "1\n";
}
if($debug){print "\n";}
#Add dummy row to train file
for my$dummy (0 .. (($range*2+1)*($numoffeatures+1)-1)){
	print OPTTRAIN "*\t";
}
print OPTTRAIN "0\n";

close OPTTRAIN;

if($debug){
	system('crf_learn '.$template.' '.$crftrain.' '.$modelfile.'crf -t');
} else {
	my$tempfile = $modelfile.'crf';
	`crf_learn $template $crftrain $tempfile -t`
}
unless($keepfiles){unlink($crftrain)}
unless($keepfiles){unlink($template)}

open(MODEL,'>'.$modelfile);
print MODEL "#CRoSSeD model, can be used in crossed_test.pl\n";
print MODEL "R:".$range." \n";
for my$feat (0 .. $numoffeatures-1){
	print MODEL "F:".($feat+1)." ";
	print MODEL "S:".$smoptima[$feat]." ";
	if (ref($binsoptima[$feat]) eq 'ARRAY'){
		print MODEL "B:".scalar(@{$binsoptima[$feat]})."\n";
		print MODEL join("\t",@{$binsoptima[$feat]});
	} else {
		print MODEL "B:0";
	}
	print MODEL "\n";
}
print MODEL "\n";

#####################
#####Subroutines#####
#####################

sub open_file {

my($filename) = @_;

open(GET_FILE_DATA, $filename) 
	or die "Cannot open file", $filename;

my(@filedata)=<GET_FILE_DATA>;

close GET_FILE_DATA;

return @filedata;

}

sub where_in_bin {
	
	my($value,@bins) = @_;
	if ($value =~ m/\+/){
		return($value);
	}
	for my$i (0 .. $#bins - 1){
		if ($value < $bins[$i]){
			return ($i);
		}
	} 
	return ($#bins);
}

sub lowess {
#Apply a lowess smoothign to a given array
	my($range,@actarray) = @_;
	#if($debug){foreach my$ele (@array){print $ele;}}
	#if($debug){print "\nLOWESS: ".$range."\n";}
	

	my$debug = 0;
	my@smoothedarray;
	
	my@array;
	my$front = 0;
	my$back = 0;
	my$switch = 0;
	for my$x (0 .. $#actarray){
		if ($actarray[$x] =~ m/\+/){
			if($switch){
				$back++;
			} else {
				$front++;
			}
		} else {
			$switch = 1;
			push @array,$actarray[$x];
		}
	}

	for my$x (0 .. $#array){
	#model is y = mx + c
	
	
	if($x - $range < 0){

		my$avey = 0;
		my$avex = 0;
		my$norm = 0;
		for my$r (0 .. 2*$range){
			my$weight = (1 - (abs(($x - $r)/(2*$range - $x +1))**3))**3;
			$avey += $weight * $array[$r];
			$avex += $weight * $r;
			$norm += $weight;
		}
		$avey = $avey / $norm;
		$avex = $avex / $norm;
		
		my$mtop = 0;
		my$mbot = 0;
		
		for my$r (0 .. 2*$range){
			my$weight = (1 - (abs(($x - $r)/(2*$range - $x +1))**3))**3;
			$mtop += $weight*(($r - $avex) * ($array[$r] - $avey));
			$mbot += $weight*(($r - $avex)**2);
		}
		
		$smoothedarray[$x] = ($mtop / $mbot) * $x + ($avey - ($mtop / $mbot) * $avex);	
		
	} elsif ($range + $x > $#array){
	
		my$avey = 0;
		my$avex = 0;
		my$norm = 0;
		for my$r ($#array - 2*$range .. $#array){
			my$weight = (1 - (abs(($x - $r)/($x - $#array - 2*$range +1))**3))**3;
			$avey += $weight * $array[$r];
			$avex += $weight * $r;
			$norm += $weight;
		}
		$avey = $avey / $norm;
		$avex = $avex / $norm;
		
		my$mtop = 0;
		my$mbot = 0;
		
		for my$r ($#array - 2*$range .. $#array){
			my$weight = (1 - (abs(($x - $r)/($x - $#array - 2*$range +1))**3))**3;
			$mtop += $weight*(($r - $avex) * ($array[$r] - $avey));
			$mbot += $weight*(($r - $avex)**2);
		}
		
		$smoothedarray[$x] = ($mtop / $mbot) * $x + ($avey - ($mtop / $mbot) * $avex)
		
	} else {
	
		my$avey = 0;
		my$avex = 0;
		my$norm = 0;
		for my$r ($x - $range .. $x + $range){
			my$weight = (1 - (abs(($x - $r)/(2*$range +1))**3))**3;
			$avey += $weight * $array[$r];
			$avex += $weight * $r;
			$norm += $weight;
		}
		$avey = $avey / $norm;
		$avex = $avex / $norm;
		
		my$mtop = 0;
		my$mbot = 0;
		
		for my$r ($x - $range .. $x + $range){
			my$weight = (1 - (abs(($x - $r)/(2*$range +1))**3))**3;
			$mtop += $weight*(($r - $avex) * ($array[$r] - $avey));
			$mbot += $weight*(($r - $avex)**2);
		}
		
		$smoothedarray[$x] = ($mtop / $mbot) * $x + ($avey - ($mtop / $mbot) * $avex)
	}
	
	}
	
	for (1 .. $front){
		unshift @smoothedarray,'+';
	}
	for (1 .. $back){
		push @smoothedarray,'+';
	}
	
	return(@smoothedarray);

}

sub value_in {
#A subroutine to check if a variable is already contained in an array
    my($value, @array) = @_;
    foreach my $element (@array)
    {
        return 1 if $value eq $element;
    }
    return 0;
}

sub cv3fold{
#Runs a 3-fold cross validation on the given negative and postive sets 
#If extra arrays are given, the subroutine will asume these arrays are supposed to replace the normal values of the feature $feat
	my@posdata = @{shift(@_)};
	my@negdata = @{shift(@_)};
	my$numoffeatures = shift @_;
	my$range = shift @_;
	my$svar = shift @_;
	my@hitsonvar;

	
	my($feat,@disposdata,@disnegdata);
	if(@_ > 0){
		$feat = shift @_;
		@disposdata = @{shift(@_)};
		@disnegdata = @{shift(@_)};
	} else {
		$feat = -1;
	}

	my$posref = 0;
	my@scores;
	my$opttrain = 'opttrain'.$$;
	my$opttest = 'opttest'.$$;
	my$optdata = 'optdata'.$$;		
	my$optmodel = 'optmodel'.$$;
	for my$ele (0 .. 2 ){
	open(OPTTRAIN, '>'.$opttrain);
	open(OPTTEST, '>'.$opttest);
	for my$x (0 .. $#posdata){
		if (($x - 2)%($svar*2+1) == 0){
		unless ($ele == floor($x/(scalar(@posdata)/3))){
			#Only add the central position for training
			for my$getfeat (0 .. $numoffeatures){
				unless ($getfeat == $feat){
					print OPTTRAIN join("\t",@{$posdata[$x][$getfeat]}),"\t";
				} else {
					print OPTTRAIN join("\t",@{$disposdata[$x]}),"\t";
				}
			}
			print OPTTRAIN "1\n";
		} else {
			#Add all five positions for testing so they can later be selected
			for my$var (-$svar .. $svar){
			for my$getfeat (0 .. $numoffeatures){
				unless ($getfeat == $feat){
					print OPTTEST join("\t",@{$posdata[$x+$var][$getfeat]}),"\t";
					
				} else {
					print OPTTEST join("\t",@{$disposdata[$x+$var]}),"\t";
				}
			}
			print OPTTEST "1\n";
			}
		}
		}
	}
	for my$x (0 .. $#negdata){
		for my$getfeat (0 .. $numoffeatures){
			unless ($getfeat == $feat){
				print OPTTEST join("\t",@{$negdata[$x][$getfeat]}),"\t";
			} else {
				print OPTTEST join("\t",@{$disnegdata[$x]}),"\t";
			}
		}
		print OPTTEST "0\n";
	}
					
	#Add dummy row to train file
	for my$dummy (0 .. (($range*2+1)*($numoffeatures+1)-1)){
		print OPTTRAIN "*\t";
	}
	print OPTTRAIN "0\n";
					
	close OPTTEST;
	close OPTTRAIN;
	
	`crf_learn $template $opttrain $optmodel`;
		
	system('crf_test -m '.$optmodel.' -v1 '.$opttest.' > '.$optdata);
	
	open (OPTITDATA, $optdata) or die "Error: Cannot open CRF++ output file";

	my(@posscores,@negscores);
	
	while (<OPTITDATA>){
		if ($_ =~ m/([0-9])\t([0-9])\/(-?[0-9]+.[0-9]+)/){
			my$truelabel = $1;
			my$predlabel = $2;
			my$pvalue = $3;
			if ($predlabel == 1){
				$pvalue = 1 - $pvalue;
			}
			if ($truelabel == 1){
				push @posscores, $pvalue;
			} else {
				push @negscores, $pvalue;
			}
		}
	}
	if ($#posscores == -1){
		die "Error: Invalid CRF++ output file"
	}
	
	unless (@posscores % ($svar*2+1) == 0){
		die "Number of positive samples not divisable by 5 at it ".$ele;
	}
	#Note scores here denote the probability of a negative hit, so the lower the better (thanks to historical reasons))
	my@lowestposscores;
	for my$i (0 .. @posscores/($svar*2+1)-1){
		#print @posscores[$i*5 .. $i*5+4],"\n";
		my$lowestpos = findlowestpos(@posscores[$i*($svar*2+1) .. $i*($svar*2+1)+4]);
		push @lowestposscores, $posscores[$i*($svar*2+1) + $lowestpos];
		$hitsonvar[$posref][$lowestpos] = 1;
		$posref++;
		#print $lowestpos,"\n";
	}
	
	my@sortedscores = sort{ $a <=> $b } @negscores;
	
	for my$pos_ele (@lowestposscores){
		my$found = 0;
		my$count = 0;
		while ($found == 0){
			if (defined ($sortedscores[$count])){
				if ($pos_ele < $sortedscores[$count]){
					$found = 1;
				} else {
					$count++;
				}
			} else {
				$found = 1;
			}
		}
		#print $count."\t";
		push @scores, $count;
	}
	}
	#print "\n";
	
	@scores = sort{$a <=> $b} @scores;
	my$score=0;
	for my$i (0 .. $#scores-1){
		$score += $scores[$i];
	}
	$score += $scores[$#scores]/2;
	# unlink($optmodel);
	# unlink($opttrain);
	# unlink($opttest);
	# unlink($optdata);
	
	return ($score,[@hitsonvar]);
}

sub fold_array {
#A subroutine to split an array into a $fold number of groups
	my($fold,@array) = @_;
	my@totalarray;
	my@outputarray;
	my$debug = 0;
	for my$i (1 .. $fold-1){
		if($debug){print $i.":\n";}
		for my$x (1 .. ceil($#array/$fold)){
			my$new = 0;
			while (not($new)){
				my$randpos = int(rand($#array+1));
				unless (value_in($randpos,@totalarray)){
					$new = 1;
					push @totalarray, $randpos;
					push @{$outputarray[$i-1]}, $array[$randpos];
					if($debug){print $randpos." ";}
				}
			}
		}
		if($debug){print "\n";}
	}
	if($debug){print $fold.":\n";}
	for my$i (0 .. $#array){
		unless (value_in($i,@totalarray)){
			push @{$outputarray[$fold-1]}, $array[$i];
			if($debug){print $i." ";}
		}
	}
	if($debug){print "\n";}
	
	return @outputarray;
}

sub findhighestpos {
	my@array = @_;
	my$highestscore = 0;
	my$highestpos;
	for my$i (0 .. $#array){
	if(defined($array[$i])){
		if($array[$i]>$highestscore){
			$highestscore = $array[$i];
			$highestpos = $i;
		}
	}
	}
	return ($highestpos);
	
}

sub findlowestpos {
	my@array = @_;
	my$lowestscore = $array[0];
	my$lowestpos = 0;
	for my$i (1 .. $#array){
		if($array[$i]<$lowestscore){
			$lowestscore = $array[$i];
			$lowestpos = $i;
		}
	}
	return ($lowestpos);
	
}

sub expand {
	my$svar = shift @_;
	my@posdata = @_;
	my@newposdata;
	foreach my$posarray (@posdata){
		for my$var (-$svar .. $svar){
			if($var == 0){
				push @newposdata,$posarray;
			} elsif ($var > 0){
				my@newposarray;
				for my$feat (0 .. $#{$posarray}){
					for (1 .. $var){
						push @{$newposarray[$feat]}, '-';
					}
					push @{$newposarray[$feat]}, @{${$posarray}[$feat]}[0 .. $#{${$posarray}[$feat]}-$var];
				}
				push @newposdata,[@newposarray];
			} else {
				my@newposarray;
				for my$feat (0 .. $#{$posarray}){
					push @{$newposarray[$feat]}, @{${$posarray}[$feat]}[-$var .. $#{${$posarray}[$feat]}];
					for (1 .. -$var){
						push @{$newposarray[$feat]}, '-';
					}
					
				}
				push @newposdata,[@newposarray];
			}
		}
	}
	return @newposdata;
}

sub expand_single {
	my$svar = shift @_;
	my@posdata = @_;
	my@newposdata;
	foreach my$posarray (@posdata){
		for my$var (-$svar .. $svar){
			if($var == 0){
				push @newposdata,$posarray;
			} elsif ($var > 0){
				my@newposarray;
				for (1 .. $var){
					push @newposarray, '-';
				}
				push @newposarray, @{$posarray}[0 .. $#{$posarray}-$var];
				push @newposdata,[@newposarray];
			} else {
				my@newposarray;
				push @newposarray, @{$posarray}[-$var .. $#{$posarray}];
				for (1 .. -$var){
					push @newposarray, '-';
				}
				push @newposdata,[@newposarray];
			}
		}
	}
	return @newposdata;
}