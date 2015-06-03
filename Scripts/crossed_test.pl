#!/usr/bin/perl -w

#CRoSSeD v1.1

use warnings;
use strict;
use POSIX qw(floor ceil);

#Configure parameters for training
my$testfile = 'testset.txt';
my$modelfile = 'crossedmodel';
my$outputfile = 'outputfile';
my$debug = 0;
my$correction = 1;
my$screen = 0;
my$var = 2;
for (my$i = 0;$i <= $#ARGV;$i++){
	if(($ARGV[$i] eq '-tf') and (defined($ARGV[$i+1]))){
		#testfile with nucleotide sequence + structural proprties (same format as used by crossed_train!)
		$testfile = $ARGV[$i+1];
		$i++;
	} elsif(($ARGV[$i] eq '-mf') and (defined($ARGV[$i+1]))){
		#a crossed model file and crf file (named $modelfile.crf)
		$modelfile = $ARGV[$i+1];
		$i++;
	}  elsif($ARGV[$i] eq '-d') {
		#debug mode, prints almost everything
		$debug = 1;
	} elsif(($ARGV[$i] eq '-out') and (defined($ARGV[$i+1]))){
		#output file with scores
		$outputfile = $ARGV[$i+1];
		$i++;
	} elsif(($ARGV[$i] eq '-corr') and ($ARGV[$i+1] eq '0')){
		#correction extension off
		$correction = 0;
	} elsif($ARGV[$i] eq '-screen') {
		#screening mode, every position will be scored instead of only the centre
		$screen = 1;
		$correction = 0;
	} elsif(($ARGV[$i] eq '-svar') and (defined($ARGV[$i+1]))){
		#maximum deviation considered for the correction extension (default is 2)
		$var = $ARGV[$i+1];
		$i++;
	}
}
my$keepfiles = 0;

#Parse test data
open (INPUT,$testfile) or die "Cannot open file with test data: ".$testfile;
my$numoffeatures;
my$numofsamples = 0;
my@testdata;
my$temprange = 0;
my$range;
my$input1 = 0;
my$input2 = 0;
my@labels;
while(<INPUT>){
	chomp $_;
	if ((not($input1))and(not($input2))){
		if ($_ =~ m/^[A-Z]\t[A-Z]/){
			$input2 = 1;
		} elsif (($_ =~ m/^[0-9]/) or ($_ =~ m/^[A-Z]\t[0-9]/)){
			$input1 = 1;
		}
	}
	if($input1){
		if(($_ =~ m/^[A-Z]\t/) or($_ =~ m/^[A-Z]$/)){
			my@rowdata = split("\t",$_);
			unless (defined($numoffeatures)){
				$numoffeatures = $#rowdata;
			}
			for my$feat (0 .. $numoffeatures){
				push @{$testdata[$numofsamples][$feat]}, $rowdata[$feat];
			}
			unless(defined($range)){
			$temprange++;
			}
			
		}  elsif ($_ =~ m/^(\d)/){
			$labels[$numofsamples] = $1;
		} else {
			$numofsamples++;
			unless(defined($range)){
				if($temprange > 0){
				$range = ($temprange-1)/2;
				}
			}
		}
	}
	if($input2){
		if ($_ =~ m/^[A-Z]\t[A-Z]/){
			
			#Read data
			my@rowdata = split("\t",$_);
			unless(defined($range)){
				my$motiflength = -1;
				while($rowdata[$motiflength+1] =~ m/[A-Z]/){
					$motiflength++;
				}
				$range = ($motiflength )/2;
				if((($motiflength % 2) > 0) and ($screen == '0')){
					die "Motif length (".$motiflength.") must be odd while screening is off\n";
				}
			}
			
			#Check data length
			my$tmpnumoffeatures = $#rowdata/($range*2 + 1) - 1; #Nucleotide data is always feature 0, so we count one too many
			if(($#rowdata % ($range*2 + 1)) > 0){
				$tmpnumoffeatures = ($#rowdata-1)/($range*2 + 1) - 1;
				if((($#rowdata-1)%($range*2 + 1))>0){
				die "Invalid motif length in one (or more) of the features\nTotal length (excluding label):".$#rowdata."\nNumber of features?:".floor($tmpnumoffeatures + 1)."\n";
				} 
			} else {
				$labels[$numofsamples] = $rowdata[$#rowdata];
			}
			unless (defined($numoffeatures)){
				$numoffeatures = $tmpnumoffeatures;
			}
			
			for my$feat (0 .. $numoffeatures){
				push @{$testdata[$numofsamples]},[@rowdata[$feat*($range*2 + 1) .. ($feat+1)*($range*2 + 1)-1]];
			}
			$numofsamples++;
		}
	}
}

#Output some data to assure the user that everything is running smoothly
print "Test file: ".$testfile."\n";
print "Model file: ".$modelfile."\n";
print "Number of samples found: ".$numofsamples."\n";
print "Number of properties found (excl nucl): ".$numoffeatures."\n";
print "Motif length: ".($range*2+1)."\n";
if($correction){
	print "Compensation for correction extension in training is ON\n";
}

my@modeldata = open_file($modelfile);
my(@smoptima,@binsoptima);
my$predefrange;
for(my$i = 0;$i <= $#modeldata;$i++){
	if($modeldata[$i] =~ m/F:([0-9]+)\sS:([0-9]+)\sB:([0-9]+)/){
		if($debug){print $modeldata[$i];}
		$smoptima[$1] = $2;
		if ($3 > 0){
			$binsoptima[$1] = [split("\t",$modeldata[$i+1])];
		} else {
			$binsoptima[$1] = 0;
		}
	} elsif ($modeldata[$i] =~ m/R:([0-9]+)/){
		$predefrange = $1;
	}
}

#Apply tranformations to entire validation set

if($debug){print "Applying smoothing and discretization\n"}
my$count = 0;
my@modtestdata;
foreach my$featdata (@testdata){
	$modtestdata[$count][0] = ${$featdata}[0];
	for my$feat (1 .. $numoffeatures){
		unless(defined($binsoptima[$feat])){
			die "Undefined bins for feature ".$feat."\n";
		}
		unless ($binsoptima[$feat] == 0){
			my@tempsmootharray;
			unless ($smoptima[$feat] == 0){
				@tempsmootharray = lowess($smoptima[$feat],@{${$featdata}[$feat]});
			} else {
				@tempsmootharray = @{${$featdata}[$feat]}
			}
			
			my@tempdisdata;
			for my$y (0 .. $#tempsmootharray){
				$tempdisdata[$y] = where_in_bin($tempsmootharray[$y],@{$binsoptima[$feat]});
			}
			$modtestdata[$count][$feat] = [@tempdisdata];
		} else {
			$modtestdata[$count][$feat] = [@{${$featdata}[$feat]}]
		}
	}
	$count++;
}

#Divide the data set if screening is on

if ($screen){
	if($debug){print "Redividing test data for screening\n"}
	unless(defined($predefrange)){
		$predefrange = 20;
		print "WARNING, range set to default of 20\n";
	}
	@modtestdata = redivide($predefrange,@modtestdata);
	
	$#labels = -1; #For screening, labels become ambiguous
}


if($correction){
	if($debug){print "Expanding test data\n"}
	@modtestdata = expand($var,@modtestdata);
	
	if($#labels > -1){
		my@temp;
		foreach my$i (0 .. $#labels){
			for (-$var .. $var){
				if(defined($labels[$i])){
					push (@temp, $labels[$i]);
				} else {
					push (@temp, 0);
				}
			}
		}
		@labels = @temp;
	}
}

if($debug){print "Running Model\n"}
my$crftest = 'crftest'.$$;
open(TEST, '>'.$crftest);
foreach my$i (0 .. $#modtestdata){
	unless ($modtestdata[$i] eq 'ENDOFSEQ'){
	foreach my$eledata (@{$modtestdata[$i]}){
		print TEST join("|\t", @{$eledata})."|\t"; #add | for easier parsing of the crf model
	}
	if(defined($labels[$i])){
		print TEST $labels[$i];
	} else {
		print TEST "0";
	}
	
	}
	print TEST "\n";
}
close TEST;
system('crf_test -m '.$modelfile.'crf -v1 '.$crftest.' > '.$outputfile);
unless($keepfiles){unlink($crftest)}


if($correction){
	if($debug){print "Shrinking test data\n"}
	my@outputdata = open_file($outputfile);
	open(OUTPUT,">".$outputfile);
	my$i = 1;
	while ($i<$#outputdata){
		my@weights;
		foreach my$x (0 .. $var*2){
			if($outputdata[$i+$x] =~ m/([0-9])\/([0-9]*\.?[0-9]+)/){
				if ($1 == 0){
					push @weights,(1 - $2);
				} else {
					push @weights,$2;
				}
			}
		}
		my$high = findhighestpos(@weights);
		print OUTPUT $outputdata[$i+$high];
		$i = $i + ($var*2+1);
	}
}

exit;
#####################
#####Subroutines#####
#####################

sub open_file {

my($filename) = @_;

open(GET_FILE_DATA, $filename) or die "Cannot open file ", $filename;

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
	my$svar = shift(@_);
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
						push @{$newposarray[$feat]}, '+';
					}
					push @{$newposarray[$feat]}, @{${$posarray}[$feat]}[0 .. $#{${$posarray}[$feat]}-$var];
				}
				push @newposdata,[@newposarray];
			} else {
				my@newposarray;
				for my$feat (0 .. $#{$posarray}){
					push @{$newposarray[$feat]}, @{${$posarray}[$feat]}[-$var .. $#{${$posarray}[$feat]}];
					for (1 .. -$var){
						push @{$newposarray[$feat]}, '+';
					}
					
				}
				push @newposdata,[@newposarray];
			}
		}
	}
	return @newposdata;
}

sub redivide {
	my($range,@posdata) = @_;
	my@newposdata;
	foreach my$posarray (@posdata){
		for my$pos (0 .. $#{${$posarray}[0]}){
			my@newposarray;
			if ($range > $pos){
				for my$feat (0 .. $#{$posarray}){
					for (1 .. $range - $pos){
						push  @{$newposarray[$feat]}, '+';
					}
					push @{$newposarray[$feat]}, @{${$posarray}[$feat]}[0 .. $pos+$range];
				}
			} elsif (($range+$pos) > $#{${$posarray}[0]}){
				for my$feat (0 .. $#{$posarray}){
					push @{$newposarray[$feat]}, @{${$posarray}[$feat]}[$pos-$range .. $#{${$posarray}[0]}];
					for ($#{${$posarray}[0]}+1 .. $range+$pos){
						push  @{$newposarray[$feat]}, '+';
					}
				}
			} else {
				for my$feat (0 .. $#{$posarray}){
					push @{$newposarray[$feat]}, @{${$posarray}[$feat]}[$pos-$range .. $pos+$range];
				}
			}
			push @newposdata,[@newposarray];
		}
		push @newposdata, 'ENDOFSEQ';
	}
	return @newposdata;
}
