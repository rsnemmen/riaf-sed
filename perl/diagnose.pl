#!/usr/bin/perl -w
#
# For computing derivatives. Download the required library from 
# http://search.cpan.org/~jarw/Math-Derivative-0.01/Derivative.pm
# and follow the readme instructions to install it.
use Math::Derivative qw(Derivative1 Derivative2); 

# Error handling
$args=@ARGV;
if ($args == 0) {
  print "Usage: diagnose.pl <filename> \n";
  exit;
}

$diag=$ARGV[0]; # takes name of output file from command-line argument
#$diag="out_bad04.dat";

# General idea: opens the output file from the dynamics code, creates vectors
# storing the Mach number and radius, computes the derivatives of the Mach number
# using these vectors. Runs a series of diagnostics on the solution to check
# if it has a consistent behavior.

# Opens the output file from the dynamics code
open (INFILE, $diag) || 
	die "Can't open $diag !";

# Goes through the log file reading the first and second columns, 
# stores them as vectors. Dismisses lines that contain "#"s.
$i=0; # counter
$weirdam=0; # 1 if the specific angular momentum gets negative (weird behavior)
$nan=0; # 1 if the results contain NaN
$failed=0; # 1 if the results contain "FAILED!"
# The "if" just below avoids comments and empty lines in the log file
while (<INFILE>) {
  if ($_ !~ /#/ && $_ ne " ") {
    $x[$i]=substr $_, 3, 21; #16; # 1st column - radius
    $y[$i]=substr $_, 26, 22; #17; # 2nd column - Mach number
    $i++;
  }
  
# Tests if the string "ssll < 0" is in the log file, then sets the "boolean"
# variable $weirdam (stands for weird angular momentum).
  if ($_ =~ /ssll \< 0/) { $weirdam=1; } 
  
# Checks for the presence of "NaN" (bad!) in the output
  if ($_ =~ /NaN/) { $nan=1; }
  
# Checks for the presence of "FAILED!" (bad!) in the output
  if ($_ =~ /FAILED/) { $failed=1; }  
}

close(INFILE);

# For testing: introduce a discontinuity (jump) in the Mach number at R=10.
# The size of the jump is given by $jump. Set $jump=0. to disable this testing
# feature. 
$i=0;
$jump=0.; # 0 to disable this testing feature! This is the size of the vertical jump
$diff[0]=0; #test: stores the difference in y between consecutive points
foreach (@x) {
  if ($_<10.) { $y[$i]=$y[$i]+$jump; }
  if ($i>0) { $diff[$i]=$y[$i]-$y[$i-1]; } #test
  $i++;
}

# Computes the first and second derivatives of the Mach number (radial velocity vs. 
# radius)
@dydx=Derivative1(\@x,\@y);
@d2ydx2=Derivative2(\@x,\@y);
@dydxdiff=Derivative1(\@x,\@diff); #test: derivative of $diff




# Diagnostics: is the solution physical?
# ========================================
# Performs two tests:
# 1. checks if the function is always increasing looking at dydx
# 2. checks if d2ydx2 has no negative values (jumps in the solution)

$i=0; # counter
$discont=0; # "boolean" variable, 1 if "discontinuous"
$increase=1; # 1 if it is an always increasing function
# This variable controls if the sonic point will be looked for or not. This is
# useful because the sonic radius corresponds only to the first occurrence of Mach>1.
$testsonic=1; 
$sonic="Problem!"; # initial "value" of sonic radius
$largest=$y[0]; # will store the largest element of the Mach number array
$largestR=$x[0]; # stores the radius corresponding to the largest Mach number

foreach (@x) {
# 1st test, we want $increase=1. Dismisses numbers too small.
  if ($dydx[$i]>0. && abs($dydx[$i])>0.01) { $increase=0; } 
    
# 2nd test (trickier), we want $discont=0. 
# Considers only R > 3 R_S, dismisses numbers too small with |Mach''|<0.01
# and Mach''>=-0.09.
  if ($d2ydx2[$i]<=-0.09 && abs($d2ydx2[$i])>0.01 && $_>3.) { 
    $discont=1; 
  } 

# Additionally, gets the location of the sonic point
  if ($y[$i]>1. && $testsonic==1) { 
    $sonic=$x[$i-1]; 
    $testsonic=0; # will not look again for a sonic point after this
  }

# Gets the largest value of the Mach number, useful to inspect if the solution
# is approaching the transonic behavior.  
  if ($y[$i]>$largest) { 
    $largest=$y[$i]; 
    $largestR=$_;
  }
  
  $i++;
}
# End of diagnostics ======================================





# Outputs the derivatives to STDOUT
$i=0; # counter
foreach (@x) {
  print "$_	$y[$i]	$dydx[$i]	$d2ydx2[$i]	$diff[$i]	$dydxdiff[$i] \n";
  $i++;
}

# Results of diagnostics
print "\n\# Increasing? $increase \"Jumps\"? $discont R_sonic~$sonic Weird A. M.? $weirdam NaN? $nan \n";
print "\# \"FAILED!\"? $failed \n";
if ($increase==1 && $discont==0 && $weirdam==0 && $nan==0 && $sonic !~ /Problem!/ && $failed==0) {
  print "\n\# Yeah, this solution is consistent!! \n";
} else {
  print "\n\# Oops, non-physical solution. \n";
}

print "\n\# Mach_max=$largest R_max=$largestR \n";
