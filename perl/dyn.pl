#!/usr/bin/perl -w
#
# Passes many arguments to the ADAF fortran codes sent by Feng Yuan.
# You need first to compile these codes (of course).
#
# Given an initial range of eigenvalues, computes solutions, checks if 
# each solution is physical and bracket the right eigenvalue automatically. 
# This is a mix of adaf_family_manyfiles and diagnose.pl.
# This script requires you to press ctrl-C when a solution hangs.
#
# Algorithm:
# - search eigenvalue space
# - set a limit time for computations (like 20-30 sec), then switch to the next value
#     - wait if program runs for too long
#     - if so, kill it
#     - how to identify the child?
# - test if Mach > 1 in the inner regions and Mach decreases outwards
# - test for smoothness nearby the sonic point

# For computing derivatives. Download the required library from 
# http://search.cpan.org/~jarw/Math-Derivative-0.01/Derivative.pm
# and follow the readme instructions to install it.
use Math::Derivative qw(Derivative1 Derivative2); 

# Module needed to benchmark the execution time of the code
use Benchmark; # see http://perldoc.perl.org/Benchmark.html
$bench0 = new Benchmark;

# Gets values of parameters from external parameter file
&readParam;

# Determines the outer boundary conditions
&findBCs;

# Initializes gnuplot
use Chart::Gnuplot;

$sl0=$sl0i;
# Calculates the increment in sl0 given the desired number of models
$d_sl0=($sl0f-$sl0i)/$nmodels; # increment

# Stores number of iterations until nice solution is found
$iterat=1;
# Stores number of "bracketing conditions" found, i.e. how many times the 
# intervals of eigenvalues are divided before reaching the physical solution.
$brackets=0;

# Main loop
while ($sl0<=$sl0f+$d_sl0) {
print "\n \nEigenvalue = $sl0 \nIteration $iterat \nBrackets $brackets \n \n";

# Creates header of each log file
# The output file keeps being rewritten until the final run.
&header;

# Calls adaf Fortran code and computes dynamical solution
&dynamics; 

# Diagnose if the computed solution is physical or not
&diagnose;

# Plots radius x v_r/c_s for the solution
#&plot;

# Decides if the solution is OK or not, and what should be done next
if ($increase==1 && $discont==0 && $weirdam==0 && $sonic !~ /Problem!/ && $nan==0 && $nooutput==0) {
  print "\nNice solution!! (iteration $iterat, brackets $brackets) \n";
  $ops=0; # signal that solution was found, used after end of loop
  last; #  ****EXIT LOOP

  } elsif ( ($increase==0) && ($sonic =~ /Problem!/) ) {
  $lastok=$sl0; # stores the last eigenvalue computed with no stops or jumps
  $sl0=$sl0+$d_sl0;
  print "\nSubsonic solution \n";
  $ops=1;
  $iterat++;

#  } elsif ($nooutput==1) { 
#  print "\nNo global solution. Change the OBCs! \n";
#  $ops=1;
#  last;
  
  } elsif ( $nooutput==1 || ($increase==1 && $sonic =~ /Problem!/) || ($increase==1 && $discont==1) ||  $weirdam==1 || ($discont==1 && $sonic !~ /Problem!/) ) { #  || $nan==1 || ($increase==0 && $sonic !~ /Problem!/)

  if ($nooutput==1) {
      print "\n**************************** \n";
      print "\nNO GLOBAL SOLUTION. Verify the OBCs! \n";
      print "**************************** \n\n";
  }

  if ($iterat==1) {
      print "\n**************************** \n";
      print "Bad eigenvalue at 1st iteration! Decrease the lower limit. \n";
      print "**************************** \n";
      $ops=1;
      last;  }

  $bad=$sl0; # stores the "bad" eigenvalue which caused a jump or hang
# BRACKET the solution between the last OK eigenvalue and the 
# current BAD one
  $sl0f=$sl0;
  $d_sl0=($bad-$lastok)/$nmodels;
  $sl0=$lastok+$d_sl0; 
  print "\nBracketing condition found \n";
  $iterat++;
  $brackets++;
    
  } else {
  print "\nSomething weird happened! Check conditions! \n";
  $ops=1;
  last;
  
}

# Prints results of diagnostics
print "\nDiagnostics: \n";
print "Increasing? $increase \"Jumps\"? $discont R_sonic~$sonic Weird A. M.? $weirdam NaN? $nan \"FAILED\"? $failed \n";
print "No output? $nooutput Mach_max=$largest R_max=$largestR \n";
} 

if ($ops == 1) {
# To print the diagnostics for the bad solution
  print "\nDiagnostics: \n";
  print "Increasing? $increase \"Jumps\"? $discont R_sonic~$sonic Weird A. M.? $weirdam NaN? $nan \"FAILED\"? $failed \n";
  print "\nEigenvalue = $sl0   ||   No output? $nooutput Mach_max=$largest R_max=$largestR \n";
  print "\nNo solution found within the eigenvalue interval. Try changing the OBCs, \nor range/step of eigenvalues. \n"
} else {
# To print the diagnostics for the nice solution
  print "\nDiagnostics: \n";
  print "Increasing? $increase \"Jumps\"? $discont R_sonic~$sonic Weird A. M.? $weirdam NaN? $nan \"FAILED\"? $failed \n";
  print "\nEigenvalue = $sl0   ||   No output? $nooutput Mach_max=$largest R_max=$largestR \n";
  print "Number of shells = $linesout (number of lines in the output file) \n";
}

# For benchmarking the execution time
$bench1 = new Benchmark;
$dbench = timediff($bench1, $bench0);
print "\nThe code took: ", timestr($dbench),"\n";










# Run several diagnostics to determine if the solution is physical or not.
# Outputs these variables to the main code: weirdam, discont, increase, sonic.
# They will be used to determine if the solution is OK.
sub diagnose {
# Opens the output file from the dynamics code
open (INFILE, $diag) || 
  die "Can't open $diag !";
  
# Resets all arrays
@x=( );
@y=( );
@dydx=( );
@d2ydx2=( );

# Goes through the log file reading the first and second columns, 
# stores them as vectors. Dismisses lines that contain "#"s.
$i=0; # counter
$weirdam=0; # 1 if the specific angular momentum gets negative (weird behavior)
$nan=0; # 1 if the results contain NaN
$failed=0; # 1 if the results contain "FAILED!"
# The "if" just below avoids comments and empty lines in the log file
while (<INFILE>) {
  if ($_ !~ /#/  && $_ ne " ") {
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

# linesout gets the number of lines in the output file from the
# dynamics code. This number corresponds to the number of shells in the
# ADAF and will be important to compute the SED.
$linesout=$i;
#print "$linesout shells (number of lines in the output) \n";
close(INFILE);

# Computes the first and second derivatives of the Mach number 
# (radial velocity vs. radius). Before that, checks if the arrays are 
# empty. This means that the Fortran code hanged before producing any output.
if ($#x!=-1) {
  @dydx=Derivative1(\@x,\@y);
  @d2ydx2=Derivative2(\@x,\@y);
  $nooutput=0;
} else {
  $nooutput=1; # BAD!
}

# Diagnostics: is the solution physical?
# ========================================
# Performs two tests:
# 1. checks if the function is always increasing looking at dydx
# 2. checks if d2ydx2 has no negative values (jumps in the solution)

$i=0; # counter
$discont=0; # "boolean" variable, 1 if "discontinuous"
$increase=1; # 1 if it is an always increasing function
# This variable controls if the sonic point will be looked for or not. 
# This is useful because the sonic radius corresponds only to the first 
# occurrence of Mach>1.
$testsonic=1; 
$sonic="Problem!"; # initial "value" of sonic radius

# If no global solution was found at all, there is no output from the
# dynamics code (arrays $x, $y and derivatives are empty). In this case, 
# skip the tests below.
if ($nooutput==1) {
  $largest="Problem!"; 
  $largestR="Problem!";
} else { 
  $largest=$y[0]; # will store the largest element of the Mach number array
  $largestR=$x[0]; # stores the radius corresponding to the largest Mach number

foreach (@x) {
# 1st test, we want $increase=1. Dismisses numbers too small and 
# tests only for R<90.
  if ($dydx[$i]>0. && abs($dydx[$i])>0.01 && $_<=90.) { $increase=0; } 
    
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
} # this curly bracket closes the if(nooutput==1) condition
# End of diagnostics ======================================
}







# Calls adaf Fortran code and computes dynamical solution
sub dynamics {
# Opens pipe to ADAF dynamics code
open(DYN,"|~/work/projects/adafjet/adaf/fortran/dynamics_new | tee -a $diag") || 
  die "Can't open program dynamics_new! \n";

# Passes arguments to the fortran code
# Adiabatic index gamma
print DYN "$gamai \n";
# Black hole mass (in Solar masses)
print DYN "$m \n";
# ratio of gas to total pressure
print DYN "$beta \n";
# alpha viscosity
print DYN "$alfa \n";
# Fraction of turbulent dissipation that directly heats electrons
print DYN "$delta \n";
# Mdot_out (Eddington units)
print DYN "$dotm0 \n";
# R_out (units of R_S)
print DYN "$rout \n";
# p_wind ("strength of wind")
print DYN "$pp0 \n";

# Boundary conditions *******************
# T_i (ion temperature)
print DYN "$ti \n";
# T_e (electron temperature)
print DYN "$te \n";
# v_R/c_s (radial velocity/sound speed)
print DYN "$vcs \n";
# eigenvalue of the problem (the "shooting" parameter)
print DYN "$sl0 \n";

close(DYN);
}








# Prints header of each solution
sub header {
# Some auxiliary calculations (needed only for the header of the log file)
#
# Virial temperature at the outer boundary
#$tvir=0.5444091492e13*(convDbl($gamai)-1.)/convDbl($rout);

# Prints header of log file
open (LOGFILE, ">$diag");
print LOGFILE "# Input parameters for ADAF model: \n";
print LOGFILE "# $gamai - gamma - adiabatic index \n";
print LOGFILE "# $m - m - black hole mass (in 10^6 Solar masses) \n";
print LOGFILE "# $beta - beta - ratio of gas to total pressure \n";
print LOGFILE "# $alfa - alpha viscosity \n";
print LOGFILE "# $delta - delta - fraction of viscous energy that directly heats electrons \n";
print LOGFILE "# $dotm0 - mdot_out (Eddington units) \n";
print LOGFILE "# $rout - R_out (units of R_S) \n";
print LOGFILE "# $pp0 - p_wind - \"strength of wind\" \n";
print LOGFILE "# BOUNDARY CONDITIONS *********** \n";
print LOGFILE "# " . (convDbl($ti)/$tvir) . " - T_i/T_vir - ion temperature \n";
print LOGFILE "# " . (convDbl($te)/$tvir) . " - T_e/T_vir - electron temperature \n";
print LOGFILE "# $vcs - v_R/c_s - radial velocity/sound speed \n";
print LOGFILE "# $sl0 - eigenvalue of the problem (\"shooting\" parameter) \n";
print LOGFILE "# \n";
print LOGFILE "# Auxiliary values: \n";
print LOGFILE "# " . ($tvir/1e9) . "e+9 - T_vir - Virial temperature at the outer boundary \n";
print LOGFILE "# \n";
print LOGFILE "# Meaning of columns below: \n";
print LOGFILE "# 1. radius  8. H \n";
print LOGFILE "# 2. v_R/c_s 9. q_rad (cooling rate/volume) \n";
print LOGFILE "# 3. log(T_e)  10. tau (optical depth) \n";
print LOGFILE "# 4. log(T_i)  11. log(l_k) (spec. ang. mom. Kepl.) \n";
print LOGFILE "# 5. q_advi/q_vis  12. log(l) (spec. ang. mom.) \n";
print LOGFILE "# 6. (q_advi+q_adve)/q_vis 13. magnetic field \n";
print LOGFILE "# 7. c_s/c 14. log(rho) \n";
print LOGFILE "# \n";
close(LOGFILE);
}








# Function that converts from the double numeric format (1.d0) to real (1.e0)
sub convDbl {
   $temp=$_[0];
   $temp =~ s/d/e/;
   return $temp;
}











# Subroutine that determines self-consistently the boundary conditions 
# at the outer boundary of the accretion flow.
# WARNING: the formulas below are valid in the range R=100-10000 R_S,
# keep that in mind.
sub findBCs {

# Error handling: stops if R<100 R_S or R>10000 R_S
#if (convDbl($rout)<100. || convDbl($rout)>10000.) {
#  print "Error: ADAF outer radius must be 100 <= Rout <= 10000 R_S \n";
#  exit;
#}

# Virial temperature at the outer boundary
#$tvir=0.5444091492e13*(convDbl($gamai)-1.)/convDbl($rout);
$tvir=3.6e12/convDbl($rout);

# T_i (ion temperature)
#$ti="12.d9";
$ti=$ti*$tvir . "d0";

# T_e (electron temperature)
#$te="8.d9";
$te=$te*$tvir . "d0";

# v_R/c_s (radial velocity/sound speed) - again the formula below was guessed from Renyi's
# suggestions of BCs.
#$vcs="0.5d0";
#$vcs=0.503-3.03e-5*convDbl($rout);
#$vcs=$vcs . "d0";
}













# Subroutine that reads a file containing the model parameters. Gets the
# values of the parameters from this file.
sub readParam {

# Parameter file
$input="in.dat";

open (PARFILE, $input) || 
  die "Can't open $input !";

# The field separator is "=". It is important that the input values in the
# parameter file are in the strict format "var=value" (no quotes).
while (<PARFILE>) {
  if ($_ !~ /#/  && $_ ne " ") {
      @fields=split /=/, $_;
      chomp @fields;

      if ($fields[0] =~ /^gamai$/) {$gamai=$fields[1];}
      if ($fields[0] =~ /^m$/) {$m=$fields[1];}
      if ($fields[0] =~ /^beta$/) {$beta=$fields[1];}
      if ($fields[0] =~ /^alfa$/) {$alfa=$fields[1];}
      if ($fields[0] =~ /^delta$/) {$delta=$fields[1];}
      if ($fields[0] =~ /^dotm0$/) {$dotm0=$fields[1];}
      if ($fields[0] =~ /^rout$/) {$rout=$fields[1];}
      if ($fields[0] =~ /^pp0$/) {$pp0=$fields[1];}
      if ($fields[0] =~ /^sl0i$/) {$sl0i=$fields[1];}
      if ($fields[0] =~ /^sl0f$/) {$sl0f=$fields[1];}
      if ($fields[0] =~ /^nmodels$/) {$nmodels=$fields[1];}
      if ($fields[0] =~ /^ti$/) {$ti=$fields[1];}
      if ($fields[0] =~ /^te$/) {$te=$fields[1];}
      if ($fields[0] =~ /^vcs$/) {$vcs=$fields[1];}
      if ($fields[0] =~ /^diag$/) {$diag=$fields[1];}
  }  
}

close PARFILE;
}






sub plot{
my $multiChart = Chart::Gnuplot->new(
    terminal => 'aqua'
);

#----------------------------------------
# Top left chart
my @charts = ();
$charts[0][0] = Chart::Gnuplot->new(
    title => "v_r/c_s",
    xrange  => "[:60]",
);
my $dataSet = Chart::Gnuplot::DataSet->new(
    xdata => \@x,
    ydata => \@y,
    style     => "linespoints",
);
$charts[0][0]->add2d($dataSet);
#----------------------------------------

#----------------------------------------
# Top right chart
$charts[0][1] = Chart::Gnuplot->new(
    title => "d/dr(v_r/c_s)",
    xrange  => "[:60]"
);
$dataSet = Chart::Gnuplot::DataSet->new(
    xdata => \@x,
    ydata => \@dydx,
    style     => "linespoints",
);
$charts[0][1]->add2d($dataSet);
#----------------------------------------

#----------------------------------------
# Bottom left chart
$charts[1][0] = Chart::Gnuplot->new(
    title => "d2/dr2 (v_r/c_s)",
    xrange  => "[:60]"
);
$dataSet = Chart::Gnuplot::DataSet->new(
    xdata => \@x,
    ydata => \@d2ydx2,
    style     => "linespoints",
);
$charts[1][0]->add2d($dataSet);
#----------------------------------------

#----------------------------------------
# Bottom right chart
$charts[1][1] = Chart::Gnuplot->new(
    title => "",
);
#----------------------------------------

# Plot the multplot chart
$multiChart->multiplot(\@charts);
}
