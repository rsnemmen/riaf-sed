#!/usr/bin/perl -w
#
# Passes many arguments to the ADAF fortran codes sent by Feng Yuan.
# You need first to compile these codes (of course).

# Needed so that I can plot with Gnuplot
use FileHandle; # see http://perl.plover.com/FAQs/Buffering.html

# Module needed to benchmark the execution time of the code
use Benchmark; # see http://perldoc.perl.org/Benchmark.html
$bench0 = new Benchmark;

# Error handling
$args=@ARGV;
if ($args == 0) {
  print "Usage: adaf.pl <eigenvalue> \n";
  exit;
}

# Input parameters
# ==================
# Adiabatic index gamma
$gamai="1.3333333333d0";  #"1.4666666666667d0";
# Black hole mass (in 10^6 Solar masses)
$m="7.76d2";
# ratio of gas to total pressure
$beta="0.9d0";
# alpha viscosity
$alfa="0.3d0";
# Fraction of turbulent dissipation that directly heats electrons
$delta="0.3d0";
# Mdot_out (Eddington units)
$dotm0="2.8d-2";
# R_out (units of R_S)
$rout="5.d2";
# p_wind ("strength of wind")
$pp0="0.25d0";

# Outer boundary conditions (OBCs)
# T_i (ion temperature) in units of the Virial temperature
$ti=0.2;
# T_e (electron temperature) 
$te=0.19;
# Mach number=v_R/c_s (radial velocity/sound speed)
$vcs="0.2d0";

# eigenvalue of the problem (the "shooting" parameter)
#print "Eigenvalue (shooting par.)? ";
#$sl0=<STDIN>; # takes eigenvalue from std. input
#chomp $sl0;
$sl0=$ARGV[0]; # takes eigenvalue from command-line argument

# Name of log file
$diag="out.dat";
# Spectrum filename
#$spec="test.dat";
# end of input parameters ==================







# Some auxiliary calculations
#
# Virial temperature at the outer boundary
$tvir=3.6e12/convDbl($rout);
$ti=$ti*$tvir . "d0";
$te=$te*$tvir . "d0";

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
print LOGFILE "# 1. radius	8. H \n";
print LOGFILE "# 2. v_R/c_s	9. q_rad (cooling rate/volume) \n";
print LOGFILE "# 3. log(T_e)	10. tau (optical depth) \n";
print LOGFILE "# 4. log(T_i)	11. log(l_k) (spec. ang. mom. Kepl.) \n";
print LOGFILE "# 5. q_advi/q_vis	12. log(l) (spec. ang. mom.) \n";
print LOGFILE "# 6. (q_advi+q_adve)/q_vis	13. magnetic field \n";
print LOGFILE "# 7. c_s/c	14. log(rho) \n";
print LOGFILE "# \n";
close(LOGFILE);

# Opens pipe to ADAF dynamics code
open(DYN,"|~/Work/doutorado/adaf_code/dynamics_new | tee -a $diag");

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




# Opens Gnuplot for plotting the behavior of v_R/c_s.
# This is useful to check if the solution you got is physically 
# tenable, i.e. it has a sonic point and decreases as R increases.
#       References: Manmoto+97, Narayan+97, see also the description
#       of the Bondi problem by Frank+02, especially Fig. 2.1.

open(GP,"|gnuplot - > /dev/null"); # pipes commands to gnuplot
# Turn filehandle "hot", uses library FileHandle
GP->autoflush(1); 

# Gnuplot input commands
print GP "set terminal x11 \n";
print GP "set title \"v_R/c_s vs. R\" \n";
print GP "set logscale x \n";
print GP "plot [:60]\"$diag\" using 1:2 notitle with linespoints, 1. \n";

# For benchmarking the execution time
$bench1 = new Benchmark;
$dbench = timediff($bench1, $bench0);
print "The code took:", timestr($dbench),"\n";

# Note that if I disaple the autoflush above, I can't plot properly.
# The reason is that Perl then will not wait for the data to be processed
# before plotting.
<STDIN>;
close(GP);










# Function that converts from the double numeric format (1.d0) to real (1.e0)
sub convDbl {
   $temp=$_[0];
   $temp =~ s/d/e/;
   return $temp;
}
