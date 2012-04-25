#!/usr/bin/perl -w
#
# Computes a family of dynamical solutions using the ADAF code.
# Passes many arguments to the ADAF fortran codes sent by Feng Yuan.
# You need first to compile these codes (of course).

# Needed so that I can plot with Gnuplot
use FileHandle; # see http://perl.plover.com/FAQs/Buffering.html

# Input parameters
# ==================
# Adiabatic index gamma
$gamai="1.5d0";  #"1.4666666666667d0";
# Black hole mass (in 10^6 Solar masses)
$m="7.d2";
# ratio of gas to total pressure
$beta="0.9d0";
# alpha viscosity
$alfa="0.3d0";
# Fraction of turbulent dissipation that directly heats electrons
$delta="0.3d0";
# Mdot_out (Eddington units)
$dotm0="3.d-3";
# R_out (units of R_S)
$rout="150.d0";
# p_wind ("strength of wind")
$pp0="0.4d0";

# Boundary conditions at the outer boundary of the flow
# T_i (ion temperature)
$ti="1.814d9";
# T_e (electron temperature)
$te="1.71d9";
# v_R/c_s (radial velocity/sound speed)
$vcs="0.4984d0";
# range of eigenvalues of the problem (the "shooting" parameter)
#print "Eigenvalue (shooting par.)? ";
$sl0i=2.7; # initial value
$sl0f=2.9; # final value
$nmodels=10; # number of models to be computed

# Name of log file
$diag="outputs.dat";
# Spectrum filename
#$spec="test.dat";
# end of input parameters ==================





# Creates new log file
open (LOGFILE, ">$diag");
close(LOGFILE);

# Calculates the increment in sl0 given the desired number of models
$d_sl0=($sl0f-$sl0i)/$nmodels; # increment

# Runs lots of models, disable log file generation
for ($sl0 = $sl0i; $sl0 <= $sl0f+$d_sl0; $sl0=$sl0+$d_sl0) {
print "\n \n Eigenvalue = $sl0 \n \n";

open (LOGFILE, ">>$diag");
print LOGFILE "# $sl0 - eigenvalue of the problem (\"shooting\" parameter) \n";
close(LOGFILE);

# Opens pipe to ADAF dynamics code
open(DYN,"|~/doutorado/adaf_code/dynamics_new | tee -a $diag") || 
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


open (LOGFILE, ">>$diag");
print LOGFILE "\n\n";
close(LOGFILE);
}








# Function that converts from the double numeric format (1.d0) to real (1.e0)
sub convDbl {
   $temp=$_[0];
   $temp =~ s/d/e/;
   return $temp;
}
