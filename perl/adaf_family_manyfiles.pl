#!/usr/bin/perl -w
#
# Passes many arguments to the ADAF fortran codes sent by Feng Yuan.
# You need first to compile these codes (of course).
#
# Generates a family of solutions with different eigenvalues. Each solution
# is dumped in a different file.

# Needed so that I can plot with Gnuplot
use FileHandle; # see http://perl.plover.com/FAQs/Buffering.html

# Input parameters
# ==================
# Adiabatic index gamma
$gamai="1.5d0";  #"1.4666666666667d0";
# Black hole mass (in 10^6 Solar masses)
$m="2.5d0";
# ratio of gas to total pressure
$beta="0.5d0";
# alpha viscosity
$alfa="0.3d0";
# Fraction of turbulent dissipation that directly heats electrons
$delta="0.001d0";
# Mdot_out (Eddington units)
$dotm0="1.3d-4";
# R_out (units of R_S)
$rout="1.d4";
# p_wind ("strength of wind")
$pp0="0.d0";

# Boundary conditions at the outer boundary of the flow
# T_i (ion temperature)
$ti="27220457.46d0";
# T_e (electron temperature)
$te="26403843.7362d0";
# v_R/c_s (radial velocity/sound speed)
$vcs="0.2d0";
# range of eigenvalues of the problem (the "shooting" parameter)
#print "Eigenvalue (shooting par.)? ";
$sl0i=1.0; # initial value
$sl0f=2.0; # final value
$nmodels=10; # number of models to be computed

# Basename of output files. Files corresponding to each solution will be
# created following the rule $diag-<eigenvalue>.dat (e.g. basename-1.01.dat, 
# basename-2.dat etc).
$diag="output";
# Spectrum filename
#$spec="test.dat";
# end of input parameters ==================

# Keeps the original basename untouched
$basename=$diag; 

# Calculates the increment in sl0 given the desired number of models
$d_sl0=($sl0f-$sl0i)/$nmodels; # increment




# Runs lots of models, creates for each solution an output file
for ($sl0 = $sl0i; $sl0 <= $sl0f; $sl0=$sl0+$d_sl0) {

# Creates the separate filenames: basename-<eigenvalue>
$diag="$basename-$sl0.dat";

print "\n \nEigenvalue = $sl0 \n";
print "Filename: $diag \n \n";

# Some auxiliary calculations (needed only for the header of the log file)
#
# Virial temperature at the outer boundary
$tvir=0.5444091492e13*(convDbl($gamai)-1.)/convDbl($rout);

# Prints header of log file
open (LOGFILE, ">$diag");
print LOGFILE "# Input parameters for ADAF model: \n";
print LOGFILE "# $gamai - gamma - adiabatic index \n";
print LOGFILE "# $m - m - black hole mass (in Solar masses) \n";
print LOGFILE "# $beta - beta - ratio of gas to total pressure \n";
print LOGFILE "# $alfa - alpha viscosity \n";
print LOGFILE "# $delta - delta - fraction of viscous energy that directly heats electrons \n";
print LOGFILE "# $dotm0 - mdot_out (Eddington units) \n";
print LOGFILE "# $rout - R_out (units of R_S) \n";
print LOGFILE "# $pp0 - p_wind - \"strength of wind\" \n";
print LOGFILE "# BOUNDARY CONDITIONS *********** \n";
print LOGFILE "# $ti - T_i - ion temperature \n";
print LOGFILE "# $te - T_e - electron temperature \n";
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


}









# Function that converts from the double numeric format (1.d0) to real (1.e0)
sub convDbl {
   $temp=$_[0];
   $temp =~ s/d/e/;
   return $temp;
}
