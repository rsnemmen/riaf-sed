#!/usr/bin/perl -w
#
# Passes many arguments to the ADAF fortran codes sent by Feng Yuan.
# You need first to compile these codes (of course).

# Needed so that I can plot with Gnuplot
use FileHandle; # see http://perl.plover.com/FAQs/Buffering.html

# Module needed to benchmark the execution time of the code
use Benchmark; # see http://perldoc.perl.org/Benchmark.html
$bench0 = new Benchmark;

# Gets values of parameters from external parameter file
&readParam;

# Opens pipe to ADAF dynamics code
open(SPEC,"|~/Work/doutorado/adaf_code/spectrum_new");

print SPEC "$beta \n";
print SPEC "$m \n";
print SPEC "$distance \n";
print SPEC "$alfa \n";
print SPEC "$dotm0 \n";
print SPEC "$rout \n";
print SPEC "$y1 \n";
print SPEC "$y2 \n";
print SPEC "$qbreset \n";
print SPEC "$compton \n";
print SPEC "$nlines \n";

close(SPEC);




# Opens Gnuplot for plotting the behavior of v_R/c_s.
# This is useful to check if the solution you got is physically 
# tenable, i.e. it has a sonic point and decreases as R increases.
#       References: Manmoto+97, Narayan+97, see also the description
#       of the Bondi problem by Frank+02, especially Fig. 2.1.

#open(GP,"|gnuplot - > /dev/null"); # pipes commands to gnuplot
# Turn filehandle "hot", uses library FileHandle
#GP->autoflush(1); 

# Gnuplot input commands
#print GP "set title \"v_R/c_s vs. R\" \n";
#print GP "plot [0:60]\"output.dat\" using 1:2 notitle with linespoints, 1. \n";

# Note that if I disaple the autoflush above, I can't plot properly.
#<STDIN>;
#close(GP);

# For benchmarking the execution time
$bench1 = new Benchmark;
$dbench = timediff($bench1, $bench0);
print "\nThe code took: ", timestr($dbench),"\n";






# Function that converts from the double numeric format (1.d0) to real (1.e0)
sub convDbl {
   $temp=$_[0];
   $temp =~ s/d/e/;
   return $temp;
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

      if ($fields[0] =~ /^distance$/) {$distance=$fields[1];}
      if ($fields[0] =~ /^m$/) {$m=$fields[1];}
      if ($fields[0] =~ /^beta$/) {$beta=$fields[1];}
      if ($fields[0] =~ /^alfa$/) {$alfa=$fields[1];}
      if ($fields[0] =~ /^y1$/) {$y1=$fields[1];}
      if ($fields[0] =~ /^dotm0$/) {$dotm0=$fields[1];}
      if ($fields[0] =~ /^rout$/) {$rout=$fields[1];}
      if ($fields[0] =~ /^y2$/) {$y2=$fields[1];}
      if ($fields[0] =~ /^qbreset$/) {$qbreset=$fields[1];}
      if ($fields[0] =~ /^compton$/) {$compton=$fields[1];}
      if ($fields[0] =~ /^nlines$/) {$nlines=$fields[1];}
  }  
}

close PARFILE;
}
