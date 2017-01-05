#!/usr/bin/perl -w
#
# Passes many arguments to the fortran code sent by Feng Yuan that
# calculates the spectrum of the standard thin disk. You need first to 
# compile these codes (of course).

# Path to ADAF spectrum executable
$specbin="~/work/projects/adafjet/adaf/fortran/ssd_alone";


# Needed so that I can plot with Gnuplot
use FileHandle; # see http://perl.plover.com/FAQs/Buffering.html

# Module needed to benchmark the execution time of the code
use Benchmark; # see http://perldoc.perl.org/Benchmark.html
$bench0 = new Benchmark;

# Module needed for copying files without using an external program
use File::Copy;

# Gets values of parameters from external parameter file
&readParam;

# Opens pipe to ADAF dynamics code
open(SSD,"| $specbin");

#print SSD "$distance \n";
print SSD "$m \n";
print SSD "$dotm0 \n";
print SSD "$rin \n";
print SSD "$theta \n";

close(SSD);


# The code dumps the spectrum at <filename>_ssd (filename is specified 
# by the user).
move("ssd_alone.dat",$outfile . "_ssd")  or 
    die "Additional spectrum file cannot be copied. \n";
print "\nCreated file: ${outfile}_ssd \n";

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
      if ($fields[0] =~ /^dotm0$/) {$dotm0=$fields[1];}
      if ($fields[0] =~ /^rout$/) {$rin=$fields[1];}
      if ($fields[0] =~ /^theta$/) {$theta=$fields[1];}
      if ($fields[0] =~ /^spec$/) {$outfile=$fields[1];}
  }  
}

close PARFILE;
}
