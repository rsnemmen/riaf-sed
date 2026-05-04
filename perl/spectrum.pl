#!/usr/bin/perl -w
#
# Passes many arguments to the ADAF fortran codes sent by Feng Yuan.
# You need first to compile these codes (of course).

# For computing derivatives. Download the required library from 
# http://search.cpan.org/~jarw/Math-Derivative-0.01/Derivative.pm
# and follow the readme instructions to install it.
use FindBin qw($Bin);
use lib "$Bin/lib";
use Math::Derivative qw(Derivative1 Derivative2);
use ADAF::Paths qw(fortran_binary parameter_file_from_args stage_data_files);
use Cwd qw(getcwd);
use IPC::Open3;
use IO::Select;
use Symbol qw(gensym);

# Needed so that I can plot with Gnuplot
use FileHandle; # see http://perl.plover.com/FAQs/Buffering.html

# Module needed to benchmark the execution time of the code
use Benchmark; # see http://perldoc.perl.org/Benchmark.html
$bench0 = new Benchmark;

# Module needed for copying files without using an external program
use File::Copy;

# Parameter file
$input=parameter_file_from_args($0, @ARGV);

# Gets values of parameters from external parameter file
&readParam;

# Path to ADAF spectrum executable
$specbin=fortran_binary($Bin, "spectrum");

# Stage IC lookup tables from data/ into working directory
stage_data_files($Bin, getcwd(),
    qw(aomi-1.dat aomi-2.dat aomi-3.dat
       aomi2-1.dat aomi2-2.dat aomi2-3.dat
       romi-1.dat romi-2.dat romi-3.dat));

# Gets the number of lines of data in the output file from the dynamics
# code.
$i=0;
open (LINES, "x.dat") || 
	die "Can't open output from dynamics!";
while (<LINES>) { $i++; }
close LINES;
$nlines=$i;
$progress_total=$nlines;
$progress_done=0;
$progress_interactive=-t STDOUT;

# Generates spectra without Comptonization and calculates the range
# of frequencies of the synchrotron peak
&withoutCompt;

# Generates spectrum with comptonization enabled and the using the 
# appropriate range of synchrotron frequencies
&withCompt;
&finishProgress;

# For benchmarking the execution time
$bench1 = new Benchmark;
$dbench = timediff($bench1, $bench0);
print "\nThe code took: ", timestr($dbench),"\n";

# Useful user information
$lognui=log(convDbl($nui))/log(10.);
$lognuf=log(convDbl($nuf))/log(10.);
print "\nRange of frequencies of comptonization (log10): \nfrom $lognui to $lognuf \n";

# The code dumps two copies of the SED, one at spectrum.dat and other with
# the filename specified by the user at the parameter file. In addition, 
# it creates two copies of the SED without comptonization, one at ses5.dat
# and other at <filename>_ (filename is specified by the user).
if ($outfile !~ /spectrum.dat/) {
    copy("spectrum.dat", $outfile) or 
	die "Additional spectrum file cannot be copied. \n";
    copy("ses5.dat",$outfile . "_");
    print "\nCreated files: spectrum.dat, ses5.dat, $outfile, ${outfile}_ \n";
} else {
    print "\nCreated files: spectrum.dat, ses5.dat \n";
}














# Function that converts from the double numeric format (1.d0) to real (1.e0)
sub convDbl {
   $temp=$_[0];
   $temp =~ s/d/e/;
   return $temp;
}













# Displays progress for the Comptonized spectrum pass.
sub renderProgress {
   my ($stage)=@_;
   return unless $progress_interactive;

   my $width=40;
   my $done=$progress_done;
   $done=$progress_total if $done>$progress_total;

   my $fraction=0;
   if ($progress_total>0) {
      $fraction=$done/$progress_total;
   }

   my $filled=int($fraction*$width+0.5);
   $filled=$width if $filled>$width;
   my $bar=("#" x $filled) . ("-" x ($width-$filled));
   my $percent=int($fraction*100+0.5);

   printf "\rSpectrum generation %-10s [%s] %3d%% (%d/%d)",
      $stage, $bar, $percent, $done, $progress_total;
   STDOUT->flush();
}

sub clearProgressLine {
   return unless $progress_interactive;
   print "\r" . (" " x 90) . "\r";
   STDOUT->flush();
}

sub finishProgress {
   $progress_done=$progress_total;
   &renderProgress("done");
   print "\n" if $progress_interactive;
}

sub handleSpectrumOutputLine {
   my ($line,$stage,$count_progress)=@_;
   chomp $line;

   if ($line =~ /^\s*\d+\s*$/) {
      if ($count_progress) {
         $progress_done++;
         &renderProgress($stage);
      }
      return;
   }

   &clearProgressLine if $count_progress;
   print "$line\n";
   &renderProgress($stage) if $count_progress;
}

sub drainSpectrumOutput {
   my ($selector,$stdout,$stderr,$stage,$count_progress)=@_;
   my %buffers;
   my $stdout_fileno=fileno($stdout);
   $buffers{$stdout}='';
   $buffers{$stderr}='';

   while (my @ready=$selector->can_read) {
      foreach my $fh (@ready) {
         my $bytes=sysread($fh, my $buffer, 4096);
         my $fh_fileno=fileno($fh);
         my $is_stdout=defined($fh_fileno) && $fh_fileno==$stdout_fileno;
         if ($bytes) {
            $buffers{$fh}.=$buffer;
            while ($buffers{$fh} =~ s/^(.*?\n)//) {
               my $line=$1;
               if ($is_stdout) {
                  &handleSpectrumOutputLine($line,$stage,$count_progress);
               } else {
                  &clearProgressLine if $count_progress;
                  print STDERR $line;
                  &renderProgress($stage) if $count_progress;
               }
            }
         } else {
            if (length $buffers{$fh}) {
               if ($is_stdout) {
                  &handleSpectrumOutputLine($buffers{$fh},$stage,$count_progress);
               } else {
                  &clearProgressLine if $count_progress;
                  print STDERR $buffers{$fh};
                  &renderProgress($stage) if $count_progress;
               }
            }
            $selector->remove($fh);
            close($fh);
         }
      }
   }
}

sub runSpectrumPass {
   my ($stage,$count_progress,@input_lines)=@_;
   my $stdout;
   my $stderr=gensym;
   my $pid=open3(\*SPEC,$stdout,$stderr,$specbin);
   &renderProgress($stage) if $count_progress;

   foreach my $line (@input_lines) {
      print SPEC "$line \n";
   }
   close(SPEC);

   my $selector=IO::Select->new($stdout,$stderr);
   &drainSpectrumOutput($selector,$stdout,$stderr,$stage,$count_progress);

   waitpid($pid,0);
   my $exit_status=$?;
   die "$specbin failed during $stage spectrum generation.\n"
      if $exit_status != 0;
}

# Subroutine that reads a file containing the model parameters. Gets the
# values of the parameters from this file.
sub readParam {

open (PARFILE, $input) || 
	die "Can't open $input: $!\n";

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
#      if ($fields[0] =~ /^y1$/) {$y1=$fields[1];}
      if ($fields[0] =~ /^dotm0$/) {$dotm0=$fields[1];}
      if ($fields[0] =~ /^rout$/) {$rout=$fields[1];}
#      if ($fields[0] =~ /^y2$/) {$y2=$fields[1];}
#      if ($fields[0] =~ /^qbreset$/) {$qbreset=$fields[1];}
#      if ($fields[0] =~ /^compton$/) {$compton=$fields[1];}
#      if ($fields[0] =~ /^nlines$/) {$nlines=$fields[1];}
      if ($fields[0] =~ /^spec$/) {$outfile=$fields[1];}
  }  
}

close PARFILE;
}






# Generates spectra without Comptonization and calculates the range
# of frequencies of the synchrotron peak
sub withoutCompt {
# Runs first without Comptonization!
&runSpectrumPass(
    "no-Compton",
    0,
    $beta,
    $m,
    $distance,
    $alfa,
    $dotm0,
    $rout,
    "1d10", # print any value, does not matter
    "1d14", # any value
    "0", # qbreset, not read from file anymore
    "1", # disable comptonization
    $nlines
);

# Now gets the frequency range of the synchrotron peak, which will set
# the range of seed photons for the comptonization in the next run of
# the spectrum code.
open (SPEC01, "ses5.dat") || 
	die "Can't open the spectrum generated by the first run!";
	
# Resets all arrays
@x=( ); # log10 of frequency
@y=( ); # log10 of nu*Lnu
@dydx=( );
@d2ydx2=( );

$i=0;
while (<SPEC01>) {
    if ($_ =~ /NAN/) { 
	print "NaN in the spectrum: be careful! \n"; 
    } else {
	$x[$i]=substr $_, 3, 12; #16; # 1st column - nu
	$y[$i]=substr $_, 18, 12; #17; # 2nd column - nu*Lnu
	$i++;
    }
}
  
close(SPEC01);

# Computes the first and second derivatives. Before that, checks if 
# the arrays are empty. 
if ($#x!=-1) {
    @dydx=Derivative1(\@x,\@y);
    @d2ydx2=Derivative2(\@x,\@y);
}

# Tests to locate the range of frequencies of the synchrotron peak 
$i=0; # counter
foreach (@x) {
    if ($i==0){ $nui=10.**$_; } # initial nu

    if ($d2ydx2[$i]>10.) {
	$nuf=10.**$_; # final nu
	last;
    }
    $i++;
}

# Such that the fortran code understand that these numbers have double
# precision
$nui=$nui . "d0";
$nuf=$nuf . "d0";
}










# Second run of the calculation of the SED. This time enables comptonization
# of synchrotron photons using the appropriate range of frequencies 
# calculated from the first run.
sub withCompt {
# 2nd run with Comptonization enabled
    &runSpectrumPass(
        "Compton",
        1,
        $beta,
        $m,
        $distance,
        $alfa,
        $dotm0,
        $rout,
        $nui,
        $nuf,
        "0",
        "0", # enables comptonization
        $nlines
    );
}
