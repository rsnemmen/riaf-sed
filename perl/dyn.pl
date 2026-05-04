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
use FindBin qw($Bin);
use lib "$Bin/lib";
use Cwd qw(getcwd);
use File::Basename qw(basename);
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use IO::Select;
use IPC::Open3;
use POSIX qw(_exit);
use Storable qw(retrieve store);
use Symbol qw(gensym);
use Math::Derivative qw(Derivative1 Derivative2); 
use ADAF::Diagnostics qw(classify_solution);
use ADAF::Paths qw(fortran_binary);

# Module needed to benchmark the execution time of the code
use Benchmark; # see http://perldoc.perl.org/Benchmark.html
$bench0 = new Benchmark;

# Path to ADAF dynamics executable
$dynbinary=fortran_binary($Bin, "dynamics");

# Gets values of parameters from external parameter file
&readParam;

# Determines the outer boundary conditions
&findBCs;

# Initializes gnuplot
use Chart::Gnuplot;

$search_method = defined $search_method ? lc($search_method) : 'adaptive';
$eig_tol = defined $eig_tol ? convDbl($eig_tol) : 1e-4;
$dyn_timeout = defined $dyn_timeout ? convDbl($dyn_timeout) : 30;
$max_workers = defined $max_workers ? int($max_workers) : 1;
$max_workers = 1 if $max_workers < 1;

if ($search_method eq 'legacy') {
  &legacy_search;
} else {
  &adaptive_search;
}

if ($ops == 1) {
# To print the diagnostics for the bad solution
  print "\nNo solution found within the eigenvalue interval. Try changing the OBCs, or range/step of eigenvalues. \n";
  print "Final eigenvalue=$sl0 status=failed no_output=$nooutput sonic=$sonic mach_max=$largest r_mach_max=$largestR \n";
} else {
# To print the diagnostics for the nice solution
  print "\nSolution found. \n";
  print "Final eigenvalue=$sl0 status=ok no_output=$nooutput sonic=$sonic mach_max=$largest r_mach_max=$largestR shells=$linesout \n";
}
print "Output file: $diag \n";

# For benchmarking the execution time
$bench1 = new Benchmark;
$dbench = timediff($bench1, $bench0);
print "Runtime: ", timestr($dbench),"\n";



sub adaptive_search {
  $iterat=1;
  $brackets=0;
  $ops=1;

  my @candidates = coarse_candidates();
  my @coarse_results = $max_workers > 1 ? run_trials_parallel(@candidates) : ();
  my %coarse_by_value = map { $_->{eigenvalue} => $_ } @coarse_results;
  my $last_subsonic;

  foreach my $candidate (@candidates) {
    my $result = $max_workers > 1 ? $coarse_by_value{$candidate} : run_trial($candidate, getcwd(), File::Spec->rel2abs($diag), 0);
    die "Missing result for eigenvalue $candidate from parallel search.\n" unless defined $result;

    if ($result->{status} eq 'nice') {
      report_trial($result);
      if (defined $last_subsonic) {
        refine_bracket($last_subsonic->{eigenvalue}, $result->{eigenvalue}, $result);
        return;
      }
      rerun_final_solution($result->{eigenvalue});
      $ops=0;
      return;
    }

    if ($result->{status} eq 'subsonic') {
      report_trial($result);
      $last_subsonic = $result;
      next;
    }

    if (is_bad_for_bracket($result)) {
      if (!defined $last_subsonic) {
        print "iter=$iterat brackets=$brackets eigenvalue=$result->{eigenvalue} status=bad eigenvalue at 1st iteration; decrease the lower limit\n";
        load_result_globals($result);
        $sl0=$result->{eigenvalue};
        return;
      }

      report_trial($result, 'bracketing');
      $brackets++;
      refine_bracket($last_subsonic->{eigenvalue}, $result->{eigenvalue});
      return;
    }

    report_trial($result, 'weird result; check conditions');
    load_result_globals($result);
    $sl0=$result->{eigenvalue};
    return;
  }
}


sub legacy_search {
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
  $current_sl0=$sl0;

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
    print "iter=$iterat brackets=$brackets eigenvalue=$current_sl0 status=nice solution\n";
    $ops=0; # signal that solution was found, used after end of loop
    last; #  ****EXIT LOOP

    } elsif ( ($increase==0) && ($sonic =~ /Problem!/) ) {
    $lastok=$sl0; # stores the last eigenvalue computed with no stops or jumps
    $sl0=$sl0+$d_sl0;
    print "iter=$iterat brackets=$brackets eigenvalue=$current_sl0 status=subsonic\n";
    $ops=1;
    $iterat++;

  #  } elsif ($nooutput==1) {
  #  print "\nNo global solution. Change the OBCs! \n";
  #  $ops=1;
  #  last;

    } elsif ( $nooutput==1 || ($increase==1 && $sonic =~ /Problem!/) || ($increase==1 && $discont==1) ||  $weirdam==1 || ($discont==1 && $sonic !~ /Problem!/) ) { #  || $nan==1 || ($increase==0 && $sonic !~ /Problem!/)

    if ($nooutput==1) {
        $status="no global solution; verify the OBCs";
    } else {
        $status="bracketing";
    }

    if ($iterat==1) {
        print "iter=$iterat brackets=$brackets eigenvalue=$current_sl0 status=bad eigenvalue at 1st iteration; decrease the lower limit\n";
        $ops=1;
        last;  }

    $bad=$sl0; # stores the "bad" eigenvalue which caused a jump or hang
  # BRACKET the solution between the last OK eigenvalue and the
  # current BAD one
    $sl0f=$sl0;
    $d_sl0=($bad-$lastok)/$nmodels;
    $sl0=$lastok+$d_sl0;
    print "iter=$iterat brackets=$brackets eigenvalue=$current_sl0 status=$status\n";
    $iterat++;
    $brackets++;

    } else {
    print "iter=$iterat brackets=$brackets eigenvalue=$current_sl0 status=weird result; check conditions\n";
    $ops=1;
    last;

  }

  }
}


sub coarse_candidates {
  my @candidates;
  my $steps = int($nmodels);
  $steps = 1 if $steps < 1;
  my $step = ($sl0f - $sl0i) / $steps;

  for my $i (0 .. $steps) {
    push @candidates, $sl0i + $i * $step;
  }

  return @candidates;
}


sub refine_bracket {
  my ($low, $high, $best_nice) = @_;

  while (abs($high - $low) > $eig_tol) {
    my $step = ($high - $low) / $nmodels;
    last if $step <= 0;
    my $last_subsonic = $low;
    my $advanced = 0;

    for (my $candidate = $low + $step; $candidate <= $high + $step / 10.0; $candidate += $step) {
      $candidate = $high if $candidate > $high;
      my $result = run_trial($candidate, getcwd(), File::Spec->rel2abs($diag), 0);

      if ($result->{status} eq 'nice') {
        report_trial($result);
        rerun_final_solution($result->{eigenvalue});
        $ops=0;
        return;
      }

      if ($result->{status} eq 'subsonic') {
        report_trial($result);
        $last_subsonic = $result->{eigenvalue};
        $advanced = 1;
        next;
      }

      if (is_bad_for_bracket($result)) {
        report_trial($result, 'bracketing');
        $low = $last_subsonic;
        $high = $result->{eigenvalue};
        $brackets++;
        $advanced = 1;
        last;
      }

      report_trial($result, 'weird result; check conditions');
      load_result_globals($result);
      $sl0=$result->{eigenvalue};
      return;
    }

    last unless $advanced;
  }

  if (defined $best_nice) {
    rerun_final_solution($best_nice->{eigenvalue});
    $ops=0;
    return;
  }

  $sl0=$low;
  $sonic='Problem!';
  $largest='Problem!';
  $largestR='Problem!';
  $nooutput=1;
}


sub run_trials_parallel {
  my @candidates = @_;
  my $root = tempdir('dyn-search-XXXX', TMPDIR => 1, CLEANUP => 1);
  my @pending = @candidates;
  my %children;
  my @results;

  while (@pending || keys %children) {
    while (@pending && keys(%children) < $max_workers) {
      my $candidate = shift @pending;
      my $workdir = File::Spec->catdir($root, "trial-$candidate");
      make_path($workdir);
      prepare_trial_directory($workdir);
      my $result_file = File::Spec->catfile($root, "trial-$candidate.storable");
      my $pid = fork();
      die "Can't fork: $!\n" unless defined $pid;

      if ($pid == 0) {
        my $result = run_trial($candidate, $workdir, File::Spec->catfile($workdir, basename($diag)), 0);
        store($result, $result_file);
        _exit(0);
      }

      $children{$pid} = $result_file;
    }

    my $done = wait();
    last if $done == -1;
    my $result_file = delete $children{$done};
    push @results, retrieve($result_file) if defined $result_file && -e $result_file;
  }

  return sort { $a->{eigenvalue} <=> $b->{eigenvalue} } @results;
}


sub prepare_trial_directory {
  my ($workdir) = @_;
  my $cwd = getcwd();

  foreach my $support_file (qw(hot.dat)) {
    my $source = File::Spec->catfile($cwd, $support_file);
    next unless -e $source;
    my $target = File::Spec->catfile($workdir, $support_file);
    symlink($source, $target) || copy($source, $target) || die "Can't stage $source in $workdir: $!\n";
  }
}


sub run_trial {
  my ($trial_sl0, $workdir, $trial_diag, $keep_output) = @_;

  header($trial_diag, $trial_sl0);
  my ($timed_out, $exit_status) = run_dynamics($trial_sl0, $workdir, $trial_diag);
  my $result = classify_solution($trial_diag);
  $result->{eigenvalue} = $trial_sl0;
  $result->{timed_out} = $timed_out;
  $result->{exit_status} = $exit_status;
  $result->{status} = trial_status($result);

  unlink $trial_diag if !$keep_output && $trial_diag ne File::Spec->rel2abs($diag);

  return $result;
}


sub trial_status {
  my ($result) = @_;

  return 'timeout' if $result->{timed_out};
  return 'nice' if $result->{is_nice};
  return 'subsonic' if $result->{increase} == 0 && $result->{sonic} eq 'Problem!';
  return 'bad' if is_bad_for_bracket($result);
  return 'weird';
}


sub is_bad_for_bracket {
  my ($result) = @_;

  return 1 if $result->{status} && $result->{status} eq 'timeout';
  return 1 if $result->{nooutput};
  return 1 if $result->{nan};
  return 1 if $result->{failed};
  return 1 if $result->{weirdam};
  return 1 if $result->{increase} == 1 && $result->{sonic} eq 'Problem!';
  return 1 if $result->{increase} == 1 && $result->{discont} == 1;
  return 1 if $result->{discont} == 1 && $result->{sonic} ne 'Problem!';

  return 0;
}


sub report_trial {
  my ($result, $override_status) = @_;
  my $status = defined $override_status ? $override_status : $result->{status};
  $status = 'nice solution' if $status eq 'nice';
  $status = 'no global solution; verify the OBCs' if $status eq 'bad' && $result->{nooutput};
  print "iter=$iterat brackets=$brackets eigenvalue=$result->{eigenvalue} status=$status\n";
  $iterat++;
}


sub rerun_final_solution {
  my ($final_sl0) = @_;

  my $result = run_trial($final_sl0, getcwd(), File::Spec->rel2abs($diag), 1);
  load_result_globals($result);
  $sl0=$final_sl0;
}


sub load_result_globals {
  my ($result) = @_;

  $nooutput=$result->{nooutput};
  $sonic=$result->{sonic};
  $largest=$result->{largest};
  $largestR=$result->{largestR};
  $linesout=$result->{linesout};
  $increase=$result->{increase};
  $discont=$result->{discont};
  $weirdam=$result->{weirdam};
  $nan=$result->{nan};
  $failed=$result->{failed};
}


sub run_dynamics {
  my ($trial_sl0, $workdir, $trial_diag) = @_;
  my $reader;
  my $err = gensym;
  my $oldcwd = getcwd();
  my $timed_out = 0;
  my $exit_status = 0;
  my $pid;

  chdir $workdir or die "Can't chdir to $workdir: $!\n";
  $pid = open3(\*DYN, $reader, $err, $dynbinary);
  chdir $oldcwd or die "Can't chdir back to $oldcwd: $!\n";

  # Passes arguments to the fortran code
  print DYN "$gamai \n";
  print DYN "$m \n";
  print DYN "$beta \n";
  print DYN "$alfa \n";
  print DYN "$delta \n";
  print DYN "$dotm0 \n";
  print DYN "$rout \n";
  print DYN "$pp0 \n";
  print DYN "$ti \n";
  print DYN "$te \n";
  print DYN "$vcs \n";
  print DYN "$trial_sl0 \n";

  close(DYN);

  open (LOGAPPEND, ">>$trial_diag") ||
    die "Can't open $trial_diag !";

  my $ok = eval {
    local $SIG{ALRM} = sub { die "dyn_timeout\n"; };
    alarm($dyn_timeout);

    my $selector = IO::Select->new($reader, $err);
    while (my @ready = $selector->can_read) {
      foreach my $fh (@ready) {
        my $bytes = sysread($fh, my $buffer, 4096);
        if ($bytes) {
          print LOGAPPEND $buffer;
        } else {
          $selector->remove($fh);
          close($fh);
        }
      }
    }

    waitpid($pid, 0);
    $exit_status = $?;
    alarm(0);
    1;
  };

  if (!$ok) {
    alarm(0);
    $timed_out = 1 if $@ =~ /dyn_timeout/;
    kill 'TERM', $pid;
    sleep 1;
    kill 'KILL', $pid;
    waitpid($pid, 0);
    $exit_status = $?;
  }

  close(LOGAPPEND);

  return ($timed_out, $exit_status);
}








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
run_dynamics($sl0, getcwd(), File::Spec->rel2abs($diag));
}








# Prints header of each solution
sub header {
my ($header_diag, $header_sl0) = @_;
$header_diag = $diag unless defined $header_diag;
$header_sl0 = $sl0 unless defined $header_sl0;
# Some auxiliary calculations (needed only for the header of the log file)
#
# Virial temperature at the outer boundary
#$tvir=0.5444091492e13*(convDbl($gamai)-1.)/convDbl($rout);

# Prints header of log file
open (LOGFILE, ">$header_diag");
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
print LOGFILE "# $header_sl0 - eigenvalue of the problem (\"shooting\" parameter) \n";
print LOGFILE "# \n";
print LOGFILE "# Auxiliary values: \n";
print LOGFILE "# " . ($tvir/1e9) . "e+9 - T_vir - Virial temperature at the outer boundary \n";
print LOGFILE "# \n";
print LOGFILE "# Meaning of columns below: \n";
print LOGFILE "# 1. radius  \n";
print LOGFILE "# 2. v_R/c_s  \n";
print LOGFILE "# 3. log(T_e)   \n";
print LOGFILE "# 4. log(T_i)  \n";
print LOGFILE "# 5. q_advi/q_vis   \n";
print LOGFILE "# 6. (q_advi+q_adve)/q_vis \n";
print LOGFILE "# 7. c_s/c  \n";
print LOGFILE "# 8. H \n";
print LOGFILE "# 9. q_rad (cooling rate/volume) \n";
print LOGFILE "# 10. tau (optical depth) \n";
print LOGFILE "# 11. log(l_k) (spec. ang. mom. Kepl.) \n";
print LOGFILE "# 12. log(l) (spec. ang. mom.) \n";
print LOGFILE "# 13. magnetic field \n";
print LOGFILE "# 14. log(rho) [g/cm^3] \n";
print LOGFILE "# 15. log(ne) [1/cm^3] \n";
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
      if ($fields[0] =~ /^search_method$/) {$search_method=$fields[1];}
      if ($fields[0] =~ /^eig_tol$/) {$eig_tol=$fields[1];}
      if ($fields[0] =~ /^dyn_timeout$/) {$dyn_timeout=$fields[1];}
      if ($fields[0] =~ /^max_workers$/) {$max_workers=$fields[1];}
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
