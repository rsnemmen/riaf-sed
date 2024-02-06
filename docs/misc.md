# Useful advice 

If you want to inspect the behaviour of the dynamical solution,  have a look at the log file `out` which you define in the 38th line of `in.dat`. Every line corresponds to a different cylindrical radius shell of the RIAF, with the first column listing the radius (in units of M) and the others giving physical quantities.

The eigenvalue parameter is usually between 1 and 3
 
Guideline to find the correct physical solution:
When running the code I should always look for the solutions which have a smooth behavior of v_R(R), the Mach number is >1 at the inner region (have a sonic point) and for which v_R decreases as R increases. 

Apparently the code "hangs" (stays processing for very long time) after passing the eigenvalue corresponding to the physical solution. This provides a way of automatically finding the correct solution!! Bracketing between the last valid solution and the "freezing" one.

When I get lots of NaN and FAILED! results from the dynamics code for a broad range of eigenvalues given some values of the basic parameters, try changing slightly the boundary conditions. 

## Doing many model realizations "in parallel"

If you are fitting a SED, you need to explore many models and combinations of parameters. Unfortunately, we do not have yet any MCMC fitting that does this automatically for you. Here is how I usually run the ADAF models in parallel ("dumb parallelization"). I created three folders inside `adaf_code/perl` named `run01`, `run02` and `run03`. Inside each of these folders there is a parameter file `in.dat`. 

I open one terminal with three tabs corresponding to each folder (or three terminals). Then I edit the three input files corresponding to a set of models, and finally run the program at the same time in each terminal. Hence why "dumb parallelization".

The image below shows a screenshot of OS X during a typical parallel run. 

![OS X running code in parallel](./osxparallel.png =300x) 




## More details about setting the boundary conditions

For large Rout (greater than a few thousands r_g):
  
* If rout~10^4r_g, we usually set vcs=0.2 or 0.3, or even 0.1; 
* T_i is about 0.1~0.3 T_viral, with T_viral being the viral temperature.
* T_e should be smaller than T_i. As rout is great, v_R is small, the ions and electrons have almost enough time to achieve thermal equilibrium, and so the difference between T_i and Te is small.

For example: rout=10^4 r_g,  T_i~6e7K, T_e~5.9e7K,  vcs~0.2
    
For Rout~100 r_g:

* vcs should be larger,  I usually set vcs=0.5. Feng often takes rout as a few 10^2 or 10^3 r_g. He once said it would be difficult to find the global solution if rout>=10^5 r_g.
* you need to adjust them to have the global solution. 

For example: \dot{M}=0.1\dot{M_Edd}, rout=100 r_g => T_i=15e9K, T_e=8e9K, vcs=0.5.

## List of useful variables

 - slk: Keplerian specific angular momentum
 - ssll: ADAF specific ang. mom.
 - tao: vertical optical depth
 - bb: magnetic field strength
 - qn: cooling rate per unit volume 
 - omigak: Keplerian angular momentum
 - omiga: ang. mom.




# Software 

## List of codes and their functions

### Original codes (Feng's group)

`./bak/dynamics.f`, `/bak/spectrum.f`, `/bak/newsp.f`, `/bak/ssd.f`

Original codes sent by Feng.

### Core numerical routines  

**./fortran/spectrum_new.f
./fortran/dynamics_new.f
./fortran/ssd_new.f**
  My modifications with respect to dynamics.f and newsp.f. Added parameter input via standard input, descriptions, modified standard output.

```
dyn.pl    
│
└───dynamics_new.f
```

**./perl/dyn.pl**
  Given an initial range of eigenvalues, computes solutions, checks if each solution is physical and bracket the right eigenvalue automatically. This is a mix of adaf_family_manyfiles and diagnose.pl. This code tends to have the *most updated* fixes to bugs etc.

```
spectrum.pl    
│
└───spectrum_new.f
```

**./perl/spectrum.pl**
  Given a dynamical solution previously computed, calculates the spectrum. The parameters must match the dyn. solution. Currently for test purposes.

#### Thin disk


```
ssd.pl    
│
└───ssd_new.f
```
  
**./perl/ssd.pl**
  Opens a pipe to ssd_new.f and runs it.

**./fortran/ssd_alone.f**
  My modification of ssd.f where I discard the effect of the illumination of the thin disk by the ADAF. In other words, this code calculates the spectrum of a standard thin disk only.

**./perl/ssd_alone.pl**
  Pipes input into ssd_alone.f.

### Other 

**./perl/adaf.pl**  
Opens a pipe to dynamics_new and computes the dynamics for a given set of parameters. Plots Mach number vs. R as a diagnostic of the behavior of the solution.

**./perl/adaf_family.pl**
  Opens a pipe to dynamics_new and computes a family of ADAF models, varying the eigenvalue.

**./fortran/spectrum_thin.f
./fortran/dynamics_thin.f**
  Codes that enable the thin disk contribution. Feng recommended not to enable the thin disk inside the ADAF dynamics code, but use ssd.f to compute the thin disk spectrum.

./perl/diagnose.pl
  Given a log file from adaf.pl, runs a series of diagnostics on the resulting solution to check if it is consistent.
  
./perl/adaf_iterate.pl
  Calculates a series of ADAF solutions, iterating over the eigenvalue. Detects if the solution was computed successfully (but not yet whether it is physical) or whether the user interrupted the run pressing ctrl-C. Modified from adaf_family.

./perl/adaf_family_manyfiles.pl
  Calculates several ADAF solutions for a range of eigenvalues. Instead of dumping all the solutions into one file, creates individual files for each solution. I wrote this script specifically for diagnose_batch.pl.

./perl/diagnose_batch.pl
  Runs diagnose on a bunch of output files. These output files may be generated with adaf_family_manyfiles.pl for instance.

./perl/group_files.pl
  Takes the files generated by adaf_family_manyfiles.pl and concatenates all of them into one file, separates each file by two newlines. Useful for plotting all the solution together.
  
./fortran/dynamics_debug.f
./fortran/debug.f
  For debugging purposes, especially with gdb.


## Description of folders

./tests
  Several tests of the SED and dynamics.

./perl/run*
  Folders created to store results from simultaneous runs of the code. Useful for exploiting the capabilities of the multi-core processors. Each run of the code goes into only one core.