Here is how I usually run the ADAF models in parallel. I created three folders inside adaf_code/perl named run01, run02 and run03. Inside each of these folders there is a parameter file in.dat. 

The guideline for setting the ADAF outer boundary conditions is:
- For large Rout (greater than a few thousands r_g):
e.g., rout=1e4rs, T_i=0.2 Tvir, T_e=0.19Tvir, vcs=0.2
- For Rout~100Rs:
e.g., rout=100 r_g, T_i=0.6 Tvir, T_e=0.08Tvir, vcs=0.5

I open one terminal with three tabs corresponding to each folder (or three terminals). Then I edit the three input files corresponding to a set of models, and finally run the program at the same time in each terminal. The relevant codes are: 
dyn.pl - computes ADAF dynamics
spectrum.pl - generates ADAF SED once you get a good global solution 
ssd.pl - computes truncated thin disk SED

For visualizing the resulting SEDs, use the code 'work/projects/finished/liners/seds/misc/model.py'

The image readme.png shows a screenshot of my Mac OS X during a typical parallel run. This is especially useful in multi-core/multi-processor architectures since Feng's code is serial and it would take a considerable amount of effort to make it parallel.

If you are having trouble finding a global solution, try playing around with dyntype.pl. Instead of trying to find automatically the "shooting value" or eigenvalue of the boundary value problem, you input eigenvalues manually and inspect the resulting plots radius vs radial velocity.
