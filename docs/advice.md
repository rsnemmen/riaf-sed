# Potentially useful advice 

If you are having trouble finding a global solution, try playing around with `dyntype.pl`. Instead of trying to find automatically the "shooting value" or eigenvalue of the boundary value problem, you input eigenvalues manually and inspect the resulting plots radius vs radial velocity.

If you want to inspect the behaviour of the dynamical solution,  have a look at the log file `out` which you define in the 38th line of `in.dat`. Every line corresponds to a different cylindrical radius shell of the RIAF, with the first column listing the radius (in units of M) and the others giving physical quantities.

The eigenvalue parameter is usually between 1 and 3
 
Guideline to find the correct physical solution:
When running the code I should always look for the solutions which have a smooth behavior of v_R(R), the Mach number is >1 at the inner region (have a sonic point) and for which v_R decreases as R increases. 

Apparently the code "hangs" (stays processing for very long time) after passing the eigenvalue corresponding to the physical solution. This provides a way of automatically finding the correct solution!! Bracketing between the last valid solution and the "freezing" one.

When I get lots of NaN and FAILED! results from the dynamics code for a broad range of eigenvalues given some values of the basic parameters, try changing slightly the boundary conditions. 

## Boundary conditions

The guideline for setting the ADAF outer boundary conditions is:

- Large Rout (R>~ a few thousand R_g): e.g., Rout=1e4 Rs, T_i=0.2 Tvir, T_e=0.19Tvir, vcs=0.2
- Rout ~ 100Rs: e.g., Rout=100 R_g, T_i=0.6 Tvir, T_e=0.08Tvir, vcs=0.5

Please refer to the Appendix A of my [PhD thesis](http://hdl.handle.net/10183/16325) or [Yuan, Ma & Narayan 2008, ApJ, 679, 984](http://iopscience.iop.org/article/10.1086/587484/meta) for more information on the BC choices.


## Doing many model realizations "in parallel"

If you are fitting a SED, you need to explore many models and combinations of parameters. Unfortunately, we do not have yet any MCMC fitting that does this automatically for you. Here is how I usually run the ADAF models in parallel ("dumb parallelization"). I created three folders inside `adaf_code/perl` named `run01`, `run02` and `run03`. Inside each of these folders there is a parameter file `in.dat`. 

I open one terminal with three tabs corresponding to each folder (or three terminals). Then I edit the three input files corresponding to a set of models, and finally run the program at the same time in each terminal. Hence why "dumb parallelization".

The image below shows a screenshot of OS X during a typical parallel run. 

![OS X running code in parallel](../docs/osxparallel.png =300x) 




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




# References

General, succinct description of SED models: [Nemmen et al. (2014)](http://mnras.oxfordjournals.org/content/438/4/2804)

Broad review about theory and application of RIAFs: [Yuan & Narayan (2014)](https://www.annualreviews.org/doi/10.1146/annurev-astro-082812-141003) 

More details about models: [Rodrigo Nemmen's PhD thesis](http://hdl.handle.net/10183/16325) (in portuguese), [Yuan et al. (2003)](http://adsabs.harvard.edu/abs/2003ApJ...598..301Y)

Global solutions: [Manmoto et al. (1997)](http://iopscience.iop.org/article/10.1086/304817/meta); [Narayan et al. (1997)](http://iopscience.iop.org/article/10.1086/303591/meta)

Boundary conditions: Appendix A of [Nemmen's PhD thesis](http://hdl.handle.net/10183/16325) (in portuguese) or [Yuan, Ma & Narayan 2008, ApJ, 679, 984](http://iopscience.iop.org/article/10.1086/587484/meta). 