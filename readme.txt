How to use this code
======================
You have to use the Perl scripts inside the perl folder. Please refer to the readme.txt inside that folder. 

Regarding which boundary conditions to choose for each radius, please refer to the Appendix A of my PhD thesis or Yuan, Ma & Narayan 2008, ApJ, 679, 984. 

Learnings
==========
(All the text below was written before 4/18/2008) 

Eigenvalue parameter is usually between 1 and 3

Input files for spectrum.f must match the output of dynamics.f: e.g. 1.dat and hot*.dat

 - slk: Keplerian specific angular momentum
 - ssll: ADAF specific ang. mom.
 - tao: vertical optical depth
 - bb: magnetic field strength
 - qn: cooling rate per unit volume 
 - omigak: Keplerian angular momentum
 - omiga: ang. mom.
 
Guideline to find the correct physical solution:
When running the code I should always look for the solutions which have a smooth behavior of v_R(R), the Mach number is >1 at the inner region (have a sonic point) and for which v_R decreases as R increases. 

Guideline for setting the boundary conditions:
  - For large Rout (greater than a few thousands r_g):
    * If rout~10^4r_g, we usually set vcs=0.2 or 0.3, or even 0.1; 
    * T_i is about 0.1~0.3 T_viral, with T_viral being the viral temperature.
      T_e should be smaller than T_i. As rout is great, v_R is small, the
       ions and electrons have almost enough time to achieve thermal
       equilibrium, and so the difference between T_i and Te is small.
       For example:
           rout=10^4 r_g,  T_i~6e7K, T_e~5.9e7K,  vcs~0.2
    
  - For Rout~100 r_g:
    * vcs should be larger,  I usually set vcs=0.5. Feng often takes rout as a few 10^2 or 10^3 r_g. He once said it would be difficult to find the global solution if rout>=10^5 r_g.
    * you need to adjust them to have the global solution. 
    For example: \dot{M}=0.1\dot{M_Edd}, rout=100 r_g => T_i=15e9K, T_e=8e9K, vcs=0.5.

Apparently the code "hangs" (stays processing for very long time) after passing the eigenvalue corresponding to the physical solution. This provides a way of automatically finding the correct solution!! Bracketing between the last valid solution and the "freezing" one. (9/05)

When I get lots of NaN and FAILED! results from the dynamics code for a broad range of eigenvalues given some values of the basic parameters, try changing slightly the boundary conditions. (11/6)

To get the spectrum of the SSD you just need to run the dynamics code. (11/21)

Todo
======
X Write a Perl program to drive the ADAF codes.

X Understand what the code does:
look Manmoto+97 and Narayan+97, global solutions
check the proper behavior of the solution for v_r/c_s

X Ideas for making the code automatic:
- search eigenvalue space
- set a limit time for computations (like 20-30 sec), then switch to the next value
	- wait if program runs for too long
	- if so, kill it
	- how to identify the child?
- see if the program crashes
- test if Mach > 1 in the inner regions and Mach decreases outwards
- test for smoothness nearby the sonic point

Compare L from dynamics to L from spectrum.

Find out why the Mac g77 compiler does not work, while the gfortran linux compiler works.
- compare results from g77, gfortran in the two architectures
X - install gfortran on ilarion

X Understand and modify the new spectrum code.

11/21: how to set the range of Comptonization correctly?
why the mismatch between mdot in the adaf and ssd spectra?
where does the code hang?

Questions
===========
X What are the input parameters? List them and ask confirmation. 
What about y(1)=rout? Ask about the meaning of some variables.

X How to get the transonic solution?

X When I run the code I notice that I get "nice" transonic solutions for slightly different values of the shooting parameter (or eigenvalue). Illustrate with the BCs I have now for 2.3 and 2.31.
How to pick the correct solution? 

X What values of the boundary conditions should I pick? Or do I keep them the same?

X How to make the search of eigenvalues automatic?


