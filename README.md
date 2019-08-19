Radiatively inefficient accretion flow: Dynamics and spectrum
==================================================

This is a set of routines to compute the spectral energy distributions (SEDs) of radiatively inefficient accretion flows (RIAFs) around black holes. This should be useful for researchers interested in modeling the electromagnetic radiation from e.g. low-luminosity active galactic nuclei, X-ray binaries and other astrophysical applications. A RIAF consists of a geometrically thick, optically thin accretion flow filled with a very hot, two-temperature gas with the ion temperatures reaching 1E12 K. 

These routines use a semi-analytical approach to treat the radiation from the RIAF (also called sometimes advection-dominated accretion flows, ADAFs) in which the accretion flow is considered stationary assuming an α-viscosity and a pseudo-Newtonian gravity, and the radiative transfer is treated in considerable detail, taking into account synchrotron, inverse Compton scattering and bremsstrahlung processes as appropriate for hot plasmas (e.g. [Yuan et al. 2005](https://iopscience.iop.org/article/10.1086/427206); [Nemmen et al. 2006](https://iopscience.iop.org/article/10.1086/500571); [Nemmen et al. 2014](https://academic.oup.com/mnras/article/438/4/2804/2907740)).

The bottleneck of the calculations is in solving the dynamical structure of the flow and computing the inverse Compton radiation. The radiative transfer calculations take advantage of parallel architectures with OpenMP. The expected speedup is `ncores/2` compared with a serial run, where `ncores` is the number of CPU cores in your machine.

![The dashed line corresponds to the SED calculated for the RIAF around the black hole at the center of galaxy M87, taken from [Wong et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...849L..17W/abstract).](./m87sed.png) 
Figure: The dashed line corresponds to the SED calculated for the RIAF around the black hole at the center of galaxy M87, taken from [Wong et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...849L..17W/abstract).

# Requirements

- Fortran compiler (e.g., gfortran) with OpenMP support
- Perl with modules `Math::Derivative`, `Chart::Gnuplot`
- Optional: Gnuplot for diagnostic plots

To install the required Perl modules use the commands:

    cpan App::cpanminus
    sudo cpan Math::Derivative
    sudo cpan Chart::Gnuplot

# Installation

To compile the routines, please clone this repository in your machine and then issue these commands inside the repo folder:

    cd fortran
    make

The Fortran binaries will be located inside the `fortran` dir. The Perl binaries are in the `perl` dir.

# Model description

These routines use a semi-analytical approach to treat the radiation from the RIAF (also called sometimes advection-dominated accretion flows, ADAFs) in which the accretion flow is considered stationary assuming an α-viscosity and a pseudo-Newtonian gravity, and the radiative transfer is treated in considerable detail, taking into account synchrotron, inverse Compton scattering and bremsstrahlung processes as appropriate for hot plasmas (e.g. [Nemmen et al. 2006](https://iopscience.iop.org/article/10.1086/500571); [Nemmen et al. 2014](https://academic.oup.com/mnras/article/438/4/2804/2907740)).

Our model for the RIAF emission is described in [Nemmen et al. (2014)](https://academic.oup.com/mnras/article/438/4/2804/2907740). RIAFs are usually characterized by the presence of outflows or winds, which prevent a considerable fraction of the gas that is available at large radii from being accreted onto the black hole (see [Yuan & Narayan 2014](https://www.annualreviews.org/doi/10.1146/annurev-astro-082812-141003) for a review). In order to take this mass-loss into account, the radial variation of the accretion rate is parameterized as $\dot{M}(R) = \dot{M}_{\rm o} \left( R/R_{\rm o} \right)^{s}$ (or $\rho(R) \propto R^{-3/2+s}$) where $\dot{M}_{\rm o}$ is the rate measured at the outer radius $R_{\rm o}$ of the RIAF (Blandford & Begelman 1999). 

The parameters of the model are:

| Parameter | Description |
|:--|:--|
| *s* | power-law index for accretion rate (or density) radial variation |
| *\dot{M}_o* | mass accretion rate at the outer radius |
| *R_o* | outer radius |
| *M* | black hole mass |
| *alpha* | Shakura-Sunyaev viscosity parameter  |
| *beta* | ratio between the gas and total pressures |
| *delta* | fraction of energy dissipated via turbulence that directly heats electrons |
| *gamma* | adiabatic index |

The units of the parameters are described in the input parameter files included in the `examples` folder. 

# Usage

## Model setup

- edit `perl/dyn.pl`, `perl/spectrum.pl` and `perl/ssd.pl` and adjust the path to the executables (variables `$dynbinary`, `$specbin` and `$specbin`, respectively)
- cd to the directory that will contain the SED
- edit the input file `in.dat` with the desired model parameters
- include in this directory the following files: `aomi*dat`, `romi*dat`

## Compute SED

1. run `perl/dyn.pl` to compute ADAF dynamics to find physical global solution, adjusting range of eigenvalues `sl0i`,`sl0f` if required
2. once you get a good (physical) global solution in step 1, run `perl/spectrum.pl` to generate ADAF SED 
3. optional: run `perl/ssd.pl` to compute truncated thin disk SED

For visualizing the resulting SEDs, use the code `work/codes/python/model.py`

If you are having trouble finding a global solution, try playing around with `dyntype.pl`. Instead of trying to find automatically the "shooting value" or eigenvalue of the boundary value problem, you input eigenvalues manually and inspect the resulting plots radius vs radial velocity.

## Examples of models

Two examples of input parameter files are included in the `examples` folder:

- `largeR.dat`: ADAF with *R_out=1E4 Rs*
- `smallR.dat`: *R_out=500 Rs*

In order to compute the corresponding models, please rename the files to `in.dat` before running.
 


## Boundary conditions

The guideline for setting the ADAF outer boundary conditions is:

- Large Rout (R>~ a few thousand R_g): e.g., Rout=1e4 Rs, T_i=0.2 Tvir, T_e=0.19Tvir, vcs=0.2
- Rout ~ 100Rs: e.g., Rout=100 R_g, T_i=0.6 Tvir, T_e=0.08Tvir, vcs=0.5

Please refer to the Appendix A of my [PhD thesis](http://hdl.handle.net/10183/16325) or [Yuan, Ma & Narayan 2008, ApJ, 679, 984](http://iopscience.iop.org/article/10.1086/587484/meta) for more information on the BC choices.


# Citation

You are morally obligated to cite the following papers in any scientific literature that results from use of any part of this code:

1. [Yuan, F.; Cui, W. & Narayan, R. An Accretion-Jet Model for Black Hole Binaries: Interpreting the Spectral and Timing Features of XTE J1118+480. ApJ, 2005 , 620 , 905](https://iopscience.iop.org/article/10.1086/427206)
2. [Nemmen, R. S.; Storchi-Bergmann, T. & Eracleous, M.
Spectral models for low-luminosity active galactic nuclei in LINERs: the role of advection-dominated accretion and jets 
MNRAS, 2014 , 438 , 2804](http://mnras.oxfordjournals.org/content/438/4/2804)
3. [Yuan, F.; Zdziarski, A. A.; Xue, Y.; Wu, X. Modeling the Hard States of XTE J1550-564 during Its 2000 Outburst. ApJ, 2007, 659, 541](https://ui.adsabs.harvard.edu/abs/2007ApJ...659..541Y/abstract)





# References

General, succint description of SED models: [Nemmen et al. (2014)](http://mnras.oxfordjournals.org/content/438/4/2804)

Broad review about theory and application of RIAFs: [Yuan & Narayan (2014)](https://www.annualreviews.org/doi/10.1146/annurev-astro-082812-141003) 

More details about models: [Rodrigo Nemmen's PhD thesis](http://hdl.handle.net/10183/16325) (in portuguese), [Yuan et al. (2003)](http://adsabs.harvard.edu/abs/2003ApJ...598..301Y)

Global solutions: [Manmoto et al. (1997)](http://iopscience.iop.org/article/10.1086/304817/meta); [Narayan et al. (1997)](http://iopscience.iop.org/article/10.1086/303591/meta)

Boundary conditions: Appendix A of [Nemmen's PhD thesis](http://hdl.handle.net/10183/16325) (in portuguese) or [Yuan, Ma & Narayan 2008, ApJ, 679, 984](http://iopscience.iop.org/article/10.1086/587484/meta). 

# TODO 

By order of priority:

- [ ] port the core fortran code to C and better organize it
- [ ] add nonthermal emission
- [ ] OpenACC version for radiative transfer
- [ ] parallelize shooting method in `dyn.pl`


---

Copyright (c) 2019, [Rodrigo Nemmen](http://rodrigonemmen.com), [Feng Yuan](http://center.shao.ac.cn/fyuan/yuan.html).
[All rights reserved](http://opensource.org/licenses/BSD-2-Clause).


