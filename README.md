Radiatively inefficient accretion flow: Dynamics and spectrum
==================================================

This is a set of routines to compute the spectral energy distributions (SEDs) of radiatively inefficient accretion flows (RIAFs) around black holes. This should be useful for researchers interested in modeling the electromagnetic radiation from e.g. low-luminosity active galactic nuclei, X-ray binaries and other astrophysical applications. A RIAF consists of a geometrically thick, optically thin accretion flow filled with a very hot, two-temperature gas with the ion temperatures reaching $10^{12}$ K. 

These routines use a semi-analytical approach to treat the radiation from the RIAF (also called sometimes advection-dominated accretion flows, ADAFs) in which the accretion flow is considered stationary assuming an α-viscosity and a pseudo-Newtonian gravity, and the radiative transfer is treated in considerable detail, taking into account synchrotron, inverse Compton scattering and bremsstrahlung processes as appropriate for hot plasmas (e.g. [Yuan et al. 2005](https://iopscience.iop.org/article/10.1086/427206); [Nemmen et al. 2006](https://iopscience.iop.org/article/10.1086/500571); [Nemmen et al. 2014](https://academic.oup.com/mnras/article/438/4/2804/2907740)).

The bottleneck of the calculations is in solving the dynamical structure of the flow and computing the inverse Compton radiation. The radiative transfer calculations take advantage of parallel architectures with OpenMP. 

![The dashed line corresponds to the SED calculated for the RIAF around the black hole at the center of galaxy M87, taken from [Wong et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...849L..17W/abstract).](./docs/m87sed.png) 
Figure: The dashed line corresponds to the SED calculated for the RIAF around the black hole at the center of galaxy M87, taken from [Wong et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...849L..17W/abstract).

# Requirements

- Fortran compiler (e.g., gfortran) with OpenMP support
- Perl with modules `Math::Derivative`, `Chart::Gnuplot`
- Python 3 with Matplotlib for the smoke/regression comparison checks and diagnostic plots

To install the required Perl modules use the commands:

```shell
cpan App::cpanminus
sudo cpan Math::Derivative
sudo cpan Chart::Gnuplot
```

# Installation

To compile the routines, please clone this repository in your machine and then issue these commands inside the repo folder:

    make build

The compiled Fortran binaries will be placed in the `bin/` directory at the repository root. The Perl wrappers are in the `perl/` dir.

Convenience targets at the repository root:

```shell
make build # compile the Fortran executables via `fortran/Makefile`
make clean # remove the compiled Fortran binaries
make smoke # run the tests 
```

# Model description

These routines use a semi-analytical approach to treat the radiation from the RIAF (also called sometimes advection-dominated accretion flows, ADAFs) in which the accretion flow is considered stationary assuming an α-viscosity and a pseudo-Newtonian gravity, and the radiative transfer is treated in considerable detail, taking into account synchrotron, inverse Compton scattering and bremsstrahlung processes as appropriate for hot plasmas (e.g. [Nemmen et al. 2006](https://iopscience.iop.org/article/10.1086/500571); [Nemmen et al. 2014](https://academic.oup.com/mnras/article/438/4/2804/2907740)).

Our model for the RIAF emission is described in [Nemmen et al. (2014)](https://academic.oup.com/mnras/article/438/4/2804/2907740). RIAFs are usually characterized by the presence of outflows or winds, which prevent a considerable fraction of the gas that is available at large radii from being accreted onto the black hole (see [Yuan & Narayan 2014](https://www.annualreviews.org/doi/10.1146/annurev-astro-082812-141003) for a review). In order to take this mass-loss into account, the radial variation of the accretion rate is parameterized as $\dot{M}(R) = \dot{M}_{\rm o} \left( R/R_{\rm o} \right)^{s}$ 
(or $\rho(R) \propto R^{-3/2+s}$) 
where $\dot{M}_{\rm o}$ is the rate measured at the outer radius $R_{\rm o}$ of the RIAF (Blandford & Begelman 1999). 

The parameters of the model are:

| Parameter | Description |
|:--|:--|
| `s` | power-law index for accretion rate (or density) radial variation |
| $\dot{M}_o$ | mass accretion rate at the outer radius |
| `R_o` | outer radius |
| `M` | black hole mass |
| `alpha` | Shakura-Sunyaev viscosity parameter  |
| `beta` | ratio between the gas and total pressures |
| `delta` | fraction of energy dissipated via turbulence that directly heats electrons |
| `gamma` | adiabatic index |

The units of the parameters are described in the input parameter files included in the `examples` folder. 

# Usage

## Model setup

1. Build the Fortran binaries with `make build`
2. `cd` to the directory that will contain the SED
3. Prepare the parameter file `model.dat` (you can start from the one in `examples/largeR.dat`)
4. Run the Perl wrappers directly from this repository; `perl/run_model.pl` (or `perl/dyn.pl`, `perl/spectrum.pl` and `perl/ssd.pl`)

## Compute SED

For the usual workflow, use the one-command runner:

```shell
cd /path/model-run
/path/perl/run_model.pl model.dat
``` 

This will compute the dynamics, checks that the solution is physical, compute the spectrum, and write a SED plot named after the parameter file basename, such as `model.png`, in the current working directory.

### Manual runs

1. run `perl/dyn.pl` to compute ADAF dynamics to find physical global solution, adjusting range of eigenvalues `sl0i`,`sl0f` if required
2. once you get a good (physical) global solution in step 1, run `perl/spectrum.pl` to generate ADAF SED
3. optional: run `perl/ssd.pl` to compute truncated thin disk SED

You can pass a parameter file explicitly:

```shell
cd /path/model-run
perl /path/perl/dyn.pl model.dat
perl /path/perl/spectrum.pl model.dat
perl /path/perl/ssd.pl model.dat
```

With no parameter-file argument, the scripts read `in.dat` from the working directory.

The scripts above are meant to be invoked from the working directory that contains the selected parameter file and, for `spectrum.pl`, the matching `x.dat` produced by `dyn.pl`.

## Parameter files: where to start

Two examples of input parameter files are included in the `examples` folder:

- `largeR.dat`: ADAF with *R_out=1E4 Rs*
- `smallR.dat`: *R_out=500 Rs*

In order to compute the corresponding models, either copy one to `in.dat` or pass the file explicitly, for example `perl perl/run_model.pl examples/largeR.dat`.

For more information and useful references, refer to [this document](./docs/advice.md). 

# Smoke tests

Test suite can be run with:

    make smoke

This command rebuilds the Fortran binaries, checks representative "nice" and "bad" dynamical solutions from `tests/testcases_for_code/`, verifies the bundled `tests/n1097/` and `tests/m81/` example outputs are still parseable and physically plausible, and runs the `examples/largeR.dat` dynamics and spectrum workflow. It compares all 15 radial-profile columns against `tests/reference/largeR_dyn.out`, compares the generated spectrum against `tests/reference/largeR_spectrum.out`, and writes diagnostic plots to `tests/artifacts/largeR_dynamics.png` and `tests/artifacts/largeR_spectrum.png`.


# Citation

You are morally obligated to cite the following papers in any scientific literature that results from use of any part of this code:

1. [Yuan, F.; Cui, W. & Narayan, R. An Accretion-Jet Model for Black Hole Binaries: Interpreting the Spectral and Timing Features of XTE J1118+480. ApJ, 2005 , 620 , 905](https://iopscience.iop.org/article/10.1086/427206)
2. [Nemmen, R. S.; Storchi-Bergmann, T. & Eracleous, M. Spectral models for low-luminosity active galactic nuclei in LINERs: the role of advection-dominated accretion and jets MNRAS, 2014 , 438 , 2804](http://mnras.oxfordjournals.org/content/438/4/2804)
3. [Yuan, F.; Zdziarski, A. A.; Xue, Y.; Wu, X. Modeling the Hard States of XTE J1550-564 during Its 2000 Outburst. ApJ, 2007, 659, 541](https://ui.adsabs.harvard.edu/abs/2007ApJ...659..541Y/abstract)



# TODO 

- [x] include general-purpose user scripts for plotting SEDs
- [ ] add nonthermal emission
- [x] improve parallelization, radiative transfer