# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ADAF (Advection-Dominated Accretion Flow) SED calculator — computes spectral energy distributions of radiatively inefficient accretion flows (RIAFs) around black holes. Used for modeling low-luminosity AGN, X-ray binaries, and similar astrophysical sources.

**Primary authors:** Feng Yuan (Fortran solvers), Rodrigo Nemmen (Perl wrappers, eigenvalue solver, OpenMP)

## Build

```bash
make build    # builds all four executables: dynamics, spectrum, ssd_new, ssd_alone
make clean    # removes binaries
make smoke    # runs lightweight smoke/regression checks
```

The top-level `Makefile` delegates to `fortran/Makefile`, which outputs binaries to `bin/` at the repo root. Individual targets are still available via `make -C fortran`.

- Compiler: `gfortran` with `-O` optimization
- `dynamics` uses `-ffpe-trap=invalid,zero` (prevents hangs on NaN)
- `spectrum` uses `-fopenmp` (parallelizes over frequency bins, ~ncores/2 speedup)

## Typical Workflow

All runs happen in a working directory containing a TOML parameter file. With no argument, the main wrappers read `in.toml`; they can also take one positional parameter file such as `model.toml`. Build the Fortran binaries first with `make build`; the Perl wrappers resolve binaries from `bin/` and automatically symlink the IC lookup tables from `data/` into the working directory — no manual copying required.

1. Edit `in.toml` with model parameters, or prepare `model.toml`
2. Find physical global solution: `perl /path/to/perl/dyn.pl` or `perl /path/to/perl/dyn.pl model.toml`
3. Compute SED: `perl /path/to/perl/spectrum.pl` or `perl /path/to/perl/spectrum.pl model.toml`
4. Optional thin disk SED: `perl /path/to/perl/ssd.pl` or `perl /path/to/perl/ssd.pl model.toml`

Use `examples/largeR.toml` or `examples/smallR.toml` as templates, either copied to `in.toml` or passed explicitly.

## Architecture

The code uses a shooting method (boundary value problem) to find the global dynamical solution of the ADAF equations, then computes the emergent spectrum.

**Data flow:**
```
in.toml or model.toml → dyn.pl → dynamics (Fortran) → x.dat (radial structure)
                                                        ↓
                                       spectrum.pl → spectrum (Fortran) → spectrum.dat
                                       ssd.pl      → ssd_alone (Fortran) → *_ssd
```

**Key components:**

| File | Role |
|------|------|
| `fortran/dynamics.f` | Solves ADAF radial ODE structure (Runge-Kutta, 4 unknowns: R, v_R, T_e, T_i) |
| `fortran/spectrum.f` | Radiative transfer: synchrotron, inverse Compton (2D lookup tables), bremsstrahlung |
| `fortran/ssd_new.f` | Thin disk spectrum with ADAF illumination |
| `fortran/ssd_alone.f` | Standalone thin disk spectrum used by `perl/ssd.pl` |
| `perl/dyn.pl` | Eigenvalue solver — searches `sl0` range, detects physical solutions automatically |
| `perl/spectrum.pl` | Wrapper; generates two outputs (with/without Compton) |
| `perl/ssd.pl` | Wrapper for the standalone thin disk spectrum executable |
| `perl/dyntype.pl` | **Debug only** — interactive eigenvalue inspector; prompts for eigenvalues one at a time and plots v_R/c_s + derivatives via gnuplot. Requires updating the hard-coded binary path before use. |
| `perl/adaf.pl` | **Debug only** — single-shot manual runner with hard-coded parameters and gnuplot inspection of v_R/c_s. Requires updating the hard-coded binary path before use. |
| `perl/lib/ADAF/Diagnostics.pm` | Solution classifier — parses output, detects sonic point, discontinuities, NaN |
| `perl/lib/ADAF/Paths.pm` | Resolves Fortran binary paths (`bin/`) and stages data files from `data/` |

**Interpolation tables** (`data/aomi-*.dat`, `data/romi-*.dat`): Pre-computed lookup tables for inverse Compton calculations. Canonical copies live in `data/`; `spectrum.pl` symlinks them into the working directory automatically.

## Key Parameters (`in.toml`)

| Parameter | Description |
|-----------|-------------|
| `dynamics.gamma` | Adiabatic index |
| `dynamics.mass` | Black hole mass (units of 10^6 M_sun) |
| `dynamics.beta` | Gas/total pressure ratio |
| `dynamics.alpha` | Shakura-Sunyaev viscosity α |
| `dynamics.delta` | Fraction of turbulent dissipation heating electrons |
| `dynamics.mdot_out` | Accretion rate at outer radius (Eddington units) |
| `dynamics.r_out` | Outer radius (Schwarzschild radii) |
| `dynamics.p_wind` | Outflow strength (`s` parameter) |
| `boundary.ti`, `boundary.te` | Outer boundary temperatures for ions/electrons |
| `boundary.mach` | v_R/c_s ratio at outer boundary |
| `shooting.sl0_initial`, `shooting.sl0_final` | Eigenvalue search range |
| `shooting.n_models` | Number of eigenvalue samples |

## Finding Physical Solutions

The eigenvalue (`sl0`) is typically between 1 and 3. Physical solutions satisfy:
- Mach number > 1 at inner region (transonic/sonic point)
- `v_R(R)` is smooth (no discontinuities or spikes)
- `v_R` decreases as R increases (outward)

The code "hangs" after passing the correct eigenvalue — this indicates the physical solution bracket. `dyn.pl` automates this by bracketing between last valid solution and the freezing one.

If many `NaN`/`FAILED!` results appear, try slightly adjusting boundary conditions (`ti`, `te`, `vcs`).

## Boundary Conditions Guidelines

- **Large Rout (≥ few thousand R_g):** `vcs=0.2–0.3`, T_i ≈ 0.1–0.3 T_virial, T_e slightly below T_i
  Example: rout=10⁴ R_g → T_i=6e7 K, T_e=5.9e7 K, vcs=0.2

- **Rout ~ 100 R_g:** `vcs=0.5`, larger T_i and T_e
  Example: Ṁ=0.1 Ṁ_Edd, rout=100 R_g → T_i=15e9 K, T_e=8e9 K, vcs=0.5

## Parallelization Strategy

For parameter space exploration, use "dumb parallelization": create separate run directories (`run01/`, `run02/`, etc.), each with its own `in.toml` or custom parameter file, and run `dyn.pl` simultaneously in separate terminal tabs.

## Dependencies

- `gfortran` with OpenMP support
- Perl modules: `Math::Derivative`, `Chart::Gnuplot`
- Python 3 with Matplotlib for dynamics/spectrum regression comparisons and plots
- Optional: `gnuplot` for diagnostic plots

Install Perl modules: `sudo cpan Math::Derivative && sudo cpan Chart::Gnuplot`

## Test Cases

`tests/testcases_for_code/` contains reference output files:
- `out_nice*.dat` — physically valid solutions
- `out_bad*.dat` — non-physical (for comparison)
- `out_discont*.dat` — discontinuous/invalid solutions

`tests/n1097/` and `tests/m81/` contain complete example model runs.

`tests/reference/largeR_dyn.out` is the fiducial dynamics profile generated from `examples/largeR.toml`. `tests/compare_dynamics.py` compares all 15 radial-profile columns, including Mach number and `log(rho)`. `tests/plot_dynamics.py` writes a gridded inspection plot to `tests/artifacts/largeR_dynamics.png`.

`tests/reference/largeR_spectrum.out` is the fiducial spectrum generated by running `perl/dyn.pl examples/largeR.toml` followed by `perl/spectrum.pl examples/largeR.toml` in the same working directory. `tests/compare_spectrum.py` compares the log-frequency grid and `log10(nu Lnu)` values with tight log-space tolerances. `tests/plot_spectrum.py` writes `tests/artifacts/largeR_spectrum.png`, with the spectrum overlay on top and a residual panel below.

Run the lightweight smoke/regression checks with `make smoke`. This rebuilds the Fortran binaries, checks representative good/bad dynamical fixtures, validates bundled example spectra, runs the `largeR` dynamics and spectrum workflow, compares both outputs with their fiducial solutions, and refreshes the dynamics and spectrum plot artifacts.

## Useful Output Variables (in `x.dat` / log file)

`slk`: Keplerian specific angular momentum | `ssll`: ADAF specific angular momentum | `tao`: vertical optical depth | `bb`: magnetic field | `qn`: cooling rate per unit volume
