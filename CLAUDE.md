# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ADAF (Advection-Dominated Accretion Flow) SED calculator — computes spectral energy distributions of radiatively inefficient accretion flows (RIAFs) around black holes. Used for modeling low-luminosity AGN, X-ray binaries, and similar astrophysical sources.

**Primary authors:** Feng Yuan (Fortran solvers), Rodrigo Nemmen (Perl wrappers, eigenvalue solver, OpenMP)

## Build

```bash
cd fortran
make          # builds all four executables: dynamics, spectrum, ssd_new, ssd_alone
make clean    # removes binaries
```

Individual targets: `make dynamics`, `make spectrum`, `make ssd`, `make ssd_alone`

- Compiler: `gfortran` with `-O` optimization
- `dynamics` uses `-ffpe-trap=invalid,zero` (prevents hangs on NaN)
- `spectrum` uses `-fopenmp` (parallelizes over frequency bins, ~ncores/2 speedup)

## Typical Workflow

All runs happen in a working directory containing `in.dat` and the interpolation tables (`aomi-*.dat`, `romi-*.dat`). Adjust paths to executables in `dyn.pl`, `spectrum.pl`, and `ssd.pl` first.

1. Edit `in.dat` with model parameters
2. Find physical global solution: `perl /path/to/perl/dyn.pl`
3. Compute SED: `perl /path/to/perl/spectrum.pl`
4. Optional thin disk SED: `perl /path/to/perl/ssd.pl`

Use `examples/largeR.dat` or `examples/smallR.dat` as `in.dat` templates.

## Architecture

The code uses a shooting method (boundary value problem) to find the global dynamical solution of the ADAF equations, then computes the emergent spectrum.

**Data flow:**
```
in.dat → dyn.pl → dynamics (Fortran) → x.dat (radial structure)
                                              ↓
                             spectrum.pl → spectrum (Fortran) → spectrum.dat
                             ssd.pl      → ssd_new  (Fortran) → *_ssd
```

**Key components:**

| File | Role |
|------|------|
| `fortran/dynamics.f` | Solves ADAF radial ODE structure (Runge-Kutta, 4 unknowns: R, v_R, T_e, T_i) |
| `fortran/spectrum.f` | Radiative transfer: synchrotron, inverse Compton (2D lookup tables), bremsstrahlung |
| `fortran/ssd_new.f` | Thin disk spectrum with ADAF illumination |
| `perl/dyn.pl` | Eigenvalue solver — searches `sl0` range, detects physical solutions automatically |
| `perl/spectrum.pl` | Wrapper; generates two outputs (with/without Compton) |
| `perl/dyntype.pl` | Manual eigenvalue inspector with gnuplot diagnostics |
| `perl/adaf.pl` | Single-model runner with diagnostic plots |

**Interpolation tables** (`aomi-*.dat`, `romi-*.dat`): Pre-computed lookup tables for inverse Compton calculations. Must be present in working directory alongside `in.dat`.

## Key Parameters (`in.dat`)

| Parameter | Description |
|-----------|-------------|
| `gamma` | Adiabatic index |
| `m` | Black hole mass (units of 10^6 M_sun) |
| `beta` | Gas/total pressure ratio |
| `alfa` | Shakura-Sunyaev viscosity α |
| `delta` | Fraction of turbulent dissipation heating electrons |
| `mdot` | Accretion rate at outer radius (Eddington units) |
| `rout` | Outer radius (Schwarzschild radii) |
| `pp0` | Outflow strength (`s` parameter) |
| `ti`, `te` | Outer boundary temperatures for ions/electrons (K) |
| `vcs` | v_R/c_s ratio at outer boundary |
| `sl0i`, `sl0f` | Eigenvalue search range |
| `nmodels` | Number of eigenvalue samples |

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

For parameter space exploration, use "dumb parallelization": create separate run directories (`run01/`, `run02/`, etc.), each with its own `in.dat`, and run `dyn.pl` simultaneously in separate terminal tabs.

## Dependencies

- `gfortran` with OpenMP support
- Perl modules: `Math::Derivative`, `Chart::Gnuplot`
- Optional: `gnuplot` for diagnostic plots

Install Perl modules: `sudo cpan Math::Derivative && sudo cpan Chart::Gnuplot`

## Test Cases

`tests/testcases_for_code/` contains reference output files:
- `out_nice*.dat` — physically valid solutions
- `out_bad*.dat` — non-physical (for comparison)
- `out_discont*.dat` — discontinuous/invalid solutions

`tests/n1097/` and `tests/m81/` contain complete example model runs.

## Useful Output Variables (in `x.dat` / log file)

`slk`: Keplerian specific angular momentum | `ssll`: ADAF specific angular momentum | `tao`: vertical optical depth | `bb`: magnetic field | `qn`: cooling rate per unit volume
