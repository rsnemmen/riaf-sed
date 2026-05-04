# Repository Guidelines

## Project Structure & Module Organization

This repository computes ADAF/RIAF dynamics and spectra using Fortran solvers with Perl wrappers.

- `fortran/`: numerical solvers. `dynamics.f` computes radial flow structure; `spectrum.f` computes emission; `ssd_new.f` and `ssd_alone.f` handle thin-disk spectra.
- `perl/`: user-facing wrappers and workflow scripts. Shared Perl modules live in `perl/lib/ADAF/`.
- `data/`: canonical inverse-Compton lookup tables staged automatically by `spectrum.pl`.
- `examples/`: tracked parameter files, including `largeR.dat` and `smallR.dat`.
- `tests/`: smoke tests, reference dynamics profiles, saved example outputs, and comparison utilities.
- `bin/`: generated Fortran binaries from `make build`; do not commit generated binaries.

## Build, Test, and Development Commands

- `make build`: compile all Fortran executables into `bin/`.
- `make smoke`: build binaries and run lightweight smoke/regression checks, including the `largeR` dynamics profile comparison.
- `make clean`: remove compiled binaries and Fortran build artifacts.
- `make clean-data`: remove top-level generated data products from test/model runs.
- Example run: `perl perl/dyn.pl examples/largeR.dat`.

Run Perl wrappers from a model working directory when using local `in.dat`; otherwise pass an explicit parameter file path.

## Coding Style & Naming Conventions

Keep changes consistent with the existing language style. Fortran sources use legacy fixed-form style; avoid broad formatting churn. Perl scripts use `strict`/`warnings` in newer code and shared helpers under `perl/lib/ADAF/`. Python test utilities should be standard-library only unless a dependency is clearly justified.

Use descriptive filenames matching existing patterns: `out_*` for dynamics outputs, `spec_*` for spectra, and `*_ssd` for thin-disk products.

## Testing Guidelines

The primary test entry point is `make smoke`. It validates good/bad dynamics fixtures, bundled example spectra, and compares `perl/dyn.pl examples/largeR.dat` against `tests/reference/largeR_dyn.out` using `tests/compare_dynamics.py`.

When changing dynamics behavior, regenerate the fiducial only if the numerical change is intentional and scientifically reviewed. Preserve diagnostic failures that catch nonphysical solutions.

## Commit & Pull Request Guidelines

Recent commits use concise subjects such as `Add dynamics profile regression test` or `Allow custom ADAF parameter files`. Prefer short imperative messages and keep each commit focused.

Pull requests should describe the scientific or workflow motivation, list commands run (`make smoke` at minimum), and call out changes to reference outputs or numerical tolerances. Include plots or profile summaries when behavior changes affect spectra or radial structure.

## Configuration Notes

Required tools are `gfortran`, Perl modules `Math::Derivative` and `Chart::Gnuplot`, and Python 3 for regression comparisons. Gnuplot is optional for diagnostic plotting.
